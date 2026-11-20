"""Exploratory modified-to-unmodified peptide MS area fractions.

These are signal fractions, NOT absolute molar site occupancy. The denominator
contains only unmodified forms and explicitly chosen PTM types; unrelated
adducts, gas-phase fragments, and nonspecific cleavage remain excluded.
"""
from __future__ import annotations
import pandas as pd

# Conservative default targets. The UI may let users choose another literal
# modification, but the output must still be described as an *MS area fraction*.
DEFAULT_PTMS = ('Deamidation', 'Oxidation', 'Double Oxidation')


def peptide_backbone(identification: str) -> str:
    """Position-resolved backbone (left side of BF identification)."""
    return str(identification).partition('=')[0].strip()


def modification_fractions(frame: pd.DataFrame, target_mods: list[str] | tuple[str, ...] = DEFAULT_PTMS) -> pd.DataFrame:
    """Per-run and per-backbone MS area fraction for each selected PTM.

    Missing unmodified reference, missing PTM, or 0 denominator are explicit:
    no occupancy is inferred when an unmodified form was not identified.
    """
    targets=[str(x).strip() for x in target_mods if str(x).strip() and x not in ('None','Not annotated')]
    cols=['experiment','sample','sample_group','raw_file','preparation','backbone','modification','unmodified_area','modified_area','denominator_area','modified_ms_area_pct','reference_present']
    if not targets:
        return pd.DataFrame(columns=cols)
    identified = frame[frame['identified']].copy()
    identified['backbone'] = identified['identification'].map(peptide_backbone)
    identified = identified[identified['backbone'].ne('')]
    identified['modification'] = identified['modification'].fillna('').astype(str).str.strip()
    # A multitype annotation cannot be assigned unambiguously to just one target.
    eligible=identified[identified['modification'].isin(['None','']+targets)].copy()
    if eligible.empty:
        return pd.DataFrame(columns=cols)
    keys=['experiment','sample','sample_group','raw_file','preparation','backbone']
    agg=eligible.groupby(keys+['modification'], dropna=False, as_index=False)['ms_area'].sum()
    unmod=agg[agg['modification'].isin(['None',''])].groupby(keys, dropna=False, as_index=False)['ms_area'].sum().rename(columns={'ms_area':'unmodified_area'})
    modified=agg[agg['modification'].isin(targets)].rename(columns={'ms_area':'modified_area'})
    if modified.empty:
        return pd.DataFrame(columns=cols)
    total_modified=modified.groupby(keys,dropna=False,as_index=False)['modified_area'].sum().rename(columns={'modified_area':'total_selected_ptm_area'})
    result=modified.merge(unmod,on=keys,how='left',validate='many_to_one').merge(total_modified,on=keys,how='left',validate='many_to_one')
    result['reference_present'] = result['unmodified_area'].notna()
    result['denominator_area'] = result['unmodified_area'] + result['total_selected_ptm_area']
    result['modified_ms_area_pct'] = (100*result['modified_area'] / result['denominator_area'].where(result['denominator_area']>0)).where(result['reference_present'])
    return result[cols].sort_values(keys+['modification']).reset_index(drop=True)
