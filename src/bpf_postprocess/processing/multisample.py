"""Multi-run comparisons and replicate summaries for peptide MS areas."""
from __future__ import annotations
import numpy as np
import pandas as pd
from bpf_postprocess.processing.comparison import long_peptide_areas, compare_runs
from bpf_postprocess.processing.id_types import combined_id_types


def reference_comparison(frame: pd.DataFrame, reference_run: str) -> pd.DataFrame:
    """All RAW runs against one reference. Nonpositive/missing baseline -> NA fold.

    This is raw MS signal comparison, NOT normalization of injection/digestion.
    """
    matrix = compare_runs(frame)
    if matrix.empty:
        return pd.DataFrame(columns=['feature', 'comparison_run', 'reference_run', 'reference_area', 'comparison_area', 'delta_area', 'fold_vs_reference', 'percent_change'])
    if reference_run not in matrix.columns:
        raise ValueError('Reference RAW not present in comparison matrix')
    run_cols = [c for c in matrix if c not in ('feature', reference_run)]
    rows=[]
    for run in run_cols:
        part = matrix[['feature', reference_run, run]].copy()
        part.columns = ['feature', 'reference_area', 'comparison_area']
        part.insert(1, 'comparison_run', run)
        part.insert(1, 'reference_run', reference_run)
        ref = pd.to_numeric(part['reference_area'], errors='coerce')
        comp = pd.to_numeric(part['comparison_area'], errors='coerce')
        part['delta_area'] = comp - ref
        part['fold_vs_reference'] = comp / ref.where(ref>0)
        part['percent_change'] = (part['fold_vs_reference']-1)*100
        rows.append(part)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['feature','reference_run','comparison_run','reference_area','comparison_area','delta_area','fold_vs_reference','percent_change'])


def replicate_statistics(frame: pd.DataFrame, by: str = 'sample') -> pd.DataFrame:
    """Mean/SD/CV of one aggregated peptide area per RAW injection.

    Replicates are deliberately not inferred from the number of components.
    Sample/group labels must be explicitly assigned by the user.
    """
    if by not in ('sample', 'sample_group'):
        raise ValueError('Group by sample or sample_group')
    areas = long_peptide_areas(frame)
    cols = [by, 'preparation', 'feature']
    if areas.empty:
        return pd.DataFrame(columns=cols + ['n_runs', 'mean_ms_area', 'sd_ms_area', 'rsd_pct', 'id_types'])
    if by == 'sample_group':
        # The group is a run-level label; avoid assuming it is a sample property.
        labels=frame[['experiment','raw_file','sample_group']].drop_duplicates()
        areas=areas.merge(labels, on=['experiment','raw_file'], how='left', validate='many_to_one')
    values=areas.groupby(cols, dropna=False, as_index=False).agg(
        n_runs=('ms_area', 'count'), mean_ms_area=('ms_area', 'mean'), sd_ms_area=('ms_area','std'),
        id_types=('id_types', lambda x: '; '.join(sorted(set('; '.join(x.fillna('')).split('; ')) - {''})))
    )
    values['rsd_pct'] = 100*values['sd_ms_area']/values['mean_ms_area'].where(values['mean_ms_area']>0)
    return values.sort_values(cols).reset_index(drop=True)
