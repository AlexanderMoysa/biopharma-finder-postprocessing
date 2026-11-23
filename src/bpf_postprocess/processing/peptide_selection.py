"""Shared selection of BF peptide identifications across analytical tabs.

A peptide is keyed by the *full identification*, matching the existing
cross-RAW comparison logic. Deselecting a peptide excludes that identification
in every RAW run, rather than one arbitrary dataframe row or component number.
The original component records stay intact for later re-selection.
"""
from __future__ import annotations

import pandas as pd
from bpf_postprocess.processing.id_types import combined_id_types


ANNOTATIONS = ("peptide_sequence", "modification", "site", "positions", "protein")


def peptide_ids(frame: pd.DataFrame) -> pd.Series:
    """Normalized full BioPharma Finder identification per input row."""
    return frame["identification"].fillna("").astype(str).str.strip()


def apply_peptide_exclusions(frame: pd.DataFrame, excluded: set[str]) -> pd.DataFrame:
    """Remove excluded identified peptide IDs, keep unidentified signals."""
    if not excluded or frame.empty:
        return frame.copy()
    mask = frame["identified"] & peptide_ids(frame).isin(excluded)
    return frame.loc[~mask].copy()


def filter_peptide_table(table: pd.DataFrame, excluded: set[str]) -> pd.DataFrame:
    """Filter a derived table with a full-BF-ID ``feature`` column.

    The *candidate* table used by a checkbox editor is retained separately,
    including deselected features, so any edit can be reversed.
    """
    if not excluded or table.empty or "feature" not in table:
        return table.copy()
    return table.loc[~table["feature"].fillna("").astype(str).str.strip().isin(excluded)].copy()


def update_peptide_exclusions(
    features: pd.Series, checked: pd.Series, excluded: set[str]
) -> set[str]:
    """Update only visible peptide IDs; preserve hidden/existing selections.

    Each displayed peptide must have one checkbox. This avoids ambiguous edits
    to duplicate feature rows (e.g. one feature × several reference comparisons).
    """
    if len(features) != len(checked):
        raise ValueError("Peptide selection length mismatch")
    labels = features.fillna("").astype(str).str.strip().tolist()
    if any(not label for label in labels):
        raise ValueError("Peptide selection cannot use blank identification")
    if len(labels) != len(set(labels)):
        raise ValueError("Peptide selector requires exactly one row per identification")
    if checked.isna().any():
        raise ValueError("Peptide selection cannot contain missing checkboxes")
    result = set(excluded)
    for feature, is_checked in zip(labels, checked.tolist()):
        if bool(is_checked):
            result.discard(feature)
        else:
            result.add(feature)
    return result


def peptide_inventory(frame: pd.DataFrame) -> pd.DataFrame:
    """Per-peptide metadata and source-component totals for checkbox editors.

    This table is assembled from component-selected data *before* peptide-level
    exclusions, so previously deselected peptides remain visible and reversible.
    No cross-peptide occupancy claim is made from MS Area sums.
    """
    columns = ["feature", *ANNOTATIONS, "id_types", "raw_runs", "components", "total_ms_area"]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    rows = frame.loc[frame["identified"]].copy()
    rows["feature"] = peptide_ids(rows)
    rows = rows.loc[rows["feature"].ne("")].copy()
    if rows.empty:
        return pd.DataFrame(columns=columns)
    rows["run_key"] = rows["experiment"].astype(str) + " / " + rows["raw_file"].astype(str)
    def joined(values: pd.Series) -> str:
        return "; ".join(dict.fromkeys(
            value for value in values.fillna("").astype(str).str.strip() if value
        ))
    result = rows.groupby("feature", as_index=False, sort=True).agg(
        **{column: (column, joined) for column in ANNOTATIONS},
        id_types=("id_type", combined_id_types),
        raw_runs=("run_key", "nunique"),
        components=("component_no", "size"),
        total_ms_area=("ms_area", lambda s: s.sum(min_count=1)),
    )
    return result[columns]
