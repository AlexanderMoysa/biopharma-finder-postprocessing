"""Display-only BF-style compaction of *already selected* peptide segments.

The compact map uses ONE representative bar for each (chain, first AA, last AA)
region, regardless of modified forms, adducts and charge states. It retains
all selected components in the numerical summary and detailed CSV/Excel.

Bar colors approximate BioPharma Finder's visual bins with an independently
specified MS-area proxy: the total selected area for a region divided by the
strongest selected region of THAT SAME CHAIN. This is not BF Peptide Recovery.
"""
from __future__ import annotations

import pandas as pd

BASE_COLUMNS = [
    "run", "experiment", "raw_file", "preparation", "protein", "chain_number",
    "start", "end", "representative_feature", "representative_modification",
    "rt_min", "regional_ms_area", "relative_signal_pct", "component_count",
    "identification_count", "metric_definition",
]


def compact_coverage_segments(segments: pd.DataFrame) -> pd.DataFrame:
    """BF-like one-bar-per-position view without reintroducing deselected rows.

    All input records MUST come from calculate_coverage(selected_rows, one_raw).
    Only display segments are merged. RT is taken from the highest-area selected
    identification in that exact interval. Missing areas remain excluded upstream.
    """
    if segments.empty:
        return pd.DataFrame(columns=BASE_COLUMNS)
    required = {"run", "experiment", "raw_file", "preparation", "protein",
                "chain_number", "start", "end", "feature", "modification",
                "rt_min", "ms_area", "component_count"}
    missing = required - set(segments.columns)
    if missing:
        raise ValueError(f"Coverage segments missing required fields: {sorted(missing)}")
    unique_runs = segments[["run", "experiment", "raw_file"]].drop_duplicates()
    if len(unique_runs) != 1:
        raise ValueError("Compact mapping requires exactly one RAW run")
    group_cols = ["run", "experiment", "raw_file", "preparation", "protein",
                  "chain_number", "start", "end"]
    frame = segments.copy()
    frame["ms_area"] = pd.to_numeric(frame["ms_area"], errors="coerce")
    frame["component_count"] = pd.to_numeric(frame["component_count"], errors="coerce")
    if (frame["ms_area"] < 0).any():
        raise ValueError("Selected peptide areas must be non-negative")
    sums = frame.groupby(group_cols, dropna=False, as_index=False).agg(
        regional_ms_area=("ms_area", lambda x: x.sum(min_count=1)),
        component_count=("component_count", "sum"),
        identification_count=("feature", "nunique"),
    )
    dominant = (frame.sort_values(["ms_area", "feature"],
                                  ascending=[False, True], na_position="last", kind="stable")
                .drop_duplicates(group_cols, keep="first"))
    representative = dominant[group_cols + ["feature", "modification", "rt_min"]].rename(
        columns={"feature": "representative_feature",
                 "modification": "representative_modification"})
    data = sums.merge(representative, on=group_cols, how="left", validate="one_to_one")
    for chain, indices in data.groupby("chain_number").groups.items():
        top = data.loc[indices, "regional_ms_area"].max()
        data.loc[indices, "relative_signal_pct"] = (
            data.loc[indices, "regional_ms_area"] / top * 100 if pd.notna(top) and top > 0 else 0.
        )
    data["metric_definition"] = (
        "Selected MS Area in this chain/interval ÷ strongest selected interval in SAME chain ×100; "
        "NOT Thermo peptide recovery"
    )
    return data[BASE_COLUMNS].sort_values(
        ["chain_number", "start", "end"], kind="stable").reset_index(drop=True)
