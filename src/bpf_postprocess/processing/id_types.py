"""BioPharma Finder RAW-specific identification types, with missing data distinguished.

Original BF CSV has a global ID Type (Both/Full/MS2/NoID) and run-scoped
ID Type:<raw>. The run-scoped value takes precedence for per-RAW analysis.
Native PMF ID Type is not decoded by this project: do not infer from confidence.
"""
from __future__ import annotations

import pandas as pd

UNAVAILABLE_ID_TYPE = "Unavailable (not decoded)"


def displayed_id_types(values: pd.Series) -> pd.Series:
    """Return reportable ID Types without silently conflating missing and NoID."""
    clean = values.fillna("").astype(str).str.strip()
    return clean.where(clean.ne(""), UNAVAILABLE_ID_TYPE)


def combined_id_types(values: pd.Series) -> str:
    """Preserve every observed category in an aggregate (no forced 'majority' ID)."""
    return "; ".join(sorted(set(displayed_id_types(values).tolist())))


def id_type_summary(frame: pd.DataFrame) -> pd.DataFrame:
    """RAW-run level distribution after all filters/component/peptide selection."""
    columns = ["experiment", "sample", "raw_file", "preparation", "id_type",
               "components", "identified_components", "ms_area_sum"]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    data = frame.copy()
    data["id_type"] = displayed_id_types(data["id_type"])
    result = data.groupby(columns[:5], dropna=False, as_index=False).agg(
        components=("component_no", "size"),
        identified_components=("identified", "sum"),
        ms_area_sum=("ms_area", lambda s: s.sum(min_count=1)),
    )
    return result[columns].sort_values(columns[:5]).reset_index(drop=True)
