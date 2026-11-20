"""Pure, reusable component filters; UI-independent."""
from __future__ import annotations
import pandas as pd


def filter_components(
    frame: pd.DataFrame, *, samples: list[str] | None = None,
    raw_files: list[str] | None = None, preparations: list[str] | None = None,
    identified: str = "All", peptide_query: str = "", modification_query: str = "",
    min_area: float = 0, min_confidence_pct: float = 0,
) -> pd.DataFrame:
    data = frame
    if samples is not None:
        data = data[data["sample"].isin(samples)]
    if raw_files is not None:
        data = data[data["raw_file"].isin(raw_files)]
    if preparations is not None:
        data = data[data["preparation"].isin(preparations)]
    if identified == "Identified":
        data = data[data["identified"]]
    if identified == "Unidentified":
        data = data[~data["identified"]]
    if peptide_query.strip():
        query = peptide_query.strip()
        found = (
            data["identification"].str.contains(query, case=False, regex=False, na=False) |
            data["peptide_sequence"].str.contains(query, case=False, regex=False, na=False)
        )
        data = data[found]
    if modification_query.strip():
        query = modification_query.strip()
        found = (
            data["modification"].str.contains(query, case=False, regex=False, na=False) |
            data["identification"].str.contains(query, case=False, regex=False, na=False)
        )
        data = data[found]
    if min_confidence_pct > 0:
        # NaN (unidentified/missing) never passes a positive confidence filter.
        data = data[data["confidence_pct"] >= min_confidence_pct]
    if min_area > 0:
        data = data[data["ms_area"] >= min_area]
    return data.copy()
