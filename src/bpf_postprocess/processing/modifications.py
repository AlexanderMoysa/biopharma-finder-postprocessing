"""Exploratory modification signal summaries, not chemically validated occupancy."""
from __future__ import annotations
import pandas as pd


def modification_areas(data: pd.DataFrame) -> pd.DataFrame:
    """Sum area by literal BioPharma Finder modification annotation and RAW run.

    A blank modification is shown as 'Not annotated'. Native PMF currently only
    retains free-text Identification, so Mod is not available unless CSV supplied.
    Do NOT interpret as per-site occupancy or percentage of modified protein.
    """
    rows = data[data["identified"]].copy()
    rows["modification"] = rows["modification"].fillna("").str.strip().replace("", "Not annotated")
    return rows.groupby(["experiment", "sample", "raw_file", "preparation", "modification"],
                        dropna=False, as_index=False).agg(
        ms_area=("ms_area", "sum"), components=("component_no", "size")
    )
