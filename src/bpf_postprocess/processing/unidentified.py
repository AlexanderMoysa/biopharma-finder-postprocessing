"""Unidentified components: filtering and runwise summary only.

Future cross-run mass/RT alignment must be validated before adding it.
"""
from __future__ import annotations
import pandas as pd
from bpf_postprocess.processing.id_types import combined_id_types


def unidentified_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = data[~data["identified"]].copy()
    return rows.groupby(["experiment", "sample", "raw_file", "preparation"],
                        dropna=False, as_index=False).agg(
        components=("component_no", "size"), nonzero_peaks=("ms_area", lambda s: int((s > 0).sum())),
        area_sum=("ms_area", "sum"), id_types=("id_type", combined_id_types),
    )
