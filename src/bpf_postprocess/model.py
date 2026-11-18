"""One normalized row = one component in one RAW run.

Do not use Component No. as a cross-experiment identification key.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import pandas as pd

COLS = [
    "experiment", "sample", "raw_file", "condition", "preparation",
    "replicate", "component_no", "identified", "identification",
    "peptide_sequence", "modification", "site", "positions", "protein",
    "rt_min", "mz", "charge", "mono_mass_da", "ms_area",
    "confidence_pct", "delta_ppm", "id_type", "integration_type",
    "source_file", "source_type", "import_note",
]
NUMERIC_COLS = ["component_no", "rt_min", "mz", "charge", "mono_mass_da", "ms_area", "confidence_pct", "delta_ppm"]

@dataclass
class ImportResult:
    data: pd.DataFrame
    experiments: pd.DataFrame
    runs: pd.DataFrame
    warnings: list[str]


def as_table(records: list[dict] | pd.DataFrame) -> pd.DataFrame:
    table = pd.DataFrame(records).reindex(columns=COLS)
    for col in NUMERIC_COLS:
        table[col] = pd.to_numeric(table[col], errors="coerce")
    table["identified"] = table["identified"].fillna(False).astype(bool)
    for col in ("identification", "peptide_sequence", "modification", "site", "positions", "protein"):
        table[col] = table[col].fillna("").astype(str)
    return table


def derive_sample_name(experiment: str) -> str:
    for suffix in ("_combined", "-combined", " combined"):
        if experiment.lower().endswith(suffix):
            return experiment[:-len(suffix)]
    return experiment


def guess_preparation(raw_file: str) -> str:
    """Only explicit name evidence; user can edit labels in the app."""
    import re
    name = Path(raw_file).name.lower()
    if re.search(r"non[-_ ]?reduced|non[-_ ]?red|unreduced", name):
        return "Non-reduced"
    if re.search(r"with[-_ ]?dtt|(?:^|[-_ ])dtt(?:[-_.]|$)|(?<!non)reduced", name):
        return "Reduced"
    return "Unspecified"


def make_runs(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["experiment", "sample", "raw_file", "condition", "preparation", "replicate", "source_file", "source_type"]
    if df.empty:
        return pd.DataFrame(columns=cols + ["components", "identified_components"])
    grouped = df.groupby(cols, dropna=False, sort=False).agg(
        components=("component_no", "size"),
        identified_components=("identified", "sum"),
    ).reset_index()
    return grouped
