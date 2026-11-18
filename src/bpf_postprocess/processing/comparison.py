"""Area comparisons. No unknown-component matching by component number."""
from __future__ import annotations
import numpy as np
import pandas as pd


def feature_key(df: pd.DataFrame) -> pd.Series:
    """Match on the full BioPharma Finder identification annotation.

    The same ID text is available in both CSV and PMF. Detailed CSV-only
    sequence/modification fields are retained for filters but NOT used to
    build incompatible cross-format comparison keys. This assumes the same
    protein and analysis method: the caller must verify comparability.
    """
    return df["identification"].fillna("").astype(str).str.strip()


def long_peptide_areas(df: pd.DataFrame) -> pd.DataFrame:
    """Sum multiple component peaks for the same peptide per RAW run.

    Explicitly named aggregation: it does NOT estimate modification occupancy.
    """
    identified = df[df["identified"]].copy()
    identified["feature"] = feature_key(identified)
    identified = identified[identified["feature"].str.strip().ne("")]
    grouped = identified.groupby(
        ["experiment", "sample", "raw_file", "preparation", "feature"],
        as_index=False, dropna=False,
    ).agg(ms_area=("ms_area", "sum"), component_count=("component_no", "size"))
    return grouped


def compare_runs(df: pd.DataFrame, run_ids: list[str] | None = None) -> pd.DataFrame:
    """Matrix: peptides x RAW runs (keyed by experiment and RAW), no auto-ratios."""
    data = long_peptide_areas(df)
    data["run"] = data["experiment"] + " / " + data["raw_file"]
    if run_ids is not None:
        data = data[data["run"].isin(run_ids)]
    if data.empty:
        return pd.DataFrame(columns=["feature"])
    matrix = data.pivot_table(index="feature", columns="run", values="ms_area", aggfunc="sum", observed=True)
    return matrix.sort_index().reset_index()


def ratio_table(matrix: pd.DataFrame, control: str, treatment: str) -> pd.DataFrame:
    """Add treatment-control and fold ratio, null for zero/missing reference."""
    if control not in matrix.columns or treatment not in matrix.columns:
        raise ValueError("Both run labels must be columns in the comparison matrix")
    result = matrix[["feature", control, treatment]].copy()
    numerator = pd.to_numeric(result[treatment], errors="coerce")
    denominator = pd.to_numeric(result[control], errors="coerce")
    result["delta_area"] = numerator - denominator
    safe_denominator = denominator.where(denominator > 0)
    result["fold_vs_reference"] = numerator.div(safe_denominator)
    result["log2_fold"] = np.log2(result["fold_vs_reference"].where(result["fold_vs_reference"] > 0))
    return result
