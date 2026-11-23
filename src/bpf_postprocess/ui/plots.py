"""Plotly chart builders for the Streamlit frontend.

Keep plotting inputs explicit: tables from another tab must never be reused
by accident. Calculations remain in the processing package.
"""
from __future__ import annotations

import pandas as pd
import plotly.express as px


def modification_fraction_figure(fractions: pd.DataFrame, modification: str):
    """Build a chart of quantified modified-peptide MS-area percentages.

    Only rows with a defined area fraction are plotted. Missing unmodified
    reference data yield NaN, not an inferred zero or 100% fraction.
    Returns None when there are no quantifiable rows for this modification.
    """
    required = {
        "raw_file", "sample", "preparation", "backbone",
        "modification", "modified_ms_area_pct",
    }
    missing = required.difference(fractions.columns)
    if missing:
        raise ValueError(f"Modified-fraction plot missing columns: {sorted(missing)}")

    rows = fractions.loc[fractions["modification"] == modification].copy()
    rows["modified_ms_area_pct"] = pd.to_numeric(
        rows["modified_ms_area_pct"], errors="coerce"
    )
    rows = rows[rows["modified_ms_area_pct"].notna()]
    if rows.empty:
        return None

    return px.bar(
        rows,
        x="raw_file",
        y="modified_ms_area_pct",
        color="backbone",
        barmode="group",
        hover_data=["sample", "preparation"],
        labels={"modified_ms_area_pct": "Modified MS Area (%)"},
    )
