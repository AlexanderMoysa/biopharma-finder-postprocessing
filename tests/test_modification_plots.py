"""Regression tests for plotting the PTM percentages (not peptide areas)."""
import pandas as pd
import pytest

from bpf_postprocess.ui.plots import modification_fraction_figure


def test_plot_uses_modified_percentages():
    fractions = pd.DataFrame({
        "raw_file": ["control.raw", "stress.raw", "stress.raw"],
        "sample": ["Control", "Stress", "Stress"],
        "preparation": ["Reduced", "Reduced", "Reduced"],
        "backbone": ["1:N18-G21"] * 3,
        "modification": ["Deamidation", "Deamidation", "Oxidation"],
        "modified_ms_area_pct": [1.5, 20.0, 70.0],
    })
    fig = modification_fraction_figure(fractions, "Deamidation")
    assert fig is not None
    assert sorted(float(v) for trace in fig.data for v in trace.y) == [1.5, 20.0]
    assert all(trace.name == "1:N18-G21" for trace in fig.data)
    assert fig.layout.yaxis.title.text == "Modified MS Area (%)"


def test_plot_skips_missing_reference_without_inventing_values():
    fractions = pd.DataFrame({
        "raw_file": ["a.raw"], "sample": ["A"], "preparation": ["Non-reduced"],
        "backbone": ["1:N18-G21"], "modification": ["Deamidation"],
        "modified_ms_area_pct": [float("nan")],
    })
    assert modification_fraction_figure(fractions, "Deamidation") is None
    assert modification_fraction_figure(fractions, "Oxidation") is None


def test_plot_rejects_wrong_dataframe_schema():
    peptide_areas = pd.DataFrame({
        "experiment": ["4197_combined"], "raw_file": ["a.raw"],
        "feature": ["1:N18-G21"], "ms_area": [1000.0],
    })
    with pytest.raises(ValueError, match="modified_ms_area_pct"):
        modification_fraction_figure(peptide_areas, "Deamidation")
