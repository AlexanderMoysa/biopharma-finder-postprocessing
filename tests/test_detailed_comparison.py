"""Regression coverage for explicit, preparation-safe run comparisons."""
import math

import pandas as pd
import pytest

from bpf_postprocess.model import as_table
from bpf_postprocess.processing.detailed_comparison import (
    annotated_area_matrix,
    compatible_run_ids,
    detailed_component_records,
    detailed_reference_comparison,
    run_catalog,
    select_comparison_runs,
)


def _sample_data():
    rows = []
    # Nonreduced reference has two charge states for one peptide: 100+50 area.
    for experiment, raw, prep, areas in (
        ("4197_combined", "4197_NR.raw", "Non-reduced", [(100, 5, 2, 10.1), (50, 6, 3, 10.6)]),
        ("4197_combined", "4197_DTT.raw", "Reduced", [(900, 5, 2, 10.2)]),
        ("4204_combined", "4204_NR.raw", "Non-reduced", [(120, 5, 2, 10.4)]),
        ("4204_combined", "4204_DTT.raw", "Reduced", [(400, 5, 2, 10.9)]),
        ("4300_combined", "4300_NR.raw", "Non-reduced", [(160, 5, 2, 10.5)]),
    ):
        for idx, (area, mz, z, rt) in enumerate(areas, start=1):
            rows.append(dict(experiment=experiment, raw_file=raw,
                             sample=experiment.split("_")[0], sample_group="DS",
                             preparation=prep, component_no=idx, identified=True,
                             identification="1:N18-G21 = 455m(N18+Deamidation)",
                             peptide_sequence="NYCG", modification="Deamidation", site="N18",
                             positions="18-21", protein="A", rt_min=rt, mz=mz, charge=z,
                             mono_mass_da=455.1, confidence_pct=99.0, ms_area=area,
                             delta_ppm=0.25, id_type="MS2", integration_type="Auto"))
        # same peptide from another run is not the same component number
    rows.append(dict(experiment="4197_combined", raw_file="4197_NR.raw", sample="4197",
                     preparation="Non-reduced", component_no=51, identified=False, ms_area=9999.0))
    return as_table(rows)


def test_run_picker_rejects_mixed_preparation_and_unselected_runs():
    data = _sample_data()
    catalog = run_catalog(data)
    reference = "4197_combined / 4197_NR.raw"
    allowed = compatible_run_ids(catalog, reference)
    assert allowed == ["4204_combined / 4204_NR.raw", "4300_combined / 4300_NR.raw"]
    chosen = select_comparison_runs(data, reference, [allowed[0]])
    assert set(chosen.raw_file) == {"4197_NR.raw", "4204_NR.raw"}
    assert "4197_DTT.raw" not in set(chosen.raw_file)
    with pytest.raises(ValueError, match="Reduced and non-reduced"):
        select_comparison_runs(data, reference, ["4204_combined / 4204_DTT.raw"])
    with pytest.raises(ValueError, match="its own comparator"):
        select_comparison_runs(data, reference, [reference])


def test_rich_reference_table_includes_annotations_and_run_specific_rt():
    data = _sample_data()
    ref = "4197_combined / 4197_NR.raw"
    cmp = "4204_combined / 4204_NR.raw"
    result = detailed_reference_comparison(data, ref, [cmp])
    assert len(result) == 1
    row = result.iloc[0]
    assert row.peptide_sequence == "NYCG"
    assert row.modification == "Deamidation"
    assert row.site == "N18"
    assert row.positions == "18-21"
    assert row.protein == "A"
    assert row.reference_area == 150.0
    assert row.comparison_area == 120.0
    assert math.isclose(row.fold_vs_reference, 0.8)
    assert row.percent_change == pytest.approx(-20)
    # Dom. peak is area 100 @ 10.1 min, not an average or a made-up 10.35 min.
    assert row.reference_rt_min == pytest.approx(10.1)
    assert row.comparison_rt_min == pytest.approx(10.4)
    assert row.reference_component_count == 2
    assert row.comparison_component_count == 1
    assert row.reference_charge == 2
    assert row.comparison_confidence_pct == 99


def test_matrix_preserves_each_run_rt_and_full_component_records():
    data = _sample_data()
    ref = "4197_combined / 4197_NR.raw"
    cmp = "4204_combined / 4204_NR.raw"
    matrix = annotated_area_matrix(data, ref, [cmp])
    assert len(matrix) == 1
    assert matrix.iloc[0][f"MS Area | {ref}"] == 150
    assert matrix.iloc[0][f"RT dominant (min) | {cmp}"] == pytest.approx(10.4)
    assert not any("DTT" in c for c in matrix.columns)
    records = detailed_component_records(data, ref, [cmp])
    assert len(records) == 3
    assert sorted(records.loc[records.raw_file == "4197_NR.raw", "rt_min"].tolist()) == [10.1, 10.6]
    assert set(records.role) == {"Reference", "Comparison"}
    assert "confidence_pct" in records.columns
    assert "mono_mass_da" in records.columns
    assert "id_type" in records.columns


def test_missing_reference_peptide_does_not_infer_zero():
    data = _sample_data()
    data = data.copy()
    new = as_table([dict(experiment="4204_combined", sample="4204", raw_file="4204_NR.raw",
                         preparation="Non-reduced", component_no=88, identified=True,
                         identification="2:R22-R32 = 1427.8m", peptide_sequence="RGFFYTPKTRR",
                         modification="None", ms_area=33., rt_min=5.3)])
    data = pd.concat([data, new], ignore_index=True)
    out = detailed_reference_comparison(data, "4197_combined / 4197_NR.raw", ["4204_combined / 4204_NR.raw"])
    other = out.loc[out.feature == "2:R22-R32 = 1427.8m"].iloc[0]
    assert math.isnan(other.reference_area)
    assert math.isnan(other.fold_vs_reference)
    assert other.comparison_area == 33.


def test_no_comparison_yields_no_fold_rows():
    data = _sample_data()
    ref = "4197_combined / 4197_NR.raw"
    out = detailed_reference_comparison(data, ref, [])
    assert out.empty
    assert len(annotated_area_matrix(data, ref, [])) == 1


def test_zero_reference_stays_nonquantifiable():
    data = _sample_data()
    data.loc[(data.raw_file == "4197_NR.raw") & (data.identified), "ms_area"] = 0
    out = detailed_reference_comparison(data, "4197_combined / 4197_NR.raw", ["4204_combined / 4204_NR.raw"])
    assert out.reference_area.iloc[0] == 0
    assert pd.isna(out.fold_vs_reference.iloc[0])


def test_two_explicit_comparators_both_same_prep():
    data = _sample_data()
    ref = "4197_combined / 4197_NR.raw"
    targets = ["4204_combined / 4204_NR.raw", "4300_combined / 4300_NR.raw"]
    out = detailed_reference_comparison(data, ref, targets)
    assert len(out) == 2
    assert set(out.comparison_run) == set(targets)
    assert out.loc[out.comparison_run == targets[1], "comparison_area"].iloc[0] == 160
    matrix = annotated_area_matrix(data, ref, targets)
    assert f"MS Area | {targets[0]}" in matrix.columns
    assert f"MS Area | {targets[1]}" in matrix.columns
    assert all("DTT" not in name for name in matrix.columns)


def test_no_detectable_area_keeps_original_rt_only_in_component_table():
    data = _sample_data()
    data.loc[(data.raw_file == "4197_NR.raw") & (data.identified), "ms_area"] = 0
    ref = "4197_combined / 4197_NR.raw"
    comp = "4204_combined / 4204_NR.raw"
    matrix = annotated_area_matrix(data, ref, [comp])
    assert pd.isna(matrix.iloc[0][f"RT dominant (min) | {ref}"])
    detail = detailed_component_records(data, ref, [comp])
    original = detail.loc[detail.raw_file == "4197_NR.raw", "rt_min"]
    assert set(original) == {10.1, 10.6}
