"""Regression tests for shared peptide/component selection across tabs."""
import pandas as pd
import pytest

from bpf_postprocess.model import as_table
from bpf_postprocess.processing.comparison import compare_runs
from bpf_postprocess.processing.modifications import modification_areas
from bpf_postprocess.processing.peptide_selection import (
    apply_peptide_exclusions,
    filter_peptide_table,
    peptide_inventory,
    update_peptide_exclusions,
)
from bpf_postprocess.processing.unidentified import unidentified_summary


def fixture_rows():
    rows = []
    for experiment, raw, area_a, area_b in (
        ("4197", "nr.raw", 100, 40),
        ("4204", "nr.raw", 200, 60),
    ):
        rows.extend([
            dict(experiment=experiment, raw_file=raw, component_no=1, identified=True,
                 identification="1:G1-E4 = 416m", peptide_sequence="GIVE",
                 modification="None", ms_area=area_a),
            dict(experiment=experiment, raw_file=raw, component_no=2, identified=True,
                 identification="2:R22-R32 = 1427m", peptide_sequence="RGFFYTPKTRR",
                 modification="Oxidation", ms_area=area_b),
            dict(experiment=experiment, raw_file=raw, component_no=3, identified=False,
                 identification="", ms_area=3),
        ])
    return as_table(rows)


def test_initially_every_peptide_is_selected_and_identifiers_are_unique():
    raw = fixture_rows()
    inventory = peptide_inventory(raw)
    assert len(inventory) == 2
    assert inventory["raw_runs"].tolist() == [2, 2]
    assert inventory["components"].tolist() == [2, 2]
    assert len(apply_peptide_exclusions(raw, set())) == 6


def test_exclude_one_peptide_from_all_runs_without_losing_unknown_components():
    raw = fixture_rows()
    inventory = peptide_inventory(raw)
    excluded = update_peptide_exclusions(
        inventory["feature"], pd.Series([False, True]), set()
    )
    processed = apply_peptide_exclusions(raw, excluded)
    assert len(processed) == 4
    assert set(processed.loc[processed["identified"], "identification"]) == {
        inventory.iloc[1]["feature"]
    }
    assert (~processed["identified"]).sum() == 2
    matrix = compare_runs(processed)
    assert matrix.shape[0] == 1
    mods = modification_areas(processed)
    assert set(mods["modification"]) == {"Oxidation"}
    assert unidentified_summary(processed)["components"].sum() == 2
    # Inventory is generated from pre-exclusion data: missing features are recoverable.
    assert len(peptide_inventory(raw)) == 2


def test_reselect_a_peptide_without_resetting_other_hidden_exclusions():
    original = {"hidden", "1:G1-E4 = 416m"}
    restored = update_peptide_exclusions(
        pd.Series(["1:G1-E4 = 416m"]), pd.Series([True]), original
    )
    assert restored == {"hidden"}
    unchanged = update_peptide_exclusions(
        pd.Series(["2:R22-R32 = 1427m"]), pd.Series([False]), restored
    )
    assert unchanged == {"hidden", "2:R22-R32 = 1427m"}


def test_duplicate_peptide_labels_refused_in_selection_editor():
    with pytest.raises(ValueError, match="one row per identification"):
        update_peptide_exclusions(
            pd.Series(["A", "A"]), pd.Series([True, False]), set()
        )


def test_empty_inventory_and_filter_do_not_crash():
    rows = fixture_rows().iloc[:0]
    inventory = peptide_inventory(rows)
    assert inventory.empty
    assert "feature" in inventory
    assert apply_peptide_exclusions(rows, {"A"}).empty


def test_derived_comparison_tables_preserve_reselectable_candidates():
    """All-peptides-off must not make the reference/comparator UI disappear."""
    from bpf_postprocess.processing.detailed_comparison import (
        run_catalog, select_comparison_runs, annotated_area_matrix,
        run_labeled_reference_matrix, detailed_reference_comparison,
    )
    raw = fixture_rows()
    # Fill preparation to make the two NR runs comparable.
    raw["preparation"] = "Non-reduced"
    names = run_catalog(raw)["run"].tolist()
    ref, cmp = names
    candidate_rows = select_comparison_runs(raw,ref,[cmp])
    assert len(candidate_rows) == 6
    all_features = set(peptide_inventory(raw)["feature"])
    filtered_rows = apply_peptide_exclusions(candidate_rows, all_features)
    assert len(filtered_rows) == 2  # unknowns, not removed by peptide filtering
    wide_candidates = annotated_area_matrix(raw,ref,[cmp])
    long_candidates = detailed_reference_comparison(raw,ref,[cmp])
    labeled_candidates = run_labeled_reference_matrix(raw,ref,[cmp])
    assert len(wide_candidates) == len(labeled_candidates) == 2
    assert len(long_candidates) == 2
    for candidate in (wide_candidates, long_candidates, labeled_candidates):
        visible = filter_peptide_table(candidate, all_features)
        assert visible.empty
        restored = filter_peptide_table(candidate, set())
        assert len(restored) == 2


def test_filter_preserves_non_feature_table_and_hidden_peptides():
    x = pd.DataFrame({"feature": ["A", "B"], "ms_area": [1.0, 3.0]})
    assert filter_peptide_table(x,{"A"})["feature"].tolist() == ["B"]
    assert filter_peptide_table(x,{"unknown"}).equals(x)
    assert filter_peptide_table(x[["ms_area"]],{"A"}).equals(x[["ms_area"]])
