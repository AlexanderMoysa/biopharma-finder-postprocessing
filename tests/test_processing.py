import pandas as pd
import pytest
from bpf_postprocess.model import as_table
from bpf_postprocess.processing.filters import filter_components
from bpf_postprocess.processing.comparison import long_peptide_areas, compare_runs, ratio_table
from bpf_postprocess.processing.unidentified import unidentified_summary


def _rows():
    items=[]
    for raw, areas in [("a.raw", [10, 5, 100]), ("b.raw", [20, 0, 200])]:
        for no, value in enumerate(areas, start=1):
            items.append(dict(experiment="E",sample="S",raw_file=raw,preparation="Unspecified",
                              component_no=no, identified=no < 3, identification="P1" if no < 3 else "",
                              ms_area=value))
    return as_table(items)


def test_compare_aggregates_same_id_within_each_run():
    matrix=compare_runs(_rows())
    assert matrix.shape == (1,3)
    assert matrix.iloc[0]["E / a.raw"] == 15
    assert matrix.iloc[0]["E / b.raw"] == 20
    ratio=ratio_table(matrix,"E / a.raw","E / b.raw")
    assert ratio.iloc[0]["fold_vs_reference"] == pytest.approx(4/3)


def test_zero_reference_does_not_use_9999_ratio():
    matrix=pd.DataFrame({"feature":["P1"], "R":[0.],"T":[4.]})
    out=ratio_table(matrix,"R","T")
    assert pd.isna(out.iloc[0]["fold_vs_reference"])


def test_filters_and_unknowns():
    rows=_rows()
    result=filter_components(rows,identified="Unidentified")
    assert len(result)==2
    summary=unidentified_summary(rows)
    assert summary["components"].tolist()==[1,1]
