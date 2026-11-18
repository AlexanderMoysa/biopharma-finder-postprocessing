import pandas as pd
import pytest
from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.importers.loader import load_files

CSV_BYTES = ("No.,Identification,Peptide Sequence,Mod,Site,Charge St.,Conf. Score,"
             "MS Area:a.raw,RT:a.raw,M/Z:a.raw,Mono Mass Exp.:a.raw,"
             "MS Area:b_with_DTT.raw,RT:b_with_DTT.raw,M/Z:b_with_DTT.raw,Mono Mass Exp.:b_with_DTT.raw\n"
             "1,P1,AA,Oxidation,M1,2,99.0,123,5.2,600.1,1198.2,456,5.3,600.2,1198.3\n"
             "2,,,,,1,,0,2.0,220.1,219.0,77,2.1,220.2,219.1\n").encode("utf-8")


def test_csv_two_raws_and_unknowns():
    data, meta, warnings = read_csv(CSV_BYTES, name="test_combined.csv")
    assert len(data) == 4
    assert data["raw_file"].nunique() == 2
    assert data["identified"].sum() == 2
    assert data["ms_area"].sum() == 656
    assert data.loc[data.raw_file == "b_with_DTT.raw", "preparation"].eq("Reduced").all()
    assert data.loc[data.raw_file == "a.raw", "preparation"].eq("Unspecified").all()
    assert meta["experiment"] == "test_combined"
    assert not warnings


def test_invalid_csv():
    with pytest.raises(ValueError, match="No."):
        read_csv(b"wrong,header\nabc,def\n", name="bad.csv")
