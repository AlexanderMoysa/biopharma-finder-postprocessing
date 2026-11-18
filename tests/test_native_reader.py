from pathlib import Path
import os
import pytest
from bpf_postprocess.importers.results_importer import read_resultsbpf


def test_native_reader_with_optional_external_fixture():
    path_str=os.getenv("BPF_TEST_CORE_ZIP")
    if not path_str:
        pytest.skip("Set BPF_TEST_CORE_ZIP to an independently validated local PMF core ZIP")
    data, metadata, warnings=read_resultsbpf(Path(path_str))
    assert len(data)>0
    assert data["raw_file"].nunique()>=1
    assert data["identified"].any()
    assert "audit" in metadata
    assert warnings
