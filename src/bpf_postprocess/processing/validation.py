"""Independent row-level agreement of experimental PMF extraction to BF CSV."""
from __future__ import annotations
import pandas as pd
from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.importers.results_importer import read_resultsbpf

TOLERANCES = {
    "rt_min": 0.006,
    "mz": 0.00055,
    "charge": 0.0,
    "mono_mass_da": 0.000055,
    "ms_area": 0.51,
    "confidence_pct": 0.055,
}


def validate_native(native_path, csv_path) -> dict:
    """Compare one PMF experiment against its own BF CSV export.

    Reports every mismatch, including non-identified components. CSV precision
    determines numeric tolerances. Validation is only for fields present in both.
    """
    native, meta, _ = read_resultsbpf(native_path)
    csvdata, csvmeta, _ = read_csv(csv_path)
    if native["raw_file"].nunique() != csvdata["raw_file"].nunique():
        # Continue and report missing runs as missing records.
        pass
    keys = ["raw_file", "component_no"]
    merged = native.merge(csvdata, on=keys, how="outer", suffixes=("_pmf", "_bf"), indicator=True,
                          validate="one_to_one")
    failures = []
    missing_pmf = int((merged["_merge"] == "right_only").sum())
    missing_bf = int((merged["_merge"] == "left_only").sum())
    both = merged[merged["_merge"] == "both"]
    id_pmf = both["identification_pmf"].fillna("").astype(str)
    id_bf = both["identification_bf"].fillna("").astype(str)
    id_diffs = both[id_pmf.ne(id_bf)]
    for _, line in id_diffs.iterrows():
        failures.append({"raw_file": line["raw_file"], "component_no": int(line["component_no"]),
                         "field": "identification", "pmf": line["identification_pmf"], "csv": line["identification_bf"]})
    for field, tol in TOLERANCES.items():
        a, b = f"{field}_pmf", f"{field}_bf"
        if a not in both or b not in both:
            continue
        left = pd.to_numeric(both[a], errors="coerce")
        right = pd.to_numeric(both[b], errors="coerce")
        differs = (left.sub(right).abs() > tol) | (left.isna() != right.isna())
        for _, line in both[differs].iterrows():
            failures.append({"raw_file": line["raw_file"], "component_no": int(line["component_no"]),
                             "field": field, "pmf": line[a], "csv": line[b]})
    return {
        "experiment": meta["experiment"],
        "bf_csv_experiment": csvmeta["experiment"],
        "native_rows": len(native),
        "bf_csv_rows": len(csvdata),
        "matched_rows": len(both),
        "missing_in_pmf": missing_pmf,
        "missing_in_csv": missing_bf,
        "field_mismatches": len(failures),
        "pass": (missing_pmf == 0 and missing_bf == 0 and not failures),
        "examples": failures[:30],
    }
