"""Source dispatch, lossless BF CSV retention, validated PMF+CSV enrichment."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

from bpf_postprocess.model import ImportResult, make_runs
from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.importers.results_importer import read_resultsbpf


MATCH_KEYS = ["experiment", "raw_file", "component_no"]
CSV_METADATA_FIELDS = [
    "peptide_sequence", "modification", "site", "positions", "protein",
    "delta_ppm", "delta_ppm_overall", "id_type", "id_type_overall", "integration_type",
    "level", "best_asr", "avg_mass_da", "mono_mass_theo_da",
    "missed_cleavages", "bf_ratio_label", "bf_ratio", "bf_max_condition",
    "bf_min_condition", "bf_max_ms_area", "bf_avg_ms_area_condition",
    "bf_cv_pct_condition", "comment", "bf_csv_source_file",
]
# Comparisons are only to establish that metadata belongs to the SAME BF
# processing result. PMF float32 values are deliberately kept after matching.
VALIDATION_TOLERANCES = {
    "ms_area": 0.51, "rt_min": 0.0051, "mz": 0.00051,
    "mono_mass_da": 0.000051, "confidence_pct": 0.051, "charge": 0.0,
}


def load_file(path: str | Path, *, display_name: str | None = None):
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == ".csv":
        return read_csv(source, name=display_name)
    if suffix in (".resultsbpf", ".zip"):
        return read_resultsbpf(source, name=display_name)
    raise ValueError(f"Unsupported input type: {suffix}")


def _checked_index(frame: pd.DataFrame, kind: str) -> pd.DataFrame:
    if frame.duplicated(MATCH_KEYS).any():
        raise ValueError(f"Duplicate {kind} experiment/RAW/component keys; unsafe to enrich")
    return frame.set_index(MATCH_KEYS, verify_integrity=True)


def enrich_pmf_with_csv(pmf: pd.DataFrame, csv: pd.DataFrame) -> pd.DataFrame:
    """Use PMF numerical precision and independently verified BF CSV annotations.

    Abort if component sets or measured values disagree. No fuzzy matching,
    no silent overwrite of unrelated experiments, no doubling of records.
    """
    left = _checked_index(pmf, "PMF")
    right = _checked_index(csv, "CSV")
    missing = left.index.difference(right.index)
    extra = right.index.difference(left.index)
    if len(missing) or len(extra):
        raise ValueError(
            "PMF/CSV records do not align: "
            f"{len(missing)} PMF keys absent from CSV; {len(extra)} CSV keys absent from PMF. "
            f"Example PMF-only: {list(missing[:2])}; CSV-only: {list(extra[:2])}"
        )
    right = right.reindex(left.index)
    id_conflict = left["identification"].fillna("").astype(str).str.strip().ne(
        right["identification"].fillna("").astype(str).str.strip()
    )
    if id_conflict.any():
        index = left.index[id_conflict][0]
        raise ValueError(f"PMF/CSV identification mismatch at {index}; use matching BioPharma Finder exports")
    for col, tolerance in VALIDATION_TOLERANCES.items():
        a = pd.to_numeric(left[col], errors="coerce")
        b = pd.to_numeric(right[col], errors="coerce")
        # One format may lack a confidence score; only compare when both exist.
        invalid = a.notna() & b.notna() & (a.sub(b).abs() > tolerance + 1e-10)
        if invalid.any():
            index = left.index[invalid][0]
            raise ValueError(f"PMF/CSV {col} mismatch at {index}: "
                             f"PMF={a.loc[index]}, CSV={b.loc[index]}; "
                             "do not combine different processing results")
    enriched = left.copy()
    for col in CSV_METADATA_FIELDS:
        enriched[col] = right[col].to_numpy()
    # Keep the PMF's higher-precision float32 measurements, which were matched
    # above; imported identification, RT, area etc. remain PMF-based.
    enriched["source_type"] = "resultsbpf+csv"
    enriched["annotation_origin"] = "bf_csv_verified_against_pmf"
    enriched["import_note"] = "PMF numeric values; BF CSV annotations (verified matching keys and measurements)"
    enriched["condition"] = right["condition"].to_numpy()
    return enriched.reset_index()


def load_files(sources: list[str | Path]) -> ImportResult:
    if not sources:
        raise ValueError("Upload at least one BioPharma Finder file")
    groups: dict[str, dict[str, tuple[pd.DataFrame, dict, list[str], str]]] = defaultdict(dict)
    notes: list[str] = []
    original_csv_tables: dict[str, pd.DataFrame] = {}
    for file in sources:
        path = Path(file)
        tab, meta, warnings = load_file(path, display_name=path.name)
        experiment = str(meta["experiment"])
        kind = "csv" if path.suffix.lower() == ".csv" else "pmf"
        if kind in groups[experiment]:
            raise ValueError(
                f"Duplicate {kind.upper()} source for {experiment}. "
                "Use at most one CSV and one PMF/resultsbpf per experiment."
            )
        raw_csv = meta.pop("original_csv_rows", None)
        if raw_csv is not None:
            original_csv_tables[experiment] = raw_csv
        groups[experiment][kind] = (tab, meta, warnings, path.name)

    merged_tables, exps = [], []
    for experiment, variants in groups.items():
        if "csv" in variants and "pmf" in variants:
            pmf, pmf_meta, pmf_notes, pmf_file = variants["pmf"]
            csv, csv_meta, csv_notes, csv_file = variants["csv"]
            merged_tables.append(enrich_pmf_with_csv(pmf, csv))
            # PMF metadata retains complete protein sequence and method.
            exps.append(pmf_meta)
            notes.append(
                f"{experiment}: PMF + BF CSV verified and enriched: "
                f"{len(pmf):,} component/RAW records, 0 duplicate rows. "
                f"Source files: {pmf_file}, {csv_file}."
            )
            # PMF ID-Type/site warnings are no longer applicable after enrichment.
            notes.extend(f"{pmf_file}: {w}" for w in pmf_notes
                         if "ID Type" not in w and "disulfide sites/connectivity" not in w)
            notes.extend(f"{csv_file}: {w}" for w in csv_notes)
        else:
            kind = next(iter(variants))
            table, metadata, warnings, source_name = variants[kind]
            merged_tables.append(table)
            exps.append(metadata)
            notes.extend(f"{source_name}: {w}" for w in warnings)
    combined = pd.concat(merged_tables, ignore_index=True)
    if combined.duplicated(MATCH_KEYS).any():
        raise ValueError("Overlapping experiment/RAW/component records remain; refusing to double count")
    return ImportResult(combined, pd.DataFrame(exps), make_runs(combined),
                        notes, source_csv_tables=original_csv_tables)
