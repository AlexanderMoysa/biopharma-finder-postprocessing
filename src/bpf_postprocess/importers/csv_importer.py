"""Import BioPharma Finder CSV without dropping original metadata.

Each original BF component row contains results from multiple RAW injections.
This reader makes one normalized component/RAW row, and keeps an independent
copy of the entire original wide CSV (all columns with original names).
"""
from __future__ import annotations

import csv
import io
import re
from pathlib import Path

import pandas as pd

from bpf_postprocess.model import as_table, derive_sample_name, guess_preparation
from bpf_postprocess.processing.id_types import UNAVAILABLE_ID_TYPE


def _read_text(blob: bytes) -> str:
    for encoding in ("utf-8-sig", "utf-16", "cp1252"):
        try:
            candidate = blob.decode(encoding)
            if "No." in candidate[:10000]:
                return candidate
        except UnicodeError:
            continue
    raise ValueError("Unrecognized text encoding or missing No. column in CSV")


def _condition_labels(header: list[str], raws: list[str]) -> dict[str, str]:
    """Link adjacent Avg MS Area/%%CV columns to their following MS Area RAW.

    Mapping by header position, rather than guessing from RAW filenames,
    handles arbitrary condition names such as 05, _04, and DTT.
    """
    found: dict[str, str] = {}
    current = None
    for heading in header:
        if heading.startswith("Avg MS Area: "):
            current = heading.removeprefix("Avg MS Area: ").strip()
        elif heading.startswith("MS Area:"):
            raw = heading.removeprefix("MS Area:")
            if raw in raws and current is not None:
                found[raw] = current
    return found


def read_csv(source: str | Path | bytes, name: str | None = None) -> tuple[pd.DataFrame, dict, list[str]]:
    if isinstance(source, bytes):
        blob = source
        source_name = name or "uploaded.csv"
    else:
        path = Path(source)
        blob = path.read_bytes()
        source_name = name or path.name
    reader = csv.DictReader(io.StringIO(_read_text(blob)))
    if not reader.fieldnames or "No." not in reader.fieldnames:
        raise ValueError("Not a BioPharma Finder component table: missing 'No.'")
    if len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("Duplicate BF CSV column names; cannot retain columns losslessly")
    headings = list(reader.fieldnames)
    area_names = [k.removeprefix("MS Area:") for k in headings
                  if k.lower().startswith("ms area:") and k.lower().endswith(".raw")]
    if not area_names:
        raise ValueError("No MS Area:<RAW filename> columns found")
    conditions = _condition_labels(headings, area_names)
    ratios = [k for k in headings if k.startswith("Ratio (")]
    if len(ratios) > 1:
        raise ValueError("Multiple BF Ratio fields; source cannot be mapped unambiguously")
    ratio_name = ratios[0] if ratios else ""
    # Windows appends (1), (2), etc. to downloaded copies of a CSV.
    # The experiment stored in a PMF is still the unsuffixed name.
    experiment = re.sub(r"(?i)\(\d+\)$", "", re.sub(r"(?i)\.csv$", "", source_name))
    sample = derive_sample_name(experiment)
    original = []
    rows = []
    def raw_or_overall(rec: dict, raw_header: str, overall_header: str):
        # If BF provided a RAW-specific column, preserve a blank cell as blank.
        # Do not fall back to an experiment-level value from another condition.
        return rec.get(raw_header) if raw_header in headings else rec.get(overall_header)

    for rec in reader:
        if not str(rec.get("No.") or "").strip():
            continue
        if None in rec:
            raise ValueError("Malformed CSV: component row has more values than header fields")
        original.append({k: rec.get(k, "") for k in headings})
        ident = str(rec.get("Identification") or "").strip()
        for raw in area_names:
            condition = conditions.get(raw, "")
            data = {
                "experiment": experiment, "sample": sample, "sample_group": "", "raw_file": raw,
                "condition": condition, "preparation": guess_preparation(raw), "replicate": "",
                "component_no": rec.get("No."), "identified": bool(ident),
                "identification": ident, "peptide_sequence": rec.get("Peptide Sequence") or "",
                "modification": rec.get("Mod") or "", "site": rec.get("Site") or "",
                "positions": rec.get("Positions") or "", "protein": rec.get("Protein") or "",
                "rt_min": raw_or_overall(rec, "RT:" + raw, "RT"),
                "mz": raw_or_overall(rec, "M/Z:" + raw, "M/Z"),
                "charge": rec.get("Charge St."),
                "mono_mass_da": raw_or_overall(rec, "Mono Mass Exp.:" + raw, "Mono Mass Exp."),
                "ms_area": rec.get("MS Area:" + raw),
                "confidence_pct": rec.get("Conf. Score"),
                "delta_ppm": raw_or_overall(rec, "∆ ppm:" + raw, "∆ ppm"),
                "delta_ppm_overall": rec.get("∆ ppm"),
                # 'Both' is an overall BF assignment, not a per-RAW ID type.
                "id_type": (rec.get("ID Type:" + raw) or "") if ("ID Type:" + raw) in headings
                           else ((rec.get("ID Type") or "") if len(area_names) == 1 else UNAVAILABLE_ID_TYPE),
                "id_type_overall": rec.get("ID Type") or "",
                "integration_type": rec.get("Integration Type") or "",
                "level": rec.get("Level") or "",
                "best_asr": rec.get("Best ASR"),
                "avg_mass_da": raw_or_overall(rec, "Avg Mass Exp.:" + raw, "Avg Mass Exp."),
                "mono_mass_theo_da": rec.get("Mono Mass Theo."),
                "missed_cleavages": rec.get("Missed Cleavages"),
                "bf_ratio_label": ratio_name,
                "bf_ratio": rec.get(ratio_name) if ratio_name else "",
                "bf_max_condition": rec.get("Max Cond.") or "",
                "bf_min_condition": rec.get("Min Cond.") or "",
                "bf_max_ms_area": rec.get("Max MS Area"),
                "bf_avg_ms_area_condition": rec.get("Avg MS Area: " + condition) if condition else "",
                "bf_cv_pct_condition": rec.get("%CV: " + condition) if condition else "",
                "comment": rec.get("Comment") or "",
                "bf_csv_source_file": source_name,
                "source_file": source_name, "source_type": "csv",
                "annotation_origin": "bf_csv", "import_note": "",
            }
            rows.append(data)
    if not rows:
        raise ValueError("No component records in the CSV")
    result = as_table(rows)
    if result.duplicated(["experiment", "raw_file", "component_no"]).any():
        raise ValueError("Duplicate component No. in CSV; cannot match BF records safely")
    meta = {
        "experiment": experiment,
        "method": "Not included in CSV export",
        "protein_chains": [],
        "raw_runs": [{"raw_file": raw, "condition": conditions.get(raw, "")} for raw in area_names],
        "source_file": source_name,
        # This sidecar preserves EVERY original BF column and the exact cell text,
        # including arbitrary condition/RAW-dependent fields.
        "original_csv_rows": pd.DataFrame(original, columns=headings),
    }
    return result, meta, []
