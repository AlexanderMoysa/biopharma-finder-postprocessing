"""Parse original BioPharma Finder wide CSV into the normalized long model."""
from __future__ import annotations
import csv
import io
import re
from pathlib import Path
import pandas as pd
from bpf_postprocess.model import as_table, derive_sample_name, guess_preparation


def read_csv(source: str | Path | bytes, name: str | None = None) -> tuple[pd.DataFrame, dict, list[str]]:
    if isinstance(source, bytes):
        blob = source
        source_name = name or "uploaded.csv"
    else:
        path = Path(source)
        blob = path.read_bytes()
        source_name = name or path.name
    text = None
    for encoding in ("utf-8-sig", "utf-16", "cp1252"):
        try:
            trial = blob.decode(encoding)
            if "No." in trial[:10000]:
                text = trial
                break
        except UnicodeError:
            pass
    if text is None:
        raise ValueError("Unrecognized text encoding or missing No. column in CSV")
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames or "No." not in reader.fieldnames:
        raise ValueError("Not a BioPharma Finder component table: missing 'No.'")
    area_names = [k[len("MS Area:"):] for k in reader.fieldnames if k.lower().startswith("ms area:") and k.lower().endswith(".raw")]
    if not area_names:
        raise ValueError("No MS Area:<RAW filename> columns found")
    experiment = re.sub(r"(?i)\.csv$", "", source_name)
    sample = derive_sample_name(experiment)
    rows = []
    for rec in reader:
        if not rec.get("No.", "").strip():
            continue
        for raw in area_names:
            ident = (rec.get("Identification") or "").strip()
            data = {
                "experiment": experiment, "sample": sample, "sample_group": "", "raw_file": raw,
                "condition": "", "preparation": guess_preparation(raw), "replicate": "",
                "component_no": rec.get("No."), "identified": bool(ident),
                "identification": ident, "peptide_sequence": rec.get("Peptide Sequence") or "",
                "modification": rec.get("Mod") or "", "site": rec.get("Site") or "",
                "positions": rec.get("Positions") or "", "protein": rec.get("Protein") or "",
                "rt_min": rec.get("RT:" + raw) or rec.get("RT"),
                "mz": rec.get("M/Z:" + raw) or rec.get("M/Z"),
                "charge": rec.get("Charge St."),
                "mono_mass_da": rec.get("Mono Mass Exp.:" + raw) or rec.get("Mono Mass Exp."),
                "ms_area": rec.get("MS Area:" + raw),
                "confidence_pct": rec.get("Conf. Score"),
                "delta_ppm": rec.get("∆ ppm:" + raw) or rec.get("∆ ppm"),
                "id_type": rec.get("ID Type:" + raw) or rec.get("ID Type"),
                "integration_type": rec.get("Integration Type"),
                "source_file": source_name, "source_type": "csv", "annotation_origin": "bf_csv", "import_note": "",
            }
            rows.append(data)
    if not rows:
        raise ValueError("No component records in the CSV")
    result = as_table(rows)
    warnings = []
    if result.duplicated(["experiment", "raw_file", "component_no"]).any():
        warnings.append("Duplicate component No. in CSV: inspect original export")
    metadata = {
        "experiment": experiment,
        "method": "Not included in CSV export",
        "protein_chains": [],
        "raw_runs": [{"raw_file": raw} for raw in area_names],
        "source_file": source_name,
    }
    return result, metadata, warnings
