"""Thermo BioPharma Finder peptide .resultsbpf reader (experimental, v2).

Reads ALL component entries, including unidentified components, directly from
PMF stored inside a native .resultsbpf or a compact *_core_for_analysis.zip.
It ignores the embedded Thermo RAW data and never modifies input files.

Validated against a CSV on one BioPharma Finder 5.0.68.8 experiment.
Cross-check against BioPharma Finder's CSV before using on other experiments,
software versions or for regulated / quantitative conclusions.

Use the public functions load_core(), extract_components(), compare_csv().

Note: RT, m/z, MS area are *run-specific*. A blank Identification in output
means BioPharma Finder did not assign an ID, not that it is a free peptide.
For cross-experiment matching, do NOT use component numbers as universal IDs.
"""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
import struct
import tempfile
import zipfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

MARKER = b"\xff\xfe\xff"
# PMF uses these TWO distinct markers before the 0-based component number.
# The earlier prototype handled only FF FF, causing unidentified records
# with the 00 00 variant to be missed.
COMPONENT_ID_PREFIXES = (b"\xff\xff", b"\x00\x00")
MAX_PMF_UNCOMPRESSED = 500 * 1024**2
FIELDNAMES = [
    "Experiment", "RAW File", "Condition", "Component No.",
    "Identified", "Identification", "RT (min)", "M/Z", "Charge",
    "MS Area", "Mono Mass Exp. (Da)", "Conf. Score (%)",
]


@contextmanager
def open_experiment_archive(source: Path):
    """Yield the inner experiment archive, never reading RawFiles.zip content."""
    with zipfile.ZipFile(source) as outer:
        names = outer.namelist()
        if "ExperimentSummary.json" in names and any(n.lower().endswith(".pmf") for n in names):
            yield outer  # Compact core ZIP
            return
        inner_name = next((n for n in names if Path(n).name.lower() == "experiment.zip"), None)
        if inner_name is None:
            raise ValueError("Expected Experiment.zip or ExperimentSummary.json + .pmf in input")
        # A ZIP nested inside another ZIP needs a seekable stream; tempfile is
        # disk-backed and automatically closed. RawFiles.zip itself is NOT read.
        with tempfile.TemporaryFile() as temp:
            with outer.open(inner_name) as fp:
                shutil.copyfileobj(fp, temp, length=1024**2)
            temp.seek(0)
            with zipfile.ZipFile(temp) as inner:
                yield inner


def load_core(source: Path) -> tuple[dict, bytes, str]:
    with open_experiment_archive(source) as z:
        pmfs = [name for name in z.namelist() if name.lower().endswith(".pmf")]
        if len(pmfs) != 1 or "ExperimentSummary.json" not in z.namelist():
            raise ValueError("Expected exactly one PMF and ExperimentSummary.json")
        if z.getinfo(pmfs[0]).file_size > MAX_PMF_UNCOMPRESSED:
            raise ValueError("PMF exceeds the reader's 500 MiB safety limit")
        metadata = json.loads(z.read("ExperimentSummary.json"))
        pmf_data = z.read(pmfs[0])
        return metadata, pmf_data, pmfs[0]


def locate_runs(pmf: bytes, info: dict) -> list[dict]:
    paths = info.get("FullPathFileNames") or []
    if not paths:
        raise ValueError("No FullPathFileNames in metadata; cannot map PMF to RAW runs")
    conditions = {
        str(x["k"]).lower(): str(x["v"])
        for x in info.get("Rawfile2ConditionDic") or []
        if isinstance(x, dict) and "k" in x and "v" in x
    }
    runs = []
    for raw_path in paths:
        raw_name = raw_path.replace("/", "\\").rsplit("\\", 1)[-1]
        offset = pmf.find(raw_path.encode("utf-16le"))
        if offset < 0:
            raise ValueError(f"Cannot identify PMF run block for {raw_name}")
        runs.append({"raw_file": raw_name, "condition": conditions.get(raw_name.lower(), ""),
                     "path": raw_path, "start_offset": offset})
    runs.sort(key=lambda r: r["start_offset"])
    if len(set(r["start_offset"] for r in runs)) != len(runs):
        raise ValueError("Cannot distinguish RAW run blocks (duplicate PMF offsets)")
    for j, run in enumerate(runs):
        run["end_offset"] = runs[j + 1]["start_offset"] if j + 1 < len(runs) else len(pmf)
    return runs


def extract_components(pmf: bytes, info: dict) -> tuple[list[dict], list[dict], dict]:
    """Scan observed PMF record structure; both identified & unidentified.

    Record:
      FF FE FF <one-byte UTF-16LE length> <UTF-16LE identification>
      <FF FF or 00 00> <uint32 zero-based component ID>
    Run-specific floats (RT, m/z, MS Area, mass) are stored at fixed offsets
    before identification start. First record for each (run, component No.)
    corresponds to an individual component; later records are aggregation
    duplicates. These rules were empirically verified on 4197 only.
    """
    runs = locate_runs(pmf, info)
    unique = [dict() for _ in runs]
    duplicates = [0] * len(runs)
    candidates = [0] * len(runs)
    accepted_kinds = [Counter() for _ in runs]
    pos = 0
    while True:
        idx = pmf.find(MARKER, pos)
        if idx < 0:
            break
        pos = idx + len(MARKER)
        if idx + 4 > len(pmf):
            break
        count = pmf[idx + 3]
        if count > 250:
            continue
        start = idx + 4
        end = start + 2 * count
        if start < 46 or end + 6 > len(pmf):
            continue
        prefix = pmf[end:end + 2]
        if prefix not in COMPONENT_ID_PREFIXES:
            continue
        try:
            pieces = pmf[start:end].decode("utf-16le").split("\t", 1)[0].splitlines()
            annotation = pieces[0].strip() if pieces else ""
            charge = struct.unpack_from("<H", pmf, start - 46)[0]
            rt = struct.unpack_from("<f", pmf, start - 44)[0]
            mz = struct.unpack_from("<f", pmf, start - 32)[0]
            area = struct.unpack_from("<f", pmf, start - 24)[0]
            monomass = struct.unpack_from("<f", pmf, start - 12)[0]
            confidence = struct.unpack_from("<f", pmf, start - 8)[0]
            comp_no = struct.unpack_from("<I", pmf, end + 2)[0] + 1
        except (ValueError, IndexError, struct.error, UnicodeError):
            continue
        # Reject false marker hits in the RAW-run blocks.
        if not (1 <= charge <= 50 and math.isfinite(rt) and 0 < rt < 150
                and math.isfinite(mz) and 30 < mz < 20000
                and math.isfinite(area) and 0 <= area < 1e13
                and math.isfinite(monomass) and 0 <= monomass < 200000
                and math.isfinite(confidence) and 0 <= confidence <= 1
                and 1 <= comp_no < 10_000_000):
            continue
        k = next((j for j, run in enumerate(runs)
                  if run["start_offset"] <= start < run["end_offset"]), None)
        if k is None:
            continue
        candidates[k] += 1
        if comp_no in unique[k]:
            duplicates[k] += 1
            continue
        accepted_kinds[k][prefix.hex()] += 1
        unique[k][comp_no] = {
            "Experiment": info.get("ExperimentName", ""),
            "RAW File": runs[k]["raw_file"],
            "Condition": runs[k]["condition"],
            "Component No.": comp_no,
            "Identified": bool(annotation),
            "Identification": annotation,
            "RT (min)": rt,
            "M/Z": mz,
            "Charge": charge,
            "MS Area": area,
            "Mono Mass Exp. (Da)": monomass,
            # Blank ID means no confidence available for reporting.
            "Conf. Score (%)": confidence * 100 if annotation else "",
        }
    rows = [item for run_rows in unique for _, item in sorted(run_rows.items())]
    audit = {
        "candidates_per_run": candidates,
        "unique_components_per_run": [len(x) for x in unique],
        "identified_components_per_run": [sum(bool(r["Identified"]) for r in v.values()) for v in unique],
        "unidentified_components_per_run": [sum(not r["Identified"] for r in v.values()) for v in unique],
        "duplicate_or_aggregate_records_skipped": duplicates,
        "component_prefix_counts": [dict(c) for c in accepted_kinds],
    }
    return rows, runs, audit


def compare_csv(rows: list[dict], csv_path: Path) -> tuple[list[str], int]:
    """Compare all records and numerical fields, using CSV's export precision."""
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        csvrows = list(csv.DictReader(f))
    reference = {int(item["No."]): item for item in csvrows}
    byrun: dict[str, dict[int, dict]] = {}
    for row in rows:
        byrun.setdefault(str(row["RAW File"]), {})[int(row["Component No."])] = row
    errors = 0
    lines = [f"Reference CSV: {csv_path}", f"Reference component rows: {len(csvrows)}"]
    for raw_name, run_rows in byrun.items():
        missing = sorted(set(reference) - set(run_rows))
        extra = sorted(set(run_rows) - set(reference))
        mistakes: list[str] = []
        for no in missing[:15]:
            mistakes.append(f"No. {no}: missing in PMF extraction")
        for no in extra[:15]:
            mistakes.append(f"No. {no}: not present in reference CSV")
        numeric = [
            ("RT (min)", f"RT:{raw_name}", 0.006),
            ("M/Z", f"M/Z:{raw_name}", 0.00055),
            ("MS Area", f"MS Area:{raw_name}", 0.51),
            ("Mono Mass Exp. (Da)", f"Mono Mass Exp.:{raw_name}", 0.000055),
            ("Charge", "Charge St.", 0),
            ("Conf. Score (%)", "Conf. Score", 0.055),
        ]
        field_checks = Counter()
        for no in sorted(set(run_rows) & set(reference)):
            got, want = run_rows[no], reference[no]
            if got["Identification"] != want.get("Identification", ""):
                mistakes.append(f"No. {no}: identification differs")
            field_checks["Identification"] += 1
            for key, csv_key, tol in numeric:
                target = want.get(csv_key)
                actual = got[key]
                if target in (None, ""):
                    continue
                if actual in (None, ""):
                    mistakes.append(f"No. {no} {key}: PMF value blank; reference {target}")
                    continue
                try:
                    a, b = float(actual), float(target)
                except (ValueError, TypeError):
                    mistakes.append(f"No. {no} {key}: non-numeric value")
                    continue
                field_checks[key] += 1
                if abs(a - b) > tol:
                    mistakes.append(f"No. {no} {key}: PMF {a:.9g}, CSV {b:.9g}")
        total_mismatches = len(missing) + len(extra) + max(0, len(mistakes) - min(15, len(missing)) - min(15, len(extra)))
        errors += total_mismatches
        lines.append(f"\nRAW run: {raw_name}")
        lines.append(f"Extracted {len(run_rows)}/{len(reference)} components")
        lines.append(f"IDs: {sum(bool(x['Identification']) for x in run_rows.values())} identified, "
                     f"{sum(not bool(x['Identification']) for x in run_rows.values())} unidentified")
        lines.append("Compared fields: " + ", ".join(f"{key}={val}" for key, val in field_checks.items()))
        lines.append(f"Missing={len(missing)} Extra={len(extra)} Mismatches={total_mismatches}")
        lines.extend("  " + item for item in mistakes[:30])
    return lines, errors

