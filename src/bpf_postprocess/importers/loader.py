"""Input dispatch, metadata collection, safe multi-file concatenation."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from bpf_postprocess.model import ImportResult, make_runs
from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.importers.results_importer import read_resultsbpf


def load_file(path: str | Path, *, display_name: str | None = None):
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == ".csv":
        return read_csv(source, name=display_name)
    if suffix in (".resultsbpf", ".zip"):
        return read_resultsbpf(source, name=display_name)
    raise ValueError(f"Unsupported input type: {suffix}")


def load_files(sources: list[str | Path]) -> ImportResult:
    tables, exps, notes = [], [], []
    for file in sources:
        tab, meta, warnings = load_file(file)
        tables.append(tab)
        exps.append(meta)
        notes.extend(f"{Path(file).name}: {m}" for m in warnings)
    if not tables:
        raise ValueError("Upload at least one file")
    combined = pd.concat(tables, ignore_index=True)
    # Importing PMF and its CSV for the same experiment can double count results.
    overlaps = combined.duplicated(["experiment", "raw_file", "component_no"], keep=False)
    if overlaps.any():
        raise ValueError(
            "Duplicated experiment/RAW/component rows. You may have uploaded both CSV and "
            "resultsbpf from the same experiment or repeated an input. Choose one source per experiment."
        )
    return ImportResult(combined, pd.DataFrame(exps), make_runs(combined), notes)
