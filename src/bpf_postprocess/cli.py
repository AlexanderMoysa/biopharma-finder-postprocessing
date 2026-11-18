"""CLI for reproducible import of BF files to normalized CSV/metadata."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from bpf_postprocess.importers.loader import load_files
from bpf_postprocess.reports.exports import csv_bytes


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert BioPharma Finder exports to normalized long table")
    parser.add_argument("files", nargs="+", type=Path, help=".resultsbpf, compact .zip, or BF .csv files")
    parser.add_argument("--out", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    result = load_files(args.files)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "components.csv").write_bytes(csv_bytes(result.data))
    (args.out / "runs.csv").write_bytes(csv_bytes(result.runs))
    (args.out / "experiments.json").write_text(json.dumps(result.experiments.to_dict("records"), ensure_ascii=False, indent=2), encoding="utf-8")
    (args.out / "warnings.txt").write_text("\n".join(result.warnings), encoding="utf-8")
    print(f"Imported {len(result.experiments)} experiments, {len(result.runs)} RAW runs, {len(result.data):,} component-run rows")
    print(f"Exports: {args.out.resolve()}")
    for note in result.warnings:
        print("WARNING:", note)

if __name__ == "__main__":
    main()
