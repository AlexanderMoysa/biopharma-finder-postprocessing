"""Validate a native experiment's PMF against original BF CSV export."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from bpf_postprocess.processing.validation import validate_native


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native", type=Path, help=".resultsbpf or compact PMF .zip")
    parser.add_argument("csv", type=Path, help="BF-exported CSV from the same experiment")
    parser.add_argument("--out", type=Path, help="Optional JSON report path")
    args = parser.parse_args()
    report = validate_native(args.native, args.csv)
    data = json.dumps(report, ensure_ascii=False, indent=2)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(data, encoding="utf-8")
        print("Saved:", args.out)
    print(data)
    if not report["pass"]:
        raise SystemExit(2)

if __name__ == "__main__":
    main()
