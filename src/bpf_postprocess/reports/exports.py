"""Programmatic table and plot export functions."""
from __future__ import annotations
import io
import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False).encode("utf-8-sig")


def excel_bytes(tables: dict[str, pd.DataFrame]) -> bytes:
    out = io.BytesIO()
    with pd.ExcelWriter(out, engine="openpyxl") as writer:
        used = set()
        for name, data in tables.items():
            sheet = re.sub(r"[\\/*?:\[\]]", "_", name)[:31] or "Sheet"
            while sheet in used:
                sheet = (sheet[:28] + f"_{len(used)}")[:31]
            used.add(sheet)
            data.to_excel(writer, sheet_name=sheet, index=False)
            writer.sheets[sheet].freeze_panes = "A2"
            writer.sheets[sheet].auto_filter.ref = writer.sheets[sheet].dimensions
    return out.getvalue()


def static_bar_figure_bytes(frame: pd.DataFrame, feature: str, ext: str = "png") -> bytes:
    """Export an actual static figure without dependence on Chrome/Kaleido."""
    if ext not in ("png", "svg"):
        raise ValueError("Only png/svg supported")
    sub = frame[(frame["feature"] == feature)].copy()
    if sub.empty:
        raise ValueError("No rows for selected feature")
    labels = (sub["sample"] + " / " + sub["raw_file"]).tolist()
    areas = sub["ms_area"].fillna(0).tolist()
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.2), 4.5))
    ax.bar(range(len(labels)), areas)
    ax.set_ylabel("MS Area (reported by BioPharma Finder)")
    ax.set_title(feature[:90])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    fig.tight_layout()
    out = io.BytesIO()
    fig.savefig(out, format=ext, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out.getvalue()
