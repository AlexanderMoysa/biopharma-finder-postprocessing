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


def heatmap_figure_bytes(matrix: pd.DataFrame, ext: str = 'png') -> bytes:
    """Static comparison heatmap, log10(MS Area+1) with NA shown as missing."""
    import numpy as np
    if ext not in ('png','svg'):
        raise ValueError('Only png/svg supported')
    if matrix.empty or 'feature' not in matrix:
        raise ValueError('No heatmap data')
    cols=[c for c in matrix if c!='feature']
    data=matrix[cols].apply(pd.to_numeric,errors='coerce').to_numpy(dtype=float)
    if not cols:
        raise ValueError('No RAW runs in heatmap')
    log_data=np.ma.masked_invalid(np.log10(np.maximum(data,0)+1))
    fig,ax=plt.subplots(figsize=(max(8,min(22,len(cols)*1.2)),max(5,min(28,len(matrix)*0.25))))
    plot=ax.imshow(log_data,aspect='auto',interpolation='nearest')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols,rotation=60,ha='right',fontsize=8)
    if len(matrix)<=120:
        ax.set_yticks(range(len(matrix)))
        ax.set_yticklabels(matrix['feature'],fontsize=5)
    else:
        ax.set_ylabel('Peptides (too many labels to display)')
    fig.colorbar(plot,ax=ax,label='log10(MS Area + 1)')
    fig.tight_layout()
    out=io.BytesIO()
    fig.savefig(out,format=ext,dpi=180,bbox_inches='tight')
    plt.close(fig)
    return out.getvalue()


def simple_bar_figure_bytes(frame: pd.DataFrame, x: str, y: str, title: str, ext: str='png') -> bytes:
    """Export a reproducible static chart without a browser/Kaleido dependency."""
    if ext not in ('png','svg'):
        raise ValueError('Only png/svg supported')
    sub=frame.dropna(subset=[x,y]).copy()
    fig,ax=plt.subplots(figsize=(max(8,min(23,len(sub)*0.55)),5))
    ax.bar(range(len(sub)),pd.to_numeric(sub[y],errors='coerce').fillna(0).tolist())
    ax.set_xticks(range(len(sub)))
    ax.set_xticklabels(sub[x].astype(str).tolist(),rotation=55,ha='right',fontsize=8)
    ax.set_ylabel(y)
    ax.set_title(title)
    fig.tight_layout()
    out=io.BytesIO()
    fig.savefig(out,format=ext,dpi=180,bbox_inches='tight')
    plt.close(fig)
    return out.getvalue()
