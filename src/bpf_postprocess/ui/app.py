"""Local Streamlit interface; all analyses consume the selected component rows."""
from __future__ import annotations
import hashlib
import sys
import tempfile
from pathlib import Path

SRC = Path(__file__).resolve().parents[2]
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pandas as pd
import plotly.express as px
import streamlit as st

from bpf_postprocess.model import make_runs
from bpf_postprocess.importers.loader import load_files
from bpf_postprocess.processing.filters import filter_components
from bpf_postprocess.processing.selection import component_keys, selected_rows, update_exclusions
from bpf_postprocess.processing.comparison import compare_runs, long_peptide_areas
from bpf_postprocess.processing.multisample import reference_comparison, replicate_statistics
from bpf_postprocess.processing.modifications import modification_areas
from bpf_postprocess.processing.relative_modifications import DEFAULT_PTMS, modification_fractions
from bpf_postprocess.processing.unidentified import unidentified_summary
from bpf_postprocess.reports.exports import (
    csv_bytes, excel_bytes, static_bar_figure_bytes, heatmap_figure_bytes, simple_bar_figure_bytes,
)

st.set_page_config(page_title="BioPharma Finder Post-processing", layout="wide", page_icon="🧪")
st.title("BioPharma Finder · Peptide Mapping")
st.caption("Local analysis · multiple experiments and RAW runs · experimental PMF reader")


def downloads(table: pd.DataFrame, base: str, *, title: str | None = None) -> None:
    """Explicit download controls beside every table (not browser-grid toolbar)."""
    a,b = st.columns([1,1])
    with a:
        st.download_button("⬇ Download CSV" if title is None else f"⬇ {title} CSV",
                           data=csv_bytes(table), file_name=f"{base}.csv", mime="text/csv",
                           key=f"csv_{base}")
    with b:
        st.download_button("⬇ Download Excel" if title is None else f"⬇ {title} Excel",
                           data=excel_bytes({base[:31]: table}), file_name=f"{base}.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                           key=f"xlsx_{base}")


with st.expander("Import experiments", expanded="dataset" not in st.session_state):
    uploaded = st.file_uploader("BioPharma Finder CSV, .resultsbpf or compact PMF ZIP",
                                type=["csv", "resultsbpf", "zip"], accept_multiple_files=True,
                                help="Choose one input representation per experiment.")
    st.info("Files are processed on this computer by Streamlit. They are not uploaded to GitHub.")
    if st.button("Import selected files", type="primary", disabled=not uploaded):
        with st.spinner("Reading components and RAW runs..."):
            try:
                with tempfile.TemporaryDirectory(prefix="bf_app_") as tmp:
                    paths=[]
                    for i,obj in enumerate(uploaded):
                        path=Path(tmp) / f"upload_{i}" / Path(obj.name).name
                        path.parent.mkdir(parents=True,exist_ok=True)
                        with path.open("wb") as fp:
                            while chunk := obj.read(1024*1024):
                                fp.write(chunk)
                        obj.seek(0)
                        paths.append(path)
                    dataset=load_files(paths)
                st.session_state["dataset"]=dataset
                st.session_state["run_labels"]=dataset.runs[[
                    "experiment","raw_file","sample","sample_group","preparation","replicate"
                ]].copy().fillna("")
                st.session_state["excluded_component_keys"]=set()
                st.session_state["select_version"]=0
                st.success(f"Imported {len(dataset.experiments)} experiment(s), {len(dataset.runs)} RAW run(s)")
            except Exception as exc:
                st.error(f"Import failed: {exc}")

if "dataset" not in st.session_state:
    st.write("Upload one or more experiments above. A BioPharma Finder file may contain several RAW runs.")
    st.stop()

dataset=st.session_state["dataset"]
for warning in dataset.warnings:
    st.warning(warning)

with st.expander("1 · Experiment / RAW run labels", expanded=True):
    st.write("Assign **Sample**, **Group** (e.g. DS, DP, stressed, control), **Preparation** and **Replicate**. "
             "Use the same Sample ID for injections you intend to treat as technical replicates.")
    edited_labels=st.data_editor(
        st.session_state["run_labels"], hide_index=True, use_container_width=True,
        disabled=["experiment","raw_file"], key="run_editor",
        column_config={
            "experiment":st.column_config.TextColumn("Experiment"),
            "raw_file":st.column_config.TextColumn("RAW file"),
            "sample":st.column_config.TextColumn("Sample ID",required=True),
            "sample_group":st.column_config.TextColumn("Group (optional)"),
            "preparation":st.column_config.SelectboxColumn("Preparation",options=["Unspecified","Non-reduced","Reduced","Other"],required=True),
            "replicate":st.column_config.TextColumn("Replicate (optional)")
        }
    )
    if st.button("Apply RAW run labels"):
        st.session_state["run_labels"]=edited_labels.copy()
        st.rerun()

# Reflect uncommitted labels in the visible session for immediate preview;
# applied labels become persistent on button click.
raw=dataset.data.drop(columns=["sample","sample_group","preparation","replicate"],errors="ignore").merge(
    edited_labels, on=["experiment","raw_file"], how="left", validate="many_to_one"
)
st.dataframe(make_runs(raw),hide_index=True,use_container_width=True)
with st.expander("Download RAW run overview"):
    downloads(make_runs(raw),"run_overview")

st.sidebar.header("Component filters")
samples=sorted(raw["sample"].dropna().astype(str).unique().tolist())
selected_samples=st.sidebar.multiselect("Samples",samples,default=samples)
preps=sorted(raw["preparation"].dropna().astype(str).unique().tolist())
selected_preps=st.sidebar.multiselect("Preparation",preps,default=preps)
run_pairs=raw[["experiment","raw_file"]].drop_duplicates().copy()
run_pairs["label"]=run_pairs["experiment"]+" / "+run_pairs["raw_file"]
run_choices=run_pairs["label"].tolist()
chosen_runs=st.sidebar.multiselect("RAW runs",run_choices,default=run_choices)
raw=raw[(raw["experiment"]+" / "+raw["raw_file"]).isin(chosen_runs)].copy()
which_id=st.sidebar.radio("Identification",["All","Identified","Unidentified"])
peptide_q=st.sidebar.text_input("Peptide sequence / identification contains")
mod_q=st.sidebar.text_input("Modification contains")
min_conf=st.sidebar.slider("Minimum confidence (%)",min_value=0.0,max_value=100.0,value=0.0,step=0.5,
                           help="At 0 no confidence filter. Above 0, unidentified and missing-confidence records are excluded.")
min_area=st.sidebar.number_input("Minimum MS Area",min_value=0.0,value=0.0)
filtered=filter_components(raw,samples=selected_samples,preparations=selected_preps, identified=which_id,
                           peptide_query=peptide_q,modification_query=mod_q,min_area=min_area,
                           min_confidence_pct=min_conf).reset_index(drop=True)

st.subheader("2 · Component explorer")
excluded=st.session_state.setdefault("excluded_component_keys",set())
selectable=filtered.copy()
selectable.insert(0,"Selected",~component_keys(selectable).isin(excluded))
selected_count=int(selectable["Selected"].sum())
a,b,c,d=st.columns(4)
a.metric("Experiments", filtered["experiment"].nunique())
b.metric("RAW runs", filtered[["experiment","raw_file"]].drop_duplicates().shape[0])
c.metric("Visible components",f"{len(filtered):,}")
d.metric("Selected components",f"{selected_count:,}")
st.caption("All visible components start selected. Uncheck rows you want to exclude, then click **Apply selection**. "
           "The selection is remembered when you change filters. Unidentified components have no peptide sequence or confidence score.")

# The version and visible keys change the editor key to prevent cached checkboxes
# from overwriting a deliberate Select/Deselect All or Apply action.
signature=hashlib.sha1((str(st.session_state.get('select_version',0))+"|"+
                         "|".join(component_keys(filtered).tolist())).encode()).hexdigest()[:16]
preferred_cols = [
    "Selected", "experiment", "sample", "raw_file", "preparation", "component_no",
    "identified", "confidence_pct", "peptide_sequence", "modification", "site",
    "identification", "ms_area", "rt_min", "mz", "charge", "positions", "protein"
]
selectable=selectable[[c for c in preferred_cols if c in selectable.columns] +
                      [c for c in selectable.columns if c not in preferred_cols]]
editor=st.data_editor(selectable,hide_index=True,use_container_width=True,height=390,
                      disabled=[c for c in selectable.columns if c!="Selected"],
                      column_config={"Selected":st.column_config.CheckboxColumn("Use in analysis",required=True)},
                      key=f"components_{signature}")
act_a,act_b,act_c=st.columns(3)
if act_a.button("Apply selection",type="primary"):
    st.session_state["excluded_component_keys"]=update_exclusions(filtered,editor["Selected"],excluded)
    st.session_state["select_version"]=st.session_state.get("select_version",0)+1
    st.rerun()
if act_b.button("Select all visible"):
    excluded.difference_update(component_keys(filtered).tolist())
    st.session_state["select_version"]=st.session_state.get("select_version",0)+1
    st.rerun()
if act_c.button("Deselect all visible"):
    excluded.update(component_keys(filtered).tolist())
    st.session_state["select_version"]=st.session_state.get("select_version",0)+1
    st.rerun()

selected=selected_rows(filtered,excluded)
st.caption(f"**Applied selection:** {len(selected):,} of {len(filtered):,} visible component-run rows. "
           "Tables and charts below use the applied selection only.")
with st.expander("Download Component Explorer tables",expanded=True):
    downloads(filtered,"component_explorer_visible",title="Visible components")
    downloads(selected,"component_explorer_selected",title="Selected components")

areas=long_peptide_areas(selected)
matrix=compare_runs(selected)
mods=modification_areas(selected)
unknown=unidentified_summary(selected)

peptide_tab,comparison_tab,mod_tab,unknown_tab,report_tab=st.tabs([
    "Peptide comparison","Multi-sample & replicates","Modifications","Unidentified","Export report"
])
with peptide_tab:
    st.write("Integrated MS areas, grouped by full BioPharma Finder identification. "
             "No normalization or occupancy inference; multiple components with the same ID are summed per run.")
    st.dataframe(matrix,use_container_width=True,hide_index=True,height=320)
    downloads(matrix,"peptide_comparison",title="Peptide matrix")
    if not areas.empty:
        features=sorted(areas["feature"].dropna().unique().tolist())
        feature=st.selectbox("Peptide / identification",features)
        plot_data=areas[areas["feature"]==feature]
        fig=px.bar(plot_data,x="raw_file",y="ms_area",color="sample",
                   hover_data=["experiment","preparation","component_count"],title=feature,
                   labels={"ms_area":"MS Area","raw_file":"RAW run"})
        st.plotly_chart(fig,use_container_width=True)
        downloads(plot_data,"selected_peptide_areas")
        st.download_button("⬇ Figure PNG",static_bar_figure_bytes(plot_data,feature,"png"),"peptide_comparison.png","image/png")
        st.download_button("⬇ Figure SVG",static_bar_figure_bytes(plot_data,feature,"svg"),"peptide_comparison.svg","image/svg+xml")
        st.download_button("⬇ Interactive HTML",fig.to_html(include_plotlyjs=True).encode("utf-8"),
                           "peptide_comparison.html","text/html")

with comparison_tab:
    st.write("Compare **all** selected RAW runs against one reference, plus replicate mean / SD / RSD%. "
             "Compare like-for-like preparations. Raw MS area differences are not concentration changes.")
    columns=[c for c in matrix.columns if c!="feature"]
    if len(columns)>=2:
        ref=st.selectbox("Reference RAW run",columns)
        ref_table=reference_comparison(selected,ref)
        st.dataframe(ref_table,use_container_width=True,hide_index=True,height=300)
        downloads(ref_table,"all_vs_reference",title="All vs reference")
    else:
        st.info("Select at least two RAW runs to compare to one reference.")
    if not matrix.empty and len(columns)>=2:
        st.write("**Heatmap** — log10(MS Area + 1); zeros are retained, missing peptides are not inferred as zero.")
        heatmap=matrix.set_index("feature")[columns]
        heatmap=heatmap.loc[heatmap.notna().sum(axis=1).sort_values(ascending=False).index]
        import numpy as np
        heatmap_log=np.log10(heatmap.clip(lower=0)+1)
        heat=px.imshow(heatmap_log,aspect="auto",color_continuous_scale="Viridis",
                       labels={"x":"RAW run","y":"Identification","color":"log10(Area+1)"})
        heat.update_layout(height=min(1100,max(350,22*len(heatmap_log))))
        st.plotly_chart(heat,use_container_width=True)
        st.download_button("⬇ Heatmap HTML",heat.to_html(include_plotlyjs=True).encode("utf-8"),
                           "area_heatmap.html","text/html")
        ha,hb=st.columns(2)
        with ha:
            st.download_button("⬇ Heatmap PNG",heatmap_figure_bytes(matrix,"png"),
                               "area_heatmap.png","image/png")
        with hb:
            st.download_button("⬇ Heatmap SVG",heatmap_figure_bytes(matrix,"svg"),
                               "area_heatmap.svg","image/svg+xml")
        downloads(matrix,"heatmap_source_areas",title="Heatmap data")
    st.write("**Replicate statistics** — one peptide area per RAW injection; n=1 produces undefined SD and RSD%.")
    group_by=st.radio("Group replicates by",["sample","sample_group"],
                      format_func=lambda x:"Sample ID" if x=="sample" else "Group (DS/DP/stressed/control)",horizontal=True)
    stats=replicate_statistics(selected,by=group_by)
    st.dataframe(stats,use_container_width=True,hide_index=True,height=300)
    downloads(stats,"replicate_stats",title="Replicate statistics")
    if not areas.empty:
        chosen_feature=st.selectbox("Peptide for grouped comparison",sorted(areas["feature"].unique()),key="group_feature")
        sub=areas[areas["feature"]==chosen_feature]
        groupfig=px.bar(sub,x="sample",y="ms_area",color="preparation",barmode="group",
                        hover_data=["raw_file","experiment"],title=chosen_feature)
        st.plotly_chart(groupfig,use_container_width=True)
        st.download_button("⬇ Grouped plot PNG",simple_bar_figure_bytes(
            sub, "raw_file", "ms_area", "Grouped peptide areas", "png"),
            "grouped_peptide_areas.png", "image/png")
        st.download_button("⬇ Grouped plot HTML",groupfig.to_html(include_plotlyjs=True).encode("utf-8"),
                           "grouped_peptide_areas.html","text/html")

with mod_tab:
    st.write("BF modification labels and MS Area per RAW run. The native PMF reader derives most labels and "
             "peptide sequences from identification text, but does **not** reconstruct disulfide sites or connectivity.")
    if selected["source_type"].eq("resultsbpf").any():
        st.warning("PMF modification fields are inferred, not decoded from dedicated PMF metadata. "
                   "Site-specific disulfide assignments require the BF CSV export or further parser validation.")
    st.dataframe(mods,use_container_width=True,hide_index=True,height=270)
    downloads(mods,"modification_areas",title="Modification areas")
    st.divider()
    st.write("**Relative modified-peptide MS Area (%)** — comparison of selected PTMs against unmodified "
             "forms of the **same sequence/position-resolved peptide backbone**. This is not molar occupancy. "
             "Gas-phase fragments, salts/adducts and unrelated species are excluded unless explicitly selected.")
    possible=sorted(set(selected.loc[selected["identified"],"modification"].fillna("").astype(str))-{ "", "None" })
    targets=st.multiselect("PTM labels to include in denominator and show as modified",
                           possible,default=[x for x in DEFAULT_PTMS if x in possible])
    fractions=modification_fractions(selected,targets)
    st.dataframe(fractions,use_container_width=True,hide_index=True,height=280)
    downloads(fractions,"relative_modification_ms_area",title="Relative modification")
    st.caption("If no unmodified reference identification exists for a backbone/run, the percentage is blank; "
               "a missing peptide is not assigned 0% or 100% occupancy.")
    if not fractions.empty:
        ptm_to_plot=st.selectbox("Modification chart",sorted(fractions["modification"].unique()))
        subset=fractions[fractions["modification"]==ptm_to_plot]
        st.plotly_chart(px.bar(sub,x="raw_file",y="modified_ms_area_pct",color="backbone",
                               barmode="group",hover_data=["sample","preparation"],
                               labels={"modified_ms_area_pct":"Modified MS Area (%)"}),use_container_width=True)

with unknown_tab:
    st.write("Unknown components are not identified peptides. Never match unknowns across experiments by component No. alone.")
    st.dataframe(unknown,use_container_width=True,hide_index=True)
    downloads(unknown,"unidentified_summary",title="Unknown summary")
    unknown_rows=selected[~selected["identified"]].copy()
    st.dataframe(unknown_rows,use_container_width=True,hide_index=True,height=300)
    downloads(unknown_rows,"unidentified_components",title="Unknown components")

with report_tab:
    st.write("Export selected RAW component data and derived tables as one Excel workbook. No Thermo RAW payload included.")
    tables={"SelectedComponents":selected,"VisibleComponents":filtered,
            "PeptideAreas":areas,"PeptideMatrix":matrix,"Modifications":mods,
            "Unidentified":unknown,"Runs":make_runs(selected)}
    st.download_button("⬇ Download full Excel report",excel_bytes(tables),
                       file_name="biopharma_postprocessing_report.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    downloads(selected,"report_selected_components",title="Report components")

st.caption("Experimental research software. Validate native PMF and quantitative metrics against BioPharma Finder CSV exports.")
