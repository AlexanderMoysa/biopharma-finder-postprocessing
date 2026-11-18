"""Run with: streamlit run src/bpf_postprocess/ui/app.py"""
from __future__ import annotations
import sys
from pathlib import Path

# Support running the uninstalled source checkout directly in Streamlit.
SRC = Path(__file__).resolve().parents[2]
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import tempfile
import pandas as pd
import plotly.express as px
import streamlit as st
from bpf_postprocess.model import make_runs
from bpf_postprocess.importers.loader import load_files
from bpf_postprocess.processing.filters import filter_components
from bpf_postprocess.processing.comparison import long_peptide_areas, compare_runs, ratio_table
from bpf_postprocess.processing.modifications import modification_areas
from bpf_postprocess.processing.unidentified import unidentified_summary
from bpf_postprocess.reports.exports import csv_bytes, excel_bytes, static_bar_figure_bytes

st.set_page_config(page_title="BioPharma Finder Post-processing", layout="wide", page_icon="🧪")
st.title("BioPharma Finder · Peptide Mapping")
st.caption("Local post-processing | multiple experiments and RAW runs | experimental native PMF import")

with st.expander("Import experiments", expanded="dataset" not in st.session_state):
    uploaded = st.file_uploader(
        "Upload BioPharma Finder CSV, native .resultsbpf, or compact PMF ZIP",
        type=["csv", "resultsbpf", "zip"], accept_multiple_files=True,
        help="Choose either CSV or resultsbpf for the same experiment, not both."
    )
    st.info("Files are processed by this local Streamlit process. It does not upload them to GitHub or any external service.")
    if st.button("Import selected files", type="primary", disabled=not uploaded):
        with st.spinner("Reading files and combining RAW runs..."):
            try:
                with tempfile.TemporaryDirectory(prefix="bf_app_") as tmp:
                    paths = []
                    for i, obj in enumerate(uploaded):
                        path = Path(tmp) / f"upload_{i}" / Path(obj.name).name
                        path.parent.mkdir(parents=True, exist_ok=True)
                        with path.open("wb") as f:
                            while True:
                                chunk = obj.read(1024 * 1024)
                                if not chunk:
                                    break
                                f.write(chunk)
                        obj.seek(0)
                        paths.append(path)
                    dataset = load_files(paths)
                    # Path labels from temporary folders are not retained; source names are.
                    st.session_state["dataset"] = dataset
                    st.session_state["run_labels"] = dataset.runs[[
                        "experiment", "raw_file", "sample", "preparation", "replicate"
                    ]].copy()
                    st.success(f"Imported {len(dataset.experiments)} experiments and {len(dataset.runs)} RAW runs")
            except Exception as err:
                st.error(f"Import failed: {err}")

if "dataset" not in st.session_state:
    st.markdown("Upload one or more files above to start. Each RAW file becomes a separate analysis run.")
    st.stop()

loaded = st.session_state["dataset"]
for note in loaded.warnings:
    st.warning(note)

with st.expander("1 · Experiments and RAW run labels", expanded=True):
    st.write("One experiment file may contain multiple RAW runs. Edit sample, preparation and replicate labels here.")
    edited = st.data_editor(
        st.session_state["run_labels"], hide_index=True, use_container_width=True,
        disabled=["experiment", "raw_file"], key="run_editor",
        column_config={
            "preparation": st.column_config.SelectboxColumn(
                "Preparation", options=["Unspecified", "Non-reduced", "Reduced", "Other"], required=True
            ),
            "sample": st.column_config.TextColumn("Sample ID", required=True),
            "replicate": st.column_config.TextColumn("Replicate (optional)"),
        }
    )
    if st.button("Apply run labels"):
        st.session_state["run_labels"] = edited.copy()
        st.rerun()
    raw = loaded.data.copy()
    raw = raw.drop(columns=["sample", "preparation", "replicate"]).merge(
        edited, on=["experiment", "raw_file"], how="left", validate="many_to_one"
    )
    st.dataframe(make_runs(raw), use_container_width=True, hide_index=True)

st.sidebar.header("Component filters")
samples = st.sidebar.multiselect("Samples", sorted(raw["sample"].dropna().unique().tolist()), default=sorted(raw["sample"].dropna().unique().tolist()))
preps = st.sidebar.multiselect("Preparation", sorted(raw["preparation"].dropna().unique().tolist()), default=sorted(raw["preparation"].dropna().unique().tolist()))
run_pairs = raw[["experiment", "raw_file"]].drop_duplicates()
run_pairs["label"] = run_pairs["experiment"] + " / " + run_pairs["raw_file"]
run_choices = run_pairs["label"].tolist()
chosen_run_labels = st.sidebar.multiselect("LC-MS runs", run_choices, default=run_choices)
selected_pairs = set(chosen_run_labels)
keep = (raw["experiment"] + " / " + raw["raw_file"]).isin(selected_pairs)
raw = raw[keep].copy()
id_type = st.sidebar.radio("Identification", ["All", "Identified", "Unidentified"], horizontal=False)
peptide_q = st.sidebar.text_input("Peptide sequence / identification contains")
mod_q = st.sidebar.text_input("Modification contains (literal text)")
min_area = st.sidebar.number_input("Minimum MS Area", min_value=0.0, value=0.0)
filtered = filter_components(raw, samples=samples, preparations=preps,
                             identified=id_type, peptide_query=peptide_q,
                             modification_query=mod_q, min_area=min_area)

st.subheader("2 · Component explorer")
a, b, c, d = st.columns(4)
a.metric("Experiments", filtered["experiment"].nunique())
b.metric("RAW runs", filtered[["experiment", "raw_file"]].drop_duplicates().shape[0])
c.metric("Component-run records", f"{len(filtered):,}")
d.metric("Identified records", f"{int(filtered['identified'].sum()):,}")
st.dataframe(filtered, use_container_width=True, hide_index=True, height=300)
st.download_button("Download filtered CSV", csv_bytes(filtered), file_name="filtered_components.csv", mime="text/csv")

areas = long_peptide_areas(filtered)
matrix = compare_runs(filtered)
mods = modification_areas(filtered)
unknown = unidentified_summary(filtered)

tab1, tab2, tab3, tab4 = st.tabs(["Peptide comparison", "Modifications", "Unidentified", "Export report"])
with tab1:
    st.write("Compare peptide identification and integrated MS areas across RAW runs. No automatic area normalization or occupancy calculation.")
    if matrix.empty:
        st.info("No identified peptide records with the current filters.")
    else:
        st.dataframe(matrix, use_container_width=True, hide_index=True, height=300)
        choices = [c for c in matrix.columns if c != "feature"]
        if len(choices) > 1:
            left, right = st.columns(2)
            reference = left.selectbox("Reference RAW", choices, index=0)
            treatment = right.selectbox("Comparison RAW", choices, index=min(1, len(choices)-1))
            if reference != treatment:
                ratios = ratio_table(matrix, reference, treatment)
                st.dataframe(ratios, use_container_width=True, hide_index=True, height=200)
                st.download_button("Download fold-change table CSV", csv_bytes(ratios), file_name="peptide_fold_change.csv", mime="text/csv")
        features = areas["feature"].unique().tolist()
        feature = st.selectbox("Peptide / identification to plot", features)
        plot_data = areas[areas["feature"] == feature]
        if not plot_data.empty:
            fig = px.bar(plot_data, x="raw_file", y="ms_area", color="sample", hover_data=["experiment", "preparation", "component_count"],
                         title=feature, labels={"ms_area": "MS Area (BF)", "raw_file": "LC-MS run"})
            st.plotly_chart(fig, use_container_width=True)
            st.download_button("Download interactive HTML figure", fig.to_html(include_plotlyjs=True).encode("utf-8"), file_name="peptide_comparison.html", mime="text/html")
            st.download_button("Download PNG", static_bar_figure_bytes(plot_data, feature, "png"), file_name="peptide_comparison.png", mime="image/png")
            st.download_button("Download SVG", static_bar_figure_bytes(plot_data, feature, "svg"), file_name="peptide_comparison.svg", mime="image/svg+xml")
with tab2:
    st.write("Literal modification labels from CSV, not a validated site occupancy measurement.")
    if filtered["source_type"].eq("resultsbpf").any():
        st.warning("The native PMF reader currently does not decode the separate Mod and Site fields. Use a BF CSV export for accurate modification category summaries.")
    st.dataframe(mods, use_container_width=True, hide_index=True)
    if not mods.empty:
        st.plotly_chart(px.bar(mods, x="modification", y="ms_area", color="sample", barmode="group", hover_data=["raw_file", "preparation"]), use_container_width=True)
with tab3:
    st.write("Unidentified components are kept as distinct signals. They are NOT automatically matched across experiments by component number.")
    st.dataframe(unknown, use_container_width=True, hide_index=True)
    st.dataframe(filtered[~filtered["identified"]], use_container_width=True, hide_index=True)
with tab4:
    st.write("Export all tables in one Excel workbook. Raw MS data are never included.")
    tables = {"Components": filtered, "PeptideAreas": areas, "Comparison": matrix,
              "Modifications": mods, "Unidentified": unknown, "Runs": make_runs(filtered)}
    st.download_button("Download Excel workbook", excel_bytes(tables), file_name="bf_postprocessing.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    st.download_button("Download comparison matrix CSV", csv_bytes(matrix), file_name="peptide_comparison.csv", mime="text/csv")

st.caption("Experimental research software. Verify quantitative conclusions against BioPharma Finder exports and appropriate method validation.")
