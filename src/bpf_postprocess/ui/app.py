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
from bpf_postprocess.processing.detailed_comparison import (
    run_catalog, compatible_run_ids, select_comparison_runs,
    annotated_area_matrix, detailed_reference_comparison, detailed_component_records,
)
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
                st.session_state["comparison_session_epoch"] = st.session_state.get("comparison_session_epoch",0)+1
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
        st.session_state["run_labels"], hide_index=True, width="stretch",
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
st.dataframe(make_runs(raw),hide_index=True,width="stretch")
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
editor=st.data_editor(selectable,hide_index=True,width="stretch",height=390,
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

# Comparison-specific RAW selection never changes the explorer or PTM filters.
# In particular, reduced and non-reduced runs must not share an area comparison.
comparison_catalog = run_catalog(selected)
comparison_reference = None
comparison_targets = []
comparison_data = selected.iloc[0:0].copy()
comparison_matrix = pd.DataFrame(columns=["feature"])
comparison_wide = pd.DataFrame(columns=["feature", "peptide_sequence", "modification", "site", "positions", "protein"])
comparison_details = pd.DataFrame()
comparison_summary = pd.DataFrame()
comparison_areas = pd.DataFrame()

with st.expander("3 · Choose RAW files for peptide comparison", expanded=True):
    st.write("Select **one reference RAW** and any number of **comparison RAW files**. "
             "Only runs with the same preparation are eligible. "
             "This selection controls the two comparison tabs and their exports, "
             "without removing runs from Component Explorer or Modification Analysis.")
    if comparison_catalog.empty:
        st.info("No selected components are available. Check your filters or component selection.")
    else:
        run_options = comparison_catalog["run"].tolist()
        comparison_reference = st.selectbox("Reference RAW run", run_options, key="comparison_reference_run")
        ref_prep = comparison_catalog.loc[
            comparison_catalog["run"] == comparison_reference, "preparation"
        ].iloc[0]
        possible_comparisons = compatible_run_ids(comparison_catalog, comparison_reference)
        # Key depends on reference/available runs; switching the reference resets
        # the picker rather than keeping stale, incompatible RAW selections.
        picker_signature = hashlib.sha1((
            str(st.session_state.get("comparison_session_epoch", 0)) + "|" +
            comparison_reference + "|" + "|".join(possible_comparisons)
        ).encode()).hexdigest()[:12]
        comparison_targets = st.multiselect(
            "Comparison RAW runs (select or deselect manually)",
            options=possible_comparisons, default=possible_comparisons,
            key=f"comparison_targets_{picker_signature}",
        )
        st.caption(f"Reference preparation: **{ref_prep}**. "
                   "Only RAW files with exactly this preparation are shown in the comparator picker.")
        if ref_prep == "Unspecified":
            st.warning("Reference preparation is Unspecified. Assign Reduced or Non-reduced "
                       "in section 1 to reliably avoid comparing different preparation types.")
        if not comparison_targets:
            st.info("Choose at least one comparison RAW file for fold-change tables and heatmaps.")
        selected_run_ids = [comparison_reference, *comparison_targets]
        st.dataframe(comparison_catalog.loc[comparison_catalog["run"].isin(selected_run_ids)],
                     width="stretch", hide_index=True)
        comparison_data = select_comparison_runs(selected, comparison_reference, comparison_targets)
        comparison_matrix = compare_runs(comparison_data, run_ids=selected_run_ids)
        if not comparison_matrix.empty:
            existing_columns = [c for c in selected_run_ids if c in comparison_matrix.columns]
            comparison_matrix = comparison_matrix[["feature", *existing_columns]]
        comparison_wide = annotated_area_matrix(selected, comparison_reference, comparison_targets)
        comparison_details = detailed_component_records(selected, comparison_reference, comparison_targets)
        comparison_summary = detailed_reference_comparison(selected, comparison_reference, comparison_targets)
        comparison_areas = long_peptide_areas(comparison_data)

areas=long_peptide_areas(selected)
mods=modification_areas(selected)
unknown=unidentified_summary(selected)

peptide_tab,comparison_tab,mod_tab,unknown_tab,report_tab=st.tabs([
    "Peptide comparison","Multi-sample & replicates","Modifications","Unidentified","Export report"
])
with peptide_tab:
    st.write("**Only your chosen RAW runs** are shown below. Rows are matched using the full "
             "BF peptide identification. MS areas are summed for duplicate components/charge "
             "states; displayed RT is from the largest-area component per run. "
             "All individual RT, m/z and charge measurements are preserved in the detailed table.")
    if comparison_wide.empty:
        st.info("No identified peptides in the selected RAW runs.")
    else:
        st.dataframe(comparison_wide, width="stretch", hide_index=True, height=370)
    downloads(comparison_wide, "comparison_annotated_matrix", title="Annotated comparison")
    st.write("**Every individual identified component**, including RT, sequence, "
             "modification, site, confidence, charge, m/z and source run.")
    st.dataframe(comparison_details, width="stretch", hide_index=True, height=320)
    downloads(comparison_details, "comparison_all_peptide_components", title="Individual components")
    if not comparison_areas.empty:
        features = sorted(comparison_areas["feature"].dropna().unique().tolist())
        feature = st.selectbox("Peptide / identification to plot", features)
        plot_data = comparison_areas.loc[comparison_areas["feature"] == feature].copy()
        fig = px.bar(plot_data, x="raw_file", y="ms_area", color="sample",
                     hover_data=["experiment", "preparation", "component_count"],
                     title=feature, labels={"ms_area": "MS Area", "raw_file": "RAW run"})
        st.plotly_chart(fig, width="stretch")
        downloads(plot_data, "selected_peptide_areas")
        st.download_button("⬇ Figure PNG", static_bar_figure_bytes(plot_data, feature, "png"),
                           "peptide_comparison.png", "image/png")
        st.download_button("⬇ Figure SVG", static_bar_figure_bytes(plot_data, feature, "svg"),
                           "peptide_comparison.svg", "image/svg+xml")
        st.download_button("⬇ Interactive HTML", fig.to_html(include_plotlyjs=True).encode("utf-8"),
                           "peptide_comparison.html", "text/html")

with comparison_tab:
    st.write("**Reference and comparison RAW files are selected in section 3 above.** "
             "The table includes peptide sequence, modification, site, positions, protein, "
             "reference/comparison areas, RT, m/z, charge, confidence and fold change. "
             "RT/mz/charge/confidence refer to the dominant-area component in each run; "
             "the detailed per-component table is downloadable above. "
             "Fold change is not normalized for injection or digestion efficiency.")
    if comparison_reference is None or not comparison_targets:
        st.info("Select a reference and one or more comparison RAW runs in section 3.")
    elif comparison_summary.empty:
        st.info("No identified peptide comparisons could be calculated for these RAW files.")
    else:
        st.dataframe(comparison_summary, width="stretch", hide_index=True, height=390)
    downloads(comparison_summary, "selected_reference_comparison", title="Reference comparison")
    columns = [c for c in comparison_matrix.columns if c != "feature"]
    if not comparison_matrix.empty and len(columns) >= 2:
        st.write("**Heatmap** — log10(MS Area + 1) for the selected RAW runs only. "
                 "Missing identifications remain missing, not zero.")
        heatmap = comparison_matrix.set_index("feature")[columns]
        heatmap = heatmap.loc[heatmap.notna().sum(axis=1).sort_values(ascending=False).index]
        import numpy as np
        heatmap_log = np.log10(heatmap.clip(lower=0) + 1)
        heat = px.imshow(heatmap_log, aspect="auto", color_continuous_scale="Viridis",
                         labels={"x": "RAW run", "y": "Identification", "color": "log10(Area+1)"})
        heat.update_layout(height=min(1100, max(350, 22 * len(heatmap_log))))
        st.plotly_chart(heat, width="stretch")
        st.download_button("⬇ Heatmap HTML", heat.to_html(include_plotlyjs=True).encode("utf-8"),
                           "area_heatmap.html", "text/html")
        ha, hb = st.columns(2)
        with ha:
            st.download_button("⬇ Heatmap PNG", heatmap_figure_bytes(comparison_matrix, "png"),
                               "area_heatmap.png", "image/png")
        with hb:
            st.download_button("⬇ Heatmap SVG", heatmap_figure_bytes(comparison_matrix, "svg"),
                               "area_heatmap.svg", "image/svg+xml")
        downloads(comparison_matrix, "selected_heatmap_source_areas", title="Heatmap data")
    st.write("**Replicate statistics** — based only on your selected RAW runs; "
             "n=1 has undefined SD and RSD%.")
    group_by = st.radio("Group replicates by", ["sample", "sample_group"],
                        format_func=lambda x: "Sample ID" if x == "sample" else "Group (DS/DP/stressed/control)",
                        horizontal=True)
    stats = replicate_statistics(comparison_data, by=group_by)
    st.dataframe(stats, width="stretch", hide_index=True, height=300)
    downloads(stats, "selected_replicate_stats", title="Replicate statistics")
    if not comparison_areas.empty:
        chosen_feature = st.selectbox("Peptide for grouped comparison",
                                      sorted(comparison_areas["feature"].unique()), key="group_feature")
        group_data = comparison_areas.loc[comparison_areas["feature"] == chosen_feature].copy()
        groupfig = px.bar(group_data, x="sample", y="ms_area", color="preparation", barmode="group",
                          hover_data=["raw_file", "experiment"], title=chosen_feature)
        st.plotly_chart(groupfig, width="stretch")
        st.download_button("⬇ Grouped plot PNG", simple_bar_figure_bytes(
            group_data, "raw_file", "ms_area", "Grouped peptide areas", "png"),
            "grouped_peptide_areas.png", "image/png")
        st.download_button("⬇ Grouped plot HTML", groupfig.to_html(include_plotlyjs=True).encode("utf-8"),
                           "grouped_peptide_areas.html", "text/html")

with mod_tab:
    st.write("BF modification labels and MS Area per RAW run. The native PMF reader derives most labels and "
             "peptide sequences from identification text, but does **not** reconstruct disulfide sites or connectivity.")
    if selected["source_type"].eq("resultsbpf").any():
        st.warning("PMF modification fields are inferred, not decoded from dedicated PMF metadata. "
                   "Site-specific disulfide assignments require the BF CSV export or further parser validation.")
    st.dataframe(mods,width="stretch",hide_index=True,height=270)
    downloads(mods,"modification_areas",title="Modification areas")
    st.divider()
    st.write("**Relative modified-peptide MS Area (%)** — comparison of selected PTMs against unmodified "
             "forms of the **same sequence/position-resolved peptide backbone**. This is not molar occupancy. "
             "Gas-phase fragments, salts/adducts and unrelated species are excluded unless explicitly selected.")
    possible=sorted(set(selected.loc[selected["identified"],"modification"].fillna("").astype(str))-{ "", "None" })
    targets=st.multiselect("PTM labels to include in denominator and show as modified",
                           possible,default=[x for x in DEFAULT_PTMS if x in possible])
    fractions=modification_fractions(selected,targets)
    st.dataframe(fractions,width="stretch",hide_index=True,height=280)
    downloads(fractions,"relative_modification_ms_area",title="Relative modification")
    st.caption("If no unmodified reference identification exists for a backbone/run, the percentage is blank; "
               "a missing peptide is not assigned 0% or 100% occupancy.")
    if not fractions.empty:
        ptm_to_plot = st.selectbox(
            "Modification chart", sorted(fractions["modification"].dropna().unique())
        )
        # Use the *modification fraction* table here, not the peptide-area
        # table from the previous tab. Undefined ratios must not become 0%.
        fraction_plot_data = fractions.loc[
            fractions["modification"] == ptm_to_plot
        ].copy()
        fraction_plot_data["modified_ms_area_pct"] = pd.to_numeric(
            fraction_plot_data["modified_ms_area_pct"], errors="coerce"
        )
        fraction_plot_data = fraction_plot_data.dropna(
            subset=["modified_ms_area_pct"]
        )
        if fraction_plot_data.empty:
            st.info(
                "No quantifiable MS-area percentages for this modification. "
                "The matching unmodified peptide may be absent or its denominator may be zero."
            )
        else:
            fig_mod = px.bar(
                fraction_plot_data,
                x="raw_file",
                y="modified_ms_area_pct",
                color="backbone",
                barmode="group",
                hover_data=["sample", "preparation"],
                labels={"modified_ms_area_pct": "Modified MS Area (%)"},
            )
            st.plotly_chart(fig_mod, width="stretch")

with unknown_tab:
    st.write("Unknown components are not identified peptides. Never match unknowns across experiments by component No. alone.")
    st.dataframe(unknown,width="stretch",hide_index=True)
    downloads(unknown,"unidentified_summary",title="Unknown summary")
    unknown_rows=selected[~selected["identified"]].copy()
    st.dataframe(unknown_rows,width="stretch",hide_index=True,height=300)
    downloads(unknown_rows,"unidentified_components",title="Unknown components")

with report_tab:
    st.write("Export selected RAW component data and derived tables as one Excel workbook. No Thermo RAW payload included.")
    tables={"SelectedComponents":selected,"VisibleComponents":filtered,
            "ComparisonRuns":comparison_catalog.loc[comparison_catalog["run"].isin(
                [comparison_reference, *comparison_targets]
            )] if comparison_reference is not None else comparison_catalog,
            "DetailedPeptides":comparison_details,
            "AnnotatedMatrix":comparison_wide,"ReferenceComparison":comparison_summary,
            "PeptideAreas":comparison_areas,"PeptideMatrix":comparison_matrix,
            "Modifications":mods,"Unidentified":unknown,"Runs":make_runs(selected)}
    st.download_button("⬇ Download full Excel report",excel_bytes(tables),
                       file_name="biopharma_postprocessing_report.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    downloads(selected,"report_selected_components",title="Report components")

st.caption("Experimental research software. Validate native PMF and quantitative metrics against BioPharma Finder CSV exports.")
