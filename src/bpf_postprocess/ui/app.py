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
from bpf_postprocess.processing.id_types import displayed_id_types, id_type_summary, UNAVAILABLE_ID_TYPE
from bpf_postprocess.processing.selection import component_keys, selected_rows, update_exclusions
from bpf_postprocess.processing.peptide_selection import (
    apply_peptide_exclusions, filter_peptide_table, peptide_inventory, update_peptide_exclusions,
)
from bpf_postprocess.processing.comparison import compare_runs, long_peptide_areas
from bpf_postprocess.processing.multisample import reference_comparison, replicate_statistics
from bpf_postprocess.processing.detailed_comparison import (
    run_catalog, compatible_run_ids, select_comparison_runs,
    annotated_area_matrix, detailed_reference_comparison, detailed_component_records,
    run_labeled_reference_matrix,
)
from bpf_postprocess.processing.modifications import modification_areas
from bpf_postprocess.processing.relative_modifications import DEFAULT_PTMS, modification_fractions
from bpf_postprocess.processing.unidentified import unidentified_summary
from bpf_postprocess.ui.mapping_dialog import open_mapping_dialog
from bpf_postprocess.reports.exports import (
    csv_bytes, excel_bytes, static_bar_figure_bytes, heatmap_figure_bytes, simple_bar_figure_bytes,
)

st.set_page_config(page_title="BioPharma Finder Post-processing", layout="wide", page_icon="🧪")
st.title("BioPharma Finder · Peptide Mapping")
st.caption("Local analysis · multiple experiments and RAW runs · experimental PMF reader")
# Stable Streamlit-key CSS classes (not brittle internal DOM selectors).
st.markdown("""<style>
.st-key-peptide_comparison_panel, .st-key-multiplesample_comparison_panel {
    border: 3px solid #28a653 !important;
    border-radius: 17px !important;
    padding: 18px 18px 12px 18px !important;
}
</style>""", unsafe_allow_html=True)


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


def peptide_table_editor(table: pd.DataFrame, editor_id: str, *, height: int = 350) -> None:
    """Display synchronized peptide checkboxes without dropping deselected rows.

    All editor instances use a single exclusion set. Every peptide is a full BF
    identification; deselecting it affects *all* RAW runs and all report tabs.
    Each editor has its own widget key, refreshed after an Apply/Select All action.
    """
    if table.empty or "feature" not in table:
        st.info("No identified peptides are available in the current selection.")
        return
    if table["feature"].duplicated().any():
        raise ValueError("Peptide selection table must contain one row per identification")
    excluded_features = st.session_state.setdefault("excluded_peptide_features", set())
    labels = table["feature"].astype(str).str.strip()
    visible = table.copy().reset_index(drop=True)
    visible.insert(0, "Selected", ~labels.reset_index(drop=True).isin(excluded_features))
    signature = hashlib.sha1((
        str(st.session_state.get("comparison_session_epoch",0)) + "|" +
        str(st.session_state.get("peptide_select_version",0)) + "|" +
        "\x1f".join(labels.tolist())
    ).encode()).hexdigest()[:16]
    editor = st.data_editor(
        visible, hide_index=True, width="stretch", height=height,
        disabled=[column for column in visible if column != "Selected"],
        column_config={"Selected": st.column_config.CheckboxColumn(
            "Use in analysis", required=True,
            help="Affects the full peptide identification across all RAW runs and tabs")},
        key=f"peptide_editor_{editor_id}_{signature}",
    )
    st.caption(f"Selected {int(visible['Selected'].sum())} of {len(visible)} peptides. "
               "Uncheck or recheck rows and click Apply selection. "
               "Changes are shared with all other tabs, plots and report exports.")
    a,b,c = st.columns(3)
    if a.button("Apply peptide selection", type="primary", key=f"apply_peptides_{editor_id}"):
        st.session_state["excluded_peptide_features"] = update_peptide_exclusions(
            visible["feature"], editor["Selected"], excluded_features
        )
        st.session_state["peptide_select_version"] = st.session_state.get("peptide_select_version",0)+1
        st.rerun()
    if b.button("Select all visible", key=f"all_peptides_{editor_id}"):
        excluded_features.difference_update(labels.tolist())
        st.session_state["peptide_select_version"] = st.session_state.get("peptide_select_version",0)+1
        st.rerun()
    if c.button("Deselect all visible", key=f"no_peptides_{editor_id}"):
        excluded_features.update(labels.tolist())
        st.session_state["peptide_select_version"] = st.session_state.get("peptide_select_version",0)+1
        st.rerun()


def unidentified_component_editor(source: pd.DataFrame) -> None:
    """Unidentified signals cannot be selected as peptides, only as components."""
    candidates = source.loc[~source["identified"]].copy().reset_index(drop=True)
    if candidates.empty:
        st.info("No unidentified components match the current filters.")
        return
    excluded_components = st.session_state.setdefault("excluded_component_keys",set())
    keys = component_keys(candidates)
    grid = candidates.copy()
    grid.insert(0,"Selected",~keys.isin(excluded_components))
    front = ["Selected","experiment","sample","raw_file","component_no","id_type","ms_area",
             "rt_min","mz","charge","mono_mass_da","confidence_pct"]
    grid = grid[[c for c in front if c in grid] + [c for c in grid if c not in front]]
    signature = hashlib.sha1((
        str(st.session_state.get("comparison_session_epoch",0)) + "|" +
        str(st.session_state.get("select_version",0)) + "|" +
        "\x1f".join(keys.tolist())
    ).encode()).hexdigest()[:16]
    edit = st.data_editor(
        grid, hide_index=True, width="stretch", height=340,
        disabled=[c for c in grid if c!="Selected"],
        column_config={"Selected":st.column_config.CheckboxColumn("Use in analysis",required=True)},
        key=f"unknown_component_editor_{signature}",
    )
    st.caption("Unknown signals are not identified peptides. Selection is per RAW/component; "
               "changes update Component Explorer and the exported unknown summary.")
    a,b,c=st.columns(3)
    if a.button("Apply unknown selection",type="primary",key="unknown_apply"):
        st.session_state["excluded_component_keys"] = update_exclusions(
            candidates,edit["Selected"],excluded_components
        )
        st.session_state["select_version"] = st.session_state.get("select_version",0)+1
        st.rerun()
    if b.button("Select all visible",key="unknown_all"):
        excluded_components.difference_update(keys.tolist())
        st.session_state["select_version"] = st.session_state.get("select_version",0)+1
        st.rerun()
    if c.button("Deselect all visible",key="unknown_none"):
        excluded_components.update(keys.tolist())
        st.session_state["select_version"] = st.session_state.get("select_version",0)+1
        st.rerun()


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
                st.session_state["excluded_peptide_features"]=set()
                st.session_state["peptide_select_version"]=0
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
# Normalize older Streamlit sessions and PMF imports: missing ID Type does not
# mean NoID. The original BF CSV supplies validated run-specific categories.
raw["id_type"] = displayed_id_types(raw["id_type"])
if raw["id_type"].eq(UNAVAILABLE_ID_TYPE).any():
    unavailable = int(raw["id_type"].eq(UNAVAILABLE_ID_TYPE).sum())
    st.warning(
        f"ID Type is unavailable for {unavailable:,} component-run records. "
        "Native PMF/.resultsbpf ID Type is not yet decoded. For MS2/Full/NoID "
        "classification, import the original BioPharma Finder CSV export. "
        "Unavailable is not the same as NoID."
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
raw_id_types=displayed_id_types(raw["id_type"])
options=sorted(raw_id_types.unique().tolist())
selected_id_types=st.sidebar.multiselect(
    "ID Type (per RAW)", options, default=options,
    key=f"id_type_filter_{st.session_state.get('comparison_session_epoch',0)}",
    help="The per-RAW BioPharma Finder ID Type (MS2 / Full / NoID). "
         "The global 'Both' flag is retained separately. PMF reader cannot yet decode this field; "
         "Unavailable is not NoID. All types are included by default.",
)
if UNAVAILABLE_ID_TYPE in options:
    st.sidebar.caption("Some imported PMF records have ID Type unavailable. "
                       "They are not automatically classified as NoID.")
filtered=filter_components(raw,samples=selected_samples,preparations=selected_preps, identified=which_id,
                           peptide_query=peptide_q,modification_query=mod_q,min_area=min_area,
                           min_confidence_pct=min_conf,id_types=selected_id_types).reset_index(drop=True)

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
    "identified", "id_type", "id_type_overall", "confidence_pct", "peptide_sequence", "modification", "site",
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
st.caption(f"**Applied component selection:** {len(selected):,} of {len(filtered):,} visible rows. "
           "You can also deselect whole peptide identifications within each analysis tab; "
           "those exclusions are shared across tabs and report downloads.")
with st.expander("Download Component Explorer tables",expanded=True):
    downloads(filtered,"component_explorer_visible",title="Visible components")
    downloads(selected,"component_explorer_selected",title="Selected components")

# Component Explorer selection and peptide-level exclusion are separate.
# Peptide selection is one shared set across all analysis tabs and RAW runs.
excluded_peptides = st.session_state.setdefault("excluded_peptide_features",set())
analysis_selected = apply_peptide_exclusions(selected,excluded_peptides)
peptide_candidates = peptide_inventory(selected)

# RAW selection remains based on component-selected rows so that deselected
# peptides can be reselected without losing their RAW runs from the picker.
comparison_catalog = run_catalog(selected)
comparison_reference = None
comparison_targets = []
selected_run_ids = []
comparison_data = selected.iloc[0:0].copy()
comparison_matrix = pd.DataFrame(columns=["feature"])
comparison_wide = pd.DataFrame(columns=["feature", "peptide_sequence", "modification", "site", "positions", "protein"])
comparison_wide_candidates = comparison_wide.copy()
comparison_labeled_candidates = comparison_wide.copy()
comparison_details = pd.DataFrame()
comparison_summary = pd.DataFrame()
comparison_labeled = pd.DataFrame()
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
        # Select run scope BEFORE peptide exclusion. Even if all identified
        # features in a RAW are unchecked, keep the RAW available for reselecting.
        comparison_data_candidates = select_comparison_runs(selected, comparison_reference, comparison_targets)
        comparison_data = apply_peptide_exclusions(comparison_data_candidates, excluded_peptides)
        comparison_matrix = compare_runs(comparison_data, run_ids=selected_run_ids)
        if not comparison_matrix.empty:
            existing_columns = [c for c in selected_run_ids if c in comparison_matrix.columns]
            comparison_matrix = comparison_matrix[["feature", *existing_columns]]
        comparison_wide_candidates = annotated_area_matrix(selected, comparison_reference, comparison_targets)
        comparison_labeled_candidates = run_labeled_reference_matrix(selected, comparison_reference, comparison_targets)
        comparison_details_candidates = detailed_component_records(selected, comparison_reference, comparison_targets)
        comparison_summary_candidates = detailed_reference_comparison(selected, comparison_reference, comparison_targets)
        # Derived tables use the same BF identification exclusion key, without
        # deleting candidate rows needed to re-select peptides in any tab.
        comparison_wide = filter_peptide_table(comparison_wide_candidates, excluded_peptides)
        comparison_labeled = filter_peptide_table(comparison_labeled_candidates, excluded_peptides)
        comparison_details = filter_peptide_table(comparison_details_candidates, excluded_peptides)
        comparison_summary = filter_peptide_table(comparison_summary_candidates, excluded_peptides)
        comparison_areas = long_peptide_areas(comparison_data)

areas=long_peptide_areas(analysis_selected)
mods=modification_areas(analysis_selected)
unknown=unidentified_summary(analysis_selected)

peptide_tab,comparison_tab,mod_tab,unknown_tab,report_tab=st.tabs([
    "Peptide comparison","Multi-sample & replicates","Modifications","Unidentified","Export report"
])
with peptide_tab:
    st.write("**Only your chosen RAW runs** are shown below. Rows are matched using the full "
             "BF peptide identification. MS areas are summed for duplicate components/charge "
             "states; displayed RT is from the largest-area component per run. "
             "All individual RT, m/z and charge measurements are preserved in the detailed table.")
    with st.container(border=True, key="peptide_comparison_panel"):
        st.write("**Select peptides directly in this table** (all checked initially). "
                 "Selection applies to every RAW and is shared across analysis tabs.")
        peptide_table_editor(comparison_wide_candidates,"peptide_comparison",height=370)
        downloads(comparison_wide, "comparison_annotated_matrix", title="Annotated comparison")
        st.caption("Mapping reads only the **applied** peptide checkboxes, component selection and "
                   "filters. Click Apply peptide selection first if you have just changed any checkboxes.")
        if st.button("🧬 Mapping · selected peptides", key="mapping_open_peptide_comparison",
                     disabled=comparison_catalog.empty):
            open_mapping_dialog(analysis_selected, dataset.data, dataset.experiments, comparison_catalog["run"].tolist())
    st.write("**Every individual identified component**, including ID Type (RAW and overall), RT, sequence, "
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
    st.write("**Wide comparison — one row per peptide.** Every MS Area, RT, "
             "Delta, Fold and Change (%) column is labeled with its RAW file; "
             "the reference is marked **[REFERENCE]**. "
             "Each additional RAW adds its own columns without mixing samples.")
    if comparison_reference is None or not comparison_targets:
        st.info("Select a reference and one or more comparison RAW runs in section 3.")
    with st.container(border=True, key="multiplesample_comparison_panel"):
        st.write("**Choose peptides to include**; all checked initially. "
                 "Excluded peptides can be rechecked here at any time.")
        peptide_table_editor(comparison_labeled_candidates,"multi_sample",height=390)
        downloads(comparison_labeled, "reference_comparison_by_raw", title="RAW-labeled comparison")
        st.caption("Maps use **only applied selections** and build one protein coverage map "
                   "per chosen RAW. No unselected peptides contribute.")
        if st.button("🧬 Mapping · selected peptides", key="mapping_open_multi_sample",
                     disabled=comparison_catalog.empty):
            open_mapping_dialog(analysis_selected, dataset.data, dataset.experiments, comparison_catalog["run"].tolist())
    with st.expander("Detailed long-format table (one peptide × comparison RAW)"):
        st.caption("The reference_run and comparison_run columns explicitly identify both files "
                   "for every row; useful for sorting, statistics and software import.")
        st.dataframe(comparison_summary, width="stretch", hide_index=True, height=320)
        downloads(comparison_summary, "selected_reference_comparison", title="Long reference comparison")
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
    st.write("**Choose peptide identifications for modification analysis**. "
             "A checkbox removes that peptide from the modification-area and "
             "relative-fraction calculations across all RAW runs.")
    peptide_table_editor(peptide_candidates,"modifications",height=330)
    st.write("**Modification MS Area summary (selected peptides)** — ID types are aggregated per RAW.")
    st.dataframe(mods,width="stretch",hide_index=True,height=270)
    downloads(mods,"modification_areas",title="Modification areas")
    st.divider()
    st.write("**Relative modified-peptide MS Area (%)** — comparison of selected PTMs against unmodified "
             "forms of the **same sequence/position-resolved peptide backbone**. This is not molar occupancy. "
             "Gas-phase fragments, salts/adducts and unrelated species are excluded unless explicitly selected.")
    possible=sorted(set(analysis_selected.loc[analysis_selected["identified"],"modification"].fillna("").astype(str))-{ "", "None" })
    targets=st.multiselect("PTM labels to include in denominator and show as modified",
                           possible,default=[x for x in DEFAULT_PTMS if x in possible])
    fractions=modification_fractions(analysis_selected,targets)
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
    st.write("Unknown components have no peptide identification, so they are selected "
             "**individually** (experiment + RAW + component number), not by peptide ID.")
    unidentified_component_editor(filtered)
    st.write("**Selected unidentified component summary** — NoID for CSV, Unavailable for undecoded PMF.")
    st.dataframe(unknown,width="stretch",hide_index=True)
    downloads(unknown,"unidentified_summary",title="Unknown summary")
    unknown_rows=analysis_selected.loc[~analysis_selected["identified"]].copy()
    downloads(unknown_rows,"unidentified_components",title="Selected unknown components")

with report_tab:
    st.write("**Review peptide selection before export.** All peptides are selected by default; "
             "changes made here apply to all analytical tabs and the report.")
    peptide_table_editor(peptide_candidates,"export_report",height=310)
    st.write("**ID Type distribution** (selected component records; filterable with Confidence and MS Area).")
    st.dataframe(id_type_summary(analysis_selected), width="stretch", hide_index=True, height=220)
    downloads(id_type_summary(analysis_selected), "id_type_summary", title="ID Type summary")
    st.write("Export selected component data and derived tables as one Excel workbook. "
             "Raw MS files are never included. The analytical sheets contain only selected "
             "peptides and components. **VisibleComponents** and **PeptideSelection** are "
             "audit sheets that deliberately retain excluded rows or their checkbox state "
             "so the selection can be reviewed.")
    export_selection=peptide_candidates.copy()
    if not export_selection.empty:
        export_selection.insert(0,"Selected",~export_selection["feature"].isin(excluded_peptides))
    tables={"SelectedComponents":analysis_selected,"VisibleComponents":filtered,
            "IDTypeSummary":id_type_summary(analysis_selected),
            "PeptideSelection":export_selection,
            "ComparisonRuns":comparison_catalog.loc[comparison_catalog["run"].isin(
                [comparison_reference, *comparison_targets]
            )] if comparison_reference is not None else comparison_catalog,
            "DetailedPeptides":comparison_details,
            "AnnotatedMatrix":comparison_wide,"LabeledComparison":comparison_labeled,
            "ReferenceComparison":comparison_summary,
            "PeptideAreas":comparison_areas,"PeptideMatrix":comparison_matrix,
            "Modifications":mods,"Unidentified":unknown,"Runs":make_runs(analysis_selected)}
    st.download_button("⬇ Download full Excel report",excel_bytes(tables),
                       file_name="biopharma_postprocessing_report.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    downloads(analysis_selected,"report_selected_components",title="Report components")

st.caption("Experimental research software. Validate native PMF and quantitative metrics against BioPharma Finder CSV exports.")
