"""Streamlit modal for selected-only, per-RAW sequence coverage maps."""
from __future__ import annotations

import hashlib
import re

import pandas as pd
import streamlit as st

from bpf_postprocess.processing.coverage import ProteinChain, reference_chains, calculate_coverage
from bpf_postprocess.reports.coverage_map import coverage_figure_bytes, report_tables, mapping_package
from bpf_postprocess.reports.exports import csv_bytes, excel_bytes


@st.dialog("Sequence Coverage Map · selected peptides only", width="large")
def open_mapping_dialog(
    selected_rows: pd.DataFrame,
    all_imported_rows: pd.DataFrame,
    experiment_metadata: pd.DataFrame,
    selected_runs: list[str],
) -> None:
    """Open map for ONE filtered RAW, with optional batch export of filtered RAWs.

    selected_rows is the already component- and peptide-filtered dataframe, not
    the original dataset. all_imported_rows are allowed to seed the *reference
    sequence only*, never the mapped peptides or coverage numerator.
    """
    st.caption("Each map uses only the **applied peptide selection** and filters, "
               "restricted to the chosen RAW file. Other detected peptides are **not** added. "
               "If you changed checkboxes, click Apply peptide selection before Mapping.")
    # The RAW selector is independent of reference/comparator selection.
    # selected_runs comes from COMPONENT-scoped, globally filtered run metadata;
    # selected_rows is the stricter already-applied peptide/component selection.
    # Therefore a RAW with zero selected peptides remains selectable, with ZERO
    # coverage (never an automatic fallback to all imported peptides).
    available=list(dict.fromkeys(selected_runs))
    scoped=selected_rows.copy()
    if not available:
        st.warning("No RAW files remain under the current global filters/component selection.")
        return
    epoch=st.session_state.get('comparison_session_epoch',0)
    selector_key=hashlib.sha1('|'.join(available).encode()).hexdigest()[:12]
    run=st.selectbox("Choose RAW file (Reduced and Non-reduced mapped independently)",available,
                     key=f'mapping_run_{epoch}_{selector_key}')
    exp=run.split(' / ',1)[0]
    this_run=scoped.loc[(scoped['experiment'].astype(str)+' / '+scoped['raw_file'].astype(str))==run].copy()
    st.caption(f"Applied, identified component rows for this RAW: "
               f"**{int(this_run['identified'].sum())}**. No other RAW data contributes to coverage. "
               "RAW selection is independent of the reference/comparison picker.")

    reference=reference_chains(exp,all_imported_rows,experiment_metadata)
    if not reference:
        st.error("No protein reference could be recovered. Import a native PMF containing protein "
                 "sequences, or a BioPharma Finder CSV with positional peptide identifications.")
        return
    st.markdown("**Protein reference sequences**")
    st.caption("Native PMF: full reference from ExperimentSummary. CSV: reconstructed from "
               "identification positions across imported files, **not** used as coverage evidence. "
               "For CSV, confirm complete chain termini before interpreting percent coverage.")
    editable_chains=[]
    for chain in reference:
        input_key=f"map_seq_{epoch}_{exp}_{chain.number}"
        seq=st.text_input(f"{chain.number}:{chain.name} — sequence", value=chain.sequence,
                          key=input_key,help="Enter full sequence (one-letter AA codes) if the inferred reference is incomplete")
        valid=bool(re.fullmatch('[A-Za-z]+',re.sub(r'\s+','',seq)))
        confirmed=st.checkbox(f"Full-length sequence verified for {chain.name}",
                              value=chain.verified,key=f"map_confirm_{epoch}_{exp}_{chain.number}",
                              help="PMF sequences are already verified against experiment metadata; "
                                   "CSV sequences are inferred and require manual confirmation")
        if not valid:
            st.warning(f"Chain {chain.name} has missing (?) or invalid residues. "
                       "Map uses its current reference, but percent coverage is provisional.")
        editable_chains.append(ProteinChain(chain.number,chain.name,
            re.sub(r'\s+','',seq).upper(),
            source=chain.source if not confirmed else 'User-confirmed full-length reference',
            verified=confirmed and valid))
    st.markdown("**Mapping parameters**")
    pa,pb,pc=st.columns(3)
    min_area=pa.number_input("Minimum selected-peptide MS Area",min_value=0.0,value=0.0,
                             help="Uses summed MS Area of selected component/charge-state records in this RAW")
    min_conf=pb.number_input("Minimum confidence (%)",min_value=0.0,max_value=100.0,value=0.0,step=1.0,
                             help="0 includes MS-only and unknown-confidence identifications; upper confidence filter only if positive")
    max_mass=pc.number_input("Maximum peptide mass (Da)",min_value=0.0,value=0.0,step=500.0,
                             help="0 disables an additional mapping mass filter")
    pd_,pe,pf=st.columns(3)
    min_signal=pd_.number_input("Minimum relative MS signal (%)",min_value=0.0,max_value=100.0,value=0.0,step=0.5,
                                help="Our signal-based approximation, not Thermo's Recovery")
    min_overlap=pe.number_input("Minimum overlap signal (%)",min_value=0.0,max_value=100.0,value=0.0,step=0.5,
                                help="Our local overlap-area proxy, not Thermo's recovery algorithm")
    residues_per_row=int(pf.number_input("Residues per row",min_value=10,max_value=200,value=50,step=5))
    settings=dict(minimum_ms_area=min_area,minimum_confidence_pct=min_conf,
                  maximum_peptide_mass_da=max_mass,minimum_relative_signal_pct=min_signal,
                  minimum_overlap_signal_pct=min_overlap)
    try:
        report=calculate_coverage(scoped,run,editable_chains,**settings)
    except ValueError as exc:
        st.error(f"Cannot generate this sequence coverage map: {exc}")
        return
    if report.peptides.empty:
        st.warning("No **selected**, positive-area peptide passed the current mapping thresholds. "
                   "No excluded/unselected peptides have been substituted.")
    for warning in dict.fromkeys(report.warnings):
        st.warning(warning)
    st.markdown("**Coverage summary — current RAW**")
    st.dataframe(report.summary,width="stretch",hide_index=True)
    st.caption("MS peak area is an explicitly allocated fraction of selected peptide signal "
               "(a disulfide-linked peak shared across chains is divided equally). "
               "Molar abundance is **not** calculated from MS Area. "
               "Colors denote selected-peptide signal relative to the strongest selected feature in this RAW.")
    st.image(coverage_figure_bytes(report,'png',residues_per_row),width='stretch')
    st.markdown("**Peptide segments used in the map**")
    st.dataframe(report.segments,width="stretch",hide_index=True,height=220)
    image_cols=st.columns(3)
    for col,ext,mime in zip(image_cols,('png','svg','pdf'),('image/png','image/svg+xml','application/pdf')):
        with col:
            st.download_button(f"⬇ Mapping {ext.upper()}",data=coverage_figure_bytes(report,ext,residues_per_row),
                               file_name=f"sequence_coverage_{report.raw_file}.{ext}",mime=mime,
                               key=f'map_{ext}_{epoch}_{run}')
    dl1,dl2=st.columns(2)
    with dl1:
        st.download_button("⬇ Mapped segments CSV",data=csv_bytes(report.segments),
                           file_name=f"coverage_segments_{report.raw_file}.csv",mime="text/csv")
    with dl2:
        st.download_button("⬇ Complete Mapping Excel",data=excel_bytes(report_tables(report)),
                           file_name=f"sequence_coverage_{report.raw_file}.xlsx",
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    with st.expander("Batch export: RAW files available under global filters",expanded=False):
        selection_fields=[c for c in ('experiment','raw_file','component_no','identification','ms_area','rt_min')
                          if c in scoped.columns]
        selection_signature=hashlib.sha1(("|".join(available)+"|"+
            scoped[selection_fields].to_csv(index=False)+"|"+
            repr([(c.number,c.sequence,c.verified) for c in editable_chains])+"|"+
            repr(settings)+"|"+str(residues_per_row)+"|"+str(epoch)).encode()).hexdigest()
        if st.session_state.get("prepared_mapping_signature") != selection_signature:
            st.session_state.pop("prepared_selected_mapping_zip",None)
            st.session_state["prepared_mapping_signature"] = selection_signature
        st.caption("The ZIP contains one PNG, SVG, PDF, CSV tables and Excel workbook **per filtered RAW**. "
                   "Coverage uses only applied peptide selections and these mapping parameters, "
                   "not all imported peptides.")
        if st.button("Prepare ZIP for selected RAW files",key="prepare_mapping_zip"):
            batch=[]
            for candidate in available:
                candidate_exp=candidate.split(' / ',1)[0]
                refs=editable_chains if candidate_exp==exp else reference_chains(candidate_exp,all_imported_rows,experiment_metadata)
                if not refs:
                    st.warning(f"Skipped {candidate}: missing reference")
                    continue
                try:
                    batch.append(calculate_coverage(scoped,candidate,refs,**settings))
                except ValueError as exc:
                    st.warning(f"Skipped {candidate}: {exc}")
            if batch:
                st.session_state["prepared_selected_mapping_zip"] = mapping_package(
                    batch,residues_per_row=residues_per_row)
        if st.session_state.get("prepared_selected_mapping_zip"):
            st.download_button("⬇ Download selected RAW maps ZIP",
                               data=st.session_state["prepared_selected_mapping_zip"],
                               file_name="selected_raw_sequence_coverage_maps.zip",mime="application/zip",
                               key="download_mapping_batch_zip")
    st.caption("Thermo-inspired visualization, not a validated Thermo-equivalent recovery or abundance calculation. "
               "All displayed peptides were selected and mapped to the explicitly chosen RAW.")
