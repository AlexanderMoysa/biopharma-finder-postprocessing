"""Regression tests: user-selected peptides are the *only* coverage evidence."""
from __future__ import annotations

import io
from pathlib import Path
import zipfile

import pandas as pd
import pytest

from bpf_postprocess.model import as_table
from bpf_postprocess.importers.loader import load_files
from bpf_postprocess.processing.coverage import (
    ProteinChain, calculate_coverage, parse_peptide_regions, reference_chains,
)
from bpf_postprocess.processing.peptide_selection import apply_peptide_exclusions
from bpf_postprocess.processing.selection import component_keys, selected_rows
from bpf_postprocess.reports.coverage_map import coverage_figure_bytes, mapping_package, report_tables


def make_selected():
    rows=[]
    for raw, prep in [('NR.raw','Non-reduced'), ('DTT.raw','Reduced')]:
        rows.extend([
            dict(experiment='exp',raw_file=raw,preparation=prep,sample='S',identified=True,
                 component_no=1, identification='1:A1-C3 = 301m',peptide_sequence='ABC',
                 ms_area=100,rt_min=4.5,confidence_pct=99,mono_mass_da=301,id_type='MS2'),
            dict(experiment='exp',raw_file=raw,preparation=prep,sample='S',identified=True,
                 component_no=2, identification='2:D1-F3 = 333m',peptide_sequence='DEF',
                 ms_area=30,rt_min=5.5,confidence_pct=95,mono_mass_da=333,id_type='Full'),
            dict(experiment='exp',raw_file=raw,preparation=prep,sample='S',identified=False,
                 component_no=3,identification='',ms_area=500,id_type='NoID'),
        ])
    return as_table(rows)


def chains():
    return [ProteinChain(1,'A','ABCDE'),ProteinChain(2,'B','DEFGH')]


def test_only_applied_peptides_contribute_and_reference_inference_does_not_leak():
    frame=make_selected()
    selected=apply_peptide_exclusions(frame, {'2:D1-F3 = 333m'})
    r=calculate_coverage(selected,'exp / NR.raw',chains())
    assert set(r.peptides.feature)=={'1:A1-C3 = 301m'}
    assert r.summary.set_index('protein').loc['A','covered_residues']==3
    assert r.summary.set_index('protein').loc['B','covered_residues']==0
    assert r.summary.set_index('protein').loc['A','sequence_coverage_pct']==60
    assert r.summary.set_index('protein').loc['B','sequence_coverage_pct']==0
    assert len(r.segments)==1
    assert set(r.segments.raw_file)=={'NR.raw'}
    assert r.peptides.ms_area.tolist()==[100.0]  # NOT 100+30 or unidentified+500


def test_deselect_all_is_zero_not_full_fallback():
    frame=make_selected()
    reduced=apply_peptide_exclusions(frame, {'1:A1-C3 = 301m','2:D1-F3 = 333m'})
    r=calculate_coverage(reduced,'exp / NR.raw',chains())
    assert len(r.peptides)==0
    assert r.segments.empty
    assert r.summary.covered_residues.tolist()==[0,0]
    assert r.summary.sequence_coverage_pct.tolist()==[0.,0.]


def test_per_raw_isolation_reduced_vs_nonreduced():
    frame=make_selected()
    # One additional selected peptide only in NR must NOT show up in reduced.
    extra=as_table([dict(experiment='exp',raw_file='NR.raw',preparation='Non-reduced',sample='S',
                         component_no=99,identified=True,identification='1:D4-E5 = 200m',
                         peptide_sequence='DE',ms_area=300,rt_min=10.0,mono_mass_da=200,id_type='Full')])
    frame=pd.concat([frame,extra],ignore_index=True)
    nr=calculate_coverage(frame,'exp / NR.raw',chains())
    reduced=calculate_coverage(frame,'exp / DTT.raw',chains())
    assert nr.summary.loc[nr.summary.protein=='A','covered_residues'].iloc[0]==5
    assert reduced.summary.loc[reduced.summary.protein=='A','covered_residues'].iloc[0]==3
    assert '1:D4-E5 = 200m' not in set(reduced.peptides.feature)
    assert set(reduced.segments.preparation)=={'Reduced'}


def test_component_level_selection_applies_before_coverage():
    frame=make_selected()
    keys=component_keys(frame)
    hidden={keys.loc[(frame.raw_file=='NR.raw')&(frame.component_no==1)].iloc[0]}
    chosen=selected_rows(frame,hidden)
    r=calculate_coverage(chosen,'exp / NR.raw',chains())
    assert r.summary.loc[r.summary.protein=='A','covered_residues'].iloc[0]==0
    assert r.summary.loc[r.summary.protein=='B','covered_residues'].iloc[0]==3


def test_disulfide_complex_spans_two_chains_but_one_identification():
    frame=as_table([dict(experiment='exp',raw_file='NR.raw',preparation='Non-reduced',
                         component_no=11,identified=True,
                         identification='1:A1-C3/2:D1-F3 = 600m[1ss]',peptide_sequence='ABC/DEF',
                         ms_area=200,rt_min=6.5,mono_mass_da=600,id_type='MS2')])
    r=calculate_coverage(frame,'exp / NR.raw',chains())
    assert len(r.peptides)==1
    assert len(r.segments)==2
    assert r.summary.covered_residues.tolist()==[3,3]
    assert r.summary.allocated_ms_area.tolist()==[100.0,100.0]
    assert r.summary.selected_ms_peak_area_pct.tolist()==[50.,50.]


def test_thresholds_and_zero_area_do_not_create_fake_coverage():
    frame=make_selected()
    r=calculate_coverage(frame,'exp / NR.raw',chains(),minimum_ms_area=50)
    assert r.summary.covered_residues.tolist()==[3,0]
    r2=calculate_coverage(frame,'exp / NR.raw',chains(),minimum_confidence_pct=98)
    assert r2.summary.covered_residues.tolist()==[3,0]
    r3=calculate_coverage(frame,'exp / NR.raw',chains(),maximum_peptide_mass_da=320)
    assert r3.summary.covered_residues.tolist()==[3,0]
    frame.loc[frame.component_no==1,'ms_area']=0
    r4=calculate_coverage(frame,'exp / NR.raw',chains())
    assert r4.summary.covered_residues.tolist()==[0,3]
    with pytest.raises(ValueError,match='non-negative'):
        calculate_coverage(frame,'exp / NR.raw',chains(),minimum_ms_area=-1)


def test_inference_reference_is_separate_from_coverage_evidence():
    frame=make_selected()
    references=reference_chains('exp',frame,pd.DataFrame([{'experiment':'exp','protein_chains':[]}]))
    assert [p.sequence for p in references]==['ABC','DEF']
    assert all(not c.verified for c in references)
    # Infer from both chains, but mapping excludes B based on selected IDs.
    r=calculate_coverage(apply_peptide_exclusions(frame,{'2:D1-F3 = 333m'}),'exp / NR.raw',references)
    assert r.summary.covered_residues.tolist()==[3,0]
    assert not r.summary.coverage_verified.any()
    assert r.warnings


def test_unknown_reference_residues_provisional_and_conflict_rejected():
    ref=[ProteinChain(1,'A','AB?DE',verified=False)]
    frame=as_table([dict(experiment='exp',raw_file='NR.raw',preparation='Non-reduced',
                         identified=True,component_no=1,identification='1:A1-B2 = 200m',
                         peptide_sequence='AB',ms_area=10,rt_min=3.0,id_type='MS2')])
    r=calculate_coverage(frame,'exp / NR.raw',ref)
    assert not bool(r.summary.coverage_verified.iloc[0])
    wrong=[ProteinChain(1,'A','ZZZZZ')]
    assert calculate_coverage(frame,'exp / NR.raw',wrong).segments.empty


def test_map_exports_have_expected_signatures_and_no_unselected_features():
    frame=make_selected()
    selected=apply_peptide_exclusions(frame,{'2:D1-F3 = 333m'})
    r=calculate_coverage(selected,'exp / NR.raw',chains())
    png=coverage_figure_bytes(r,'png',50)
    svg=coverage_figure_bytes(r,'svg',50)
    pdf=coverage_figure_bytes(r,'pdf',50)
    assert png[:8]==b'\x89PNG\r\n\x1a\n'
    assert b'<svg' in svg
    assert pdf[:4]==b'%PDF'
    tables=report_tables(r)
    assert set(['ProteinCoverage','MappedPeptides','MappedSegments','Residues','Parameters'])<=set(tables)
    assert set(tables['MappedSegments'].feature)=={'1:A1-C3 = 301m'}
    archive=mapping_package([r],residues_per_row=50)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        filenames=z.namelist()
        assert any(x.endswith('.pdf') for x in filenames)
        assert any(x.endswith('MappedPeptides.csv') for x in filenames)
        content=z.read(next(x for x in filenames if x.endswith('MappedPeptides.csv'))).decode('utf-8-sig')
        assert '2:D1-F3' not in content


def test_real_4197_csv_and_pmf_references():
    csv=Path('/mnt/data/4197_combined.csv')
    pmf=Path('/mnt/data/4197_combined_peptide_core_for_analysis.zip')
    if not csv.exists() or not pmf.exists():
        pytest.skip('Optional local BioPharma Finder experiment fixtures unavailable')
    data=load_files([csv])
    names=reference_chains('4197_combined',data.data,data.experiments)
    assert [(n.number,len(n.sequence)) for n in names]==[(1,21),(2,32)]
    # Selected ONLY one real feature -- must not return all-RAW 100% coverage.
    rsubset=data.data.loc[(data.data.raw_file=='4197p26_05.raw') &
                 (data.data.identification=='1:G1-E4 = 416.2271m')].copy()
    r=calculate_coverage(rsubset,'4197_combined / 4197p26_05.raw',names)
    assert r.summary.covered_residues.tolist()==[4,0]
    assert r.summary.sequence_coverage_pct.iloc[0]==pytest.approx(4/21*100)
    assert len(r.peptides)==1
    native=load_files([pmf])
    pnames=reference_chains('4197_combined',native.data,native.experiments)
    assert all(c.verified for c in pnames)
    assert [c.sequence for c in pnames]==[c.sequence for c in names]


def test_parse_regions_rejects_missing_annotations():
    assert parse_peptide_regions('1:A1-C3/2:D1-F3 = 600m[1ss]') == [(1,1,3,'A','C'),(2,1,3,'D','F')]
    assert not parse_peptide_regions('unidentified')


def test_no_applied_components_never_falls_back_to_raw():
    original=make_selected()
    empty=original.iloc[0:0].copy()
    report=calculate_coverage(empty,'exp / DTT.raw',chains())
    assert report.peptides.empty and report.segments.empty
    assert report.summary.covered_residues.tolist()==[0,0]
    assert report.preparation=='No selected components'
