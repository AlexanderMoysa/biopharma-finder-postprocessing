"""ID Type fidelity, selection and display regression tests."""
from __future__ import annotations

import math
import pandas as pd
import pytest

from bpf_postprocess.model import as_table
from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.processing.filters import filter_components
from bpf_postprocess.processing.id_types import (
    UNAVAILABLE_ID_TYPE, combined_id_types, displayed_id_types, id_type_summary
)
from bpf_postprocess.processing.detailed_comparison import (
    annotated_area_matrix, detailed_reference_comparison, run_labeled_reference_matrix
)
from bpf_postprocess.processing.modifications import modification_areas
from bpf_postprocess.processing.unidentified import unidentified_summary
from bpf_postprocess.processing.peptide_selection import peptide_inventory
from bpf_postprocess.processing.relative_modifications import modification_fractions


def test_run_specific_ids_override_overall_both():
    source = (
        "No.,Identification,Peptide Sequence,Mod,ID Type,ID Type:ref.raw,ID Type:treated.raw,"
        "MS Area:ref.raw,MS Area:treated.raw,RT:ref.raw,RT:treated.raw,Conf. Score\n"
        "1,1:G1-E4 = 416.2m,GIVE,None,Both,MS2,Full,120,150,4.31,4.32,99.4\n"
        "2,,,,NoID,NoID,NoID,55,66,4.22,4.21,\n"
    ).encode()
    rows, _, _ = read_csv(source, name='experiment.csv')
    assert rows.shape[0] == 4
    p = rows[rows.component_no == 1].set_index('raw_file')
    assert p.loc['ref.raw','id_type'] == 'MS2'
    assert p.loc['treated.raw','id_type'] == 'Full'
    assert set(p.id_type_overall) == {'Both'}
    assert rows[rows.component_no == 2].id_type.eq('NoID').all()
    assert set(displayed_id_types(rows.id_type)) == {'MS2', 'Full', 'NoID'}
    assert filter_components(rows,id_types=['Full']).shape[0] == 1
    assert filter_components(rows,id_types=['NoID']).shape[0] == 2
    assert filter_components(rows,id_types=['MS2','Full','NoID']).shape[0] == len(rows)
    assert filter_components(rows,id_types=[]).empty
    assert len(filter_components(rows,id_types=['MS2'],min_confidence_pct=90,min_area=100))==1


def test_unavailable_remains_distinct_from_noid():
    df=as_table([
        dict(experiment='E',raw_file='pmf.raw',component_no=1,ms_area=500,
             identified=True,identification='P',id_type='',sample='S',preparation='Non-reduced'),
        dict(experiment='E',raw_file='csv.raw',component_no=2,ms_area=200,
             identified=False,id_type='NoID',sample='S',preparation='Non-reduced'),
    ])
    assert displayed_id_types(df.id_type).tolist()==[UNAVAILABLE_ID_TYPE,'NoID']
    assert len(filter_components(df,id_types=['NoID']))==1
    assert len(filter_components(df,id_types=[UNAVAILABLE_ID_TYPE]))==1
    assert combined_id_types(df.id_type) == f'NoID; {UNAVAILABLE_ID_TYPE}'
    summary=id_type_summary(df)
    assert set(summary.id_type)=={'NoID',UNAVAILABLE_ID_TYPE}
    assert summary.components.sum()==2
    assert summary.ms_area_sum.sum()==700


def _reference_dataframe():
    r=[]
    for raw,prep,rows in [
        ('reference.raw','Non-reduced',[(1,'Full',120.,'None'),(2,'MS2',30.,'None'),(3,'NoID',100.,'None')]),
        ('comparison.raw','Non-reduced',[(1,'MS2',90.,'None'),(2,'Full',80.,'Oxidation'),(3,'NoID',80.,'None')]),
        ('with_DTT.raw','Reduced',[(1,'MS2',500.,'None')]),
    ]:
        for no,typ,area,mod in rows:
            r.append(dict(experiment='E',sample='S',sample_group='Control',raw_file=raw,
                          preparation=prep,component_no=no,identified=typ!='NoID',
                          identification='1:G1-E4 = 416.2m' if no==1 else ('1:G1-E4 = 416.2m(Ox)' if no==2 else ''),
                          peptide_sequence='GIVE' if no!=3 else '',modification=mod,
                          id_type=typ,ms_area=area,rt_min=4.2,confidence_pct=99.))
    return as_table(r)


def test_id_types_in_peptide_comparison_and_report_tables():
    df=_reference_dataframe()
    ref='E / reference.raw';other='E / comparison.raw'
    annotated=annotated_area_matrix(df,ref,[other])
    assert f'ID Type | {ref}' in annotated
    assert annotated.loc[annotated.feature=='1:G1-E4 = 416.2m',f'ID Type | {ref}'].iloc[0]=='Full'
    labeled=run_labeled_reference_matrix(df,ref,[other])
    assert f'ID Type | {ref} [REFERENCE]' in labeled
    assert f'ID Type | {other}' in labeled
    assert 'DTT' not in str(labeled.columns.tolist())
    long=detailed_reference_comparison(df,ref,[other])
    assert 'reference_id_type' in long and 'comparison_id_type' in long
    assert 'reference_id_types' in long and 'comparison_id_types' in long
    inv=peptide_inventory(df)
    assert 'id_types' in inv
    assert 'Full; MS2' in set(inv['id_types'])
    mods=modification_areas(df)
    assert 'id_types' in mods
    unknown=unidentified_summary(df)
    assert 'id_types' in unknown
    assert set(unknown['id_types'])=={'NoID'}


def test_relative_modification_reports_id_types_of_both_species():
    data=as_table([
        dict(experiment='E',sample='S',sample_group='C',raw_file='a.raw',preparation='Reduced',
             component_no=1,identified=True,identification='1:G1-E4 = 416m',
             peptide_sequence='GIVE',modification='None',id_type='Full',ms_area=80.),
        dict(experiment='E',sample='S',sample_group='C',raw_file='a.raw',preparation='Reduced',
             component_no=2,identified=True,identification='1:G1-E4 = 416m(M1+Oxidation)',
             peptide_sequence='GIVE',modification='Oxidation',id_type='MS2',ms_area=20.),
    ])
    out=modification_fractions(data,['Oxidation'])
    assert len(out)==1
    assert out.iloc[0].modified_id_types=='MS2'
    assert out.iloc[0].unmodified_id_types=='Full'
    assert out.iloc[0].modified_ms_area_pct==20


def test_empty_type_table_and_empty_filter_safe():
    empty=as_table([])
    assert len(id_type_summary(empty))==0
    assert len(filter_components(empty,id_types=['MS2']))==0
