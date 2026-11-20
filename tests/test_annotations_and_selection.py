import pandas as pd
import pytest
from bpf_postprocess.importers.annotation_parser import parse_identification
from bpf_postprocess.processing.filters import filter_components
from bpf_postprocess.processing.selection import component_keys, selected_rows, update_exclusions
from bpf_postprocess.model import as_table

CHAINS = [
    {'number':1,'name':'A','sequence':'GIVEQCCTSICSLYQLENYCG'},
    {'number':2,'name':'B','sequence':'FVNQHLCGSHLVEALYLVCGERGFFYTPKTRR'},
]

def test_pmf_annotations_recover_sequences_and_mods():
    d=parse_identification('1:N18-G21/2:F1-E13 = 1934.8404m[1ss]',CHAINS)
    assert d['peptide_sequence']=='NYCG/FVNQHLCGSHLVE'
    assert d['positions']=='18-21 | 1-13'
    assert d['modification']=='1ss'
    assert d['protein']=='A'
    assert d['site']==''  # disulfide pairing cannot be assumed
    e=parse_identification('1:N18-G21 = 455.1475m(N18+Deamidation)',CHAINS)
    assert e['modification']=='Deamidation'
    assert e['site']=='N18'
    f=parse_identification('2:F1-E13 = 1481.7085m(~Q4+Deamidation)',CHAINS)
    assert f['site']=='~Q4'
    g=parse_identification('2:F1-E13 = 1481.7085m(Fe3+)',CHAINS)
    assert g['modification']=='Fe3+'


def test_invalid_peptide_is_not_guessed():
    assert not parse_identification('1:X18-G21 = 455m',CHAINS)['peptide_sequence']
    assert not parse_identification('1:N18-G29 = 455m',CHAINS)['peptide_sequence']


def test_confidence_filter_and_persistent_selection():
    rows=[]
    for experiment, raw, comp, conf in [('A','A.raw',1,99.5),('A','A.raw',2,60.0),('B','B.raw',1,90.0),('B','B.raw',2,None)]:
        rows.append({'experiment':experiment,'sample':experiment,'raw_file':raw,
                     'component_no':comp,'identified':conf is not None,'confidence_pct':conf,'ms_area':100,
                     'identification':f'peptide {comp}' if conf is not None else ''})
    data=as_table(rows)
    assert len(filter_components(data,min_confidence_pct=95))==1
    assert len(filter_components(data,min_confidence_pct=0))==4
    excluded=update_exclusions(data,pd.Series([True,False,True,True]),set())
    assert len(selected_rows(data,excluded))==3
    subset=data[data['experiment']=='B']
    assert len(selected_rows(subset,excluded))==2
    # The same component number 1 in another experiment is not excluded.
    assert component_keys(data).nunique()==4
    again=update_exclusions(subset,pd.Series([False,True]),excluded)
    assert len(selected_rows(data,again))==2
