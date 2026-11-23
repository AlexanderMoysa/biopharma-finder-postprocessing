"""No field loss and strict BF CSV / native-PMF reconciliation."""
from __future__ import annotations

import csv
import io
from pathlib import Path

import pandas as pd
import pytest

from bpf_postprocess.importers.csv_importer import read_csv
from bpf_postprocess.importers.loader import enrich_pmf_with_csv, load_files
from bpf_postprocess.processing.detailed_comparison import annotated_area_matrix, run_labeled_reference_matrix


def bf_export() -> bytes:
    """Two conditions, global Both vs per-RAW MS2/Full; dynamic BF headers."""
    columns = [
        "Level", "No.", "Identification", "Positions", "Peptide Sequence", "Mod", "Site",
        "∆ ppm", "Conf. Score", "Best ASR", "ID Type", "Integration Type",
        "RT", "M/Z", "Charge St.", "Mono Mass Exp.", "Avg Mass Exp.",
        "Mono Mass Theo.", "Missed Cleavages", "Ratio (07/05)", "Max Cond.",
        "Min Cond.", "Max MS Area", "Avg MS Area: 05", "%CV: 05",
        "MS Area:05.raw", "∆ ppm:05.raw", "ID Type:05.raw", "RT:05.raw", "M/Z:05.raw",
        "Mono Mass Exp.:05.raw", "Avg Mass Exp.:05.raw",
        "Avg MS Area: 07", "%CV: 07", "MS Area:with_DTT_07.raw",
        "∆ ppm:with_DTT_07.raw", "ID Type:with_DTT_07.raw", "RT:with_DTT_07.raw",
        "M/Z:with_DTT_07.raw", "Mono Mass Exp.:with_DTT_07.raw",
        "Avg Mass Exp.:with_DTT_07.raw", "Protein", "Comment",
    ]
    fields = {
        "Level":"Component", "No.":"7", "Identification":"1:G1-E4 = 416.2271m",
        "Positions":"1-4", "Peptide Sequence":"GIVE", "Mod":"None", "Site":"G1",
        "∆ ppm":"-1.5", "Conf. Score":"99.9", "Best ASR":"1.2", "ID Type":"Both",
        "Integration Type":"Automatic", "RT":"4.32", "M/Z":"417.233",
        "Charge St.":"1", "Mono Mass Exp.":"416.2269", "Avg Mass Exp.":"416.23",
        "Mono Mass Theo.":"416.2271", "Missed Cleavages":"2",
        "Ratio (07/05)":"2.0", "Max Cond.":"07", "Min Cond.":"05",
        "Max MS Area":"200", "Avg MS Area: 05":"100", "%CV: 05":"0",
        "MS Area:05.raw":"100", "∆ ppm:05.raw":"-1.5", "ID Type:05.raw":"MS2",
        "RT:05.raw":"4.32", "M/Z:05.raw":"417.233", "Mono Mass Exp.:05.raw":"416.2269",
        "Avg Mass Exp.:05.raw":"416.23", "Avg MS Area: 07":"200", "%CV: 07":"0",
        "MS Area:with_DTT_07.raw":"200", "∆ ppm:with_DTT_07.raw":"-0.8",
        "ID Type:with_DTT_07.raw":"Full", "RT:with_DTT_07.raw":"4.31",
        "M/Z:with_DTT_07.raw":"417.232", "Mono Mass Exp.:with_DTT_07.raw":"416.2270",
        "Avg Mass Exp.:with_DTT_07.raw":"416.24", "Protein":"A", "Comment":"verified",
    }
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=columns)
    writer.writeheader()
    writer.writerow(fields)
    return out.getvalue().encode('utf-8-sig')


def test_all_source_fields_and_condition_mapping_survive():
    data, meta, warnings = read_csv(bf_export(), name="test_combined(1).csv")
    assert meta['experiment']=='test_combined'
    original = meta['original_csv_rows']
    assert original.shape==(1,43)
    assert original.loc[0,'Ratio (07/05)']=='2.0'
    assert original.loc[0,'Comment']=='verified'
    assert original.loc[0,'ID Type']=='Both'
    assert len(data)==2
    d=data.set_index('raw_file')
    assert d.loc['05.raw','id_type']=='MS2'
    assert d.loc['with_DTT_07.raw','id_type']=='Full'
    assert d.id_type_overall.eq('Both').all()
    assert d.loc['05.raw','delta_ppm']==-1.5
    assert d.loc['with_DTT_07.raw','delta_ppm']==-0.8
    assert d.loc['05.raw','bf_avg_ms_area_condition']==100
    assert d.loc['with_DTT_07.raw','bf_avg_ms_area_condition']==200
    assert d.bf_cv_pct_condition.tolist()==[0,0]
    assert d.loc['with_DTT_07.raw','avg_mass_da']==pytest.approx(416.24)
    assert d.loc['05.raw','best_asr']==pytest.approx(1.2)
    assert d.loc['05.raw','mono_mass_theo_da']==pytest.approx(416.2271)
    assert d.loc['05.raw','missed_cleavages']==2
    assert d.loc['05.raw','bf_ratio']==pytest.approx(2.0)
    assert d.loc['05.raw','bf_ratio_label']=='Ratio (07/05)'
    assert d.comment.eq('verified').all()
    assert d.loc['05.raw','condition']=='05'
    assert d.loc['with_DTT_07.raw','condition']=='07'
    assert d.bf_csv_source_file.eq('test_combined(1).csv').all()
    assert not warnings


def test_verified_merge_keeps_precise_pmf_and_csv_metadata():
    csv_data, _, _ = read_csv(bf_export(), name="test_combined.csv")
    pmf=csv_data.copy()
    pmf['ms_area']=[100.125,200.125]
    pmf['rt_min']=[4.32001,4.31001]
    pmf['id_type']='Unavailable (not decoded)'
    pmf['id_type_overall']=''
    pmf['site']=''
    pmf['delta_ppm']=float('nan')
    pmf['source_type']='resultsbpf'
    pmf['source_file']='test.resultsbpf'
    enriched=enrich_pmf_with_csv(pmf,csv_data)
    assert len(enriched)==2
    assert enriched.ms_area.tolist()==[100.125,200.125]
    assert enriched.id_type.tolist()==['MS2','Full']
    assert enriched.delta_ppm.tolist()==[-1.5,-0.8]
    assert enriched.site.eq('G1').all()
    assert enriched.source_type.eq('resultsbpf+csv').all()
    assert enriched.source_file.eq('test.resultsbpf').all()
    assert enriched.bf_csv_source_file.eq('test_combined.csv').all()
    assert enriched.annotation_origin.eq('bf_csv_verified_against_pmf').all()


def test_wrong_experiment_or_reprocessed_csv_rejected():
    csv_data, _, _=read_csv(bf_export(), name='test_combined.csv')
    pmf=csv_data.copy()
    changed=csv_data.copy()
    changed.loc[0,'ms_area']=12345
    with pytest.raises(ValueError,match='ms_area mismatch'):
        enrich_pmf_with_csv(pmf,changed)
    changed=csv_data.copy()
    changed.loc[0,'identification']='Wrong peptide'
    with pytest.raises(ValueError,match='identification mismatch'):
        enrich_pmf_with_csv(pmf,changed)
    changed=csv_data.iloc[:1].copy()
    with pytest.raises(ValueError,match='records do not align'):
        enrich_pmf_with_csv(pmf,changed)


def test_extra_fields_are_run_labeled_in_first_peptide_table():
    data, _, _=read_csv(bf_export(),name='test_combined.csv')
    ref='test_combined / 05.raw'; cmp='test_combined / with_DTT_07.raw'
    # Here we explicitly test run-labeled extra fields with same preparation;
    # compare_runs otherwise intentionally rejects Reduced vs Non-reduced.
    data['preparation']='Unspecified'
    annotated=annotated_area_matrix(data,ref,[cmp])
    assert annotated.loc[0,f'Delta ppm dominant | {ref}']==-1.5
    assert annotated.loc[0,f'Delta ppm dominant | {cmp}']==-0.8
    assert annotated.loc[0,f'Best ASR dominant | {ref}']==1.2
    assert annotated.loc[0,f'Avg mass exp. dominant (Da) | {cmp}']==pytest.approx(416.24)
    wide=run_labeled_reference_matrix(data,ref,[cmp])
    assert wide.loc[0,f'Delta ppm dominant | {ref} [REFERENCE]']==-1.5
    assert wide.loc[0,f'Mono mass theo. dominant (Da) | {cmp}']==pytest.approx(416.2271)
    assert wide.loc[0,f'Missed cleavages dominant | {cmp}']==2


def test_real_file_pair_enriches_exactly_without_missing_components():
    csv_4197=Path('/mnt/data/4197_combined.csv')
    pmf_4197=Path('/mnt/data/4197_combined_peptide_core_for_analysis.zip')
    csv_4204=Path('/mnt/data/4204_combined.csv')
    if not (csv_4197.exists() and pmf_4197.exists() and csv_4204.exists()):
        pytest.skip('Optional original BioPharma Finder files not available')
    dataset=load_files([csv_4204,pmf_4197,csv_4197])
    assert len(dataset.data)==3322
    assert len(dataset.runs)==4
    assert len(dataset.experiments)==2
    assert {k:len(v.columns) for k,v in dataset.source_csv_tables.items()}=={
        '4197_combined':43,'4204_combined':43,
    }
    assert not dataset.data.duplicated(['experiment','raw_file','component_no']).any()
    assert dataset.data['id_type'].eq('Unavailable (not decoded)').sum()==0
    assert set(dataset.data['id_type'])=={'Full','MS2','NoID'}
    assert dataset.data['delta_ppm'].notna().sum()==442
    assert dataset.data['site'].astype(str).ne('').sum()==144
    assert dataset.data['best_asr'].notna().sum()==442
    assert dataset.data['source_type'].eq('resultsbpf+csv').sum()==1694


def test_blank_run_specific_fields_are_not_invented_from_global_values():
    text = (
        'No.,Identification,ID Type,∆ ppm,ID Type:ref.raw,ID Type:other.raw,'
        '∆ ppm:ref.raw,∆ ppm:other.raw,MS Area:ref.raw,MS Area:other.raw\n'
        '1,1:G1-E4 = 416m,Both,12.3,MS2,,,-1.2,100,200\n'
    ).encode('utf-8')
    data, _, _ = read_csv(text, name='exp.csv')
    ref = data.loc[data.raw_file == 'ref.raw'].iloc[0]
    other = data.loc[data.raw_file == 'other.raw'].iloc[0]
    assert ref.id_type == 'MS2'
    assert pd.isna(ref.delta_ppm)  # Empty RAW Δppm is NOT global 12.3.
    assert other.id_type == ''     # Empty RAW ID Type is NOT global 'Both'.
    assert other.delta_ppm == -1.2
    assert ref.delta_ppm_overall == 12.3
    assert other.id_type_overall == 'Both'


def test_source_columns_retained_after_file_loading(tmp_path):
    original = bf_export()
    source = tmp_path/'test_combined.csv'
    source.write_bytes(original)
    dataset = load_files([source])
    assert dataset.data.shape[0] == 2
    assert set(dataset.source_csv_tables) == {'test_combined'}
    assert dataset.source_csv_tables['test_combined'].columns.tolist() == list(
        csv.DictReader(io.StringIO(original.decode('utf-8-sig'))).fieldnames
    )
    assert dataset.data.delta_ppm.notna().sum() == 2
    assert dataset.data.site.eq('G1').all()
