import math
import pandas as pd
from bpf_postprocess.model import as_table
from bpf_postprocess.processing.multisample import reference_comparison, replicate_statistics
from bpf_postprocess.processing.relative_modifications import modification_fractions


def sample_data():
    rows=[]
    for sample, run, amount, mod in [
        ('Control','C1',100.,0.),('Control','C2',120.,10.),
        ('Stress','S1',80.,20.),('Other','O1',150.,0.)
    ]:
        for n, (annotation, area) in enumerate([('None',amount),('Deamidation',mod)],1):
            rows.append(dict(experiment=sample,sample=sample,sample_group='grp_'+sample,
                             raw_file=run,preparation='Reduced',component_no=n,
                             identified=True, identification='1:N18-G21 = 455m'+(('(N18+Deamidation)') if n==2 else ''),
                             peptide_sequence='NYCG',modification=annotation,ms_area=area))
    return as_table(rows)


def test_reference_multi_run():
    data=sample_data()
    out=reference_comparison(data,'Control / C1')
    assert set(out.comparison_run)=={'Control / C2','Stress / S1','Other / O1'}
    stress=out[(out.comparison_run=='Stress / S1') & (out.feature=='1:N18-G21 = 455m')].iloc[0]
    assert stress.fold_vs_reference==0.8
    assert math.isclose(stress.percent_change,-20)


def test_replicate_sd_rsd():
    data=sample_data()
    stats=replicate_statistics(data)
    control=stats[(stats['sample']=='Control') & (stats['feature']=='1:N18-G21 = 455m')].iloc[0]
    assert control.n_runs==2
    assert control.mean_ms_area==110
    assert math.isclose(control.sd_ms_area, 14.1421356237)
    single=stats[stats['sample']=='Stress'].iloc[0]
    assert math.isnan(single.sd_ms_area)


def test_relative_modification_with_and_without_reference():
    df=sample_data()
    f=modification_fractions(df,['Deamidation'])
    stress=f[f.raw_file=='S1'].iloc[0]
    assert stress.reference_present
    assert stress.modified_ms_area_pct==20
    assert f[f.raw_file=='C1'].iloc[0].modified_ms_area_pct==0
    # Missing plain form is NOT imputed with zero and NOT treated as 100%.
    missing=df[(df.raw_file=='S1') & (df.modification=='Deamidation')]
    f2=modification_fractions(missing,['Deamidation'])
    assert len(f2)==1
    assert not f2.iloc[0].reference_present
    assert math.isnan(f2.iloc[0].modified_ms_area_pct)


def test_zero_area_denominator():
    data=sample_data()
    zero=data[data.raw_file=='C1'].copy()
    zero['ms_area']=0
    f=modification_fractions(zero,['Deamidation'])
    assert pd.isna(f.iloc[0]['modified_ms_area_pct'])
