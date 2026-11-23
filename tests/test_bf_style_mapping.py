"""Independent BF-inspired compact rendering and selected-only invariants."""
import math

import pandas as pd
import pytest

from bpf_postprocess.model import as_table
from bpf_postprocess.processing.coverage import ProteinChain, calculate_coverage
from bpf_postprocess.processing.coverage_display import compact_coverage_segments
from bpf_postprocess.processing.peptide_selection import apply_peptide_exclusions
from bpf_postprocess.reports.coverage_map import (
    _compact_lanes, coverage_figure_bytes, report_tables,
)


REFERENCE = [ProteinChain(1, 'A', 'ABCDE'), ProteinChain(2, 'B', 'DEFGH')]
RUN = 'exp / dtt.raw'


def _frame():
    rows = [
        # Two modified forms of SAME region, only one displayed bar.
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=True,
             component_no=1, identification='1:A1-C3 = 300m',
             peptide_sequence='ABC', modification='None', ms_area=300, rt_min=4.3),
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=True,
             component_no=2, identification='1:A1-C3 = 300m(Na+)',
             peptide_sequence='ABC', modification='Na+', ms_area=20, rt_min=4.32),
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=True,
             component_no=3, identification='1:C3-E5 = 400m',
             peptide_sequence='CDE', modification='None', ms_area=120, rt_min=8.4),
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=True,
             component_no=4, identification='2:D1-F3 = 330m',
             peptide_sequence='DEF', modification='None', ms_area=900, rt_min=5.6),
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=True,
             component_no=5, identification='2:F3-H5 = 500m',
             peptide_sequence='FGH', modification='None', ms_area=45, rt_min=3.9),
        # Another RAW is present, must never leak into dtt map.
        dict(experiment='exp', raw_file='nr.raw', preparation='Non-reduced', identified=True,
             component_no=6, identification='1:A1-E5 = 600m',
             peptide_sequence='ABCDE', modification='None', ms_area=999999, rt_min=9.0),
        dict(experiment='exp', raw_file='dtt.raw', preparation='Reduced', identified=False,
             component_no=7, ms_area=1e9),
    ]
    return as_table(rows)


def test_compact_has_one_bar_per_peptide_region_and_includes_all_selected_areas():
    report = calculate_coverage(_frame(), RUN, REFERENCE)
    compact = compact_coverage_segments(report.segments)
    assert len(report.segments) == 5
    assert len(compact) == 4
    a = compact.loc[(compact.chain_number == 1) & (compact.start == 1)].iloc[0]
    assert a.regional_ms_area == 320
    assert a.component_count == 2
    assert a.identification_count == 2
    assert a.rt_min == pytest.approx(4.3)  # strongest annotation, not area-weighted mean
    assert a.representative_feature == '1:A1-C3 = 300m'
    assert set(compact.raw_file) == {'dtt.raw'}
    assert a.relative_signal_pct == pytest.approx(100)
    assert compact[(compact.chain_number == 2) & (compact.start == 1)].iloc[0].relative_signal_pct == 100
    assert compact[(compact.chain_number == 2) & (compact.start == 3)].iloc[0].relative_signal_pct == 5
    assert report.summary.number_of_ms_peaks.tolist() == [3, 2]


def test_no_unselected_peptide_returns_to_compact_visual_or_export():
    rows = apply_peptide_exclusions(_frame(), {'1:A1-C3 = 300m', '1:A1-C3 = 300m(Na+)',
                                                 '2:F3-H5 = 500m'})
    report = calculate_coverage(rows, RUN, REFERENCE)
    compact = compact_coverage_segments(report.segments)
    assert compact[['chain_number','start','end']].values.tolist() == [[1,3,5], [2,1,3]]
    assert report.summary.number_of_ms_peaks.tolist() == [1,1]
    assert report.summary.covered_residues.tolist() == [3,3]
    assert set(report_tables(report)['BFStyleRegions'].representative_feature) == {
        '1:C3-E5 = 400m', '2:D1-F3 = 330m'
    }


def test_deselect_all_gives_empty_compact_not_original_coverage():
    rows=_frame().iloc[0:0]
    report=calculate_coverage(rows,RUN,REFERENCE)
    assert compact_coverage_segments(report.segments).empty
    assert report.summary.covered_residues.tolist() == [0,0]
    assert coverage_figure_bytes(report,'png',view_mode='compact').startswith(b'\x89PNG')


def test_both_views_export_and_invalid_view_rejected():
    report=calculate_coverage(_frame(),RUN,REFERENCE)
    for view in ['compact','detailed']:
        assert coverage_figure_bytes(report,'pdf',view_mode=view).startswith(b'%PDF')
        assert b'<svg' in coverage_figure_bytes(report,'svg',view_mode=view)
    with pytest.raises(ValueError,match='view_mode'):
        coverage_figure_bytes(report,'png',view_mode='incorrect')


def test_lane_packing_never_overlaps_and_uses_longest_first():
    bars=[dict(start=1,end=21,regional_ms_area=10),
          dict(start=1,end=4,regional_ms_area=5),
          dict(start=5,end=17,regional_ms_area=8),
          dict(start=18,end=21,regional_ms_area=4)]
    packed=_compact_lanes(bars)
    assert packed[0][0]['end']==21 and packed[0][1]==0
    for i,(r,lane) in enumerate(packed):
        for other,other_lane in packed[i+1:]:
            assert lane!=other_lane or r['end'] < other['start'] or other['end']<r['start']


def test_mixed_raw_compaction_cannot_forge_combined_map():
    r1=calculate_coverage(_frame(),RUN,REFERENCE)
    r2=calculate_coverage(_frame(),'exp / nr.raw',REFERENCE)
    with pytest.raises(ValueError,match='exactly one RAW'):
        compact_coverage_segments(pd.concat([r1.segments,r2.segments],ignore_index=True))
