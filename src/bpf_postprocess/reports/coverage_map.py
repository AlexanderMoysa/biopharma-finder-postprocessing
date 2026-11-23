"""Static editable-free sequence-coverage figure and reproducible report exports.

The visual structure is inspired by Thermo BioPharma Finder's Coverage Map,
but its colors are explicitly based on relative MS-area signal and are NOT
Thermo's proprietary peptide-recovery score or abundance (mol).
"""
from __future__ import annotations

import io
import math
import re
import zipfile

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from bpf_postprocess.processing.coverage import CoverageResult
from bpf_postprocess.processing.coverage_display import compact_coverage_segments
from bpf_postprocess.reports.exports import csv_bytes, excel_bytes

# The requested BF-like red/yellow/green/blue intensity steps.
BINS = [
    (50, '#ff7d80'), (20, '#ffbb65'), (10, '#ffe77a'),
    (5, '#bcdd8a'), (2, '#7bdddd'), (1, '#a3ccff'),
    (0.5, '#c4b4ff'), (0.2, '#d5d5d5'), (0.1, '#efefef'), (0, '#ffffff'),
]


def signal_color(percent: float) -> str:
    for threshold, color in BINS:
        if percent >= threshold:
            return color
    return BINS[-1][1]


def _lanes(overlaps: list[dict]) -> list[tuple[dict,int]]:
    """Greedy interval packing: avoid bars covering each other."""
    occupied: list[int] = []
    result=[]
    for record in sorted(overlaps,key=lambda x:(x['start'], x['end'], -x['ms_area'])):
        start,end=int(record['start']),int(record['end'])
        index=next((i for i,last in enumerate(occupied) if last < start),len(occupied))
        if index == len(occupied): occupied.append(end)
        else: occupied[index]=end
        result.append((record,index))
    return result



def _compact_lanes(bars: list[dict]) -> list[tuple[dict,int]]:
    """Put the longest BF-style fragments on top and pack disjoint spans."""
    occupied: list[list[tuple[int,int]]] = []
    result=[]
    for item in sorted(bars, key=lambda r: (
        -(int(r['end'])-int(r['start'])+1), -float(r['regional_ms_area']),
        int(r['start']), int(r['end']))):
        begin,end=int(item['start']),int(item['end'])
        lane=next((i for i, spans in enumerate(occupied)
                   if all(end < low or begin > high for low,high in spans)), len(occupied))
        if lane==len(occupied):
            occupied.append([])
        occupied[lane].append((begin,end))
        result.append((item,lane))
    return result


def build_coverage_figure(report: CoverageResult, *, residues_per_row: int=50,
                          view_mode: str='compact'):
    """BF-style region-compact or fully detailed selected-only coverage figure.

    Neither view attempts to reproduce the proprietary BF Recovery algorithm.
    The compact view normalizes regional MS Areas independently within each chain.
    """
    if view_mode not in ('compact','detailed'):
        raise ValueError('view_mode must be compact or detailed')

    """Return a Matplotlib figure with summary, sequence letters, RT-labeled bars."""
    if not (10 <= residues_per_row <= 200):
        raise ValueError('Residues per row must be in [10, 200]')
    visual_segments=(compact_coverage_segments(report.segments) if view_mode=='compact'
                     else report.segments)
    chunks=[]
    for chain in report.chains:
        for left in range(1,len(chain.sequence)+1,residues_per_row):
            right=min(len(chain.sequence),left+residues_per_row-1)
            sub=[]
            for entry in visual_segments.loc[visual_segments['chain_number']==chain.number].to_dict('records'):
                if int(entry['start'])<=right and int(entry['end'])>=left:
                    sub.append({**entry,'start':max(left,int(entry['start']))-left+1,
                                'end':min(right,int(entry['end']))-left+1})
            bands=(_compact_lanes(sub) if view_mode=='compact' else _lanes(sub))
            count=max((i+1 for _,i in bands),default=0)
            chunks.append((chain,left,right,bands,count))
    height=sum(2.05+max(1,n)*0.38 for _,_,_,_,n in chunks)+4.6+0.38*len(report.chains)
    width=max(11, min(25,residues_per_row*0.29))
    fig,ax=plt.subplots(figsize=(width,height),facecolor='white')
    ax.set_xlim(-2.7,residues_per_row+1.2)
    ax.set_ylim(0,height)
    ax.axis('off')
    y=height-0.6
    ax.text(-2.4,y,'Sequence Coverage Map',fontsize=17,fontweight='bold',ha='left')
    y-=0.55
    ax.text(-2.4,y,f'{report.raw_file}  |  {report.experiment}  |  {report.preparation}',
            fontsize=10,ha='left')
    y-=0.43
    subtitle=('BF-style compact  |  Selected peptide regions only  |  Color: per-chain regional MS Area proxy'
              if view_mode=='compact' else
              'Detailed  |  All selected peptide identifications  |  Color: % of strongest selected peptide MS Area')
    ax.text(-2.4,y,subtitle,fontsize=9,ha='left',color='#454545')
    y-=0.43
    if not report.summary.empty:
        records=report.summary.to_dict('records')
        headers=['Proteins','Selected MS peaks*','Selected MS area (%)*','Sequence coverage','Abundance (mol)']
        cells=[]
        for s in records:
            inferred='~' if not s['coverage_verified'] else ''
            cells.append([f"{s['chain_number']}:{s['protein']}",str(s['number_of_ms_peaks']),
                          f"{s['selected_ms_peak_area_pct']:.1f}%",
                          f"{inferred}{s['sequence_coverage_pct']:.1f}% "
                          f"({s['covered_residues']}/{s['sequence_length']})",'N/A'])
        table_h=0.40*(len(cells)+1)
        table=ax.table(cellText=cells,colLabels=headers,cellLoc='left',colLoc='left',
                       colWidths=[0.12,0.19,0.23,0.26,0.20],
                       bbox=[0.005,max(.0,(y-table_h)/height),0.98,table_h/height])
        table.auto_set_font_size(False)
        table.set_fontsize(8.3)
        for (row,col),cell in table.get_celld().items():
            cell.set_linewidth(0.6)
            cell.set_edgecolor('#777777')
            if row==0:
                cell.set_facecolor('#ecf2f4')
                cell.get_text().set_weight('bold')
        y-=table_h+0.20
    ax.text(-2.4,y,'* Counts/areas include all selected components, not only displayed bars. '
            'Colors are approximate; molar abundance not inferred.',
            fontsize=7.8,ha='left',color='#555555')
    y-=0.48
    for i,(threshold,color) in enumerate(BINS):
        span=residues_per_row/len(BINS)
        x=-1.8+i*span
        ax.add_patch(Rectangle((x,y-0.12),span*0.92,0.28,
                               linewidth=.6,edgecolor='#404040',facecolor=color))
        label=f'>={threshold}%' if threshold else '>0.0%'
        ax.text(x+span*0.46,y+0.02,label,ha='center',va='center',fontsize=7.5,color='#282828')
    y-=0.75
    for chain,left,right,bands,count in chunks:
        ax.text(-2.4,y,f'Chain {chain.number}: {chain.name}  ({left}-{right})',fontsize=10,fontweight='bold',ha='left')
        y-=0.38
        for i in range(left,right+1):
            x=i-left+1
            ax.text(x,y,str(i),ha='center',fontsize=5.7,fontfamily='DejaVu Sans Mono',color='#626262')
        y-=0.25
        for i in range(left,right+1):
            x=i-left+1
            aa=chain.sequence[i-1]
            ax.text(x,y,aa,ha='center',fontsize=11.5,fontweight='bold',fontfamily='DejaVu Sans Mono')
        y-=0.47
        for entry,lane in bands:
            seg_y=y-lane*0.37
            bar_left=entry['start']-.46
            bar_w=entry['end']-entry['start']+.92
            ax.add_patch(Rectangle((bar_left,seg_y-.18),bar_w,.31,linewidth=.7,
                                   edgecolor='#282828',facecolor=signal_color(float(entry['relative_signal_pct']))))
            rt=entry['rt_min']
            if rt is not None and math.isfinite(float(rt)) and bar_w >= 2.0:
                ax.text((entry['start']+entry['end'])/2,seg_y-.01,f'{float(rt):.1f}',
                        ha='center',va='center',fontsize=8,color='#262626')
        y-=max(1,count)*.38+0.48
    y-=0.24
    note=('BF-style proxy: selected regional MS Area / strongest selected interval in the SAME chain (%)'
          if view_mode=='compact' else
          'Detailed: relative selected-peptide MS Area / strongest selected identification in this RAW (%)')
    ax.text(-2.4,max(.28,y),note+'  |  Not Thermo Peptide Recovery',fontsize=8.5,ha='left')
    fig.tight_layout(pad=1.0)
    return fig


def coverage_figure_bytes(report: CoverageResult, ext: str='png', residues_per_row: int=50,
                          view_mode: str='compact') -> bytes:
    if ext not in ('png','svg','pdf'):
        raise ValueError('Supported exports: png, svg, pdf')
    fig=build_coverage_figure(report,residues_per_row=residues_per_row,view_mode=view_mode)
    stream=io.BytesIO()
    fig.savefig(stream,format=ext,dpi=190,facecolor='white',bbox_inches='tight')
    plt.close(fig)
    return stream.getvalue()


def report_tables(report: CoverageResult) -> dict[str,object]:
    import pandas as pd
    params = {'raw_run': report.run, **report.parameters, 'coverage_rule':'Only already-selected peptide components in this RAW',
             'area_sharing_rule':'For cross-chain complexes, MS Area divided equally among distinct chains',
             'signal_color_rule':'Feature MS Area / strongest selected feature in this RAW',
             'mol_abundance':'Not estimated; BF formula not reproduced'}
    return {'ProteinCoverage':report.summary,'MappedPeptides':report.peptides,
            'MappedSegments':report.segments,
            'BFStyleRegions':compact_coverage_segments(report.segments),
            'Residues':report.residue_table,
            'Parameters':pd.DataFrame([params]), 'Warnings':pd.DataFrame({'warning':report.warnings})}


def mapping_package(reports: list[CoverageResult], *, residues_per_row: int=50) -> bytes:
    """A ZIP of plots/results for each selected RAW, never all imported RAWs."""
    out=io.BytesIO()
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as bundle:
        for i,report in enumerate(reports,1):
            directory=re.sub(r'[^A-Za-z0-9_.-]+', '_', f'raw_{i:02d}_{report.experiment}_{report.raw_file}')
            for mode in ('compact','detailed'):
                for ext in ('png','svg','pdf'):
                    bundle.writestr(f'{directory}/sequence_coverage_{mode}.{ext}',
                                    coverage_figure_bytes(report,ext,residues_per_row,view_mode=mode))
            for name,table in report_tables(report).items():
                bundle.writestr(f'{directory}/{name}.csv',csv_bytes(table))
            bundle.writestr(f'{directory}/coverage_results.xlsx',excel_bytes(report_tables(report)))
    return out.getvalue()
