"""Selected-only, RAW-scoped protein sequence coverage for peptide mapping.

A coverage residue is supported only by explicitly selected, confidently
mapped peptide records from *one* RAW injection. Selection occurs upstream;
this module NEVER consults additional/hidden or unselected components to fill
coverage. Imported data may be used to construct the protein REFERENCE alone.

Thermo BioPharma Finder's recovery, abundance (mol), and its signal-based color
calibration are proprietary/processing-method dependent; this module exposes
its own deterministic, explicitly labeled MS-area proxies instead.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import re
from collections import defaultdict
from typing import Mapping, Sequence

import pandas as pd

SEGMENT = re.compile(r"(?P<chain>\d+):(?P<start_aa>[A-Z])(?P<start>\d+)-(?P<end_aa>[A-Z])(?P<end>\d+)")
LABEL_MASS = re.compile(r"=\s*(\d+(?:\.\d+)?)m")


@dataclass(frozen=True)
class ProteinChain:
    number: int
    name: str
    sequence: str
    source: str = "PMF experiment metadata"
    verified: bool = True


@dataclass
class CoverageResult:
    run: str
    experiment: str
    raw_file: str
    preparation: str
    chains: list[ProteinChain]
    summary: pd.DataFrame
    peptides: pd.DataFrame
    segments: pd.DataFrame
    residue_table: pd.DataFrame
    parameters: dict
    warnings: list[str]


def parse_peptide_regions(feature: str) -> list[tuple[int, int, int, str, str]]:
    """Use the chain-index and peptide region in BF ID (NOT ambiguous Site)."""
    region_text = str(feature).partition("=")[0].strip()
    if not region_text:
        return []
    found = []
    for part in region_text.split("/"):
        match = SEGMENT.fullmatch(part.strip())
        if match is None:
            return []
        idx, begin, end = int(match['chain']), int(match['start']), int(match['end'])
        if begin < 1 or end < begin:
            return []
        found.append((idx, begin, end, match['start_aa'], match['end_aa']))
    return found


def _metadata_chains(experiment: str, metadata: pd.DataFrame) -> list[ProteinChain]:
    if metadata is None or metadata.empty or 'protein_chains' not in metadata:
        return []
    candidates = metadata.loc[metadata['experiment'].astype(str) == str(experiment)]
    for _, entry in candidates.iterrows():
        items = entry.get('protein_chains')
        if not isinstance(items, list) or not items:
            continue
        chains = []
        for index, item in enumerate(items, 1):
            seq = str(item.get('sequence') or '').strip().upper()
            if not seq or not re.fullmatch('[A-Z]+', seq):
                continue
            chains.append(ProteinChain(index, str(item.get('name') or index), seq))
        if chains:
            return chains
    return []


def _infer_chains(experiment: str, source: pd.DataFrame) -> list[ProteinChain]:
    """Conservatively seed references from available BF peptide sequences.

    This is REFERENCE reconstruction only. It must not change the list of
    selected peptides used to compute coverage. Missing residues are '?' and
    must be supplied manually to give meaningful full-chain percentages.
    """
    records = source.loc[
        (source['experiment'].astype(str) == str(experiment)) & source['identified']
    ]
    seen = defaultdict(dict)
    end_positions = defaultdict(int)
    names: dict[int, str] = {}
    for _, record in records.iterrows():
        regions = parse_peptide_regions(record['identification'])
        sequences = str(record.get('peptide_sequence') or '').strip().upper().split('/')
        if not regions or len(sequences) != len(regions):
            continue
        for (idx, start, end, first, last), peptide in zip(regions, sequences):
            if len(peptide) != end-start+1 or not peptide or peptide[0] != first or peptide[-1] != last:
                continue
            end_positions[idx] = max(end_positions[idx], end)
            name = str(record.get('protein') or '').strip()
            if name and len(regions) == 1:
                names.setdefault(idx, name)
            for pos, aa in enumerate(peptide, start):
                if pos in seen[idx] and seen[idx][pos] != aa:
                    seen[idx][pos] = '?'
                elif pos not in seen[idx]:
                    seen[idx][pos] = aa
    return [ProteinChain(
        number=idx, name=names.get(idx, f'Chain {idx}'),
        sequence=''.join(seen[idx].get(p, '?') for p in range(1, length+1)),
        source='Inferred from all imported identifications (verify length and residues)',
        verified=False,
    ) for idx, length in sorted(end_positions.items())]


def reference_chains(experiment: str, source: pd.DataFrame, metadata: pd.DataFrame | None = None) -> list[ProteinChain]:
    """Prefer complete native PMF reference; otherwise infer from peptide IDs.

    Independent sequence metadata may include unselected features. Those
    features are *never* used to set coverage or draw peptide bars.
    """
    native = _metadata_chains(experiment, metadata)
    inferred = _infer_chains(experiment, source)
    by_number = {c.number: c for c in inferred}
    by_number.update({c.number: c for c in native})
    return [by_number[i] for i in sorted(by_number)]


def check_sequences(chains: Sequence[ProteinChain]) -> list[ProteinChain]:
    out = []
    for c in chains:
        seq = re.sub(r'\s+', '', c.sequence).upper()
        if not seq or not re.fullmatch(r'[A-Z?]+', seq):
            raise ValueError(f'Invalid or empty reference sequence for chain {c.number}')
        if c.number < 1:
            raise ValueError('BF chain numbers start at 1')
        out.append(replace(c, sequence=seq, verified=bool(c.verified and '?' not in seq)))
    if len({c.number for c in out}) != len(out):
        raise ValueError('Duplicate protein chain number')
    return out


def calculate_coverage(
    selected_components: pd.DataFrame,
    run: str,
    chains: Sequence[ProteinChain],
    *, minimum_ms_area: float = 0.0,
    minimum_confidence_pct: float = 0.0,
    maximum_peptide_mass_da: float = 0.0,
    minimum_relative_signal_pct: float = 0.0,
    minimum_overlap_signal_pct: float = 0.0,
) -> CoverageResult:
    """Calculate residue coverage per RAW from ONLY selected component records.

    Relative signal = peptide MS Area / highest selected peptide MS Area in that
    RAW x 100. Overlap proxy = local peptide area / strongest locally-overlapping
    peptide area x 100; neither is Thermo's proprietary recovery metric.
    """
    p = dict(minimum_ms_area=minimum_ms_area,
             minimum_confidence_pct=minimum_confidence_pct,
             maximum_peptide_mass_da=maximum_peptide_mass_da,
             minimum_relative_signal_pct=minimum_relative_signal_pct,
             minimum_overlap_signal_pct=minimum_overlap_signal_pct)
    if any(float(v) < 0 for v in p.values()):
        raise ValueError('Coverage thresholds must be non-negative')
    references = check_sequences(chains)
    ids = selected_components['experiment'].astype(str) + ' / ' + selected_components['raw_file'].astype(str)
    rows = selected_components.loc[ids == run].copy()
    if rows.empty:
        experiment, separator, raw_file = run.partition(' / ')
        if not separator:
            raise ValueError('Expected RAW run ID: experiment / raw_file')
        # This is intentional: deselecting every peptide must yield ZERO coverage,
        # never fall back to using hidden or original RAW components.
        preparation = 'No selected components'
    else:
        experiment, raw_file = str(rows.iloc[0]['experiment']), str(rows.iloc[0]['raw_file'])
        preparation = str(rows.iloc[0]['preparation'])
    cmap = {c.number: c for c in references}
    warnings = []
    for chain in references:
        if not chain.verified:
            warnings.append(f'{chain.name}: reference sequence inferred/unverified; coverage percentage is provisional.')
        if '?' in chain.sequence:
            warnings.append(f'{chain.name}: incomplete reference with unknown residues; complete it manually for reliable coverage.')

    identified = rows.loc[rows['identified'] & rows['identification'].fillna('').astype(str).str.strip().ne('')].copy()
    identified['feature'] = identified['identification'].astype(str).str.strip()
    numeric_cols = ['ms_area', 'rt_min', 'confidence_pct', 'mono_mass_da']
    for field in numeric_cols:
        identified[field] = pd.to_numeric(identified[field], errors='coerce')
    identified = identified.loc[identified['ms_area'].notna() & (identified['ms_area'] > 0)]
    # NOTE: zero MS area is not a detected coverage peak.
    peptides = []
    for feature, members in identified.groupby('feature', sort=True):
        regions = parse_peptide_regions(feature)
        if not regions:
            warnings.append(f'Cannot map unrecognized identification: {feature}')
            continue
        segments = []
        for index, start, end, first, last in regions:
            chain = cmap.get(index)
            if chain is None or end > len(chain.sequence):
                warnings.append(f'No valid reference chain/positions for: {feature}')
                segments = []
                break
            expected = chain.sequence[start-1:end]
            sequence_parts = str(members.iloc[0].get('peptide_sequence') or '').upper().split('/')
            # ID anchor and sequence span must agree. '?' reference residues
            # are allowed provisionally but never filled from hidden peptides.
            if expected[0] not in ('?', first) or expected[-1] not in ('?', last):
                warnings.append(f'ID residues disagree with reference for: {feature}')
                segments = []
                break
            if len(sequence_parts) == len(regions):
                pep = sequence_parts[len(segments)].strip()
                if pep and (len(pep) != end-start+1 or any(a != '?' and a != b for a,b in zip(expected, pep))):
                    warnings.append(f'Peptide Sequence disagrees with reference for: {feature}')
                    segments = []
                    break
            segments.append((index,start,end))
        if not segments:
            continue
        total_area = float(members['ms_area'].sum())
        if total_area < minimum_ms_area:
            continue
        conf = members['confidence_pct'].max()
        if minimum_confidence_pct > 0 and (pd.isna(conf) or conf < minimum_confidence_pct):
            continue
        # BF annotations encode the theoretical unmodified species mass,
        # while the measured component mono mass accounts for metal adducts.
        masses = members['mono_mass_da'].dropna()
        match = LABEL_MASS.search(feature)
        mass = float(masses.max()) if not masses.empty else (float(match.group(1)) if match else float('nan'))
        if maximum_peptide_mass_da > 0 and (pd.isna(mass) or mass > maximum_peptide_mass_da):
            continue
        dominant = members.sort_values('ms_area', ascending=False,kind='stable').iloc[0]
        chains_for_feature = sorted({index for index,_,_ in segments})
        peptides.append({
            'experiment':experiment,'raw_file':raw_file,'run':run,'preparation':preparation,
            'feature':feature, 'peptide_sequence':str(dominant['peptide_sequence']),
            'modification':str(dominant['modification']), 'id_type':'; '.join(sorted(set(members['id_type'].fillna('').astype(str)))),
            'ms_area':total_area, 'rt_min':dominant['rt_min'], 'confidence_pct':conf,
            'mono_mass_da':mass, 'component_count':len(members),
            'chain_count':len(chains_for_feature), 'chain_numbers':','.join(map(str,chains_for_feature)),
            'regions':segments,
        })
    peak = max((rec['ms_area'] for rec in peptides),default=0.)
    for rec in peptides:
        rec['relative_signal_pct'] = 100*rec['ms_area']/peak if peak else 0.
    if minimum_relative_signal_pct > 0:
        peptides = [rec for rec in peptides if rec['relative_signal_pct']>=minimum_relative_signal_pct]
    # Optional local overlap proxy, comparing signal of overlapping peptide
    # annotations on the SAME chain in this RAW.
    if minimum_overlap_signal_pct > 0:
        supported=[]
        for p0 in peptides:
            local_max = max((p1['ms_area'] for p1 in peptides
                if any(i0==i1 and lo0<=hi1 and lo1<=hi0
                       for i0,lo0,hi0 in p0['regions'] for i1,lo1,hi1 in p1['regions'])),default=p0['ms_area'])
            if local_max and 100*p0['ms_area']/local_max >= minimum_overlap_signal_pct:
                supported.append(p0)
        peptides = supported
    segment_rows=[]
    covered={c.number:set() for c in references}
    assigned_area=defaultdict(float)
    assigned_peaks=defaultdict(int)
    for rec in peptides:
        unique_chains={index for index,_,_ in rec['regions']}
        # Allocate one detected component equally between involved chains.
        # This is a defined proxy, NOT Thermo's response-normalized mol value.
        for idx in unique_chains:
            assigned_area[idx]+=rec['ms_area']/len(unique_chains)
            assigned_peaks[idx]+=rec['component_count']
        for idx,start,end in rec['regions']:
            covered[idx].update(range(start,end+1))
            segment_rows.append(dict(run=run,experiment=experiment,raw_file=raw_file,preparation=preparation,
                protein=cmap[idx].name,chain_number=idx,start=start,end=end,feature=rec['feature'],
                peptide_sequence=rec['peptide_sequence'],modification=rec['modification'],
                id_type=rec['id_type'],ms_area=rec['ms_area'],rt_min=rec['rt_min'],
                confidence_pct=rec['confidence_pct'],relative_signal_pct=rec['relative_signal_pct'],
                component_count=rec['component_count']))
    total_mapped=sum(assigned_area.values())
    summary=[]
    residue_rows=[]
    for chain in references:
        seq=chain.sequence
        cov=covered[chain.number]
        for index, aa in enumerate(seq,1):
            residue_rows.append({'run':run,'protein':chain.name,'chain_number':chain.number,
                                 'position':index,'amino_acid':aa,'covered':index in cov})
        summary.append({'run':run,'experiment':experiment,'raw_file':raw_file,'preparation':preparation,
                        'protein':chain.name,'chain_number':chain.number,'sequence_length':len(seq),
                        'covered_residues':len(cov),
                        'sequence_coverage_pct':100.*len(cov)/len(seq),
                        'coverage_verified':chain.verified,'number_of_ms_peaks':assigned_peaks[chain.number],
                        'allocated_ms_area':assigned_area[chain.number],
                        'selected_ms_peak_area_pct':100*assigned_area[chain.number]/total_mapped if total_mapped else 0.,
                        'abundance_mol_pct':float('nan'),
                        'reference_source':chain.source})
    report_cols = ['experiment','raw_file','run','preparation','feature','peptide_sequence','modification','id_type','ms_area',
                   'rt_min','confidence_pct','mono_mass_da','component_count','chain_count','chain_numbers','relative_signal_pct']
    segcols=['run','experiment','raw_file','preparation','protein','chain_number','start','end','feature','peptide_sequence',
             'modification','id_type','ms_area','rt_min','confidence_pct','relative_signal_pct','component_count']
    return CoverageResult(run,experiment,raw_file,preparation,references,
                          pd.DataFrame(summary),
                          pd.DataFrame([{k:v for k,v in rec.items() if k!='regions'} for rec in peptides],columns=report_cols),
                          pd.DataFrame(segment_rows,columns=segcols),
                          pd.DataFrame(residue_rows),p,warnings)
