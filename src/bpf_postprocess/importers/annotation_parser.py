"""Conservatively derive *partial* peptide annotations from BioPharma Finder ID text.

Native PMF does not yet expose its separate sequence/Mod/Site fields. This
module uses protein sequences from ExperimentSummary.json and does not claim to
recover undecodable disulfide connectivity or unreported modification sites.
"""
from __future__ import annotations
import re

SEGMENT_RE = re.compile(r"(?P<chain>\d+):(?P<start_aa>[A-Z])(?P<start>\d+)-(?P<end_aa>[A-Z])(?P<end>\d+)")
MOD_RE = re.compile(r"(?:\(([^()]*)\)|\[([^\[\]]*)\])")
SITE_MOD_RE = re.compile(r"(~?[A-Z]\d+)\+([^,;]+)")


def parse_identification(identification: str, protein_chains: list[dict]) -> dict[str, str]:
    """Return BF-like columns, or blanks where the native ID is insufficient.

    Protein index is mapped by ProteinChainNumber when present, otherwise by the
    ordered list from BioPharma Finder. Bounds and residue identities are checked.
    """
    result = dict(peptide_sequence="", modification="", site="", positions="", protein="", annotation_origin="pmf_id_derived_partial")
    if not identification:
        return result
    left, separator, right = identification.partition("=")
    if not separator:
        return result
    chains = {
        int(c.get("number", i)): c
        for i, c in enumerate(protein_chains, start=1)
    }
    regions = [p.strip() for p in left.strip().split("/")]
    peptides, positions, names = [], [], []
    for region in regions:
        match = SEGMENT_RE.fullmatch(region)
        if not match:
            return result
        g = match.groupdict()
        chain = chains.get(int(g["chain"]))
        if chain is None:
            return result
        seq = str(chain.get("sequence", "")).strip().upper()
        start, end = int(g["start"]), int(g["end"])
        if start < 1 or start > end or end > len(seq):
            return result
        fragment = seq[start - 1:end]
        if fragment[0] != g["start_aa"] or fragment[-1] != g["end_aa"]:
            return result
        peptides.append(fragment)
        positions.append(f"{start}-{end}")
        names.append(str(chain.get("name", g["chain"])))
    result["peptide_sequence"] = "/".join(peptides)
    result["positions"] = " | ".join(positions)
    result["protein"] = names[0] if names else ""
    # The right hand portion is a mass plus optional modification annotations.
    # Avoid interpreting the mass, charge, or fragments in the left peptide part.
    annotations = [a.strip() for m in MOD_RE.finditer(right) for a in m.groups() if a and a.strip()]
    labels, sites = [], []
    for annotation in annotations:
        for part in annotation.split(","):
            part = part.strip()
            if not part:
                continue
            found = SITE_MOD_RE.fullmatch(part)
            if found:
                sites.append(found.group(1))
                labels.append(found.group(2))
            else:
                labels.append(part)
    result["modification"] = ",".join(labels) if labels else "None"
    result["site"] = ",".join(dict.fromkeys(sites))
    return result
