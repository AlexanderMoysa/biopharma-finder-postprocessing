"""Import native .resultsbpf and compact experiment ZIP (RAW payload ignored)."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from bpf_postprocess.model import as_table, derive_sample_name, guess_preparation
from bpf_postprocess.importers.pmf_binary import load_core, extract_components
from bpf_postprocess.importers.annotation_parser import parse_identification


def read_resultsbpf(source: str | Path, name: str | None = None) -> tuple[pd.DataFrame, dict, list[str]]:
    path = Path(source)
    info, pmf, pmf_name = load_core(path)
    raw_rows, raw_runs, audit = extract_components(pmf, info)
    if not raw_rows:
        raise ValueError("PMF reader found no valid components; unsupported layout?")
    experiment = info.get("ExperimentName") or path.stem
    sample = derive_sample_name(experiment)
    source_file = name or path.name
    protein_chains = [
        {"number": x.get("ProteinChainNumber", i),
         "name": x.get("ProteinChainName", ""),
         "sequence": x.get("ProteinChainSequence", "")}
        for i, x in enumerate(
            ((info.get("Method") or {}).get("ProteinDigestParameters") or {}).get("ListOfProteinAliasInfo") or [], start=1)
    ]
    rows = []
    for r in raw_rows:
        annotations = parse_identification(r["Identification"], protein_chains)
        rows.append({
            "experiment": experiment, "sample": sample, "sample_group": "",
            "raw_file": r["RAW File"], "condition": r["Condition"],
            "preparation": guess_preparation(r["RAW File"]), "replicate": "",
            "component_no": r["Component No."], "identified": r["Identified"],
            "identification": r["Identification"],
            "peptide_sequence": annotations["peptide_sequence"],
            "modification": annotations["modification"], "site": annotations["site"],
            "positions": annotations["positions"], "protein": annotations["protein"],
            "rt_min": r["RT (min)"], "mz": r["M/Z"], "charge": r["Charge"],
            "mono_mass_da": r["Mono Mass Exp. (Da)"], "ms_area": r["MS Area"],
            "confidence_pct": r["Conf. Score (%)"], "delta_ppm": None,
            "id_type": "", "integration_type": "",
            "source_file": source_file, "source_type": "resultsbpf",
            "annotation_origin": annotations["annotation_origin"],
            "import_note": "Sequence and Mod inferred from PMF ID; disulfide connectivity/site may be unavailable",
        })
    metadata = {
        "experiment": experiment,
        "method": (info.get("Method") or {}).get("MethodName", ""),
        "enzyme": ((info.get("Method") or {}).get("ProteinDigestParameters") or {}).get("EnzymeName", ""),
        "protein_chains": [
            {"name": x.get("ProteinChainName"), "sequence": x.get("ProteinChainSequence")}
            for x in (((info.get("Method") or {}).get("ProteinDigestParameters") or {}).get("ListOfProteinAliasInfo") or [])
        ],
        "pmf_name": pmf_name,
        "raw_runs": [{"raw_file": x["raw_file"], "condition": x["condition"]} for x in raw_runs],
        "audit": audit,
        "source_file": source_file,
    }
    warnings = [
        "EXPERIMENTAL PMF READER: independently cross-check against Thermo CSV for new proteins, "
        "BioPharma Finder versions, and quantitative/regulated use.",
        "PMF Sequence and Mod are derived from ID text and protein reference; disulfide sites/connectivity "
        "cannot be recovered reliably from these annotations. Check against CSV for method validation.",
    ]
    if len(set(audit["unique_components_per_run"])) > 1:
        warnings.append("Different component counts between runs: check PMF compatibility")
    result = as_table(rows)
    return result, metadata, warnings
