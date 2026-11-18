# Architecture

```
BioPharma Finder .resultsbpf / core ZIP / original CSV
                |
         importers/loader.py
         /                \
 importers/pmf_binary.py   importers/csv_importer.py
 importers/results_importer.py
         \                /
          model.py (long-format DataFrame)
               |
     experiment -> RAW run -> component
               |
 processing/{filters,comparison,modifications,unidentified}.py
               |
  reports/exports.py <-> ui/app.py (Streamlit)
```

## Data normalization contract

One row = one **component, in one RAW run, in one experiment**. Keep non-reduced and reduced signals separate. Component number is stable *within a given experiment's component table* but never a universal peptide ID across experiments.

- `experiment`: BF experiment name
- `sample`: editable user label, initialized from experiment name
- `raw_file`: exact RAW filename reported by BF
- `condition`: BF condition, when available in experiment JSON
- `preparation`: editable label (only explicit DTT/reduced filenames are inferred)
- `component_no`: BF component number
- `identified`, `identification`: whether and how BF assigned the component
- `ms_area`, `mz`, `rt_min`, `charge`, `mono_mass_da`: run-specific signal fields
- `peptide_sequence`, `modification`, `site`, `positions`: supported by original CSV, not yet separated from native PMF text

## Processing contract

- Every analysis is a pure function accepting normalized dataframes.
- `comparison` keys identified features by the full BF `Identification` annotation, available in BOTH CSV and PMF. Cross-experiment comparisons are meaningful only when the same target protein and compatible digestion/processing methods were used; inspect experiment metadata and run labels before comparing.
- Summation across multiple features with the same ID is explicit; no area normalization or cross-condition merging is assumed.
- Fold change with zero/missing reference remains null (not `9999.99`).
- Unidentified features remain visible; no cross-experiment mass/RT alignment yet.
- Site occupancy, recovery, absolute quantification, and statistical testing need a separately validated definition.

## PMF reverse-engineering limitations

The PMF format is undocumented and reader offsets/markers are empirically derived from one BioPharma Finder 5.0.68.8 peptide experiment. Although both unidentified and identified components were verified against its CSV export, this **does not prove general compatibility**. Safeguard new runs through CSV cross-validation, optional PMF fixture tests, audit checks and error handling. Never treat this research software as a validated regulated analytical method.

## Planned extensions

- Site-specific modification parsing and abundance normalization
- Replicate grouping with statistical treatment
- Reference protein sequence coverage and PTM localization
- Matching unidentified components by validated mass + RT tolerances
- PDF reports and reproducible provenance manifests
- Native PMF support for further versions / method types
