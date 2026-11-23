# BioPharma Finder source CSV fields, source fidelity, and PMF enrichment

## Why there were missing fields

A BioPharma Finder peptide CSV contains **one row per component across all RAW
injections**, including experiment-global annotations and run-specific columns.
The native PMF importer produces **one normalized row per component per RAW**.
It originally reconstructed peptide sequences and modifications but did not
extract BF ID Type, Site, or Delta ppm. They cannot be guessed from PMF confidence
or mass alone.

## v0.4.1 import options

| Input | ID Type / Delta ppm / Site | Numerical measurements | BF source table |
|---|---|---|---|
| Original BF CSV only | Direct CSV | CSV export precision | Preserved completely |
| PMF / resultsbpf only | Unavailable / partial | PMF float32 precision | Not available |
| **Matching PMF + BF CSV** | **Direct verified CSV** | **PMF precision** | **Preserved completely** |

To combine PMF with CSV, upload **both** files simultaneously. Matching is
strict on `(experiment, raw_file, component_no)` and exact identification text.
MS Area, RT, m/z, charge, experimental monoisotopic mass and confidence scores
are compared with export-rounding tolerances of **0.51**, **0.0051 min**,
**0.00051 m/z**, **0**, **0.000051 Da**, and **0.051 percentage points**,
respectively, whenever both source values exist. The component key sets must
match exactly. A mismatch **aborts import**, not a partial/silent merge.

The original PMF measurement values (`ms_area`, `rt_min`, `mz`,
`mono_mass_da`, `confidence_pct`) are kept after validation. The CSV supplies:
`id_type` per RAW, `id_type_overall` (`Both` is possible), `delta_ppm` per RAW,
`delta_ppm_overall`, `site`, `peptide_sequence`, `modification`,
`positions`, `protein`, `integration_type`, `best_asr`, `avg_mass_da`,
`mono_mass_theo_da`, `missed_cleavages`, `bf_ratio_label`, `bf_ratio`,
`bf_max_condition`, `bf_min_condition`, `bf_max_ms_area`,
`bf_avg_ms_area_condition`, `bf_cv_pct_condition`, and `comment`.

**Per-RAW blanks remain blank.** If a RAW-specific Delta ppm or ID Type column
exists but its cell is empty, it is not populated from a different run or
from global `ID Type: Both`. Missing RAW-specific fields are not invented.

## All original columns retained

All original BF CSV column headings and cell texts are stored in an audit table
per experiment. For the actual 4197 and 4204 CSVs, there are **43 source
columns** per experiment. Generic classes are:

- Component identification: `Level`, `No.`, `Identification`, `Positions`,
  `Peptide Sequence`, `Mod`, `Site`, `∆ ppm`, `Conf. Score`, `Best ASR`,
  `ID Type`, `Integration Type`, `Protein`, `Comment`
- Component measurements: `RT`, `M/Z`, `Charge St.`, `Mono Mass Exp.`,
  `Avg Mass Exp.`, `Mono Mass Theo.`, `Missed Cleavages`
- Experimental condition comparison: `Ratio (<conditions>)`, `Max Cond.`,
  `Min Cond.`, `Max MS Area`, `Avg MS Area: <condition>`, `%CV: <condition>`
- RAW measurements: `MS Area:<RAW>`, `∆ ppm:<RAW>`, `ID Type:<RAW>`,
  `RT:<RAW>`, `M/Z:<RAW>`, `Mono Mass Exp.:<RAW>`, `Avg Mass Exp.:<RAW>`

The last groups repeat once per condition and/or RAW. Column counts may differ
between experiments, so the audit view preserves the actual source headings
rather than assuming there are always 43.

## How to find the original unmodified data

In **Component Explorer**, expand **Original BioPharma Finder CSV · all source
columns (audit)**, select the experiment, and download CSV or Excel. The
**Export report** workbook includes an audit sheet `BF_<experiment>` for each
uploaded CSV. These BF sheets deliberately contain every original source row,
including components deselected from analytical calculations.

The **Peptide comparison** table is *not* a reproduction of the original BF
source CSV: it aggregates same-identification components **within each RAW**.
Its representative RT, Delta ppm, Best ASR, average mass, theoretical mass and
missed-cleavage fields are labeled **dominant** because they come from the
largest-MS-area component, not an average. Component-level data remain in
**Individual components**. Comparison areas are observed signals, not molar
concentrations.

`Ratio (<conditions>)`, Max/Min conditions, Avg MS Area and `%CV` describe the
**original BF experiment**. Do not substitute these for our independently
selected comparisons, especially when contrasting Reduced and Non-reduced.
A source `%CV = 0` for one injection per condition is not evidence of perfect
reproducibility.

## Limitations

- PMF parsing remains experimental and is not validated for every BF version.
- Source CSV/PMF pairing assumes consistent BF processing output. Reprocessed
  CSVs with changed IDs or values fail validation.
- The source CSV audit is lossless in column values and names, but the normalized
  data table intentionally converts numerical fields to typed Python numbers.
- The app remains a local research tool, not a GxP-validated calculation system.
