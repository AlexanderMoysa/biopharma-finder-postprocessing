# Selected-only Sequence Coverage Map (v0.4.0)

## Contract: ONLY applied selected peptides, per RAW

1. Import BF CSV or PMF/.resultsbpf and assign preparation (Reduced / Non-reduced) in section 1.
2. Use the global filters (samples, preparation, RAW files, minimum Confidence, minimum MS Area, per-RAW ID Type) and optionally deselect individual **component** rows in Component Explorer; click **Apply selection**.
3. In Peptide comparison or Multi-sample, deselect **peptide identifications** and click **Apply peptide selection**. This shared exclusion applies to all RAW injections.
4. Click **Mapping · selected peptides** under the highlighted green table. A dedicated Streamlit modal opens.
5. Choose any RAW still available under global filters. This RAW selector is **independent** of the comparison reference and therefore supports non-reduced and DTT-reduced samples separately.
6. The coverage calculation receives only the intersection of global filters, applied selected components, and applied selected peptide identifications **in this RAW**. It never restores excluded peptides and never borrows identified peaks from another RAW.
7. If all peptide identifications have been unchecked, the map is empty / 0% coverage, even though the original uploaded experiment contains peptides.

The full protein *reference sequence* can be read from `ExperimentSummary.json` in PMF. An original BF CSV lacks that reference, so it can be **provisionally inferred from the peptide IDs and sequences across the imported experiment**. This reference reconstruction does not contribute any evidence to the coverage numerator. If the reference is partial or terminal length is unknown, edit the sequence and verify its completeness manually. An unverified reference is marked with `~` for an estimated coverage percentage.

## Calculations (explicit reproducible rules)

- Map a component by BF `identification`, e.g. `1:Q5-E17/2:F1-E13 = ...[2ss]`: one selected identification contributes spans 5–17 on chain 1 and 1–13 on chain 2.
- Assign a region only when the expected amino acids and the supplied protein reference agree (or the reference has explicit unknown `?` characters). No approximate sequence matching or unreported site assignment.
- Count each amino-acid position **once** for residue coverage even if dozens of selected charge states, modifications or disulfide-linked peptides overlap it.
- Coverage = 100 × (number of distinct covered residues) / (length of the reference chain). If the full reference length is inferred rather than confirmed, the value is **provisional**.
- MS Area for a selected peptide identification is the sum of its selected positive-area component records **within one RAW**. RT is reported from the selected, largest-area component.
- Number of selected MS peaks = number of selected component records assigned to each chain (not the number of replicate injections).
- For disulfide-linked peptides spanning multiple different chains, *allocated* area is divided equally among those unique chains. This is an **exploratory area proxy**, not Thermo's exact protein MS Peak Area or molar abundance.
- The reported `Selected MS area (%)` is the allocated area for a chain / total allocated area of all mapped chains × 100. Non-mapped/unknown signals are not included in the map.
- `Abundance (mol)` remains N/A. Do **not** infer molar fractions from MS peak areas without an appropriate response/quantification model.
- Intensity colors are based on each selected peptide feature's total MS Area divided by the largest selected peptide's MS Area **in that RAW**. These are not Thermo-specific peptide recovery values.
- Additional map thresholds are **optional and default to zero**. They include minimum selected-peptide MS Area, minimum confidence, maximum peptide mass, minimum relative MS-signal percentage, and an exploratory overlap signal threshold. Setting any positive threshold can intentionally exclude some checked peptides; default zero includes all checked, detected positive-area peptides.
- Missing/zero-area peptides do not contribute coverage. Missing or ambiguous protein positions are not fabricated.

## Visual elements and downloads

- BF-inspired top summary table (`Proteins`, selected MS peaks, allocated MS Area %, residue coverage %, `Abundance (mol)` N/A).
- One chain-at-a-time number line and amino-acid sequence, with colored selected-peptide bars and RT labels, including explicit segments of interchain/intrachain disulfide peptides.
- Editable reference sequences, per-RAW selector and mapping threshold controls.
- Direct PNG, SVG, PDF image exports, mapped segments CSV, complete Excel workbook with `ProteinCoverage`, `MappedPeptides`, `MappedSegments`, `Residues`, `Parameters`, `Warnings`.
- Optional ZIP export for each currently filtered RAW, always using applied peptide selection (not all imported peptides).

## Thermo comparison / limitations

Thermo BioPharma Finder documents the map's protein coverage summary, colored peptide bars and RT labels, manual inclusion/exclusion of components, and the following parameters: maximum peptide mass, confidence, number of residues per row, minimum recovery and minimum relative recovery of overlapping peptides.

Thermo references (official):

- [BioPharma Finder 5.2 — View the Sequence Coverage Map](https://docs.thermofisher.com/r/BioPharma-Finder-5.2-Peptide-Mapping-Analysis-User-Guide/1291918731v2)
- [BioPharma Finder 5.2 — Modify Sequence Coverage Map parameters](https://docs.thermofisher.com/r/BioPharma-Finder-5.2-Peptide-Mapping-Analysis-User-Guide/1291931019v2)
- [BioPharma Finder 5.2 — Coverage Results table parameters](https://docs.thermofisher.com/r/BioPharma-Finder-5.2-Peptide-Mapping-Analysis-User-Guide/1291896459v2)

Our map is **BF-inspired, not algorithmically identical**. Thermo's exact recovery/color rules and abundance (mol) calculation are not reproduced from a native CSV or PMF. Results should be compared to BioPharma Finder for method validation before regulated use. Chain sequences inferred from CSV require verification, especially when total coverage appears 100%.
