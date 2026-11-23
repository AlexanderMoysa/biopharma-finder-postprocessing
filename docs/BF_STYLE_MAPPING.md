# Thermo-style mapping comparison — implementation notes

Based on selected-only Reduced `4197p26_with_DTT_07.raw` from BioPharma Finder.
The 80 selected component records give 22 A-chain and 58 B-chain peaks;
coverage 100%/100%; selected MS Area A 24.1%, B 75.9%.

## Evidence from supplied BioPharma Finder installation libraries

- `ProteinSequenceMap.htm`: HTML/JS sequence highlighting, no recovery calculation.
- `ThermoFisher.BioPharma.PeptideMapping.Shared.dll`: .NET user-interface resources,
  coverage map controls and parameters (`MinimumRecovery`,
  `MinimumRelativeRecoveryOfOverlappingPeptides`), adapters for obtaining
  processed coverage results.
- `ThermoFisher.BioPharma.PeptideMapping.Algorithm.Parts.dll`: sample coverage
  HTML and a client to the peptide-mapping REST service; not enough to establish
  the numerical Peptide Recovery and molar-abundance algorithm.

## Supported comparison and limitations

The compact display groups selected identifications by (chain, start, end),
chooses the most abundant identity's RT and uses a per-chain normalized regional
MS-area proxy for color classes. All included components still contribute to
counts and areas; `BFStyleRegions` documents all intervals. The detailed view
remains available. No unselected peptide is ever added to coverage.

The main visually similar BF color assignments for the selected insulin glargine
example are expected for A1-4 (red), A5-17 (red), A1-17 (green), A18-21
(orange), B1-13 (red), B14-21 (orange), B22-32 (red), B1-32 (blue).

The visual proxy and BF Recovery may diverge for complex mixtures, linked
peptides, and distinct compounds; do not interpret color as true recovery,
chemical occupancy or molar abundance.
