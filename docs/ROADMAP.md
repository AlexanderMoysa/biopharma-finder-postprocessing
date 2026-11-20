# Roadmap

## M0 — Completed scaffold
- Independent import: native `.resultsbpf`/PMF and multiple-RAW BF CSV
- Unified component-run model, including unidentified signals
- Editable sample/preparation/replicate run labels in Streamlit
- Peptide identification comparisons with conservative fold handling
- Modifications and unidentified explorers
- Data and figure exports, tests, native-vs-CSV validation CLI

## M1 — Native-reader robustness
- Validate separate 4204 PMF (not yet provided) and additional proteins / BF versions
- Capture all BF results fields (sequence, `Mod`, `Site`, mass errors) from PMF
- Add input and validation status panel, exhaustive data integrity checks
- Reduce temporary archive-copy overhead; support very large experiments

## M2 — Advanced post-processing modules
- Site-level PTM analysis with chemically justified denominators
- Sample controls, replicates and normalization policies
- Statistical comparisons and uncertainty reporting
- Cross-experiment unidentified feature alignment with mass/RT tolerance and audit
- Sequence coverage diagrams and publication-quality report templates

## M3 — Deployment and productization
- Automated integration fixtures with anonymized licensed test files
- Persisted sessions/projects, reusable metadata mapping, reproducibility manifests
- User permissions and lab deployment options


## Version 0.2.0 implemented

- PMF text-derived sequence/Mod/positions with safeguards against claiming unknown disulfide sites.
- Confidence filtering, manual component selection, more table exports.
- Multi-RAW reference comparison, heatmap, sample/group replicate statistics, exploratory modification-area fractions.

## Next

- Validate native PMF site-level annotations on additional proteins and BF versions.
- More controlled peptide/modification grouping and more explicit replicate QC.
- Add independently validated sequence coverage and aligned unidentified-feature analysis.
