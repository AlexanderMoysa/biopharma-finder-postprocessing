# BioPharma Finder Post-processing

**Version 0.4.0** — adds RAW-specific BioPharma Finder **ID Type** values (MS2 / Full / NoID) from original CSV, a default-all ID Type multiselect alongside confidence and MS Area filters, ID Type throughout comparison/modification/unidentified tables, and an exportable ID Type summary. Missing native PMF ID Type is clearly labeled **Unavailable (not decoded)**, never guessed as NoID. Earlier multi-RAW comparison and peptide selection features are retained.

A modular, **local web application** for Thermo BioPharma Finder peptide mapping results. Import multiple peptide `.resultsbpf` archives (or compact PMF ZIPs or original BF CSV exports), retain **each RAW run separately**, filter identified and unidentified components, compare peptide areas, and export tables and figures.

> **Experimental research tool, not validated for GxP / quantitative use.** Native PMF structure has only been independently validated against a CSV for one BioPharma Finder 5.0.68.8 experiment. Verify every new file type/version against its original BioPharma Finder output.

## Quick start — Windows / VS Code

Install [Python 3.11+](https://www.python.org/), then open the project folder in VS Code and run in PowerShell:

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
streamlit run src/bpf_postprocess/ui/app.py
```

If `py -3.13` is not available, use `python -m venv .venv` instead. The app will open at `http://localhost:8501`. `.streamlit/config.toml` raises the upload limit to 1 GiB so original >200 MiB `.resultsbpf` exports can be selected. For faster uploads and lower memory use, a compact experiment ZIP without `RawFiles.zip` is recommended when available.

## Inputs

- **Native `.resultsbpf`**: inner `Experiment.zip` → `ExperimentSummary.json` + `.pmf`. Embedded `RawFiles.zip` is **not** extracted/read; the native reader scans the PMF's component records and run boundaries.
- **Compact PMF ZIP**: `ExperimentSummary.json` + `.pmf`, useful for testing without Thermo RAW data.
- **Original BF `.csv`**: supports *multiple* `MS Area:<raw>.raw` columns; other run-specific fields are normalized into individual rows. CSV provides direct `Mod`, `Site`, `Positions`, `Peptide Sequence`. In v0.3.0 the PMF importer reconstructs the peptide sequence, positions, chain name, and modification labels from the ID annotation and the embedded protein sequences. These were checked against the 4197 CSV (234/234 matched for sequence and Mod). **Disulfide connectivity and all site annotations cannot be reliably recovered from the ID string** and may remain blank; prefer CSV for complete native site assignments.

**ID Type (v0.3.3):** BioPharma Finder CSV exports contain both the experiment-wide `ID Type` (which can be `Both`) and run-scoped columns `ID Type:<RAW file>`. The importer preserves the run-specific class (`Full`, `MS2`, or `NoID`) as `id_type`, and the global class separately as `id_type_overall`. When using a native PMF/`.resultsbpf` input, ID Type is **not currently decoded from the binary format**: these records show `Unavailable (not decoded)`. This does **not** imply `NoID` or a low confidence score. If you need the true identification type, **import the original BioPharma Finder CSV** for that experiment instead of its native archive; the CSV retains the component records for both RAW runs. An on-screen warning appears when type is unavailable.

Load one representation **per experiment**. Uploading the corresponding CSV *and* `.resultsbpf` together is rejected as duplicate data.

A single experiment can contain many RAW runs (for example, reduced and non-reduced injections). The interface initially infers `Reduced` only when explicitly indicated by filename; set the other preparation labels manually in the editable RAW run table.

## Modules / pages (v0.4.0)

- **ID Type filter and displays (v0.3.3)**: Sidebar multiselect `ID Type (per RAW)` defaults to **all** categories and combines with Confidence (%) and Minimum MS Area. Component Explorer shows `id_type` and `id_type_overall`; peptide and multi-sample comparisons show ID Type for every RAW, preserving multiple different types as an explicit list where peptide components are summed. Modifications, relative PTM area fractions, unknown-component summaries and Excel exports retain the corresponding types. Export Report includes an `IDTypeSummary` sheet.
- **Component Explorer**: all rows start checked; uncheck individual component/RUN records and click **Apply selection** to exclude them from downstream analysis. Selection persists across filters. Select/Deselect All visible controls are provided. Set a minimum `confidence_pct` in the sidebar. A positive threshold excludes missing-confidence/unidentified records.
- **Experiment/RAW metadata**: edit Sample ID, Group (DS, DP, stressed, control), Preparation and Replicate labels.
- **RAW comparison configuration**: choose one reference RAW run and manually select any number of comparison RAW files. By default only RAW runs with the same `Preparation` label (Reduced / Non-reduced / Other / Unspecified) are eligible. Set missing preparation labels before comparison. This comparison-specific selection does not remove runs from Component Explorer or Modifications.
- **Peptide comparison**: per-feature annotated matrix with peptide sequence, modification, site, positions, protein and **MS Area, dominant-peak RT (min), component count for each selected RAW**. Areas of multiple components with the exact same BF identification are summed. RT/charge/mz/confidence in aggregated comparisons refer to the highest-area component in that RAW, never an invented averaged peak. If the summed area is zero, dominant RT and other measurements are blank. A separate downloadable **Individual components** table preserves *all* component RT, m/z, charge, mass, confidence, site, source/annotation fields and areas. Downloadable figures (PNG/SVG/HTML) remain available.
- **Multi-sample**: an explicitly selected RAW reference compared only with the **manually chosen RAW runs of the same preparation**, fold/percent changes (undefined if reference area is 0/missing). The reference comparison table includes peptide annotation and separately labelled run-specific representative RT, m/z, charge, mass, confidence, component count and area. Heatmap PNG/SVG/HTML and grouped plots use the same chosen run subset; no automatic reduced/non-reduced mixing.
- **Replicate statistics**: one aggregated peptide area per RAW run; n, mean, sample SD and %RSD grouped by Sample ID or user-assigned Group and preparation. With n=1, SD and %RSD are undefined.
- **Modification analysis**: BF modification labels, exploratory relative MS-area percentage for a selected modification among *unmodified + selected target PTMs of the same peptide backbone*. A missing unmodified reference returns blank, not 100%. This is not site occupancy and is not corrected for ionization efficiency.
- **Unidentified components**: kept separate; no cross-experiment matching by Component No.
- **Selection in every analysis tab (v0.3.2)**: use the **Use in analysis** checkbox at the start of each peptide matrix, RAW-labeled comparison, modification peptide inventory or export selection grid; click **Apply peptide selection**. All peptides are checked by default. The exclusion is keyed to full BF peptide identification and applies to every RAW and all tabs, including plots and selected-result exports. The editor retains unchecked rows so they can be checked again. **Unidentified** has no peptide ID, so it provides independent RAW/component checkboxes. All controls also offer **Select all visible** / **Deselect all visible**. Peptide selection is in addition to Component Explorer row selection.
- **Downloads**: selected analytical tables and full workbook explicitly as CSV/Excel; charts as PNG/SVG/HTML where available. `SelectedComponents` contains only chosen peptides/components; `VisibleComponents` and `PeptideSelection` are audit sheets that intentionally include filtered-out items with selection metadata for reproducibility.

The reference and comparator controls appear in **section 3 (Choose RAW files for peptide comparison)** above the analysis tabs. The first comparison table and the heatmap reflect only those RAW files. If the run's preparation shows `Unspecified`, edit it in **section 1** and click Apply RAW run labels. For peptide comparability, also verify protein identity and digestion/BF analysis method.

For component equivalence between experiments, peptide comparisons use an exact identification key across native PMF and CSV, assuming samples share a compatible protein and BioPharma Finder processing method. **No automated deconvolution, unidentified alignment, automatic response normalization, or modification occupancy interpretation.**

## Sequence Coverage Map — selected peptides ONLY (v0.4.0)

The **Peptide comparison** and **Multi-sample** checkbox areas now have a green rounded border. Under each highlighted table, click **🧬 Mapping · selected peptides** to open a separate modal. **First click Apply peptide selection** after unchecking peptides. The map is built independently for each RAW in the global filter scope, including Reduced and Non-reduced RAWs, using only **applied** component-level and peptide-level selections. An unselected peptide is never reintroduced, even if it was present in the original BF export. An empty selection yields 0% coverage.

The map shows the complete protein/chains, numbered amino acids, color-coded overlapping peptide bars and RT labels, per-chain sequence coverage percentages, selected component peak/area summary, and exports in **PNG / SVG / PDF / CSV / XLSX**. A ZIP can export separate maps/results for multiple filtered RAWs. PMF provides full reference sequences; original CSV can seed reference sequences from positional peptide identifications, but the user must verify inferred sequence completeness. This inferred reference *does not supply* coverage evidence.

**Not Thermo-equivalent**: color thresholds and area shares are explicitly defined MS-signal proxies, not Thermo's recovery or calibrated molar abundance. Molar abundance is left as N/A; no undocumented occupancy/quantitation is guessed. See [complete mapping design, calculations, limitations and Thermo documentation](docs/SEQUENCE_COVERAGE.md).

## Updating a downloaded ZIP project

The initial folder downloaded as ZIP does not contain Git history. Download the updated ZIP from GitHub **Code → Download ZIP** and replace the older project folder, then run `.\.venv\Scripts\python.exe -m pip install -e ".[test]"` after creating a venv if needed. For subsequent updates, prefer cloning the repository with `git clone https://github.com/AlexanderMoysa/biopharma-finder-postprocessing.git` and `git pull` in the project folder.

## Command-line import

```powershell
bpf-postprocess "C:\path\experiment_1.resultsbpf" "C:\path\experiment_2.csv" --out outputs
```

For a **new** PMF dataset, validate the reader against BioPharma Finder's own CSV export before interpreting results:

```powershell
bpf-validate-pmf "C:\path\experiment.resultsbpf" "C:\path\experiment.csv" --out outputs\pmf_validation.json
```

Emits `components.csv` (one row per RAW/component), `runs.csv`, `experiments.json`, and warnings. Run in the same Python virtual environment as the app.

## Repository layout

```
src/bpf_postprocess/
  importers/    # native nested ZIP + PMF binary parser, CSV parser, file dispatch
  model.py      # normalized long-form component/run schema
  processing/   # filters, peptide comparison, modifications, unidentified
  reports/      # downloadable CSV / XLSX / PNG / SVG
  ui/app.py     # Streamlit UI (no business logic inside importers)
  cli.py        # reproducible batch export
tests/          # synthetic unit tests; optional local instrument fixture
.github/workflows/tests.yml
```

## Validation and privacy

- Run `pytest -q` (or `BPF_TEST_CORE_ZIP=/path/to/core.zip pytest -q` on Linux/macOS for optional native fixture; Windows can set `$env:BPF_TEST_CORE_ZIP="..."`).
- Use original BF CSV as reference when bringing in a new PMF format, instrument processing method or BF version.
- `.resultsbpf`, `.pmf`, `.raw`, experiment exports, and real datasets are excluded from Git by `.gitignore`. Do not commit confidential lab data.
- The Streamlit app runs locally and has no network/API integration. Dependencies need normal pip access during installation.

See [Russian quick start](docs/GETTING_STARTED_RU.md), [roadmap](docs/ROADMAP.md), and [architecture](docs/ARCHITECTURE.md) for the data contract and known limitations.

### v0.3.1 — RAW-labeled multi-sample comparison

The primary **Multi-sample & replicates** table is now one row per identified peptide.
Every quantitative column includes the exact `Experiment / RAW filename`, and the
reference metrics carry `[REFERENCE]`. Each selected comparator has its own
`MS Area`, `RT`, `component count`, `Delta MS Area`, `Fold`, and `Change (%)` columns.
This applies to 2, 3, or more RAW runs without mixing reduced and non-reduced
preparations. The original long-form detailed comparison remains available in
a collapsible section and both forms are included in Excel/CSV exports.
Missing peptide signals remain blank (NaN), not assumed to be zero.
