# BioPharma Finder Post-processing

**Version 0.2.0** — multi-run comparison, user-selected components, confidence filtering, PMF annotation reconstruction.

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
- **Original BF `.csv`**: supports *multiple* `MS Area:<raw>.raw` columns; other run-specific fields are normalized into individual rows. CSV provides direct `Mod`, `Site`, `Positions`, `Peptide Sequence`. In v0.2.0 the PMF importer reconstructs the peptide sequence, positions, chain name, and modification labels from the ID annotation and the embedded protein sequences. These were checked against the 4197 CSV (234/234 matched for sequence and Mod). **Disulfide connectivity and all site annotations cannot be reliably recovered from the ID string** and may remain blank; prefer CSV for complete native site assignments.

Load one representation **per experiment**. Uploading the corresponding CSV *and* `.resultsbpf` together is rejected as duplicate data.

A single experiment can contain many RAW runs (for example, reduced and non-reduced injections). The interface initially infers `Reduced` only when explicitly indicated by filename; set the other preparation labels manually in the editable RAW run table.

## Modules / pages (v0.2.0)

- **Component Explorer**: all rows start checked; uncheck individual component/RUN records and click **Apply selection** to exclude them from downstream analysis. Selection persists across filters. Select/Deselect All visible controls are provided. Set a minimum `confidence_pct` in the sidebar. A positive threshold excludes missing-confidence/unidentified records.
- **Experiment/RAW metadata**: edit Sample ID, Group (DS, DP, stressed, control), Preparation and Replicate labels.
- **Peptide comparison**: identified component areas summed by exact BF identification within each RAW run. Downloadable MS area matrix, figures (PNG/SVG/HTML) and a peptide chart.
- **Multi-sample**: a selected RAW reference compared with **all** other RAW runs, fold/percent changes (undefined if reference=0/missing); heatmap PNG/SVG/HTML and grouped plots.
- **Replicate statistics**: one aggregated peptide area per RAW run; n, mean, sample SD and %RSD grouped by Sample ID or user-assigned Group and preparation. With n=1, SD and %RSD are undefined.
- **Modification analysis**: BF modification labels, exploratory relative MS-area percentage for a selected modification among *unmodified + selected target PTMs of the same peptide backbone*. A missing unmodified reference returns blank, not 100%. This is not site occupancy and is not corrected for ionization efficiency.
- **Unidentified components**: kept separate; no cross-experiment matching by Component No.
- **Downloads**: visible/selected Component Explorer data, each analytical table, and full workbook explicitly as CSV/Excel; charts as PNG/SVG/HTML where available.

For component equivalence between experiments, peptide comparisons use an exact identification key across native PMF and CSV, assuming samples share a compatible protein and BioPharma Finder processing method. **No automated deconvolution, unidentified alignment, automatic response normalization, or modification occupancy interpretation.**

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
