# BioPharma Finder Post-processing

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
- **Original BF `.csv`**: supports *multiple* `MS Area:<raw>.raw` columns; other run-specific fields are normalized into individual rows. The CSV currently includes richer `Mod`, `Site`, `Positions`, `Peptide Sequence` annotations than the native PMF decoder.

Load one representation **per experiment**. Uploading the corresponding CSV *and* `.resultsbpf` together is rejected as duplicate data.

A single experiment can contain many RAW runs (for example, reduced and non-reduced injections). The interface initially infers `Reduced` only when explicitly indicated by filename; set the other preparation labels manually in the editable RAW run table.

## Pages / modules (first milestone)

- Component explorer with IDs, annotations, RT, m/z, mass, area, per-run details
- Multi-experiment peptide area comparison and fold ratios; denominator zero → undefined
- CSV modification summary with explicit qualification (NOT site occupancy)
- Unidentified component explorer, kept separate
- Export CSV, Excel workbook, and static PNG/SVG and interactive HTML charts

For component equivalence between experiments, peptide comparisons use an exact identification key or the detailed peptide/modification fields from CSV. **No automated deconvolution, unidentified alignment, automatic response normalization, or modification occupancy interpretation.**

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
