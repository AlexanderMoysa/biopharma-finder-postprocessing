"""Run-scoped peptide comparisons, preserving annotations and run-specific MS measurements.

A feature is the full BF Identification, not a component number. Run-specific
values (RT, m/z, charge and confidence) can differ between injections, so the
wide table reports the dominant-area component for those values and retains
an independently downloadable table with *every* component/charge state.

Raw area sums/fold ratios are observational MS-signal comparisons, not
chemical concentration, site occupancy or an injection-normalized result.
"""
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from bpf_postprocess.processing.multisample import reference_comparison

RUN_ID_COLS = ("experiment", "raw_file")
ANNOTATION_COLS = ("peptide_sequence", "modification", "site", "positions", "protein")
RUN_MEASURE_COLS = ("ms_area", "rt_min", "mz", "charge", "mono_mass_da", "confidence_pct", "delta_ppm")


def run_id(experiment: str, raw_file: str) -> str:
    """Human-readable, unique within a collection of BF experiments."""
    return f"{experiment} / {raw_file}"


def run_catalog(frame: pd.DataFrame) -> pd.DataFrame:
    """The distinct runs in the *currently selected* component subset."""
    cols = ["run", "experiment", "raw_file", "sample", "sample_group", "preparation"]
    if frame.empty:
        return pd.DataFrame(columns=cols)
    runs = frame[["experiment", "raw_file", "sample", "sample_group", "preparation"]].copy()
    runs["run"] = runs["experiment"].astype(str) + " / " + runs["raw_file"].astype(str)
    if runs.duplicated(["run"]).any():
        # Duplicate identical labels are routine because each component has one row;
        # inconsistent prep/sample labels within the run are not.
        conflict = runs.groupby("run", dropna=False)[["sample", "preparation"]].nunique(dropna=False)
        if (conflict > 1).any().any():
            raise ValueError("Inconsistent sample/preparation labels in the same RAW run")
    return runs[cols].drop_duplicates("run").sort_values(["experiment", "raw_file"]).reset_index(drop=True)


def compatible_run_ids(catalog: pd.DataFrame, reference: str) -> list[str]:
    """Only same-preparation runs; Unspecified is a separate, unverified group."""
    match = catalog.loc[catalog["run"] == reference]
    if len(match) != 1:
        raise ValueError(f"Reference RAW run not found: {reference}")
    preparation = match.iloc[0]["preparation"]
    return catalog.loc[(catalog["preparation"] == preparation) & (catalog["run"] != reference), "run"].tolist()


def select_comparison_runs(frame: pd.DataFrame, reference: str, comparisons: Sequence[str]) -> pd.DataFrame:
    """Select only an explicit reference and comparator list, blocking mixed preps.

    Run selection is independent from the global filters and component checkboxes,
    so excluded components are not reintroduced into a comparison.
    """
    catalog = run_catalog(frame)
    names = set(catalog["run"].tolist())
    if reference not in names:
        raise ValueError(f"Reference RAW run not available: {reference}")
    comparisons = list(comparisons)
    if reference in comparisons:
        raise ValueError("A reference cannot also be selected as its own comparator")
    if len(comparisons) != len(set(comparisons)):
        raise ValueError("Duplicate comparison RAW run")
    unknown = set(comparisons) - names
    if unknown:
        raise ValueError(f"Unknown comparison RAW run(s): {sorted(unknown)}")
    compatible = set(compatible_run_ids(catalog, reference))
    different_preparation = set(comparisons) - compatible
    if different_preparation:
        raise ValueError("Reduced and non-reduced RAW runs cannot be combined in one peptide comparison: "
                         + ", ".join(sorted(different_preparation)))
    wanted = {reference, *comparisons}
    labels = frame["experiment"].astype(str) + " / " + frame["raw_file"].astype(str)
    return frame.loc[labels.isin(wanted)].copy()


def _combine_annotations(values: pd.Series) -> str:
    """Show every distinct nonblank annotation, never silently pick one."""
    seen = []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s and s not in seen:
            seen.append(s)
    return "; ".join(seen)


def _identified(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.assign(feature=pd.Series(dtype="str"), run=pd.Series(dtype="str"))
    df = frame.loc[frame["identified"]].copy()
    df["feature"] = df["identification"].fillna("").astype(str).str.strip()
    df = df.loc[df["feature"].ne("")].copy()
    df["run"] = df["experiment"].astype(str) + " / " + df["raw_file"].astype(str)
    return df


def _feature_metadata(identified: pd.DataFrame) -> pd.DataFrame:
    cols = ["feature", *ANNOTATION_COLS]
    if identified.empty:
        return pd.DataFrame(columns=cols)
    return identified.groupby("feature", as_index=False, dropna=False)[list(ANNOTATION_COLS)].agg(_combine_annotations)[cols]


def _run_summaries(identified: pd.DataFrame) -> pd.DataFrame:
    """Sum areas, but use a REAL dominant-area component for other measurements.

    A summed peak's RT/charge/mz is not an independent measurement. These fields
    are shown from the largest-area component in that run and labeled as such.
    The full-component table is available for other peaks/charge states.
    """
    cols = ["feature", "run", "ms_area", "component_count", "rt_min", "mz", "charge", "mono_mass_da", "confidence_pct", "delta_ppm", "component_no"]
    if identified.empty:
        return pd.DataFrame(columns=cols)
    keys = ["feature", "run"]
    count = identified.groupby(keys, as_index=False, dropna=False).agg(
        ms_area=("ms_area", lambda v: v.sum(min_count=1)),
        component_count=("component_no", "size"),
    )
    # NaN areas sort last; stable ordering makes repeated outputs deterministic.
    dominant = (identified.sort_values(["ms_area", "component_no"], ascending=[False, True],
                                      na_position="last", kind="stable")
                .drop_duplicates(keys, keep="first"))
    dominant_cols = ["rt_min", "mz", "charge", "mono_mass_da", "confidence_pct", "delta_ppm", "component_no"]
    result = count.merge(dominant[keys + dominant_cols], on=keys, validate="one_to_one")[cols]
    # BioPharma Finder may carry an RT/mz/charge for a component that has a
    # reported area of zero in this injection. Such values are retained in
    # detailed_component_records(), but cannot be treated as an observed peak
    # in the aggregated comparison table.
    no_detectable_signal = result["ms_area"].isna() | (result["ms_area"] <= 0)
    if no_detectable_signal.any():
        for measure in ("rt_min", "mz", "charge", "mono_mass_da", "confidence_pct", "delta_ppm"):
            result.loc[no_detectable_signal, measure] = float("nan")
    return result


def detailed_reference_comparison(frame: pd.DataFrame, reference: str, comparisons: Sequence[str]) -> pd.DataFrame:
    """One feature × comparator, with annotation and run-specific peak details.

    Fold is only available if both areas are known and reference area > 0.
    Different charge states/components sharing a BF Identification are summed
    for area, while RT/mz/charge/confidence refer to the dominant-area component.
    """
    selected = select_comparison_runs(frame, reference, comparisons)
    identified = _identified(selected)
    base_cols = ["feature", *ANNOTATION_COLS, "reference_run", "comparison_run", "reference_area", "comparison_area", "delta_area", "fold_vs_reference", "percent_change"]
    if not comparisons or identified.empty:
        return pd.DataFrame(columns=base_cols)
    try:
        result = reference_comparison(selected, reference)
    except ValueError as exc:
        if "Reference RAW not present" in str(exc):
            return pd.DataFrame(columns=base_cols)
        raise
    if result.empty:
        return pd.DataFrame(columns=base_cols)
    annotations = _feature_metadata(identified)
    result = result.merge(annotations, on="feature", how="left", validate="many_to_one")
    observations = _run_summaries(identified)
    ref_obs = observations.loc[observations["run"] == reference].drop(columns="run")
    result = result.merge(ref_obs.rename(columns={c: "reference_" + c for c in ref_obs if c != "feature" and c != "ms_area"}).drop(columns="ms_area"),
                          on="feature", how="left", validate="many_to_one")
    comp_obs = observations.drop(columns="ms_area").rename(columns={c: "comparison_" + c for c in observations if c not in ("feature", "run", "ms_area")})
    comp_obs = comp_obs.rename(columns={"run": "comparison_run"})
    result = result.merge(comp_obs, on=["feature", "comparison_run"], how="left", validate="many_to_one")
    front = [c for c in base_cols if c in result.columns]
    other = [c for c in result if c not in front]
    return result[front + other].sort_values(["feature", "comparison_run"]).reset_index(drop=True)


def annotated_area_matrix(frame: pd.DataFrame, reference: str, comparisons: Sequence[str]) -> pd.DataFrame:
    """Wide per-feature table: annotations + MS Area and dominant RT per RAW.

    Only requested, same-preparation RAW files are present. Missing values stay
    NaN, not zero. Detailed per-component measurements remain downloadable.
    """
    selected = select_comparison_runs(frame, reference, comparisons)
    identified = _identified(selected)
    metadata = _feature_metadata(identified)
    if metadata.empty:
        return metadata
    obs = _run_summaries(identified)
    for name in [reference, *comparisons]:
        current = obs.loc[obs["run"] == name, ["feature", "ms_area", "rt_min", "component_count"]].rename(
            columns={"ms_area": f"MS Area | {name}",
                     "rt_min": f"RT dominant (min) | {name}",
                     "component_count": f"Components | {name}"})
        metadata = metadata.merge(current, on="feature", how="left", validate="one_to_one")
    return metadata.sort_values("feature").reset_index(drop=True)


def detailed_component_records(frame: pd.DataFrame, reference: str, comparisons: Sequence[str]) -> pd.DataFrame:
    """Every identified component in chosen RAWs, with all original measurements.

    Unlike the wide summary, no charge state or RT value is collapsed.
    """
    selected = select_comparison_runs(frame, reference, comparisons)
    identified = _identified(selected)
    if identified.empty:
        return identified
    identified["role"] = identified["run"].eq(reference).map({True: "Reference", False: "Comparison"})
    front = ["role", "run", "experiment", "sample", "sample_group", "raw_file", "preparation", "component_no",
             "peptide_sequence", "modification", "site", "positions", "protein", "identification", "confidence_pct",
             "rt_min", "mz", "charge", "mono_mass_da", "delta_ppm", "ms_area"]
    return identified[[c for c in front if c in identified] + [c for c in identified if c not in front]].sort_values(
        ["run", "identification", "component_no"], kind="stable").reset_index(drop=True)


def run_labeled_reference_matrix(frame: pd.DataFrame, reference: str, comparisons: Sequence[str]) -> pd.DataFrame:
    """One row per peptide, with *RAW file names in every numeric header*.

    For 2, 3, or more samples, each comparator gets its own MS Area, RT and
    fold/delta/%-change columns. No different-preparation run is accepted.
    Missing peptide areas/ratios remain NaN instead of being imputed as zero.
    The detailed long-form reference table is kept separately for auditing.
    """
    matrix = annotated_area_matrix(frame, reference, comparisons)
    if matrix.empty:
        return matrix

    base = ["feature", *ANNOTATION_COLS]
    # Mark reference measurements explicitly, including the full RAW filename.
    reference_columns = {
        f"MS Area | {reference}": f"MS Area | {reference} [REFERENCE]",
        f"RT dominant (min) | {reference}": f"RT dominant (min) | {reference} [REFERENCE]",
        f"Components | {reference}": f"Components | {reference} [REFERENCE]",
    }
    matrix = matrix.rename(columns=reference_columns)
    ordering = [*base, *reference_columns.values()]
    detailed = detailed_reference_comparison(frame, reference, comparisons)

    for comparator in comparisons:
        # Area and RT were already built from exactly this RAW, not from an
        # arbitrary 'comparison_area' belonging to another selected sample.
        ordering.extend([
            f"MS Area | {comparator}",
            f"RT dominant (min) | {comparator}",
            f"Components | {comparator}",
        ])
        suffix = f"{comparator} vs {reference}"
        names = {
            "delta_area": f"Delta MS Area | {suffix}",
            "fold_vs_reference": f"Fold | {suffix}",
            "percent_change": f"Change (%) | {suffix}",
        }
        if detailed.empty:
            for name in names.values():
                matrix[name] = float("nan")
        else:
            rows = detailed.loc[
                detailed["comparison_run"] == comparator,
                ["feature", *names],
            ].rename(columns=names)
            matrix = matrix.merge(rows, on="feature", how="left", validate="one_to_one")
        ordering.extend(names.values())

    return matrix[ordering].sort_values("feature").reset_index(drop=True)
