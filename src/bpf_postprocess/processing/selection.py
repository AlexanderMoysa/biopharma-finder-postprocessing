"""Persistent per-run component selection keys; selections never cross experiments."""
from __future__ import annotations
import pandas as pd


def component_keys(frame: pd.DataFrame) -> pd.Series:
    # Compose from immutable source identities, never only from Component No.
    return (frame['experiment'].astype(str) + '\x1f' +
            frame['raw_file'].astype(str) + '\x1f' +
            frame['component_no'].astype('Int64').astype(str))


def selected_rows(frame: pd.DataFrame, excluded: set[str]) -> pd.DataFrame:
    return frame.loc[~component_keys(frame).isin(excluded)].copy()


def update_exclusions(frame: pd.DataFrame, chosen: pd.Series, excluded: set[str]) -> set[str]:
    """Only change selection for visible rows, preserve hidden rows across filters."""
    if len(frame) != len(chosen):
        raise ValueError('Selection length mismatch')
    keys = component_keys(frame).tolist()
    result = set(excluded)
    for key, picked in zip(keys, chosen.tolist()):
        if bool(picked):
            result.discard(key)
        else:
            result.add(key)
    return result
