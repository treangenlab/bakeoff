"""Prep-layer normalization helpers (Changes 2–3): drop the unclassified token,
and (Change 3) materialize norm_abund / classified_fraction."""
from __future__ import annotations
import numpy as np
import pandas as pd

# Genuine unclassified/unassigned-reads tokens, matched by EXACT name — never by
# taxid<=0 and never as a substring. Real species that only failed name->NCBI-taxid
# mapping carry taxid 0 but keep their place (Sylph_default GTDB names,
# Sourmash_default's Faecalibacterium). Only Sourmash emits a token row
# ('unclassified', taxid 0). 'root' is a name Sourmash/ganon can carry; 'no rank'
# is a rank label, never a name, so it is deliberately excluded.
UNCLASSIFIED_TOKENS = {"unclassified", "unassigned", "unknown", "root"}


def drop_unclassified_token(df: pd.DataFrame, name_col: str = "name") -> pd.DataFrame:
    """Drop the unclassified/unassigned-reads token rows by exact name so they never
    enter the per-rank composition. No-op for tools that emit none."""
    if df is None or df.empty or name_col not in df.columns:
        return df
    mask = df[name_col].astype(str).str.strip().str.lower().isin(UNCLASSIFIED_TOKENS)
    return df[~mask].reset_index(drop=True)


def add_norm_columns(df: pd.DataFrame, group_cols=None) -> pd.DataFrame:
    """Materialize composition + completeness from the (token-dropped) `value` (% of total):
      - `norm_abund`          = value renormalized to sum to 100 per group (composition),
      - `classified_fraction` = Σ(value)/100 per group (the kept fraction).
    `group_cols` groups the (sample, rank) cells; None treats the whole frame as one group
    (a single sample+rank, as in the DYN per-report path). Call AFTER any `value` rewrite."""
    if df is None or df.empty or "value" not in df.columns:
        return df
    df = df.copy()
    v = pd.to_numeric(df["value"], errors="coerce").fillna(0.0)
    if group_cols:
        tot = df.assign(_v=v).groupby(group_cols)["_v"].transform("sum")
    else:
        tot = pd.Series(v.sum(), index=df.index)
    df["norm_abund"] = np.where(tot > 0, v * 100.0 / tot, 0.0)
    df["classified_fraction"] = tot / 100.0
    return df
