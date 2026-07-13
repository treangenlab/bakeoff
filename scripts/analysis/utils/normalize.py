"""Prep-layer normalization helpers (Changes 2–3): drop the unclassified token,
and (Change 3) materialize norm_abund / classified_fraction."""
from __future__ import annotations
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
