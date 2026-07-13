"""
Alpha-diversity primitives — generic, dataset-agnostic.

Functions
---------
shannon_index(p) / simpson_index(p)
    Indices on a probability vector p (entries that sum to 1).

compute_alpha_diversity(long_df, min_abund_percent, *, key, value_col)
    Long-form DataFrame → per-(sample, tool_db, rank) records carrying
    Shannon, Simpson, richness, classified_fraction, and kept_fraction.

The long DataFrame must contain at minimum:
    sample_id_core, tool_db, rank, {key}, {value_col}
and optionally `assay` (cohort label). `value_col` is assumed to be in
percent units (0–100); the threshold is applied AFTER per-sample
renormalization over classified mass.

DYN cohort-cache helpers
------------------------
find_latest_cohort_cache(out_root)
    Locate the most recent <out>/metadata/<ts>/dyn-prep/tables/cohorts/
    written by `dyn_prep.py`.

load_cohort_plus_illumina(cache_dir, cohort)
    Concatenate one long-read cohort and the illumina baseline from the
    cohort-cache TSVs into a single long-form DataFrame (no intersection
    filter — that's applied per-analysis downstream).

No DYN-specific assumptions in compute_alpha_diversity itself; reusable
across mock / simulated / clinical datasets by callers that prepare the
right long-form input.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# DYN cohort-cache schema constants (must match dyn_prep.py)
_DYN_TOOLS    = ["Kraken2", "Centrifuge", "Centrifuger", "sourmash", "sylph", "ganon2"]
_DYN_DB_MODES = ["default", "unified"]
_DYN_RANKS    = ["species", "genus"]
_DYN_BASELINE = "illumina"
_DYN_REQUIRED_MIN_COLUMNS = {"value", "name", "taxid"}


def shannon_index(p) -> float:
    p = np.asarray(p, float)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if p.size else 0.0


def simpson_index(p) -> float:
    p = np.asarray(p, float)
    p = p[p > 0]
    return float(1.0 - np.sum(p ** 2)) if p.size else 0.0


def compute_alpha_diversity(
    long_df: pd.DataFrame,
    min_abund_percent: float,
    *,
    key: str = "name",
    value_col: str = "norm_abund",
    cf_col: str = "classified_fraction",
) -> pd.DataFrame:
    """Per-(sample, tool_db, rank) Shannon/Simpson/richness records.

    Reads the prep-materialized composition (`norm_abund`, already summed to 100
    per sample/rank over classified taxa) and the materialized `classified_fraction`
    (Change 3, "design B") — it does NOT recompute classified_fraction, which would
    collapse to 1.0 off `norm_abund`. Threshold is applied to the composition
    (≥ min_abund_percent %), then Shannon/Simpson are computed on the survivors,
    re-renormalized to 1. Key stays `name` (intra-profile).
    """
    if key not in {"name", "taxid"}:
        raise ValueError(f"key must be 'name' or 'taxid', got {key!r}")
    if value_col not in long_df.columns:
        raise ValueError(f"value_col={value_col!r} not in columns: {long_df.columns.tolist()}")

    group_cols = ["tool_db", "rank"]
    if "assay" in long_df.columns:
        group_cols = ["assay"] + group_cols

    thr = float(min_abund_percent) / 100.0
    records: list[dict] = []

    for group_key, df_sub in long_df.groupby(group_cols):
        if "assay" in long_df.columns:
            assay, tool_db, rank = group_key
        else:
            assay, tool_db, rank = None, group_key[0], group_key[1]

        mat_raw = df_sub.pivot_table(
            index="sample_id_core",
            columns=key,
            values=value_col,
            aggfunc="sum",
            fill_value=0.0,
        ).astype(float)

        classified_mass = mat_raw.sum(axis=1)
        # classified_fraction is materialized in prep (Change 3, design B); read it per
        # sample rather than recompute (Σnorm_abund is 100 → would collapse to 1.0).
        # Fall back to the old derivation for legacy tables that lack the column.
        if cf_col in df_sub.columns:
            classified_fraction = df_sub.groupby("sample_id_core")[cf_col].first().astype(float)
        else:
            classified_fraction = classified_mass / 100.0

        denom = classified_mass.replace(0.0, np.nan)
        mat_comp = mat_raw.div(denom, axis=0).fillna(0.0)

        if "_" in tool_db:
            tool, db_mode = tool_db.split("_", 1)
        else:
            tool, db_mode = tool_db, "NA"

        for sample_id_core, row in mat_comp.iterrows():
            p_full = row.values.astype(float)
            mask = p_full >= thr
            if not mask.any():
                continue
            p_sel = p_full[mask]
            total = p_sel.sum()
            if total <= 0:
                continue
            p_norm = p_sel / total

            rec = {
                "sample_id_core":      sample_id_core,
                "tool":                tool,
                "db_mode":             db_mode,
                "tool_db":             tool_db,
                "rank":                rank,
                "shannon":             shannon_index(p_norm),
                "simpson":             simpson_index(p_norm),
                "richness":            int(mask.sum()),
                "classified_fraction": float(classified_fraction.loc[sample_id_core]),
                "kept_fraction":       float(total),
            }
            if assay is not None:
                rec["assay"] = assay
            records.append(rec)

    out = pd.DataFrame(records)
    front = (["assay"] if "assay" in out.columns else []) + [
        "sample_id_core", "tool", "db_mode", "tool_db", "rank"
    ]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]


# ----------------------------------------------------------------------
# DYN cohort-cache helpers (used by dyn_alpha_div.ipynb / dyn_heatmap.ipynb).
# Kept here rather than in a separate module because every downstream
# alpha-div consumer already imports utils.alpha_div.
# ----------------------------------------------------------------------
def find_latest_cohort_cache(out_root: Path | str) -> Path:
    """Return the most recent <out_root>/metadata/<ts>/dyn-prep/tables/cohorts/
    directory written by dyn_prep.py.  Raise FileNotFoundError if none
    exist.  The "most recent" timestamp is the lexicographically largest
    one, which matches the YYYYMMDD_HHMMSS format dyn_prep.py uses."""
    out_root = Path(out_root)
    candidates = sorted(out_root.glob("metadata/*/dyn-prep/tables/cohorts"))
    if not candidates:
        raise FileNotFoundError(
            f"No cohort cache found under {out_root}/metadata/*/dyn-prep/tables/cohorts/. "
            "Run `python scripts/analysis/dyn_prep.py` to build one."
        )
    return candidates[-1]


def _read_one_cached_dyn(cache_dir: Path, cohort: str, tool: str,
                         db_mode: str, rank: str) -> pd.DataFrame:
    """Load one (cohort, tool, db, rank) cohort-cache TSV and normalize
    its column contract.  Returns an empty DataFrame if the file is
    missing or empty."""
    p = cache_dir / f"{cohort}_{tool}_{db_mode}_{rank}.tsv"
    if not p.exists():
        return pd.DataFrame()

    df = pd.read_csv(p, sep="\t")
    if df.empty:
        return df

    missing = _DYN_REQUIRED_MIN_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"{p.name} missing columns {sorted(missing)}. "
            f"Got: {df.columns.tolist()}"
        )

    if "sample_id_core" not in df.columns:
        if "sample_id" in df.columns:
            df["sample_id_core"] = df["sample_id"].astype(str)
        else:
            raise ValueError(f"{p.name} must contain 'sample_id_core' or 'sample_id'.")

    df = df.copy()
    df["assay"]   = df.get("assay",   cohort)
    df["tool"]    = df.get("tool",    tool)
    df["db_mode"] = df.get("db_mode", db_mode)
    df["tool_db"] = df.get("tool_db", f"{tool}_{db_mode}")
    df["rank"]    = df.get("rank",    rank)

    df["value"] = pd.to_numeric(df["value"], errors="coerce").fillna(0.0).astype(float)
    df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").fillna(0).astype("Int64")
    df["name"]  = df["name"].astype(str)
    if "sample_id" not in df.columns:
        df["sample_id"] = df["sample_id_core"]
    return df


def load_cohort_plus_illumina(cache_dir: Path | str, cohort: str) -> pd.DataFrame:
    """Concatenate one cohort + the illumina baseline across every
    (tool, db_mode, rank) cell of the cohort cache into a single
    long-form DataFrame.  No intersection filter is applied here —
    that is a per-analysis concern handled downstream."""
    cache_dir = Path(cache_dir)
    parts = []
    for source_cohort in (cohort, _DYN_BASELINE):
        for tool in _DYN_TOOLS:
            for db_mode in _DYN_DB_MODES:
                for rank in _DYN_RANKS:
                    df = _read_one_cached_dyn(cache_dir, source_cohort, tool, db_mode, rank)
                    if not df.empty:
                        parts.append(df)
    if not parts:
        raise RuntimeError(
            f"No cached TSVs were loaded from {cache_dir} for cohort={cohort}. "
            "Run `python scripts/analysis/dyn_prep.py` to (re)build the cache."
        )
    return pd.concat(parts, ignore_index=True)
