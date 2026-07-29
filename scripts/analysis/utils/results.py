"""
Shared helpers for downstream figure notebooks: resolve the latest
analysis_prep / dyn_prep timestamp, gate save sites behind a SAVE knob,
and extract host / data-group names from the report trees.

Notebooks write under <results>/<ts>/<analysis>/threshold_<X>/{tables,figures}/
(or <results>/resource-benchmark/<host>/<group>/{tables,figures}/ for resource_benchmark).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Union

# `<YYYYMMDD>_<HHMMSS>` (optionally with a `_<tag>` postfix, e.g. `_mock-sim`, `_dyn`) —
# the format analysis_prep.py / dyn_prep.py use when stamping run dirs. The literal `dry-run` subdir is excluded so
# dry-run stubs never resolve as "latest".
_TS_RE = re.compile(r"^\d{8}_\d{6}(?:_[\w.-]+)?$")

# Top-level results subdir holding derived/prepared data (was "metadata", a misnomer).
RESULTS_SUBDIR = "prepared"


# Timestamp discovery
def _list_timestamps_with_suffix(results_root: Path, suffix_subpath: str) -> list[str]:
    """Return all `<ts>` strings under <results_root>/prepared/ whose
    `<ts>/<suffix_subpath>` exists, sorted ascending. The suffix should
    identify a stage's canonical artifact (e.g. ``analysis-prep`` or
    ``dyn-prep/tables/cohorts``)."""
    results_root = Path(results_root)
    parent = results_root / RESULTS_SUBDIR
    if not parent.exists():
        return []
    hits = []
    for child in parent.iterdir():
        if not child.is_dir() or not _TS_RE.match(child.name):
            continue
        if (child / suffix_subpath).exists():
            hits.append(child.name)
    return sorted(hits)


def find_latest_analysis_prep_ts(results_root: Union[Path, str]) -> str:
    """Return the most recent `<ts>` under `<results_root>/prepared/` that
    contains an analysis-prep run. Raises FileNotFoundError if none exist."""
    ts_list = _list_timestamps_with_suffix(Path(results_root), "analysis-prep")
    if not ts_list:
        raise FileNotFoundError(
            f"No analysis_prep run found under {results_root}/prepared/*/analysis-prep. "
            "Run `python scripts/analysis/analysis_prep.py` first."
        )
    return ts_list[-1]


def find_latest_dyn_prep_ts(results_root: Union[Path, str]) -> str:
    """Return the most recent `<ts>` under `<results_root>/prepared/` that
    contains a dyn-prep run (cohort cache). Raises FileNotFoundError if none."""
    ts_list = _list_timestamps_with_suffix(
        Path(results_root), "dyn-prep/tables/cohorts"
    )
    if not ts_list:
        raise FileNotFoundError(
            f"No dyn_prep run found under {results_root}/prepared/*/dyn-prep/tables/cohorts. "
            "Run `python scripts/analysis/dyn_prep.py` first."
        )
    return ts_list[-1]


# Host extraction (used by resource_benchmark.ipynb)
_STATE_LOG_RE = re.compile(r"^run_state_([A-Za-z0-9_.-]+)_\d{8}_\d{6}\.log$")


def extract_hosts_from_reports(report_roots) -> set[str]:
    """Scan each report root for `run_state_<host>_<ts>.log` files written by
    process.sh and return the set of unique host names. Returns an empty set
    if no state logs are found."""
    hosts: set[str] = set()
    for root in report_roots:
        root = Path(root)
        if not root.exists():
            continue
        for f in root.iterdir():
            m = _STATE_LOG_RE.match(f.name)
            if m:
                hosts.add(m.group(1))
    return hosts


def extract_data_groups_from_reports(report_roots) -> set[str]:
    """Return the set of top-level data-group dirs (e.g. ZymoMockD6331,
    simulated, DYN) present under any of the report roots. Hidden dirs
    and dot-files are ignored."""
    groups: set[str] = set()
    for root in report_roots:
        root = Path(root)
        if not root.exists():
            continue
        for f in root.iterdir():
            if f.is_dir() and not f.name.startswith("."):
                groups.add(f.name)
    return groups


# Path builders
def tables_dir(out_dir: Union[Path, str], threshold=None) -> Path:
    """Return `<out_dir>/[threshold_<X>/]tables/`."""
    p = Path(out_dir)
    if threshold is not None:
        p = p / f"threshold_{threshold}"
    return p / "tables"


def figures_dir(out_dir: Union[Path, str], threshold=None) -> Path:
    """Return `<out_dir>/[threshold_<X>/]figures/`."""
    p = Path(out_dir)
    if threshold is not None:
        p = p / f"threshold_{threshold}"
    return p / "figures"


# Save gates
def save_if(save: bool, path: Union[Path, str], writer: Callable[[Path], None]) -> None:
    """Generic gate: call `writer(path)` only when `save` is True. Creates
    the parent dir on demand. Use this for save sites that don't fit the
    `to_csv`/`savefig` pattern (e.g. JSON, Excel)."""
    if not save:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer(path)


def save_csv(df, path: Union[Path, str], save: bool, **to_csv_kwargs) -> None:
    """Gated `df.to_csv(path, **kwargs)`. No-op when `save` is False."""
    if not save:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, **to_csv_kwargs)


def save_fig(fig, path: Union[Path, str], save: bool, **savefig_kwargs) -> None:
    """Gated `fig.savefig(path, **kwargs)`. No-op when `save` is False.
    `plt.show()` (called elsewhere) still fires regardless, so plots
    remain visible inline."""
    if not save:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, **savefig_kwargs)
