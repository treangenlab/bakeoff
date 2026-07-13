#!/usr/bin/env python3
"""
analysis_prep.py
================
Preprocess raw per-tool taxonomic-profiling reports into the standardized
per-(sample, tool, db) CSVs that analysis_detection.ipynb and
analysis_abundance.ipynb consume. Also builds the unified ground-truth
tables used by both notebooks.

Pipeline position:
    process.sh        -> raw tool outputs
    postprocess.sh    -> kreport / .tre / .sylphmpa
    analysis_prep.py  -> metadata/ground_truth + metadata/preprocessed   (← this script)
    notebooks         -> figures + metrics

Usage (run from bakeoff/):
    python scripts/analysis/analysis_prep.py
    python scripts/analysis/analysis_prep.py --mode detection --jobs 8
    python scripts/analysis/analysis_prep.py \
        --unified-reports reports/unified-reports \
        --default-reports reports/default-reports \
        --data-groups ZymoMockD6331,simulated,DYN \
        --out results

Every run lands under a fresh timestamp dir so reruns never overwrite. The
layout mirrors dyn_prep.py — both scripts write under <--out>/metadata/<ts>/
with their own <script>-prep/ subdir holding logs + outputs:

    <--out>/metadata/<YYYYMMDD_HHMMSS>/analysis-prep/
      ├── analysis_prep.log         (stdout for this run)
      ├── analysis_prep.err         (stderr for this run)
      ├── ground_truth/
      │   ├── <label>_gt.csv        (one per registry entry)
      │   └── ...
      └── preprocessed/
          ├── detection/            (manifest, summary, totals, <Tool>_<db>/...)
          └── abundance/            (same layout)

Notebooks/scripts that need "the most recent run" should glob
<--out>/metadata/*/analysis-prep/ and pick the lexicographically largest
timestamp (matches dyn_prep.find_latest_cohort_cache).

Ground-truth source files are looked up via a registry CSV at
`<--data-root>/ground_truth/registry.csv` with columns:
    label,relpath,sep,name_col,abund_col
Each row points the truth-builder at one source file. Add new datasets by
appending rows; no code changes required.

By default, DYN is **excluded** from preprocessing (no truth → no detection
scoring). Pass `--data-groups ZymoMockD6331,simulated,DYN` to include it.

The two preprocessed/{detection,abundance} subtrees are NOT redundant — they
come from different source files for several tools (e.g. Centrifuge reads
*_kreport.tsv for detection vs *_report.tsv for abundance). One run produces
both side-by-side.
"""

from __future__ import annotations
import argparse
import os
import re
import sys
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

# Import sibling package: scripts/analysis/utils/
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from utils.parser import (                                 # noqa: E402
    build_unified_ground_truth,
    load_kreport_by_rank,
    parse_centrifuge_report_ete3,
    parse_centrifuger_report,
    load_ganon_tre,
    parse_sourmash_kreport_ete3,
    parse_sourmash_report_ete3,
    parse_sylph_mpa_ete3,
)

# CONFIG — paths to ete3 sqlite snapshots (one per (tool, db) combo that needs
# taxid resolution). Subpaths under --db-dir; override via CLI if your layout
# differs.
ETE3_SUBPATHS = {
    # ncbi_id → subpath under --db-dir
    "refseq_032025":     "refseq03032025/ete3_taxa/taxa032025.sqlite",
    "centrifuge_122016": "default_db/cf_default/ete3_taxa/taxa122016.sqlite",
    "sourmash_032022":   "default_db/sm_default/ete3_taxa/taxa032022.sqlite",
    "sylph_042024":      "default_db/sylph_default/taxa042024.sqlite",
    "sylph_gtdb_r220":   "default_db/sylph_default/ete3_taxa/taxaGTDB-r220.sqlite",
    "centrifuger_102023":"default_db/cfer_default/ete3_taxa/taxa102023.sqlite",
}

# Map (tool, db) → which ete3 instance to use. Tools not in this map don't
# need an NCBI handle (their parsers take native taxids).
NCBI_FOR_TOOL_DB = {
    # All unified-DB tools that need ete3 → the unified snapshot
    ("Sourmash",    "unified"): "refseq_032025",
    ("Sylph",       "unified"): "refseq_032025",
    ("Centrifuge",  "unified"): "refseq_032025",
    ("Centrifuger", "unified"): "refseq_032025",
    # Default-DB tools that need ete3 → the per-tool default snapshot
    ("Sourmash",    "default"): "sourmash_032022",
    ("Sylph",       "default"): "sylph_042024",
    ("Centrifuge",  "default"): "centrifuge_122016",
    ("Centrifuger", "default"): "centrifuger_102023",
}

# Ground-truth registry: loaded from <--data-root>/ground_truth/registry.csv.
# CSV columns: label,relpath,sep,name_col,abund_col
# - `relpath`  is resolved against --data-root
# - `sep`      is the field separator of the source file ("," or "\t"; default ",")
# - `name_col` is the column holding the organism name (default "name")
# - `abund_col` is the column holding abundance, or empty if no abundance
GROUND_TRUTH_REGISTRY_RELPATH = "ground_truth/registry.csv"


def load_ground_truth_registry(data_root: Path) -> list[dict]:
    """Return list of dicts (label, relpath, sep, name_col, abund_col).
    Missing optional cells get sensible defaults."""
    import csv as _csv
    path = data_root / GROUND_TRUTH_REGISTRY_RELPATH
    if not path.is_file():
        print(f"[WARN] truth registry not found: {path} — skipping ground-truth build",
              file=sys.stderr)
        return []
    rows = []
    with path.open() as f:
        for r in _csv.DictReader(f):
            sep = (r.get("sep") or ",").strip()
            if sep in (r"\t", "\\t", "tab", "TAB"):
                sep = "\t"
            rows.append({
                "label":     r["label"].strip(),
                "relpath":   r["relpath"].strip(),
                "sep":       sep,
                "name_col":  (r.get("name_col")  or "name").strip(),
                "abund_col": (r.get("abund_col") or "").strip() or None,
            })
    return rows


# Allowed data groups and technologies (filter the report walk).
# DYN is opt-in: it has no truth, so it's excluded from preprocessing by default.
DEFAULT_DATA_GROUPS = ["ZymoMockD6331", "simulated"]
DEFAULT_TECHS       = ["pacbio", "ont"]


COMPOUND_SUFFIXES = [
    ".kreport.tsv", "_kreport.tsv", ".kreport.txt",
    "-report.csv.summarized.csv", "-report.csv",
    ".report.tsv", "_report.tsv", ".report.txt",
    "_reads.tre", "_profile.tsv", "_profiling.tsv",
    ".fastq.sylphmpa",
]
SIMPLE_SUFFIXES = [".tsv", ".csv", ".txt", ".tre", ".sylphmpa"]
TRASH_DASH_SUFFIXES = ("report", "profiling", "profile", "reads", "summary", "kreport")


# TOOL REGISTRIES (one per mode)
# Each entry: (parser_callable, input_glob_patterns, exclude_globs, parser_kwargs, needs_ncbi)
# `rank_labels`: if set, the parser expects 'S'/'G' instead of 'species'/'genus'.
def _build_registries():
    detection = {
        "Kraken2": {
            "parser": load_kreport_by_rank,
            "patterns": ["*_report.tsv"], "exclude": [],
            "kwargs": {"isKraken": True}, "needs_ncbi": False,
            "rank_labels": {"species": "S", "genus": "G"},
            "subdir": "Kraken2-results",
        },
        "Centrifuge": {
            "parser": load_kreport_by_rank,
            "patterns": ["*_kreport.tsv"], "exclude": [],
            "kwargs": {"isKraken": False}, "needs_ncbi": False,
            "rank_labels": {"species": "S", "genus": "G"},
            "subdir": "Centrifuge-results",
        },
        "Centrifuger": {
            "parser": load_kreport_by_rank,
            "patterns": ["*_kreport.tsv"], "exclude": [],
            "kwargs": {"isKraken": False}, "needs_ncbi": False,
            "rank_labels": {"species": "S", "genus": "G"},
            "subdir": "Centrifuger-results",
        },
        "Ganon2": {
            "parser": load_ganon_tre,
            "patterns": ["*_reads.tre"], "exclude": [],
            "kwargs": {"reads": True}, "needs_ncbi": False,
            "subdir": "ganon2-results",
        },
        "Sourmash": {
            "parser": parse_sourmash_kreport_ete3,
            "patterns": ["*.kreport.txt"], "exclude": [],
            "kwargs": {}, "needs_ncbi": True,
            "rank_labels": {"species": "S", "genus": "G"},
            "subdir": "sourmash-results",
        },
        "Sylph": {
            "parser": parse_sylph_mpa_ete3,
            "patterns": ["*.sylphmpa"], "exclude": [],
            "kwargs": {"reads": True}, "needs_ncbi": True,
            "subdir": "sylph-results",
        },
    }
    abundance = {
        "Kraken2": {
            "parser": load_kreport_by_rank,
            "patterns": ["*_report.tsv"], "exclude": [],
            "kwargs": {"isKraken": True}, "needs_ncbi": False,
            "rank_labels": {"species": "S", "genus": "G"},
            "subdir": "Kraken2-results",
        },
        "Centrifuge": {
            "parser": parse_centrifuge_report_ete3,
            "patterns": ["*_report.tsv"], "exclude": ["*_kreport.tsv"],
            "kwargs": {}, "needs_ncbi": True,
            "subdir": "Centrifuge-results",
        },
        "Centrifuger": {
            "parser": parse_centrifuger_report,
            "patterns": ["*_report.tsv"], "exclude": ["*_kreport.tsv"],
            "kwargs": {}, "needs_ncbi": True,
            "subdir": "Centrifuger-results",
        },
        "Ganon2": {
            "parser": load_ganon_tre,
            "patterns": ["*.tre"], "exclude": ["*_reads.tre"],
            "kwargs": {"reads": False}, "needs_ncbi": False,
            "subdir": "ganon2-results",
        },
        "Sourmash": {
            "parser": parse_sourmash_report_ete3,
            "patterns": ["*-report.csv.summarized.csv"], "exclude": [],
            "kwargs": {"use_weighted": True}, "needs_ncbi": True,
            "subdir": "sourmash-results",
        },
        "Sylph": {
            "parser": parse_sylph_mpa_ete3,
            "patterns": ["*.sylphmpa"], "exclude": [],
            "kwargs": {"reads": False}, "needs_ncbi": True,
            "subdir": "sylph-results",
        },
    }
    return {"detection": detection, "abundance": abundance}


# Sample-name normalization
def _strip_known_suffixes(name: str) -> str:
    for suf in COMPOUND_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    for suf in SIMPLE_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return name


def _strip_dash_trash(name: str) -> str:
    m = re.match(r"(.+?)-(" + "|".join(TRASH_DASH_SUFFIXES) + r")(?:.*)$", name, re.I)
    if m:
        return m.group(1)
    return name


def extract_sample_name(path: Path) -> str:
    """Sample name = report filename minus known tool/file suffixes. Raw
    dataset prefixes (e.g. D6331_, sim_) are kept so downstream lookups
    can construct <project>_<technology>_<sample> deterministically."""
    name = path.name
    name = _strip_known_suffixes(name)
    name = _strip_dash_trash(name)
    name = name.rstrip("_").strip()
    return name or path.name


# NCBI registry
class NCBIRegistry:
    """Lazy ete3 NCBITaxa loader keyed by ncbi_id (see ETE3_SUBPATHS)."""

    def __init__(self, db_dir: Path):
        self.db_dir = Path(db_dir)
        self._cache: dict[str, object] = {}

    def get(self, ncbi_id: str):
        if ncbi_id in self._cache:
            return self._cache[ncbi_id]
        if ncbi_id not in ETE3_SUBPATHS:
            raise KeyError(f"Unknown ncbi_id '{ncbi_id}'. Known: {list(ETE3_SUBPATHS)}")
        sqlite = self.db_dir / ETE3_SUBPATHS[ncbi_id]
        if not sqlite.is_file():
            raise FileNotFoundError(
                f"ete3 sqlite for '{ncbi_id}' missing at {sqlite}\n"
                f"  (configure via --db-dir or edit ETE3_SUBPATHS in this script.)"
            )
        # Import here so --help works even without ete3 installed.
        from ete3 import NCBITaxa
        ncbi = NCBITaxa(dbfile=str(sqlite))
        self._cache[ncbi_id] = ncbi
        return ncbi

    def for_tool_db(self, tool: str, db_label: str):
        key = (tool, db_label.lower())
        if key not in NCBI_FOR_TOOL_DB:
            return None
        return self.get(NCBI_FOR_TOOL_DB[key])


# Per-file totals (root + unclassified for kraken-style; sum(direct_reads) for
# sourmash). Detection-side analysis_detection.ipynb needs these to normalize
# Kraken/Centrifuge/Centrifuger/Sourmash values to "fraction of total reads",
# which is what the published 0.001% threshold is applied to. Sylph + Ganon2
# already report fractions, so no totals are emitted for them.
def _total_kraken_like(kreport_path: str) -> float:
    """Sum numReads for 'root' and 'unclassified' from a kreport. Column 2."""
    df = pd.read_csv(
        kreport_path, sep="\t", header=None,
        usecols=[1, 5], names=["numReads", "name"], dtype=str,
    )
    df["numReads"] = pd.to_numeric(df["numReads"], errors="coerce").fillna(0)
    names = df["name"].astype(str).str.strip()
    return float(df.loc[names.isin(["root", "unclassified"]), "numReads"].sum())


def _total_sourmash_kreport(kreport_path: str) -> float:
    """Sum of direct_reads column (col 3) of a sourmash kreport.txt."""
    df = pd.read_csv(
        kreport_path, sep="\t", header=None,
        usecols=[2], names=["direct_reads"], dtype=str,
    )
    df["direct_reads"] = pd.to_numeric(df["direct_reads"], errors="coerce").fillna(0)
    return float(df["direct_reads"].sum())


def compute_total_for_normalization(tool: str, raw_path: str,
                                    mode: str = "detection") -> float | None:
    """Return the normalization denominator for this (tool, file), or None if
    the tool/mode doesn't need one.

    Tools that parse a kreport-style file (Kraken2 in either mode, plus
    Centrifuge/Centrifuger/Sourmash in detection mode) get a real total so the
    precise per-taxon value can be derived from raw read/kmer counts. In
    abundance mode, Centrifuge/Centrifuger parse their *native* (non-kreport)
    reports which already carry precise abundances — no total is needed and
    we return None to avoid emitting bogus 0 rows.
    """
    try:
        # Kraken2 only emits its kreport; the total is meaningful in either mode.
        if tool == "Kraken2":
            return _total_kraken_like(raw_path)
        # Other kreport-derived tools: only meaningful in detection mode.
        if mode == "detection":
            if tool in ("Centrifuge", "Centrifuger"):
                return _total_kraken_like(raw_path)
            if tool == "Sourmash":
                return _total_sourmash_kreport(raw_path)
    except Exception as e:
        print(f"[WARN] could not compute totals for {tool} {raw_path}: {e}",
              file=sys.stderr)
    return None


def _apply_precise_value(df: pd.DataFrame, total: float | None) -> pd.DataFrame:
    """Overwrite df['value'] with `count / total * 100` for kreport-derived
    tools. No-op when total is None (non-kreport tools) or the count column
    is absent (Sylph/Ganon2) — those parsers already emit precise values."""
    if total is None or total <= 0:
        return df
    if "numReads" in df.columns:
        col = "numReads"
    elif "clade_reads" in df.columns:
        col = "clade_reads"
    else:
        return df
    num = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    df["value"] = (num / total) * 100.0
    return df


# Manifest building (walk the reports tree)
def walk_reports(root: Path, db_label: str, registry: dict,
                 data_types: list[str], techs: list[str]) -> list[dict]:
    """Find every report file under <root>/<data_type>/<tech>/<Tool-results>/
    that matches the active registry. Returns one dict per (tool, sample) row."""
    rows: list[dict] = []
    if not root.exists():
        print(f"[WARN] reports root missing: {root}", file=sys.stderr)
        return rows

    for dt_dir in sorted(root.iterdir()):
        if not dt_dir.is_dir() or dt_dir.name not in data_types:
            continue
        for tech_dir in sorted(dt_dir.iterdir()):
            if not tech_dir.is_dir() or tech_dir.name not in techs:
                continue
            for tool_name, cfg in registry.items():
                tool_dir = tech_dir / cfg["subdir"]
                if not tool_dir.is_dir():
                    continue

                matched: set[Path] = set()
                for pat in cfg["patterns"]:
                    matched.update(p for p in tool_dir.glob(pat) if p.is_file())
                for ex in cfg["exclude"]:
                    matched.difference_update(tool_dir.glob(ex))

                for f in sorted(matched):
                    rows.append({
                        "tool":        tool_name,
                        "db":          db_label,
                        "tool_db":     f"{tool_name}_{db_label}",
                        "data_type":   dt_dir.name,
                        "project":     dt_dir.name,
                        "technology":  tech_dir.name,
                        "sample":      extract_sample_name(f),
                        "raw_path":    str(f),
                    })
    return rows


# Parsing
def _normalize_rank(df: pd.DataFrame, logical_rank: str) -> pd.DataFrame:
    if df is None or len(df) == 0:
        return df
    if "rank" not in df.columns:
        df = df.copy()
        df["rank"] = logical_rank
        return df
    df = df.copy()
    r = df["rank"].astype(str).str.upper()
    df.loc[r.isin(["S", "SPECIES"]), "rank"] = "species"
    df.loc[r.isin(["G", "GENUS"]),   "rank"] = "genus"
    others = ~r.isin(["S", "SPECIES", "G", "GENUS"])
    df.loc[others, "rank"] = logical_rank
    return df


def parse_one(row: dict, registry: dict, ncbi_reg: NCBIRegistry) -> pd.DataFrame:
    """Parse one manifest row at species + genus and tag with metadata."""
    tool = row["tool"]
    cfg = registry[tool]
    parser = cfg["parser"]
    base_kwargs = dict(cfg["kwargs"])

    if cfg["needs_ncbi"]:
        ncbi = ncbi_reg.for_tool_db(tool, row["db"])
        if ncbi is None:
            raise ValueError(f"No ete3 instance configured for {tool}/{row['db']}")
        base_kwargs["ncbi"] = ncbi

    rank_labels = cfg.get("rank_labels")
    frames = []
    for logical_rank in ("species", "genus"):
        raw_rank = rank_labels[logical_rank] if rank_labels else logical_rank
        df_rank = parser(str(row["raw_path"]), rank=raw_rank, **base_kwargs)
        df_rank = _normalize_rank(df_rank, logical_rank)
        if df_rank is not None and len(df_rank) > 0:
            frames.append(df_rank)

    df_all = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    # Tag with manifest metadata.
    for c in ("tool", "db", "tool_db", "data_type", "project", "technology", "sample"):
        df_all[c] = row[c]

    # Column ordering: metadata first, parser fields after.
    meta_cols = ["tool", "db", "tool_db", "data_type", "project", "technology", "sample"]
    other = [c for c in df_all.columns if c not in meta_cols]
    return df_all[meta_cols + other]


def _parse_one_worker(args):
    """Top-level wrapper for ProcessPoolExecutor."""
    row, registry, db_dir, out_dir, mode = args
    ncbi_reg = NCBIRegistry(db_dir)  # one per worker; ete3 doesn't pickle
    df = parse_one(row, registry, ncbi_reg)
    total = compute_total_for_normalization(row["tool"], row["raw_path"], mode=mode)
    df = _apply_precise_value(df, total)
    out_path = out_dir / f"{row['tool']}_{row['db']}"
    out_path.mkdir(parents=True, exist_ok=True)
    fname = f"{row['project']}_{row['technology']}_{row['sample']}.csv"
    df.to_csv(out_path / fname, index=False)
    return {
        "tool":       row["tool"],
        "db":         row["db"],
        "project":    row["project"],
        "technology": row["technology"],
        "sample":     row["sample"],
        "rows":       len(df),
        "out_path":   str(out_path / fname),
        "total_for_normalization": total,
    }


# Ground truth
def build_truth_tables(data_root: Path, out_dir: Path, ncbi_reg: NCBIRegistry):
    """Build per-label unified truth tables from the registry CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ncbi = ncbi_reg.get("refseq_032025")
    entries = load_ground_truth_registry(data_root)
    if not entries:
        return
    for e in entries:
        src = data_root / e["relpath"]
        if not src.is_file():
            print(f"[WARN] truth source missing for label '{e['label']}': {src}",
                  file=sys.stderr)
            continue
        print(f"[truth] {e['label']}  ←  {src}")
        build_unified_ground_truth(
            dataset_path=src,
            ncbi=ncbi,
            out_dir=out_dir,
            name_col=e["name_col"],
            abund_col=e["abund_col"],
            sep=e["sep"],
            label=e["label"],
        )


# Summary
def write_summary(manifest: pd.DataFrame, results: list[dict], out_csv: Path):
    rows = []
    by_key = {(r["tool"], r["db"], r["project"], r["technology"], r["sample"]): r
              for r in results}
    for _, m in manifest.iterrows():
        key = (m["tool"], m["db"], m["project"], m["technology"], m["sample"])
        r = by_key.get(key, {})
        out_path = Path(r.get("out_path", ""))
        if out_path.is_file():
            df = pd.read_csv(out_path)
            species_mask = df["rank"].astype(str).eq("species") if "rank" in df else pd.Series(dtype=bool)
            genus_mask   = df["rank"].astype(str).eq("genus")   if "rank" in df else pd.Series(dtype=bool)
            n_species = int(species_mask.sum())
            n_genus   = int(genus_mask.sum())
            no_taxid_s = int((species_mask & (pd.to_numeric(df["taxid"], errors="coerce").fillna(0) == 0)).sum()) \
                          if "taxid" in df else 0
            no_taxid_g = int((genus_mask   & (pd.to_numeric(df["taxid"], errors="coerce").fillna(0) == 0)).sum()) \
                          if "taxid" in df else 0
            val_sum_s = float(df.loc[species_mask, "value"].fillna(0).sum()) if "value" in df else float("nan")
            val_sum_g = float(df.loc[genus_mask,   "value"].fillna(0).sum()) if "value" in df else float("nan")
        else:
            n_species = n_genus = no_taxid_s = no_taxid_g = 0
            val_sum_s = val_sum_g = float("nan")

        rows.append({**m.to_dict(),
                     "rows": r.get("rows", 0),
                     "species_rows": n_species, "genus_rows": n_genus,
                     "species_no_taxid": no_taxid_s, "genus_no_taxid": no_taxid_g,
                     "species_value_sum": val_sum_s, "genus_value_sum": val_sum_g,
                     "out_path": r.get("out_path", "")})
    pd.DataFrame(rows).to_csv(out_csv, index=False)


# Orchestrator
def run_mode(mode: str, manifest: pd.DataFrame, registry: dict,
             db_dir: Path, out_root: Path, jobs: int):
    mode_dir = out_root / "preprocessed" / mode
    mode_dir.mkdir(parents=True, exist_ok=True)

    # Save manifest first so it's available even if parsing fails partway.
    manifest.to_csv(mode_dir / "manifest.csv", index=False)
    print(f"[{mode}] manifest: {len(manifest)} files -> {mode_dir / 'manifest.csv'}")

    if manifest.empty:
        print(f"[{mode}] manifest is empty, skipping parse.")
        return

    # Submit work. Each worker initializes its own NCBIRegistry (ete3 is not
    # picklable across processes).
    task_args = [(row.to_dict(), registry, db_dir, mode_dir, mode)
                 for _, row in manifest.iterrows()]
    results: list[dict] = []
    errors: list[tuple[dict, str]] = []

    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as exe:
            futs = {exe.submit(_parse_one_worker, t): t for t in task_args}
            for i, fut in enumerate(as_completed(futs), 1):
                t = futs[fut]
                row = t[0]
                try:
                    r = fut.result()
                    results.append(r)
                    print(f"[{mode}] [{i}/{len(task_args)}] OK   "
                          f"{row['tool']}/{row['db']} {row['project']}/{row['technology']}/{row['sample']} "
                          f"({r['rows']} rows)")
                except Exception as e:
                    errors.append((row, repr(e)))
                    print(f"[{mode}] [{i}/{len(task_args)}] FAIL "
                          f"{row['tool']}/{row['db']} {row['project']}/{row['technology']}/{row['sample']}: {e}",
                          file=sys.stderr)
    else:
        ncbi_reg = NCBIRegistry(db_dir)
        for i, t in enumerate(task_args, 1):
            row = t[0]
            try:
                df = parse_one(row, registry, ncbi_reg)
                total = compute_total_for_normalization(row["tool"], row["raw_path"], mode=mode)
                df = _apply_precise_value(df, total)
                out_path = mode_dir / f"{row['tool']}_{row['db']}"
                out_path.mkdir(parents=True, exist_ok=True)
                fname = f"{row['project']}_{row['technology']}_{row['sample']}.csv"
                df.to_csv(out_path / fname, index=False)
                results.append({"tool": row["tool"], "db": row["db"],
                                "project": row["project"], "technology": row["technology"],
                                "sample": row["sample"], "rows": len(df),
                                "out_path": str(out_path / fname),
                                "total_for_normalization": total})
                print(f"[{mode}] [{i}/{len(task_args)}] OK   "
                      f"{row['tool']}/{row['db']} {row['project']}/{row['technology']}/{row['sample']} "
                      f"({len(df)} rows)")
            except Exception as e:
                errors.append((row, repr(e)))
                print(f"[{mode}] [{i}/{len(task_args)}] FAIL "
                      f"{row['tool']}/{row['db']} {row['project']}/{row['technology']}/{row['sample']}: {e}",
                      file=sys.stderr)

    # Totals file (kraken-style + sourmash only; others are None and we drop them).
    totals_rows = [
        {"tool": r["tool"], "db": r["db"], "tool_db": f"{r['tool']}_{r['db']}",
         "project": r["project"], "technology": r["technology"],
         "sample": r["sample"], "total_for_normalization": r["total_for_normalization"]}
        for r in results if r.get("total_for_normalization") is not None
    ]
    totals_csv = mode_dir / "totals.csv"
    pd.DataFrame(totals_rows).to_csv(totals_csv, index=False)
    print(f"[{mode}] totals  -> {totals_csv} ({len(totals_rows)} rows)")

    # Summary file.
    summary_csv = mode_dir / "summary.csv"
    write_summary(manifest, results, summary_csv)
    print(f"[{mode}] summary -> {summary_csv}")
    print(f"[{mode}] done: {len(results)} ok, {len(errors)} fail")


class _Tee:
    """File-like wrapper that mirrors writes to multiple streams (e.g. terminal + log file)."""
    def __init__(self, *streams):
        self._streams = streams
    def write(self, s):
        for st in self._streams:
            st.write(s)
            try: st.flush()
            except Exception: pass
    def flush(self):
        for st in self._streams:
            try: st.flush()
            except Exception: pass


# Rough per-worker RAM headroom estimate. analysis_prep.py loads an
# NCBIRegistry (ete3 sqlite) per worker (~0.7 GB resident at steady state).
_PER_WORKER_GB_ESTIMATE = 0.7


def _check_resources_safety(jobs: int, label: str = "jobs") -> None:
    """Print a CPU+RAM sanity check for the requested parallelism and warn
    if it looks unsafe. Always inside dry-run; harmless if no warnings fire."""
    cpu_total = os.cpu_count() or 1
    try:
        cpu_avail = len(os.sched_getaffinity(0))
    except AttributeError:
        cpu_avail = cpu_total
    mem_total_gb = mem_avail_gb = None
    try:
        with open("/proc/meminfo") as f:
            mi = {ln.split(":")[0]: ln.split()[1] for ln in f if ":" in ln}
        mem_total_gb = int(mi.get("MemTotal", 0))     / 1024 / 1024
        mem_avail_gb = int(mi.get("MemAvailable", 0)) / 1024 / 1024
    except Exception:
        pass

    print(f"\n[resource check]")
    print(f"  requested {label}: {jobs}")
    print(f"  cpus available: {cpu_avail} (of {cpu_total} total)")
    if mem_total_gb is not None:
        print(f"  memory available: {mem_avail_gb:.0f} GB (of {mem_total_gb:.0f} GB total)")

    warns = []
    if jobs > cpu_avail:
        warns.append(f"{label}={jobs} > available CPUs ({cpu_avail}) — will oversubscribe and may slow the run")
    if mem_avail_gb is not None:
        peak_gb = jobs * _PER_WORKER_GB_ESTIMATE
        if peak_gb > mem_avail_gb * 0.85:
            warns.append(f"estimated peak ~{peak_gb:.0f} GB ({label}={jobs} × ~{_PER_WORKER_GB_ESTIMATE} GB/worker) "
                         f"approaches available memory ({mem_avail_gb:.0f} GB) — OOM risk")
    if jobs == 1 and cpu_avail >= 4:
        warns.append(f"{label}=1 on {cpu_avail} cores is sequential — pass {label} 4-8 to use more cores")
    if warns:
        for w in warns:
            print(f"  ⚠ {w}")
    else:
        print(f"  OK — within available resources")


def _dry_run(args, unified_root: Path, default_root: Path, data_root: Path,
             out_root: Path, data_groups: list[str], techs: list[str],
             modes: list[str]) -> int:
    """Walk reports + truth registry, print the plan, write dry-run_manifest.csv,
    exit without parsing. Output mirrored to terminal AND a dry-run.log.

    Lands at a fixed path under <out>/metadata/dry-run/analysis-prep/ — the
    contents are wiped on each invocation so the dir always reflects the
    most recent dry-run only. This keeps the metadata/ tree clean (no
    timestamped stubs) and the "latest" glob ignores it (find_latest_*
    matches \\d{8}_\\d{6} timestamps, not the literal 'dry-run')."""
    import datetime, shutil
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = out_root / "metadata" / "dry-run" / "analysis-prep"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    # Tee stdout to terminal + dry-run.log for live viewing + later inspection.
    log_path = run_dir / "dry-run.log"
    log_file = open(log_path, "w")
    orig_stdout = sys.stdout
    sys.stdout = _Tee(orig_stdout, log_file)

    try:
        return _dry_run_body(args, unified_root, default_root, data_root,
                             out_root, run_dir, data_groups, techs, modes,
                             ts, log_path)
    finally:
        sys.stdout = orig_stdout
        log_file.close()


def _dry_run_body(args, unified_root, default_root, data_root, out_root, run_dir,
                  data_groups, techs, modes, ts, log_path) -> int:
    _check_resources_safety(args.jobs, label="jobs")
    print()
    print(f"=== analysis_prep DRY-RUN ({ts}) ===")
    print(f"unified-reports : {unified_root}")
    print(f"default-reports : {default_root}")
    print(f"data-root       : {data_root}")
    print(f"out             : {out_root}  -> {run_dir}")
    print(f"data-groups     : {data_groups}")
    print(f"techs           : {techs}")
    print(f"modes           : {modes}")
    print()

    # Truth registry preview
    if not args.skip_truth:
        entries = load_ground_truth_registry(data_root)
        if entries:
            print(f"[truth] {len(entries)} registry entries:")
            for e in entries:
                src = data_root / e["relpath"]
                tag = "[ok]" if src.is_file() else "[MISSING]"
                print(f"  {tag:10s} {e['label']:20s}  {src}")
            print()
        else:
            print("[truth] no registry entries — nothing to build.\n")

    # Per-mode manifest plan
    if args.skip_reports:
        print("[reports] --skip-reports set; nothing to plan.")
        return 0

    registries = _build_registries()
    all_rows: list[dict] = []
    for mode in modes:
        registry = registries[mode]
        manifest_rows = []
        if unified_root.exists():
            manifest_rows.extend(walk_reports(unified_root, "unified", registry, data_groups, techs))
        if default_root.exists():
            manifest_rows.extend(walk_reports(default_root, "default", registry, data_groups, techs))
        for r in manifest_rows:
            r["mode"] = mode
        print(f"[{mode}] plan: {len(manifest_rows)} files")
        for i, row in enumerate(manifest_rows, 1):
            raw = Path(row.get("raw_path", ""))
            tag = "[ok]" if raw.is_file() else "[MISSING]"
            print(f"  [{i:3d}/{len(manifest_rows):3d}] {tag:10s} "
                  f"{row['tool']}/{row['db']:8s} "
                  f"{row['project']}/{row['technology']}/{row['sample']}")
            if tag == "[MISSING]":
                print(f"            raw: {raw}")
        all_rows.extend(manifest_rows)
        print()

    # Write the consolidated manifest for inspection (in the timestamped run dir).
    manifest_csv = run_dir / "dry-run_manifest.csv"
    if all_rows:
        pd.DataFrame(all_rows).to_csv(manifest_csv, index=False)
        print(f"manifest written: {manifest_csv} ({len(all_rows)} rows)")
    else:
        print("manifest: no rows to write.")
    print(f"log written:      {log_path}")
    return 0


def _setup_run_dir_and_logs(out_root: Path) -> tuple[Path, Path, Path]:
    """Create <out_root>/metadata/<timestamp>/analysis-prep/ and redirect
    stdout/stderr to analysis_prep.{log,err} inside it. Layout matches the
    sibling dyn_prep.py:
        <out>/metadata/<ts>/analysis-prep/
        ├── analysis_prep.log
        ├── analysis_prep.err
        ├── ground_truth/
        └── preprocessed/{detection,abundance}/
    Returns (run_dir, log_path, err_path). Prints the new file locations on
    the *original* stderr first so the user can `tail -f` them."""
    import datetime
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = out_root / "metadata" / ts / "analysis-prep"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "analysis_prep.log"
    err_path = run_dir / "analysis_prep.err"

    # Heads-up on the original stderr before we redirect.
    orig_err_fd = os.dup(2)
    with os.fdopen(orig_err_fd, "w", closefd=True) as orig:
        orig.write(f"=== analysis_prep ({ts}) ===\n")
        orig.write(f"log:  {log_path}\n")
        orig.write(f"err:  {err_path}\n")
        orig.write(f"watch: tail -f {log_path}\n\n")
        orig.flush()

    # Redirect FDs 1 and 2; this propagates to child processes via fork.
    log_fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    err_fd = os.open(str(err_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    os.dup2(log_fd, 1)
    os.dup2(err_fd, 2)
    os.close(log_fd); os.close(err_fd)
    # Replace Python-level objects to honor line buffering on the new FDs.
    sys.stdout = os.fdopen(1, "w", buffering=1)
    sys.stderr = os.fdopen(2, "w", buffering=1)

    return run_dir, log_path, err_path


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="analysis_prep",
        description="Preprocess raw per-tool reports into standardized per-(sample, tool, db) CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run from inside bakeoff/. Most defaults assume the standard bakeoff layout.",
    )
    p.add_argument("--unified-reports", default="reports/unified-reports",
                   help="Reports tree from process.sh --db unified. (default: %(default)s)")
    p.add_argument("--default-reports", default="reports/default-reports",
                   help="Reports tree from process.sh --db default. (default: %(default)s)")
    p.add_argument("--data-root", default="data",
                   help="Root containing dataset subdirs + ground_truth/registry.csv. (default: %(default)s)")
    p.add_argument("--db-dir", default="data/ref_db",
                   help="Reference-DB root used by process.sh; ETE3 sqlite snapshots are resolved relative to this. (default: %(default)s)")
    p.add_argument("--out", default="results",
                   help="Parent output dir; a fresh <YYYYMMDD_HHMMSS>/metadata/ subdir is created inside on every run. (default: %(default)s)")
    p.add_argument("--mode", choices=["detection", "abundance", "both"], default="both",
                   help="Which preprocessing mode(s) to run. (default: %(default)s)")
    p.add_argument("--data-groups", default=",".join(DEFAULT_DATA_GROUPS),
                   help="Comma-separated data-group subdirs to scan (e.g. ZymoMockD6331,simulated). "
                        "DYN is opt-in: add it explicitly to include the cohort. (default: %(default)s)")
    p.add_argument("--techs", default=",".join(DEFAULT_TECHS),
                   help=f"Comma-separated technology subdirs to scan. (default: %(default)s)")
    p.add_argument("--jobs", "-j", type=int, default=1,
                   help="Parallel workers (process-pool). (default: %(default)s)")
    p.add_argument("--skip-truth", action="store_true",
                   help="Skip ground-truth table building (faster reruns).")
    p.add_argument("--skip-reports", action="store_true",
                   help="Skip report preprocessing (e.g. only build ground truth).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the planned (mode, tool, db, sample) rows with [ok]/[MISSING] markers; "
                        "write dry-run.log + dry-run_manifest.csv at the fixed path "
                        "<--out>/metadata/dry-run/analysis-prep/ (overwrites on rerun); "
                        "do not parse, do not write any other outputs.")
    args = p.parse_args(argv)

    unified_root = Path(args.unified_reports)
    default_root = Path(args.default_reports)
    data_root    = Path(args.data_root)
    db_dir       = Path(args.db_dir)
    out_root     = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    data_groups = [s.strip() for s in args.data_groups.split(",") if s.strip()]
    techs       = [s.strip() for s in args.techs.split(",")       if s.strip()]
    modes       = ["detection", "abundance"] if args.mode == "both" else [args.mode]

    if args.dry_run:
        return _dry_run(args, unified_root, default_root, data_root, out_root,
                       data_groups, techs, modes)

    # Fresh timestamped run dir + log/err redirect.
    # run_dir = <out>/metadata/<ts>/analysis-prep/  (matches dyn_prep layout)
    run_dir, log_path, err_path = _setup_run_dir_and_logs(out_root)

    print("=== analysis_prep ===")
    print(f"unified-reports : {unified_root}")
    print(f"default-reports : {default_root}")
    print(f"data-root       : {data_root}")
    print(f"db-dir          : {db_dir}")
    print(f"out             : {out_root}  -> {run_dir}")
    print(f"mode            : {args.mode}")
    print(f"data-groups     : {data_groups}")
    print(f"techs           : {techs}")
    print(f"jobs            : {args.jobs}")
    print()

    # 1) Ground truth (uses unified ete3 only; registry-driven).
    if not args.skip_truth:
        ncbi_reg = NCBIRegistry(db_dir)
        build_truth_tables(data_root, run_dir / "ground_truth", ncbi_reg)

    # 2) Per-mode preprocessing.
    if args.skip_reports:
        return 0

    registries = _build_registries()
    for mode in modes:
        registry = registries[mode]
        manifest_rows = []
        if unified_root.exists():
            manifest_rows.extend(walk_reports(unified_root, "unified", registry, data_groups, techs))
        if default_root.exists():
            manifest_rows.extend(walk_reports(default_root, "default", registry, data_groups, techs))
        manifest = pd.DataFrame(manifest_rows)
        run_mode(mode, manifest, registry, db_dir, run_dir, args.jobs)

    return 0


if __name__ == "__main__":
    sys.exit(main())
