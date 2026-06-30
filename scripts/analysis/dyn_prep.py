#!/usr/bin/env python
"""
dyn_prep.py — DYN cohort cache builder.

Walks the raw per-tool report trees under --reports-default and
--reports-unified, applies each tool's parser via utils/parser.py, and
writes one long-form TSV per (cohort, tool, db_mode, rank) into a
timestamped run directory:

    <out>/metadata/<timestamp>/dyn-prep/
    ├── dyn_prep.log
    ├── dyn_prep.err
    └── tables/cohorts/
        └── {cohort}_{tool}_{db_mode}_{rank}.tsv   (96 files total)

Downstream notebooks (dyn_alpha_div.ipynb, dyn_heatmap.ipynb) auto-discover
the most recent <timestamp> dir under <out>/metadata/, so any successful
run just becomes the new default cache without any path updates.

Cohort assembly here is intentionally raw and un-intersected: every cohort
(illumina, pacbio, ont_qiagen, ont_zymo) emits ALL its samples. The
intersection filter that pairs each long-read cohort to the subset of
illumina samples with the same biological IDs is applied per-analysis in
the downstream notebooks, so a single cache supports several intersection
rules.

This is the DYN counterpart to analysis_prep.py (mock + simulated
datasets). Both walk reports under the same parser layer in utils/parser;
they differ only in dataset family and output schema.
"""

from __future__ import annotations

import argparse
import datetime
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

# Make utils.* importable regardless of CWD.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from ete3 import NCBITaxa  # noqa: E402

from utils.parser import (  # noqa: E402
    load_kreport_by_rank,
    parse_centrifuge_report_ete3,
    parse_centrifuger_report,
    parse_sourmash_report_ete3,
    parse_sylph_mpa_ete3,
    load_ganon_tre,
)

COHORTS_ALL = ["illumina", "pacbio", "ont_qiagen", "ont_zymo"]
DB_MODES = ["default", "unified"]

# (subdir under DYN/, ONT kit filename token). ont_qiagen/ont_zymo share the
# same ont/ subdir and are split by filename.
COHORT_TO_FOLDER_AND_TAG = {
    "illumina":   ("illumina", None),
    "pacbio":     ("pacbio",   None),
    "ont_qiagen": ("ont",      "Qiagen"),
    "ont_zymo":   ("ont",      "Zymo"),
}

TOOL_CONFIG = {
    "Kraken2":     {"pattern": "{tool}-results/*_report.tsv",                   "parser": "kreport"},
    "Centrifuge":  {"pattern": "{tool}-results/*_report.tsv",                   "parser": "centrifuge"},
    "Centrifuger": {"pattern": "{tool}-results/*_report.tsv",                   "parser": "centrifuger"},
    "sourmash":    {"pattern": "{tool}-results/*-report.csv.summarized.csv",    "parser": "sourmash"},
    "sylph":       {"pattern": "{tool}-results/*.fastq.sylphmpa",               "parser": "sylph"},
    "ganon2":      {"pattern": "{tool}-results/*.tre",                          "parser": "ganon2"},
}

RANKS = {"species": {"kreport_code": "S"}, "genus": {"kreport_code": "G"}}


# Workers each instantiate their own NCBITaxa from these paths — ete3 isn't
# picklable across processes.
def _ete3_paths(db_dir: Path) -> dict[tuple[str, str], str]:
    return {
        ("Centrifuge", "default"): str(db_dir / "default_db" / "cf_default"     / "ete3_taxa" / "taxa122016.sqlite"),
        ("Centrifuge", "unified"): str(db_dir / "refseq03032025"                / "ete3_taxa" / "taxa032025.sqlite"),
        ("sourmash",   "default"): str(db_dir / "default_db" / "sm_default"     / "ete3_taxa" / "taxa032022.sqlite"),
        ("sourmash",   "unified"): str(db_dir / "refseq03032025"                / "ete3_taxa" / "taxa032025.sqlite"),
        ("sylph",      "default"): str(db_dir / "default_db" / "sylph_default"  / "taxa042024.sqlite"),
        ("sylph",      "unified"): str(db_dir / "refseq03032025"                / "ete3_taxa" / "taxa032025.sqlite"),
    }


def _list_sample_stems(root: Path, tool: str, cohort: str) -> set[str]:
    """List file stems for one (root, tool, cohort), with ONT kit filtering."""
    pattern = TOOL_CONFIG[tool]["pattern"].format(tool=tool)
    stems = {p.stem for p in root.glob(pattern)}
    if cohort == "ont_qiagen":
        return {s for s in stems if "_Qiagen" in s}
    if cohort == "ont_zymo":
        return {s for s in stems if "_Zymo" in s}
    return stems


def _strip_tool_postfix(stem: str, tool: str) -> str:
    s = stem
    if tool in {"Kraken2", "Centrifuge", "Centrifuger"}:
        return s[:-len("_report")] if s.endswith("_report") else s
    if tool == "sourmash":
        suffix = "-report.csv.summarized"
        return s[:-len(suffix)] if s.endswith(suffix) else s
    if tool == "sylph":
        return re.sub(r"\.(fastq|fq|fasta|fa)(\.gz)?$", "", s, flags=re.IGNORECASE)
    return s


def _normalize_dataset_id(stem: str, cohort: str, tool: str) -> str:
    s = stem
    if cohort == "ont_qiagen":
        s = s.split("_Qiagen", 1)[0]
    elif cohort == "ont_zymo":
        s = s.split("_Zymo", 1)[0]
    return _strip_tool_postfix(s, tool)


def _infer_sample_id_core(stem: str, ont_kit_tag: str | None) -> str:
    """Strip the kit tag so e.g. DYN_0004_D10_Zymo.fastq.sylphmpa -> DYN_0004_D10."""
    s = str(stem)
    if ont_kit_tag:
        s = s.split(f"_{ont_kit_tag}", 1)[0]
    return s


def _cohort_root(reports_base: Path, cohort: str) -> Path:
    folder, _ = COHORT_TO_FOLDER_AND_TAG[cohort]
    return reports_base / "DYN" / folder


def check_cohort_against_illumina(cohort: str,
                                  base_default: Path, base_unified: Path) -> None:
    """For each tool, log how many cohort samples / illumina samples exist
    under both DBs and warn if any cohort sample lacks an illumina match."""
    print(f"\n=== Illumina correspondence check (cohort={cohort}) ===")
    if cohort == "illumina":
        print(f"  (skipped — cohort is the baseline)")
        return

    cohort_def = _cohort_root(base_default, cohort)
    cohort_uni = _cohort_root(base_unified, cohort)
    ilu_def    = _cohort_root(base_default, "illumina")
    ilu_uni    = _cohort_root(base_unified, "illumina")

    any_missing = False
    counts: dict[str, int] = {}
    for tool in TOOL_CONFIG:
        ts_def = _list_sample_stems(cohort_def, tool, cohort)
        ts_uni = _list_sample_stems(cohort_uni, tool, cohort)
        if not ts_def and not ts_uni:
            continue
        norm_def = {_normalize_dataset_id(s, cohort, tool) for s in ts_def}
        norm_uni = {_normalize_dataset_id(s, cohort, tool) for s in ts_uni}
        i_def = {_normalize_dataset_id(s, "illumina", tool) for s in _list_sample_stems(ilu_def, tool, "illumina")}
        i_uni = {_normalize_dataset_id(s, "illumina", tool) for s in _list_sample_stems(ilu_uni, tool, "illumina")}

        miss_def = sorted(norm_def - i_def)
        miss_uni = sorted(norm_uni - i_uni)
        print(f"\n[{tool}]")
        print(f"  cohort default : {len(norm_def):>3} | illumina default : {len(i_def):>3} | missing : {len(miss_def):>3}")
        if miss_def:
            print("    - examples (default):", ", ".join(miss_def[:10]))
        print(f"  cohort unified : {len(norm_uni):>3} | illumina unified : {len(i_uni):>3} | missing : {len(miss_uni):>3}")
        if miss_uni:
            print("    - examples (unified):", ", ".join(miss_uni[:10]))
        if miss_def or miss_uni:
            any_missing = True
        counts[tool] = len(norm_def | norm_uni)

    if any_missing:
        print(f"\n⚠ Some {cohort} datasets do NOT have a matching Illumina dataset (see above).")
    else:
        print(f"\n✔ {cohort}: every dataset has a matching Illumina dataset.")

    if counts:
        print(f"\nDataset counts per tool ({cohort}):")
        for tool, n in counts.items():
            print(f"  {tool:11s} : {n}")
        all_n = list(counts.values())
        if len(set(all_n)) == 1:
            print(f"\n✔ {cohort} NUM_DATASETS = {all_n[0]}")
        else:
            print(f"\n⚠ {cohort} tools disagree on dataset count: {counts}")


def _build_one_tooldb(cohort: str, db_mode: str, tool: str, rank: str,
                      root: str, out_dir: str,
                      ete3_dbfiles: dict[tuple[str, str], str]) -> tuple[str, str]:
    """Worker: parse one (cohort, tool, db, rank) batch → write its TSV."""
    root = Path(root)
    out_dir = Path(out_dir)

    parser_name = TOOL_CONFIG[tool]["parser"]
    pattern = TOOL_CONFIG[tool]["pattern"].format(tool=tool)
    paths = sorted(root.glob(pattern))

    _, ont_kit_tag = COHORT_TO_FOLDER_AND_TAG[cohort]
    if ont_kit_tag:
        paths = sorted(p for p in paths if f"_{ont_kit_tag}" in p.stem)

    out_path = out_dir / f"{cohort}_{tool}_{db_mode}_{rank}.tsv"
    key = f"{cohort}|{tool}|{db_mode}|{rank}"

    if not paths:
        # Empty file keeps downstream glob() loaders simple.
        pd.DataFrame(columns=[
            "taxid", "name", "value", "sample_id", "sample_id_core",
            "assay", "tool", "db_mode", "tool_db", "rank",
        ]).to_csv(out_path, sep="\t", index=False)
        return key, str(out_path)

    ncbi = None
    if parser_name in {"centrifuge", "sourmash", "sylph"}:
        dbfile = ete3_dbfiles.get((tool, db_mode))
        if not dbfile:
            raise ValueError(f"Missing NCBITaxa dbfile for {(tool, db_mode)}")
        ncbi = NCBITaxa(dbfile=dbfile)

    rank_code = RANKS[rank]["kreport_code"]

    dfs = []
    for p in paths:
        sid = p.stem
        # _normalize_dataset_id strips BOTH the ONT kit tag (e.g. _Zymo)
        # AND the tool-specific postfix (e.g. sourmash's "-report.csv.summarized",
        # kreport-style "_report"), giving the canonical biological sample ID
        # — same normalization the illumina-correspondence check uses, so the
        # downstream intersection filter can match across (tool, db) cells.
        sid_core = _normalize_dataset_id(sid, cohort, tool)

        if parser_name == "kreport":
            df = load_kreport_by_rank(p, rank=rank_code)
        elif parser_name == "centrifuge":
            df = parse_centrifuge_report_ete3(p, rank=rank, ncbi=ncbi)
        elif parser_name == "centrifuger":
            df = parse_centrifuger_report(p, rank=rank)
        elif parser_name == "sourmash":
            df = parse_sourmash_report_ete3(p, rank=rank, ncbi=ncbi)
        elif parser_name == "sylph":
            df = parse_sylph_mpa_ete3(p, rank=rank, reads=False, ncbi=ncbi)
        elif parser_name == "ganon2":
            df = load_ganon_tre(p, rank=rank)
        else:
            raise ValueError(f"Unknown parser '{parser_name}' for tool={tool}")

        df = df.copy()
        needed = {"taxid", "name", "value"}
        missing = needed - set(df.columns)
        if missing:
            raise ValueError(f"{tool}_{db_mode} missing {sorted(missing)} from parser={parser_name} file={p.name}")

        df["sample_id"]      = sid
        df["sample_id_core"] = sid_core
        df["assay"]   = cohort
        df["tool"]    = tool
        df["db_mode"] = db_mode
        df["tool_db"] = f"{tool}_{db_mode}"
        df["rank"]    = rank
        dfs.append(df)

    pd.concat(dfs, ignore_index=True).to_csv(out_path, sep="\t", index=False)
    return key, str(out_path)


def verify_cohort_tsvs(out_dir: Path, cohorts: list[str]) -> None:
    """Quick post-build sanity report: row counts, sample counts per file."""
    print(f"\n=== Post-build verification ({out_dir}) ===")
    files = sorted(out_dir.glob("*.tsv"))
    print(f"  files written: {len(files)}")
    by_cohort: dict[str, int] = {c: 0 for c in cohorts}
    empties: list[str] = []
    for p in files:
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception as e:
            print(f"  [WARN] could not read {p.name}: {e}")
            continue
        if df.empty:
            empties.append(p.name)
            continue
        cohort = p.stem.split("_")[0] if not p.stem.startswith("ont_") else "_".join(p.stem.split("_")[:2])
        if cohort in by_cohort:
            by_cohort[cohort] += 1
    print("  files-per-cohort (non-empty):")
    for c, n in by_cohort.items():
        print(f"    {c:11s} : {n}")
    if empties:
        print(f"  empty files ({len(empties)}):")
        for name in empties[:10]:
            print(f"    - {name}")


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


# Rough per-worker RAM headroom estimate. dyn_prep workers instantiate
# NCBITaxa (ete3 sqlite) at ~0.5 GB resident steady state.
_PER_WORKER_GB_ESTIMATE = 0.5


def _check_resources_safety(threads: int, label: str = "threads") -> None:
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
    print(f"  requested {label}: {threads}")
    print(f"  cpus available: {cpu_avail} (of {cpu_total} total)")
    if mem_total_gb is not None:
        print(f"  memory available: {mem_avail_gb:.0f} GB (of {mem_total_gb:.0f} GB total)")

    warns = []
    if threads > cpu_avail:
        warns.append(f"{label}={threads} > available CPUs ({cpu_avail}) — will oversubscribe and may slow the run")
    if mem_avail_gb is not None:
        peak_gb = threads * _PER_WORKER_GB_ESTIMATE
        if peak_gb > mem_avail_gb * 0.85:
            warns.append(f"estimated peak ~{peak_gb:.0f} GB ({label}={threads} × ~{_PER_WORKER_GB_ESTIMATE} GB/worker) "
                         f"approaches available memory ({mem_avail_gb:.0f} GB) — OOM risk")
    if threads == 1 and cpu_avail >= 4:
        warns.append(f"{label}=1 on {cpu_avail} cores is sequential — pass {label} 4-8 to use more cores")
    if warns:
        for w in warns:
            print(f"  ⚠ {w}")
    else:
        print(f"  OK — within available resources")


# os.dup2 so the ProcessPool workers inherit the redirected FDs.
def _setup_run_dir_and_logs(out_root: Path) -> tuple[Path, Path, Path]:
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = out_root / "metadata" / ts / "dyn-prep"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "dyn_prep.log"
    err_path = run_dir / "dyn_prep.err"

    orig_err_fd = os.dup(2)
    with os.fdopen(orig_err_fd, "w", closefd=True) as orig:
        orig.write(f"=== dyn_prep ({ts}) ===\n")
        orig.write(f"log:  {log_path}\n")
        orig.write(f"err:  {err_path}\n")
        orig.write(f"watch: tail -f {log_path}\n\n")
        orig.flush()

    log_fd = os.open(str(log_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    err_fd = os.open(str(err_path), os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    os.dup2(log_fd, 1)
    os.dup2(err_fd, 2)
    os.close(log_fd)
    os.close(err_fd)
    return run_dir, log_path, err_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reports-default", type=Path,
                   default=Path("reports/default-reports"),
                   help="root of default-DB reports (e.g. /dodo/.../default-reports)")
    p.add_argument("--reports-unified", type=Path,
                   default=Path("reports/unified-reports"),
                   help="root of unified-DB reports")
    p.add_argument("--db-dir", type=Path,
                   default=Path("data/ref_db"),
                   help="reference-DB root (for ETE3 NCBITaxa SQLite files)")
    p.add_argument("--out", type=Path,
                   default=Path("results"),
                   help="output root; cache lands at <out>/metadata/<ts>/dyn-prep/")
    p.add_argument("--cohorts", nargs="+",
                   default=COHORTS_ALL,
                   choices=COHORTS_ALL,
                   help="which cohorts to process (default: all four)")
    p.add_argument("--threads", type=int,
                   default=min(50, (os.cpu_count() or 4)),
                   help="max parallel workers (default: min(50, cpu_count))")
    p.add_argument("--no-log-redirect", action="store_true",
                   help="print to terminal instead of redirecting to log/err files")
    p.add_argument("--dry-run", action="store_true",
                   help="Print the planned (cohort, db, tool, rank) jobs with input-file "
                        "counts; write dry-run.log + dry-run_manifest.csv at the fixed path "
                        "<--out>/metadata/dry-run/dyn-prep/ (overwrites on rerun); "
                        "do not parse, do not write any cohort TSVs.")
    return p.parse_args()


def _dry_run(args) -> int:
    """Walk planned job inputs, print the plan, write dry-run_manifest.csv,
    exit without parsing. Output mirrored to terminal AND a dry-run.log.

    Lands at a fixed path under <out>/metadata/dry-run/dyn-prep/ — the
    contents are wiped on each invocation so the dir always reflects the
    most recent dry-run only. This keeps the metadata/ tree clean (no
    timestamped stubs) and find_latest_cohort_cache ignores it (it globs
    \\d{8}_\\d{6} timestamps, not the literal 'dry-run')."""
    import shutil
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = args.out / "metadata" / "dry-run" / "dyn-prep"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "dry-run.log"
    log_file = open(log_path, "w")
    orig_stdout = sys.stdout
    sys.stdout = _Tee(orig_stdout, log_file)
    try:
        return _dry_run_body(args, run_dir, ts, log_path)
    finally:
        sys.stdout = orig_stdout
        log_file.close()


def _dry_run_body(args, run_dir: Path, ts: str, log_path: Path) -> int:
    _check_resources_safety(args.threads, label="threads")
    print()
    print(f"=== dyn_prep DRY-RUN ({ts}) ===")
    print(f"reports-default : {args.reports_default}")
    print(f"reports-unified : {args.reports_unified}")
    print(f"db-dir          : {args.db_dir}")
    print(f"out             : {args.out}  -> {run_dir}")
    print(f"cohorts         : {args.cohorts}")
    print(f"threads         : {args.threads}")
    print()

    # Input-tree existence
    for tag, p in (("--reports-default", args.reports_default),
                   ("--reports-unified", args.reports_unified)):
        status = "[ok]" if p.exists() else "[MISSING]"
        print(f"  {status:10s} {tag}: {p}")
    print()

    # Build the plan: cohort × db × tool × rank
    jobs = []
    total = len(args.cohorts) * len(DB_MODES) * len(TOOL_CONFIG) * len(RANKS)
    print(f"[plan] {total} jobs ({len(args.cohorts)} cohorts × "
          f"{len(DB_MODES)} dbs × {len(TOOL_CONFIG)} tools × {len(RANKS)} ranks)")
    i = 0
    for cohort in args.cohorts:
        print(f"\n[{cohort}]")
        for db_mode in DB_MODES:
            reports_base = args.reports_default if db_mode == "default" else args.reports_unified
            root = _cohort_root(reports_base, cohort)
            for tool in TOOL_CONFIG:
                pattern = TOOL_CONFIG[tool]["pattern"].format(tool=tool)
                if root.exists():
                    stems = _list_sample_stems(root, tool, cohort)
                    n_files = len(stems)
                else:
                    n_files = 0
                for rank in RANKS:
                    i += 1
                    if not root.exists():
                        status = "[no root]"
                    elif n_files == 0:
                        status = "[empty]"
                    else:
                        status = f"[ok {n_files:3d}f]"
                    out_name = f"{cohort}_{tool}_{db_mode}_{rank}.tsv"
                    print(f"  [{i:3d}/{total:3d}] {db_mode:7s} {tool:11s} {rank:7s} "
                          f"{status:11s} -> {out_name}")
                    jobs.append({
                        "cohort":   cohort,
                        "db_mode":  db_mode,
                        "tool":     tool,
                        "rank":     rank,
                        "n_files":  n_files,
                        "root":     str(root),
                        "pattern":  pattern,
                        "out_name": out_name,
                    })

    # Write the manifest
    manifest_csv = run_dir / "dry-run_manifest.csv"
    pd.DataFrame(jobs).to_csv(manifest_csv, index=False)
    print(f"\nmanifest written: {manifest_csv} ({len(jobs)} rows)")
    print(f"log written:      {log_path}")
    return 0


def main() -> int:
    args = parse_args()

    # Resolve & validate paths
    args.reports_default = args.reports_default.resolve()
    args.reports_unified = args.reports_unified.resolve()
    args.db_dir          = args.db_dir.resolve()
    args.out             = args.out.resolve()

    # ---- dry-run: walk inputs, print plan, exit (terminal output, no redirect) ----
    if args.dry_run:
        return _dry_run(args)

    if not args.no_log_redirect:
        run_dir, log_path, err_path = _setup_run_dir_and_logs(args.out)
    else:
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        run_dir = args.out / "metadata" / ts / "dyn-prep"
        run_dir.mkdir(parents=True, exist_ok=True)
        log_path = err_path = None

    cohort_out = run_dir / "tables" / "cohorts"
    cohort_out.mkdir(parents=True, exist_ok=True)

    print("=== Config ===")
    print(f"reports-default : {args.reports_default}")
    print(f"reports-unified : {args.reports_unified}")
    print(f"db-dir          : {args.db_dir}")
    print(f"out             : {args.out}")
    print(f"run dir         : {run_dir}")
    print(f"cohorts         : {args.cohorts}")
    print(f"threads         : {args.threads}")
    print(f"output cache    : {cohort_out}")

    if not args.reports_default.exists():
        print(f"\n✖ --reports-default does not exist: {args.reports_default}", file=sys.stderr)
        return 2
    if not args.reports_unified.exists():
        print(f"\n✖ --reports-unified does not exist: {args.reports_unified}", file=sys.stderr)
        return 2

    # Sanity check each long-read cohort against the illumina baseline
    for cohort in args.cohorts:
        if cohort == "illumina":
            continue
        check_cohort_against_illumina(cohort, args.reports_default, args.reports_unified)

    # Submit (cohort × db × tool × rank) jobs
    ete3_dbfiles = _ete3_paths(args.db_dir)

    jobs: list[tuple] = []
    for cohort in args.cohorts:
        for db_mode in DB_MODES:
            reports_base = args.reports_default if db_mode == "default" else args.reports_unified
            root = _cohort_root(reports_base, cohort)
            for tool in TOOL_CONFIG:
                for rank in RANKS:
                    jobs.append((cohort, db_mode, tool, rank, str(root), str(cohort_out), ete3_dbfiles))

    print(f"\n=== Building cohort cache ===")
    print(f"submitting {len(jobs)} jobs with max_workers={args.threads}")

    results: dict[str, str] = {}
    errors: list[Exception] = []
    with ProcessPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(_build_one_tooldb, *args_tuple) for args_tuple in jobs]
        for fut in as_completed(futs):
            try:
                k, saved = fut.result()
                results[k] = saved
            except Exception as e:
                errors.append(e)

    print(f"\n✔ wrote {len(results)} tables into {cohort_out}")
    if errors:
        print(f"\n✖ {len(errors)} jobs failed:")
        for e in errors[:10]:
            print("  -", repr(e))
        return 1

    verify_cohort_tsvs(cohort_out, args.cohorts)

    print(f"\n=== Done ===")
    if log_path:
        print(f"log : {log_path}")
        print(f"err : {err_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
