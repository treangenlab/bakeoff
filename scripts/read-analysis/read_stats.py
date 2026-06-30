#!/usr/bin/env python3
"""Per-dataset read statistics via NanoStat.

Wraps NanoStat to compute per-dataset read summary statistics for one or more
FASTQ list files. Each input is a text file listing FASTQ paths (one per line,
'#' comments allowed). The script invokes NanoStat once per FASTQ (in parallel
when --workers > 1) and parses its plain-text output into a tabular CSV.

Outputs (parent dirs auto-created):
  - Per-dataset CSV: one row per FASTQ (default: data/read_stats/read_stats.csv)
  - Grouped CSV (optional via --grouped): one row per group
    (default: data/read_stats/grouped.csv)

The `group` column is derived from the input list filename — `datasets_mock_pacbio.txt`
becomes group `mock_pacbio`.

Feeds Table 1 (grouped) and `resource_benchmark.ipynb` / Fig 11 (per-dataset).

Requires: NanoStat (pip install nanostat  or  conda install -c bioconda nanostat).
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

DEFAULT_OUT = Path("data/read_stats/read_stats.csv")
DEFAULT_GROUPED = Path("data/read_stats/grouped.csv")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Q-cutoffs we'll surface as columns (NanoStat reports >Q5,Q7,Q10,Q15,Q20,Q25,Q30)
Q_CUTOFFS = (10, 15, 20, 25, 30)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("fastq_lists", type=Path, nargs="+",
                    help="One or more FASTQ list files (e.g. datasets_mock_pacbio.txt).")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT,
                    help=f"Per-dataset output CSV (default: {DEFAULT_OUT}).")
    ap.add_argument("--grouped", type=Path, nargs="?", default=None, const=DEFAULT_GROUPED,
                    help=f"Aggregated one-row-per-group CSV. Bare --grouped uses {DEFAULT_GROUPED}; "
                         "--grouped PATH customizes; omit to skip.")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing outputs. Without it, the script aborts "
                         "if any output already exists with content.")
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel NanoStat workers (default: 1 — sequential).")
    return ap.parse_args()


# --- helpers ---------------------------------------------------------------

def derive_group(list_path: Path) -> str:
    """`datasets_mock_pacbio.txt` → `mock_pacbio`; otherwise the bare stem."""
    return list_path.stem.removeprefix("datasets_")


def check_inputs_exist(paths: list[Path]) -> None:
    missing = [p for p in paths if not p.is_file()]
    if missing:
        raise SystemExit("Input file(s) not found:\n  " + "\n  ".join(str(p) for p in missing))


def check_overwrite(paths: list[Path], overwrite: bool) -> None:
    """Abort if any of `paths` exists with content, unless --overwrite is set."""
    conflicts = [p for p in paths if p.is_file() and p.stat().st_size > 0]
    if not conflicts:
        return
    if overwrite:
        for p in conflicts:
            print(f"[overwrite] will replace {p}")
        return
    msg = ["Refusing to overwrite existing outputs:"]
    msg.extend(f"  {p.resolve()}" for p in conflicts)
    msg.append("Re-run with --overwrite to replace them.")
    raise SystemExit("\n".join(msg))


def check_nanostat_present() -> None:
    try:
        subprocess.run(["NanoStat", "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise SystemExit(
            f"NanoStat not found ({exc}). Install with:\n"
            "  pip install nanostat        # or\n"
            "  conda install -c bioconda nanostat"
        )


def read_fastq_list(path: Path) -> list[Path]:
    paths: list[Path] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        paths.append(Path(line.split()[0]))
    if not paths:
        raise SystemExit(f"No FASTQ paths in {path}.")
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise SystemExit("Missing FASTQ files:\n  " + "\n  ".join(missing))
    return paths


def dataset_label(p: Path) -> str:
    name = p.name
    for sfx in FASTQ_SUFFIXES:
        if name.endswith(sfx):
            return name.removesuffix(sfx)
    return p.stem


# --- NanoStat invocation + parsing -----------------------------------------

_Q_RE = re.compile(r">Q(\d+):\s*([\d,]+)\s*\(([\d.]+)%\)")
_TOP_RE = re.compile(r"^\s*1:\s*([\d,]+)\s*\(([\d.]+)\)")

# General-section keys NanoStat emits → our column names
_GENERAL_MAP = {
    "Mean read length":   "mean_read_length",
    "Median read length": "median_read_length",
    "STDEV read length":  "stdev_read_length",
    "Read length N50":    "n50",
    "Number of reads":    "num_reads",
    "Total bases":        "total_bases",
    "Mean read quality":  "mean_q",
    "Median read quality": "median_q",
}
_INT_FIELDS = {"num_reads", "total_bases", "n50"}


def _num(s: str) -> float:
    return float(s.replace(",", ""))


def parse_nanostat(text: str) -> dict:
    """Parse NanoStat plain-text output into a flat dict matching our schema."""
    out: dict = {}
    section: str | None = None

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("General summary"):
            section = "general"
            continue
        if "above quality cutoffs" in stripped:
            section = "q_pct"
            continue
        if stripped.startswith("Top") and "longest reads" in stripped:
            section = "longest"
            continue
        if stripped.startswith("Top"):
            section = None  # other "Top N highest …" blocks we don't capture
            continue

        if section == "general" and ":" in stripped:
            key, val = stripped.split(":", 1)
            col = _GENERAL_MAP.get(key.strip())
            if col:
                try:
                    n = _num(val.strip())
                    out[col] = int(n) if col in _INT_FIELDS else n
                except ValueError:
                    pass
        elif section == "q_pct":
            m = _Q_RE.match(stripped)
            if m:
                q = int(m.group(1))
                if q in Q_CUTOFFS:
                    out[f"reads_above_Q{q}_pct"] = float(m.group(3))
        elif section == "longest":
            m = _TOP_RE.match(stripped)
            if m:
                out["longest_read_bp"] = int(_num(m.group(1)))
                out["longest_read_q"] = float(m.group(2))
                section = None  # only need rank 1

    # Fill any Q-cutoffs NanoStat didn't emit (small inputs sometimes omit some)
    for q in Q_CUTOFFS:
        out.setdefault(f"reads_above_Q{q}_pct", float("nan"))
    return out


def summarize(path: Path, group: str) -> dict:
    """Run NanoStat on `path` and return one CSV row."""
    result = subprocess.run(
        ["NanoStat", "--fastq", str(path)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"NanoStat failed for {path}:\n{result.stderr}")

    row: dict = {
        "dataset": dataset_label(path),
        "group": group,
        "file_size_bytes": path.stat().st_size,
        "file_size_gb": path.stat().st_size / (1024 ** 3),
    }
    row.update(parse_nanostat(result.stdout))
    return row


def _summarize_task(task: tuple[Path, str]) -> dict:
    """Module-level wrapper so ProcessPoolExecutor can pickle it."""
    return summarize(*task)


# --- grouped aggregation ---------------------------------------------------

def group_stats(rows: list[dict]) -> list[dict]:
    groups: dict[str, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["group"], []).append(r)

    out: list[dict] = []
    for key in sorted(groups):
        members = groups[key]
        nreads = sum(r["num_reads"] for r in members)
        tbases = sum(r["total_bases"] for r in members)
        fsize  = sum(r["file_size_bytes"] for r in members)

        def wmean(field: str) -> float:
            if nreads == 0:
                return float("nan")
            return sum(r[field] * r["num_reads"] for r in members) / nreads

        row = {
            "group": key,
            "n_datasets": len(members),
            "file_size_bytes_sum": fsize,
            "file_size_gb_sum": fsize / (1024 ** 3),
            "num_reads_sum": nreads,
            "total_bases_sum": tbases,
            "total_bases_gbp_sum": tbases / 1e9,
            "mean_read_length_wmean":   wmean("mean_read_length"),
            "median_read_length_wmean": wmean("median_read_length"),
            "stdev_read_length_wmean":  wmean("stdev_read_length"),
            "n50_max": max(r["n50"] for r in members),
            "mean_q_wmean":   wmean("mean_q"),
            "median_q_wmean": wmean("median_q"),
            "longest_read_bp_max": max(r.get("longest_read_bp", 0) for r in members),
        }
        for q in Q_CUTOFFS:
            row[f"reads_above_Q{q}_pct_wmean"] = wmean(f"reads_above_Q{q}_pct")
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# --- main ------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    check_nanostat_present()
    check_inputs_exist(args.fastq_lists)

    outputs = [args.output] + ([args.grouped] if args.grouped else [])
    check_overwrite(outputs, args.overwrite)

    tasks: list[tuple[Path, str]] = []
    for list_path in args.fastq_lists:
        group = derive_group(list_path)
        for fq in read_fastq_list(list_path):
            tasks.append((fq, group))

    print(f"[plan] {len(tasks)} FASTQ(s); workers={args.workers}", flush=True)
    for i, (fq, _) in enumerate(tasks, 1):
        size_mb = fq.stat().st_size / 1e6
        print(f"  [{i}] {fq.name}  ({size_mb:,.0f} MB)", flush=True)

    rows: list[dict] = []
    t0 = time.time()
    if args.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(_summarize_task, t): t for t in tasks}
            for done_i, fut in enumerate(as_completed(futures), 1):
                row = fut.result()
                rows.append(row)
                print(f"  [done {done_i}/{len(tasks)}, {time.time()-t0:,.0f}s] "
                      f"{row['dataset']} — {row['num_reads']:,} reads, "
                      f"{row['total_bases'] / 1e9:.2f} Gbp",
                      flush=True)
    else:
        for done_i, t in enumerate(tasks, 1):
            row = _summarize_task(t)
            rows.append(row)
            print(f"  [done {done_i}/{len(tasks)}, {time.time()-t0:,.0f}s] "
                  f"{row['dataset']} — {row['num_reads']:,} reads, "
                  f"{row['total_bases'] / 1e9:.2f} Gbp",
                  flush=True)

    write_csv(args.output, rows)
    print(f"Wrote {len(rows)} dataset rows → {args.output.resolve()}")

    if args.grouped:
        grows = group_stats(rows)
        write_csv(args.grouped, grows)
        print(f"Wrote {len(grows)} group rows → {args.grouped.resolve()}")


if __name__ == "__main__":
    main()
