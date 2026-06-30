#!/usr/bin/env python3
"""Per-read length tables for the Fig 1 read-length boxplots.

NanoStat (see read_stats.py) gives per-dataset summaries — Table 1 / Fig 11 —
but Fig 1 needs the full length distribution per dataset, which a summary can't
provide. This script parses each FASTQ with pysam and writes one read-length
column per dataset; read_analysis.ipynb loads them and draws the boxplots.

Each input is a FASTQ list file (one path per line, '#' comments allowed). The
`group` is taken from the list filename — datasets_mock_pacbio.txt → mock_pacbio
— matching read_stats.py. FASTQs are parsed in parallel (--workers).

Output (default dir: data/read_stats/read-lengths/):
  {group}__{dataset}.csv   one read_length per line, one file per FASTQ

Requires: pysam.
"""

from __future__ import annotations

import argparse
import datetime
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

try:
    from tqdm import tqdm
except ImportError:  # progress bar is optional
    tqdm = None

_TTY = sys.stderr.isatty()

DEFAULT_OUT_DIR = Path("data/read_stats/read-lengths")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


class _Tee:
    """Write to several streams at once (for --log mirroring)."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def setup_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    f = open(log_file, "w", buffering=1)
    sys.stdout = _Tee(sys.__stdout__, f)
    sys.stderr = _Tee(sys.__stderr__, f)


def make_bar(iterable=None, **kw):
    """tqdm bar routed to the real terminal (never the log); off when not a TTY."""
    if tqdm is None or not _TTY:
        return iterable
    kw.setdefault("file", sys.__stderr__)
    return tqdm(iterable, **kw)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("fastq_lists", type=Path, nargs="+",
                    help="One or more FASTQ list files (e.g. datasets_mock_pacbio.txt).")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                    help=f"Output directory (default: {DEFAULT_OUT_DIR}).")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing per-dataset CSVs. Without it, the script "
                         "aborts if any output already exists.")
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel FASTQ workers (default: 1 — sequential).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the FASTQ→CSV plan and exit without parsing.")
    ap.add_argument("--log", type=Path, default=None,
                    help="Mirror stdout/stderr to this file. Default: an auto-named "
                         "log under <out-dir>/logs/. Use --no-log to disable.")
    ap.add_argument("--no-log", action="store_true",
                    help="Do not write a log file.")
    return ap.parse_args()


def log_header(args: argparse.Namespace, n_tasks: int) -> None:
    print(f"=== read_info.py — {datetime.datetime.now().isoformat(timespec='seconds')} ===")
    print(f"cwd:  {Path.cwd()}")
    print(f"argv: {' '.join(sys.argv)}")
    print(f"lists: {', '.join(str(p) for p in args.fastq_lists)}")
    print(f"plan: {n_tasks} FASTQ(s); workers={args.workers}; out={args.out_dir.resolve()}")
    print(flush=True)


def derive_group(list_path: Path) -> str:
    """`datasets_mock_pacbio.txt` → `mock_pacbio`; otherwise the bare stem."""
    return list_path.stem.removeprefix("datasets_")


def dataset_label(path: Path) -> str:
    name = path.name
    for sfx in FASTQ_SUFFIXES:
        if name.endswith(sfx):
            return name.removesuffix(sfx)
    return path.stem


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


def collect_tasks(fastq_lists: list[Path], out_dir: Path) -> list[tuple[Path, Path]]:
    """Resolve list files to (fastq, out_csv) pairs, named {group}__{dataset}.csv."""
    tasks: list[tuple[Path, Path]] = []
    for list_path in fastq_lists:
        group = derive_group(list_path)
        for fq in read_fastq_list(list_path):
            out_csv = out_dir / f"{group}__{dataset_label(fq)}.csv"
            tasks.append((fq, out_csv))
    return tasks


def extract_lengths(fastq: Path, out_csv: Path) -> int:
    """Write one read_length per line; return the read count."""
    lengths = [len(rec.sequence) for rec in pysam.FastxFile(str(fastq))]
    pd.DataFrame({"read_length": np.asarray(lengths, dtype=np.int64)}).to_csv(
        out_csv, index=False)
    return len(lengths)


def _task(task: tuple[Path, Path]) -> tuple[str, int]:
    fastq, out_csv = task
    return out_csv.stem, extract_lengths(fastq, out_csv)


def main() -> None:
    args = parse_args()
    missing = [p for p in args.fastq_lists if not p.is_file()]
    if missing:
        raise SystemExit("List file(s) not found:\n  " + "\n  ".join(str(p) for p in missing))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tasks = collect_tasks(args.fastq_lists, args.out_dir)

    if args.dry_run:
        print(f"[dry-run] {len(tasks)} FASTQ(s) would be parsed:")
        for fq, out_csv in tasks:
            if out_csv.is_file() and out_csv.stat().st_size > 0:
                note = "OVERWRITE" if args.overwrite else "BLOCKED (exists, --overwrite required)"
            else:
                note = "parse → write"
            print(f"  {fq.name} → {out_csv.name}: {note}")
        return

    existing = [out for _, out in tasks if out.is_file() and out.stat().st_size > 0]
    if existing and not args.overwrite:
        msg = ["Refusing to overwrite existing outputs:"]
        msg.extend(f"  {p.resolve()}" for p in existing)
        msg.append("Re-run with --overwrite to replace them.")
        raise SystemExit("\n".join(msg))

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = None if args.no_log else (args.log or args.out_dir / "logs" / f"read_info_{ts}.log")
    if log_path:
        setup_logging(log_path)
        print(f"[log] {log_path.resolve()}")
    log_header(args, len(tasks))

    bar = make_bar(total=len(tasks), unit="file", desc="parsing")

    def done(name: str, n: int) -> None:
        if bar is not None:
            bar.update(1)
            bar.set_postfix_str(f"{name}: {n:,} reads")
        else:
            print(f"  {name} — {n:,} reads", flush=True)

    if args.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for fut in as_completed([ex.submit(_task, t) for t in tasks]):
                done(*fut.result())
    else:
        for t in tasks:
            done(*_task(t))

    if bar is not None:
        bar.close()
    print(f"Wrote {len(tasks)} length table(s) → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
