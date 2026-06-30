#!/usr/bin/env python3
"""Per-read alignment error rates from FASTQ reads.

Pipeline (single CLI):
  1. Build a combined reference FASTA. References can be supplied either by
     pointing at a directory of FASTAs (--ref-dir) or by handing the script
     `taxid2genome.py`'s output CSV (--ref-csv), which has one row per taxid
     and a `genome_fasta_path` column. Each contig is renamed
     `{source_stem}__{orig_id}` and recorded in a seq_to_source.tsv companion.
  2. Align each FASTQ with minimap2 (preset map-ont or map-hifi),
     sort + index with samtools, and write a flagstat report. Multiple
     --fastq-list files of different techs can be processed in one run;
     --tech is applied to all, or inferred per list from its filename when
     omitted ('pacbio'/'hifi' → hifi, 'ont' → ont).
  3. Parse the resulting BAM with pysam and emit one per-read CSV per
     dataset with the schema consumed by error_analysis.ipynb (Fig S1).

Output layout under --out-dir (default: results/error-rate/):
  reference/combined_reference.fa[.mmi][.seq_to_source.tsv]
  alignments/{dataset}.sorted.bam[.bai]   alignments/{dataset}.sorted_flagstat.txt
  error-rate-results/per_read_{dataset}.csv

Existing outputs are reused; delete them to force a rebuild.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import pysam

try:
    from tqdm import tqdm
except ImportError:  # progress bar is optional
    tqdm = None

_TTY = sys.stderr.isatty()


def make_bar(iterable=None, **kw):
    """tqdm bar routed to the real terminal (never the log); off when not a TTY."""
    if tqdm is None or not _TTY:
        return iterable
    kw.setdefault("file", sys.__stderr__)
    return tqdm(iterable, **kw)


DEFAULT_OUT_DIR = Path("results/error-rate")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTA_GLOBS = ("*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz")

# pysam CIGAR op codes (BAM spec): M=0 I=1 D=2 S=4 H=5 ==7 X=8 (N/P unused)
CIG_M, CIG_I, CIG_D, CIG_S, CIG_H, CIG_EQ, CIG_X = 0, 1, 2, 4, 5, 7, 8


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--ref-dir", type=Path,
                     help="Directory of per-source FASTA files (.fa/.fasta/.fna, optionally .gz).")
    src.add_argument("--ref-csv", type=Path,
                     help="CSV with a 'genome_fasta_path' column "
                          "(e.g. taxid2genome.py output _best.csv or _perfile.csv).")
    ap.add_argument("--fastq-list", required=True, type=Path, nargs="+", dest="fastq_lists",
                    help="One or more FASTQ list files (one path per line, '#' comments allowed). "
                         "Lists of different techs can be mixed in a single run.")
    ap.add_argument("--tech", choices=["ont", "hifi"], default=None,
                    help="minimap2 preset (map-ont/map-hifi). Applies to all lists. If omitted, "
                         "inferred per list from its filename ('pacbio'/'hifi' → hifi, 'ont' → ont).")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                    help=f"Output directory (default: {DEFAULT_OUT_DIR}).")
    ap.add_argument("--threads", type=int, default=40,
                    help="Threads for minimap2 and samtools alignment (default: 40).")
    ap.add_argument("--parse-workers", type=int, default=1,
                    help="Parse BAMs into per-read CSVs in parallel, one worker per dataset "
                         "(default: 1 — sequential, with a per-read progress bar).")
    ap.add_argument("--bigI", default="32G",
                    help="minimap2 -I batch size for a single-part index (default: 32G).")
    ap.add_argument("--min-mapq", type=int, default=None,
                    help="Optional MAPQ floor applied during BAM parsing.")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing per_read_*.csv outputs. Without it, any dataset "
                         "whose per_read CSV already exists triggers an abort. "
                         "(combined_reference.fa and BAMs are always treated as cache.)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the planned per-dataset actions and exit without doing work.")
    ap.add_argument("--log", type=Path, default=None,
                    help="Mirror stdout/stderr to this file (header records git SHA, CLI args, "
                         "input sizes). Default: an auto-named log under <out-dir>/logs/. "
                         "Use --no-log to disable.")
    ap.add_argument("--no-log", action="store_true",
                    help="Do not write a log file.")
    return ap.parse_args()


class _Tee:
    def __init__(self, *streams): self.streams = streams
    def write(self, s):
        for st in self.streams: st.write(s)
    def flush(self):
        for st in self.streams: st.flush()


def setup_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    f = open(log_file, "w", buffering=1)
    sys.stdout = _Tee(sys.__stdout__, f)
    sys.stderr = _Tee(sys.__stderr__, f)


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def _git_sha() -> str | None:
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "diff", "--quiet"]).returncode != 0
        return f"{sha}{' (dirty)' if dirty else ''}"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def log_header(args: argparse.Namespace, tasks: list[tuple[Path, str]]) -> None:
    print(f"=== error_rates.py — {datetime.datetime.now().isoformat(timespec='seconds')} ===")
    print(f"cwd:  {Path.cwd()}")
    print(f"argv: {' '.join(sys.argv)}")
    sha = _git_sha()
    if sha:
        print(f"git:  {sha}")
    print("inputs:")
    for lp in args.fastq_lists:
        print(f"  --fastq-list {lp}: {_human_size(lp.stat().st_size)}")
    ref = args.ref_dir if args.ref_dir else args.ref_csv
    reflabel = "--ref-dir" if args.ref_dir else "--ref-csv"
    if ref and ref.exists():
        print(f"  {reflabel} {ref}: {_human_size(ref.stat().st_size if ref.is_file() else 0)}")
    print(f"fastqs ({len(tasks)}):")
    for fq, tech in tasks:
        print(f"  [{tech}] {fq}: {_human_size(fq.stat().st_size)}")
    print()


def check_inputs_exist(paths: list[Path]) -> None:
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit("Input path(s) not found:\n  " + "\n  ".join(str(p) for p in missing))


def collect_ref_fastas(ref_dir: Path | None, ref_csv: Path | None) -> list[Path]:
    """Resolve --ref-dir or --ref-csv to a sorted, deduplicated list of FASTA paths."""
    if ref_dir is not None:
        fastas: list[Path] = []
        for pat in FASTA_GLOBS:
            fastas.extend(ref_dir.glob(pat))
        if not fastas:
            raise SystemExit(f"No FASTA files found in {ref_dir}.")
        return sorted(set(fastas))

    df = pd.read_csv(ref_csv)
    if "genome_fasta_path" not in df.columns:
        raise SystemExit(f"--ref-csv {ref_csv} has no 'genome_fasta_path' column.")
    paths = [Path(p) for p in df["genome_fasta_path"].dropna().astype(str).unique()]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        raise SystemExit("FASTA files referenced by --ref-csv not found:\n  "
                         + "\n  ".join(str(p) for p in missing))
    if not paths:
        raise SystemExit(f"--ref-csv {ref_csv} has no usable genome_fasta_path entries.")
    return sorted(set(paths))


def check_dependencies() -> None:
    for tool in ("minimap2", "samtools"):
        try:
            subprocess.run([tool, "--version"], capture_output=True, check=True)
        except Exception as exc:
            raise SystemExit(f"required tool not found: {tool} ({exc})")


def read_fastq_list(path: Path) -> list[Path]:
    paths: list[Path] = []
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        paths.append(Path(line.split()[0]))
    if not paths:
        raise SystemExit(f"No FASTQ paths found in {path}.")
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise SystemExit("Missing FASTQ files:\n  " + "\n  ".join(missing))
    return paths


def dataset_label(path: Path) -> str:
    name = path.name
    for sfx in FASTQ_SUFFIXES:
        if name.endswith(sfx):
            return name[: -len(sfx)]
    return path.stem


def infer_tech(list_path: Path) -> str | None:
    """Guess the minimap2 preset from a list filename; None if ambiguous."""
    name = list_path.name.lower()
    if "hifi" in name or "pacbio" in name:
        return "hifi"
    if "ont" in name or "nanopore" in name:
        return "ont"
    return None


def build_tasks(fastq_lists: list[Path], forced_tech: str | None) -> list[tuple[Path, str]]:
    """Resolve list files to [(fastq, tech), ...], picking tech per list."""
    tasks: list[tuple[Path, str]] = []
    for lp in fastq_lists:
        tech = forced_tech or infer_tech(lp)
        if tech is None:
            raise SystemExit(
                f"Cannot infer --tech from '{lp.name}'. Name the list with 'pacbio'/'hifi' "
                "or 'ont', or pass --tech ont|hifi explicitly.")
        for fq in read_fastq_list(lp):
            tasks.append((fq, tech))
    return tasks


def source_stem(path: Path) -> str:
    return Path(path.name.removesuffix(".gz")).stem


def build_combined_reference(fastas: list[Path], out_fa: Path) -> None:
    map_tsv = out_fa.with_suffix(out_fa.suffix + ".seq_to_source.tsv")
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    n_seq = 0
    with out_fa.open("w") as fout, map_tsv.open("w") as mout:
        mout.write("new_id\torig_id\tsource\n")
        for fa in fastas:
            src = source_stem(fa)
            opener = gzip.open if fa.suffix == ".gz" else open
            with opener(fa, "rt") as fin:
                for line in fin:
                    if line.startswith(">"):
                        orig_id = line[1:].split(None, 1)[0]
                        new_id = f"{src}__{orig_id}"
                        n_seq += 1
                        fout.write(f">{new_id}\n")
                        mout.write(f"{new_id}\t{orig_id}\t{src}\n")
                    else:
                        fout.write(line)
    print(f"[ref] {len(fastas)} FASTA(s), {n_seq} sequences → {out_fa}")


def ensure_mmi_index(fa: Path, bigI: str) -> None:
    mmi = fa.with_suffix(".mmi")
    if not mmi.exists():
        print(f"[ref] indexing → {mmi} (-I {bigI})")
        subprocess.run(["minimap2", "-I", bigI, "-d", str(mmi), str(fa)], check=True)


def align(ref_fa: Path, fastq: Path, bam_out: Path, *, tech: str, threads: int, bigI: str) -> None:
    preset = {"ont": "map-ont", "hifi": "map-hifi"}[tech]
    mm2 = ["minimap2", "-ax", preset, "--eqx", "--MD", "-I", bigI,
           "-t", str(threads), str(ref_fa), str(fastq)]
    sort = ["samtools", "sort", "-@", str(threads), "-o", str(bam_out)]
    print(f"  [align] {' '.join(mm2)}")
    p1 = subprocess.Popen(mm2, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(sort, stdin=p1.stdout)
    p1.stdout.close()  # let p1 receive SIGPIPE if p2 exits early
    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 != 0 or rc2 != 0:
        raise SystemExit(f"alignment failed (minimap2 rc={rc1}, samtools rc={rc2})")

    subprocess.run(["samtools", "index", str(bam_out)], check=True)
    flagstat = bam_out.with_name(bam_out.stem + "_flagstat.txt")
    with flagstat.open("w") as f:
        subprocess.run(["samtools", "flagstat", str(bam_out)], stdout=f, check=True)


def tally_cigar(aln: pysam.AlignedSegment) -> dict:
    """Return op tallies. If '='/'X' are present they define matches/mismatches;
    otherwise we derive mismatches from NM (NM - I - D)."""
    out = dict(M_total=0, I_total=0, D_total=0, S_total=0, H_total=0,
               matches_eq=0, mismatches_X=0)
    if aln.cigartuples is None:
        return out

    has_eqx = False
    for op, ln in aln.cigartuples:
        if op == CIG_M:    out["M_total"] += ln
        elif op == CIG_I:  out["I_total"] += ln
        elif op == CIG_D:  out["D_total"] += ln
        elif op == CIG_S:  out["S_total"] += ln
        elif op == CIG_H:  out["H_total"] += ln
        elif op == CIG_EQ:
            out["matches_eq"] += ln
            has_eqx = True
        elif op == CIG_X:
            out["mismatches_X"] += ln
            has_eqx = True

    if has_eqx:
        out["M_total"] = out["matches_eq"] + out["mismatches_X"]
    else:
        try:
            nm = aln.get_tag("NM")
            out["mismatches_X"] = max(0, nm - out["I_total"] - out["D_total"])
            out["matches_eq"] = max(0, out["M_total"] - out["mismatches_X"])
        except KeyError:
            out["matches_eq"] = out["M_total"]
    return out


def per_read_row(aln: pysam.AlignedSegment, dataset: str,
                 seq2src: dict[str, tuple[str, str]]) -> dict | None:
    if aln.is_unmapped:
        return None
    c = tally_cigar(aln)
    aligned_ref = c["M_total"]                    # ref-side positions (matches + mismatches)
    aligned_read = c["M_total"] + c["I_total"]    # read-side positions
    if aligned_ref == 0 or aligned_read == 0:
        return None
    errors = c["mismatches_X"] + c["I_total"] + c["D_total"]
    ref_name = aln.reference_name or ""
    orig_id, source = seq2src.get(ref_name, ("", ""))
    return {
        "dataset": dataset,
        "read_id": aln.query_name,
        "is_primary": not (aln.is_secondary or aln.is_supplementary),
        "is_secondary": bool(aln.is_secondary),
        "is_supplementary": bool(aln.is_supplementary),
        "mapq": int(aln.mapping_quality),
        "flag": int(aln.flag),
        "ref_name": ref_name,
        "orig_id": orig_id,
        "source": source,
        "read_length": int(aln.query_length or 0),
        "aligned_ref": aligned_ref,
        "aligned_read": aligned_read,
        "matches_eq": c["matches_eq"],
        "mismatches_X": c["mismatches_X"],
        "M_total": c["M_total"],
        "I_total": c["I_total"],
        "D_total": c["D_total"],
        "S_total": c["S_total"],
        "H_total": c["H_total"],
        "errors_total": errors,
        "den_ref": aligned_ref,
        "den_read": aligned_read,
        "error_rate_ref": errors / aligned_ref,
        "error_rate_read": errors / aligned_read,
        "mismatch_rate_ref": c["mismatches_X"] / aligned_ref,
        "mismatch_rate_read": c["mismatches_X"] / aligned_read,
        "ins_rate_ref": c["I_total"] / aligned_ref,
        "ins_rate_read": c["I_total"] / aligned_read,
        "del_rate_ref": c["D_total"] / aligned_ref,
        "del_rate_read": c["D_total"] / aligned_read,
    }


def parse_bam(bam: Path, dataset: str, seq2src: dict[str, tuple[str, str]],
              min_mapq: int | None, progress: bool = False) -> pd.DataFrame:
    rows: list[dict] = []
    with pysam.AlignmentFile(str(bam), "rb") as bam_in:
        alns = bam_in.fetch(until_eof=True)
        if progress:
            total = bam_in.mapped + bam_in.unmapped if bam_in.has_index() else None
            alns = make_bar(alns, total=total, unit="aln", desc=f"parse {dataset}", leave=False)
        for aln in alns:
            if aln.is_secondary or aln.is_supplementary:
                continue  # keep primary alignments only
            if min_mapq is not None and aln.mapping_quality < min_mapq:
                continue
            row = per_read_row(aln, dataset, seq2src)
            if row is not None:
                rows.append(row)
    return pd.DataFrame(rows)


def _parse_task(task: tuple) -> tuple[str, int]:
    """Module-level wrapper so ProcessPoolExecutor can pickle the parse step."""
    bam, dataset, seq2src, min_mapq, per_read_csv = task
    df = parse_bam(bam, dataset, seq2src, min_mapq)
    df.to_csv(per_read_csv, index=False)
    return dataset, len(df)


def load_seq2src(map_tsv: Path) -> dict[str, tuple[str, str]]:
    df = pd.read_csv(map_tsv, sep="\t")
    return {r.new_id: (r.orig_id, r.source) for r in df.itertuples(index=False)}


def main() -> None:
    args = parse_args()
    check_dependencies()

    # Validate input paths upfront so we don't fail deep in the pipeline.
    check_inputs_exist(list(args.fastq_lists) + ([args.ref_dir] if args.ref_dir else [args.ref_csv]))

    ref_dir = args.out_dir / "reference"
    bam_dir = args.out_dir / "alignments"
    err_dir = args.out_dir / "error-rate-results"
    for d in (ref_dir, bam_dir, err_dir):
        d.mkdir(parents=True, exist_ok=True)

    tasks = build_tasks(args.fastq_lists, args.tech)

    # --dry-run: report what each dataset would do and exit.
    if args.dry_run:
        print(f"[dry-run] {len(tasks)} dataset(s) would be processed:")
        for fq, tech in tasks:
            ds = dataset_label(fq)
            per_read = err_dir / f"per_read_{ds}.csv"
            bam = bam_dir / f"{ds}.sorted.bam"
            if per_read.is_file() and per_read.stat().st_size > 0:
                note = "OVERWRITE per_read CSV" if args.overwrite else "BLOCKED (per_read exists, --overwrite required)"
            elif bam.exists():
                note = "parse existing BAM → per_read CSV"
            else:
                note = "align FASTQ → BAM → per_read CSV"
            print(f"  [{tech}] {ds}: {note}")
        return

    # Per-dataset overwrite gate — abort up front if any per_read CSV exists with content.
    conflicts = [err_dir / f"per_read_{dataset_label(fq)}.csv" for fq, _ in tasks]
    conflicts = [p for p in conflicts if p.is_file() and p.stat().st_size > 0]
    if conflicts and not args.overwrite:
        msg = ["Refusing to overwrite existing per-read outputs:"]
        msg.extend(f"  {p.resolve()}" for p in conflicts)
        msg.append("Re-run with --overwrite to replace them, "
                   "or remove the listed datasets from the list file(s) to skip them.")
        raise SystemExit("\n".join(msg))
    if conflicts:
        for p in conflicts:
            print(f"[overwrite] will replace {p}")

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = None if args.no_log else (args.log or args.out_dir / "logs" / f"error_rates_{ts}.log")
    if log_path:
        setup_logging(log_path)
        print(f"[log] {log_path.resolve()}")
    log_header(args, tasks)

    combined_fa = ref_dir / "combined_reference.fa"
    if not combined_fa.exists():
        fastas = collect_ref_fastas(args.ref_dir, args.ref_csv)
        build_combined_reference(fastas, combined_fa)
    else:
        print(f"[ref] reusing {combined_fa}")
    ensure_mmi_index(combined_fa, args.bigI)
    seq2src = load_seq2src(combined_fa.with_suffix(combined_fa.suffix + ".seq_to_source.tsv"))

    print(f"[plan] {len(tasks)} FASTQ(s) across {len(args.fastq_lists)} list(s); "
          f"out={args.out_dir.resolve()}")

    # Phase 1 — alignment (sequential; minimap2/samtools are already multithreaded).
    bams: list[tuple[Path, str, Path]] = []  # (bam, dataset, per_read_csv)
    for i, (fq, tech) in enumerate(tasks, 1):
        ds = dataset_label(fq)
        bam = bam_dir / f"{ds}.sorted.bam"
        print(f"\n[{i}/{len(tasks)}] {ds} ({tech})")
        if not bam.exists():
            align(combined_fa, fq, bam, tech=tech, threads=args.threads, bigI=args.bigI)
        else:
            print(f"  [bam] reusing {bam}")
        bams.append((bam, ds, err_dir / f"per_read_{ds}.csv"))

    # Phase 2 — parse BAMs → per-read CSVs (parallel across datasets if requested).
    print(f"\n[parse] {len(bams)} BAM(s); parse-workers={args.parse_workers}")
    if args.parse_workers > 1 and len(bams) > 1:
        pool = [(bam, ds, seq2src, args.min_mapq, csv) for bam, ds, csv in bams]
        with ProcessPoolExecutor(max_workers=args.parse_workers) as ex:
            futures = make_bar(as_completed([ex.submit(_parse_task, t) for t in pool]),
                               total=len(pool), unit="bam", desc="parse")
            for fut in futures:
                ds, n = fut.result()
                print(f"  [done] {ds}: {n:,} primary alignments")
    else:
        for bam, ds, csv in bams:
            print(f"  [parse] {bam} → {csv}")
            df = parse_bam(bam, ds, seq2src, args.min_mapq, progress=True)
            df.to_csv(csv, index=False)
            print(f"  [done] {len(df):,} primary alignments → {csv.resolve()}")


if __name__ == "__main__":
    main()
