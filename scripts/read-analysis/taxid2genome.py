#!/usr/bin/env python3
"""Map NCBI taxids to genome FASTA files in a RefSeq library.

Emits three CSVs under --out-dir/:
  {prefix}_perfile.csv   one row per FASTA encountered (match status flagged)
  {prefix}_grouped.csv   one row per taxid, with all matched files joined by ';'
  {prefix}_best.csv      one row per taxid, picking the best assembly using
                         RefSeq > category > assembly_level > version > date > size

Optional helper for error_rates.py: the resulting _best.csv / _perfile.csv can be
passed as --ref-csv to build a RefSeq-by-taxid reference. The manuscript's D6331
Fig S1 panel did NOT use this — it aligned against a directory of per-strain
ZymoBIOMICS FASTAs via error_rates.py --ref-dir, so this script is off that path.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import pandas as pd

DEFAULT_OUT_DIR = Path("data/taxid_to_genome")
FA_EXTS = {".fa", ".fasta", ".fna"}
ACC_RE = re.compile(r"(GCF|GCA)_\d{9}\.\d+")

REFSEQ_CATEGORY_RANK = {"reference genome": 3, "representative genome": 2}
ASSEMBLY_LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}
VERSION_STATUS_RANK = {"latest": 2, "replaced": 1}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--assembly-summary", required=True, type=Path,
                    help="NCBI assembly_summary file (tab-delimited).")
    ap.add_argument("--lib-root", required=True, type=Path,
                    help="Library root containing FASTAs (recursively scanned).")
    ap.add_argument("--taxid-list", required=True, type=Path,
                    help="Text file with one NCBI taxid per line ('#' comments allowed).")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                    help=f"Output directory (default: {DEFAULT_OUT_DIR}).")
    ap.add_argument("--prefix", default="taxid_to_genome",
                    help="Output filename prefix (default: taxid_to_genome).")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing output files. Without it, the script aborts "
                         "if any of the three output CSVs already exists with content.")
    return ap.parse_args()


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


def check_inputs_exist(paths: list[Path]) -> None:
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit("Input path(s) not found:\n  " + "\n  ".join(str(p) for p in missing))


def load_taxids(path: Path) -> set[int]:
    out: set[int] = set()
    for raw in path.read_text().splitlines():
        s = raw.split("#", 1)[0].strip()
        if not s:
            continue
        try:
            out.add(int(s))
        except ValueError:
            pass
    if not out:
        raise SystemExit(f"No taxids in {path}.")
    return out


def read_assembly_summary(path: Path) -> pd.DataFrame:
    header: list[str] | None = None
    with path.open() as f:
        for line in f:
            if line.startswith("#assembly_accession"):
                header = line.lstrip("#").strip().split("\t")
                break
    if header is None:
        raise SystemExit("'#assembly_accession' header line not found.")

    df = pd.read_csv(path, sep="\t", comment="#", header=None,
                     names=header, low_memory=False)
    for c in ("taxid", "species_taxid", "genome_size", "genome_size_ungapped",
              "replicon_count", "scaffold_count", "contig_count"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    if "seq_rel_date" in df.columns:
        df["seq_rel_date_parsed"] = pd.to_datetime(
            df["seq_rel_date"], errors="coerce", format="%Y/%m/%d"
        )
    df["is_refseq"] = df["assembly_accession"].astype(str).str.startswith("GCF_")
    return df


def domain_hint(path: Path) -> str:
    p = str(path).lower()
    for d in ("archaea", "bacteria", "fungi", "viral"):
        if f"/{d}/" in p:
            return d
    return ""


def best_rank(row: pd.Series) -> tuple:
    return (
        int(bool(row.get("is_refseq", False))),
        REFSEQ_CATEGORY_RANK.get(str(row.get("refseq_category", "")).lower(), 0),
        ASSEMBLY_LEVEL_RANK.get(str(row.get("assembly_level", "")), 0),
        VERSION_STATUS_RANK.get(str(row.get("version_status", "")).lower(), 0),
        row.get("seq_rel_date_parsed") or pd.NaT,
        float(row.get("genome_size_ungapped") or row.get("genome_size") or 0.0),
        -float(row.get("contig_count") or 1e9),
    )


def scan_library(lib_root: Path, asm_by_acc: dict[str, pd.Series]) -> pd.DataFrame:
    rows: list[dict] = []
    for dirpath, _, files in os.walk(lib_root):
        for fn in files:
            if Path(fn).suffix.lower() not in FA_EXTS:
                continue
            fp = Path(dirpath) / fn
            m = ACC_RE.search(fn)
            acc = m.group(0) if m else None
            base = {
                "genome_fasta_path": str(fp.resolve()),
                "file_basename": fn,
                "domain_hint": domain_hint(fp),
                "assembly_accession": acc or "",
            }
            if acc and acc in asm_by_acc:
                rec = asm_by_acc[acc].to_dict()
                rec.update(base)
                rec["match_status"] = "ok"
                try:
                    rec["file_size_bytes"] = fp.stat().st_size
                except OSError:
                    rec["file_size_bytes"] = None
                rows.append(rec)
            elif acc:
                rows.append({**base, "match_status": "accession_not_in_selected_taxids"})
            else:
                rows.append({**base, "match_status": "no_accession_in_name"})
    return pd.DataFrame(rows)


def join_unique(values) -> str:
    return ";".join(sorted(set(str(v) for v in values if pd.notna(v))))


def main() -> None:
    args = parse_args()
    check_inputs_exist([args.assembly_summary, args.lib_root, args.taxid_list])

    perfile_csv = args.out_dir / f"{args.prefix}_perfile.csv"
    grouped_csv = args.out_dir / f"{args.prefix}_grouped.csv"
    best_csv    = args.out_dir / f"{args.prefix}_best.csv"
    check_overwrite([perfile_csv, grouped_csv, best_csv], args.overwrite)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    taxids = load_taxids(args.taxid_list)
    print(f"[info] {len(taxids):,} taxids loaded")

    asm = read_assembly_summary(args.assembly_summary)
    print(f"[info] {len(asm):,} assembly rows")

    asm_sel = asm[asm["taxid"].isin(taxids)].copy()
    asm_by_acc = {acc: r for acc, r in asm_sel.set_index("assembly_accession").iterrows()}
    print(f"[info] {len(asm_by_acc):,} assemblies match requested taxids")

    perfile = scan_library(args.lib_root, asm_by_acc)
    perfile.to_csv(perfile_csv, index=False)
    print(f"[ok] {perfile_csv.resolve()}")

    matched = perfile[perfile["match_status"] == "ok"]
    if matched.empty:
        print("[warn] no matched files; check lib-root and accessions")
        return

    grouped = (
        matched.groupby("taxid", dropna=False)
        .agg(
            species_taxid=("species_taxid", "first"),
            organism_name=("organism_name", "first"),
            assembly_accessions=("assembly_accession", join_unique),
            genome_fasta_paths=("genome_fasta_path", join_unique),
            domains=("domain_hint", join_unique),
        )
        .reset_index()
    )
    grouped.to_csv(grouped_csv, index=False)
    print(f"[ok] {grouped_csv.resolve()}")

    best_rows: list[pd.Series] = []
    for _taxid, sub in matched.groupby("taxid", dropna=False):
        ranked = sorted(sub.iterrows(), key=lambda kv: best_rank(kv[1]), reverse=True)
        best_rows.append(ranked[0][1])
    best_df = pd.DataFrame(best_rows).reset_index(drop=True)
    keep = [c for c in (
        "taxid", "species_taxid", "organism_name", "assembly_accession", "asm_name",
        "refseq_category", "assembly_level", "version_status", "seq_rel_date",
        "genome_size", "genome_size_ungapped", "gc_percent", "contig_count",
        "is_refseq", "genome_fasta_path", "file_size_bytes", "domain_hint",
    ) if c in best_df.columns]
    best_df = best_df[keep]
    best_df.to_csv(best_csv, index=False)
    print(f"[ok] {best_csv.resolve()}")


if __name__ == "__main__":
    main()
