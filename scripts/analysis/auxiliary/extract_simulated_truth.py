#!/usr/bin/env python3
"""
extract_simulated_truth.py — build source ground-truth CSVs for the simulated
PacBio and ONT datasets by counting reads per contig in sim_*.fq.gz and
aggregating to species level.

Each sim_NNNN.fq.gz was simulated by pbsim from one contig at fixed depth,
so read count per contig is proportional to genomic DNA mass. Aggregating
per species gives the relative abundance.

Outputs (one CSV per --techs entry), standard source-GT schema:

    <data-root>/simulated/<tech>/simulated_<tech>_gt.csv
        name,abundance        # one row per species, abundance in %

Usage:
    python extract_simulated_truth.py \\
        --data-root /home/Users/pacbio_bakeoff/data \\
        --joint-fasta /home/Users/yt52/mimic/Mimic/filter_for_pbsim/new_sim/joint_fasta.fasta \\
        --ete3-sqlite /home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/ete3_taxa/taxa032025.sqlite \\
        --species-source /home/Users/pacbio_bakeoff/data/simulated/pacbio/sim_pacbio_gt.tsv \\
        --bakeoff-root /home/Users/pacbio_bakeoff/bakeoff

`--species-source` is a one-column TSV/CSV listing the intended (strain-level)
organism names — it constrains ete3 resolution to the species set the
simulation was designed around. Either sim_pacbio_gt.tsv or sim_ont_gt.tsv
works (they're identical lists).

`--counts-cache` optionally points to a precomputed TSV with rows
`<contig_num>\\t<tech>\\t<reads>` to skip re-reading the 472 fastqs.
"""

from __future__ import annotations
import argparse
import csv
import gzip
import sys
from pathlib import Path

_DESC_STOPWORDS = (" chromosome", " plasmid", " complete", " whole genome",
                   " genome assembly", " contig", " scaffold", ",")


def _desc_organism(desc: str) -> str:
    cut = len(desc)
    for sw in _DESC_STOPWORDS:
        i = desc.find(sw)
        if i != -1 and i < cut:
            cut = i
    return " ".join(desc[:cut].replace("_", " ").split())


def _read_species_source(path: Path) -> list[str]:
    """Return the list of organism names from the source manifest."""
    sep = "\t" if path.suffix == ".tsv" else ","
    names = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            # accept any of these header variants
            for k in ("name", "Name", "Species"):
                if k in row and row[k]:
                    names.append(row[k].strip())
                    break
    if not names:
        raise RuntimeError(f"no names found in {path}")
    return names


def _build_valid_species(names: list[str], ete3_sqlite: Path,
                         resolve_candidates, project_to_rank) -> dict[int, str]:
    """Resolve each input name to a species-rank taxid; return {taxid: name}."""
    from ete3 import NCBITaxa
    ncbi = NCBITaxa(dbfile=str(ete3_sqlite))
    species: dict[int, str] = {}
    for organism in names:
        for cand in resolve_candidates(organism):
            try:
                hit = ncbi.get_name_translator([cand])
            except Exception:
                hit = {}
            for taxid in hit.get(cand, []):
                sp_tid = project_to_rank(int(taxid), "species", ncbi)
                if sp_tid:
                    if sp_tid not in species:
                        species[sp_tid] = ncbi.get_taxid_translator([sp_tid]).get(sp_tid, str(sp_tid))
                    break
            else:
                continue
            break
    return species


def build_accession_to_species(joint_fasta: Path, valid_species: dict[int, str],
                               ete3_sqlite: Path,
                               resolve_candidates, project_to_rank) -> dict[str, str]:
    from ete3 import NCBITaxa
    ncbi = NCBITaxa(dbfile=str(ete3_sqlite))
    cache: dict[str, str | None] = {}

    def resolve(organism: str) -> str | None:
        if organism in cache:
            return cache[organism]
        for cand in resolve_candidates(organism):
            try:
                hit = ncbi.get_name_translator([cand])
            except Exception:
                hit = {}
            for taxid in hit.get(cand, []):
                sp_tid = project_to_rank(int(taxid), "species", ncbi)
                if sp_tid and sp_tid in valid_species:
                    cache[organism] = valid_species[sp_tid]
                    return cache[organism]
        cache[organism] = None
        return None

    acc2species: dict[str, str] = {}
    unmatched: list[str] = []
    with joint_fasta.open() as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            head = ln[1:].rstrip()
            acc, _, desc = head.partition(" ")
            organism = _desc_organism(desc)
            sp = resolve(organism)
            if sp is None:
                unmatched.append(head)
            else:
                acc2species[acc] = sp
    if unmatched:
        raise RuntimeError(
            f"{len(unmatched)} contig headers could not be resolved: "
            + "; ".join(unmatched[:5])
        )
    return acc2species


def accession_from_ref(ref_path: Path) -> str:
    with ref_path.open() as f:
        return f.readline().lstrip(">").split()[0]


def count_reads_fastq_gz(fq_path: Path) -> int:
    n = 0
    with gzip.open(fq_path, "rb") as f:
        for i, _ in enumerate(f):
            if i % 4 == 0:
                n += 1
    return n


def _load_counts_cache(cache_path: Path) -> dict[tuple[str, str], int]:
    """Return {(contig_num, tech): reads}."""
    out = {}
    with cache_path.open() as f:
        for ln in f:
            parts = ln.strip().split("\t")
            if len(parts) >= 3:
                out[(parts[0], parts[1])] = int(parts[2])
    return out


def build_profile(sim_dir: Path, tech: str, acc2species: dict[str, str],
                  counts_cache: dict[tuple[str, str], int] | None
                  ) -> dict[str, int]:
    species_reads: dict[str, int] = {sp: 0 for sp in set(acc2species.values())}
    for ref in sorted((sim_dir / "sim").glob("sim_*.ref")):
        acc = accession_from_ref(ref)
        sp = acc2species.get(acc)
        if sp is None:
            raise RuntimeError(f"Accession {acc} from {ref.name} not in mapping")
        contig_num = ref.stem.split("_")[1]
        if counts_cache is not None and (contig_num, tech) in counts_cache:
            n_reads = counts_cache[(contig_num, tech)]
        else:
            n_reads = count_reads_fastq_gz(ref.with_suffix(".fq.gz"))
        species_reads[sp] += n_reads
    return species_reads


def write_source_gt(species_reads: dict[str, int], out_path: Path) -> None:
    total = sum(species_reads.values())
    # Sort alphabetically by species name for a stable, diff-friendly output
    rows = sorted(species_reads.items(), key=lambda kv: kv[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("name,abundance\n")
        for sp, n in rows:
            pct = 100.0 * n / total if total else 0.0
            f.write(f"{sp},{pct:.6f}\n")
    print(f"  wrote {out_path}  ({len(rows)} species, {total} reads)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", type=Path, required=True,
                   help="data root containing simulated/<tech>/")
    p.add_argument("--joint-fasta", type=Path, required=True,
                   help="joint FASTA containing source contig headers")
    p.add_argument("--ete3-sqlite", type=Path, required=True,
                   help="ete3 NCBITaxa sqlite snapshot")
    p.add_argument("--species-source", type=Path, required=True,
                   help="TSV/CSV listing intended organism names (one 'name' column)")
    p.add_argument("--bakeoff-root", type=Path, required=True,
                   help="bakeoff repo root (for importing scripts/analysis/utils/)")
    p.add_argument("--techs", default="pacbio,ont",
                   help="comma-list of tech subdirs (default: %(default)s)")
    p.add_argument("--counts-cache", type=Path, default=None,
                   help="optional precomputed per-(contig,tech) read counts TSV")
    args = p.parse_args()

    sys.path.insert(0, str(args.bakeoff_root / "scripts/analysis"))
    from utils.parser import _iter_truth_resolution_candidates  # noqa: E402
    from utils.util import project_taxid_to_rank                # noqa: E402

    print("[1/4] resolving source species names to taxids via ete3")
    names = _read_species_source(args.species_source)
    print(f"      {len(names)} input names from {args.species_source.name}")
    valid_species = _build_valid_species(
        names, args.ete3_sqlite, _iter_truth_resolution_candidates, project_taxid_to_rank
    )
    print(f"      resolved to {len(valid_species)} species-rank taxids")

    print("[2/4] mapping joint-FASTA contigs to species")
    acc2species = build_accession_to_species(
        args.joint_fasta, valid_species, args.ete3_sqlite,
        _iter_truth_resolution_candidates, project_taxid_to_rank
    )
    print(f"      {len(acc2species)} contig accessions mapped")

    counts_cache = _load_counts_cache(args.counts_cache) if args.counts_cache else None
    if counts_cache:
        print(f"[3/4] using cached read counts ({len(counts_cache)} entries) from {args.counts_cache}")
    else:
        print("[3/4] counting reads from sim_*.fq.gz (no cache)")

    print("[4/4] aggregating per species and writing source GTs")
    techs = [t.strip() for t in args.techs.split(",") if t.strip()]
    for tech in techs:
        sim_dir = args.data_root / "simulated" / tech
        if not (sim_dir / "sim").exists():
            print(f"  WARN: {sim_dir}/sim missing, skipping {tech}")
            continue
        profile = build_profile(sim_dir, tech, acc2species, counts_cache)
        out_path = sim_dir / f"simulated_{tech}_gt.csv"
        write_source_gt(profile, out_path)


if __name__ == "__main__":
    main()
