#!/usr/bin/env python3

"""
Build taxonomy metadata TSV for sylph-tax.

Output format (no header, tab-separated):

    genome_fasta_basename    d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species

Column 1 must match the FASTA filenames in your sylph DB (usually the basename of the files you sketched).

Edit the PATH constants below and run:

    python build_sylph_tax_metadata_table.py
"""

import os
import re
from io import StringIO

import pandas as pd
from ete3 import NCBITaxa

# ==========================
# PATHS (EDIT THESE FIRST)
# ==========================

# NCBI assembly summary file, e.g. assembly_summary_refseq.txt
ASSEMBLY_SUMMARY_PATH = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/assembly_summary_refseq.txt"

# Text file with one FASTA path per line (the genomes you sketched with sylph)
GENOME_LIST_PATH = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/file.list"

# Output taxonomy metadata file for sylph-tax (2-column TSV, no header)
OUTPUT_TAXONOMY_PATH = "./sylph_refseq030325_taxonomy_metadata.tsv"

# Local ete3 sqlite built from the same NCBI taxdump used to build the DB.
ETE3_UNIFIED = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/ete3_taxa/taxa032025.sqlite"


# =========================================================
# Helper: load assembly_summary and fix weird header lines
# =========================================================

def load_assembly_summary(path: str) -> pd.DataFrame:
    """
    Load NCBI assembly_summary_* file into a DataFrame.

    Handles:
      - leading '##' comment line(s)
      - header line starting with '#assembly_accession'
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"assembly_summary file not found: {path}")

    with open(path) as f:
        # Drop only lines starting with '##' (true comments),
        # keep the '#assembly_accession' header line.
        kept = [line for line in f if not line.startswith("##")]

    df = pd.read_csv(StringIO("".join(kept)), sep="\t", dtype=str)

    # NCBI header first column is often '#assembly_accession'
    if "#assembly_accession" in df.columns:
        df = df.rename(columns={"#assembly_accession": "assembly_accession"})

    if "assembly_accession" not in df.columns or "taxid" not in df.columns:
        raise ValueError(
            "assembly_summary is missing required columns 'assembly_accession' and/or 'taxid'.\n"
            f"Columns found: {list(df.columns)}"
        )

    return df


# =========================================================
# Helper: extract accession from FASTA filename
# =========================================================

ACC_RE = re.compile(r"^(GCF|GCA)_\d+\.\d+")

def extract_accession_from_fname(fname: str) -> str | None:
    """
    Given a FASTA basename like:
        GCF_002863645.1_ASM286364v1_genomic.fna.gz
    return:
        GCF_002863645.1

    Returns None if no accession can be found.
    """
    m = ACC_RE.match(fname)
    return m.group(0) if m else None


# =========================================================
# Main
# =========================================================

def main():
    print("=== Building sylph taxonomy metadata table ===\n")

    # 1. Load assembly_summary and build accession -> taxid map
    print(f"Loading assembly summary: {ASSEMBLY_SUMMARY_PATH}")
    asm = load_assembly_summary(ASSEMBLY_SUMMARY_PATH)

    acc_to_taxid = dict(zip(asm["assembly_accession"], asm["taxid"]))
    print(f"Loaded {len(asm)} rows from assembly summary.")
    print(f"Unique accessions in mapping: {len(acc_to_taxid)}\n")

    print(f"Initializing ete3.NCBITaxa from {ETE3_UNIFIED}")
    ncbi = NCBITaxa(ETE3_UNIFIED)
    print("NCBITaxa ready.\n")

    wanted_ranks = [
        "superkingdom",  # d__
        "phylum",        # p__
        "class",         # c__
        "order",         # o__
        "family",        # f__
        "genus",         # g__
        "species",       # s__
    ]
    rank_prefix = {
        "superkingdom": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }

    # 3. Process the genome list
    if not os.path.exists(GENOME_LIST_PATH):
        raise FileNotFoundError(f"Genome list file not found: {GENOME_LIST_PATH}")

    print(f"Reading genome list from: {GENOME_LIST_PATH}")
    rows = []
    n_total = 0
    n_used = 0
    n_missing_acc = 0
    n_missing_taxid = 0
    n_taxonomy_errors = 0

    with open(GENOME_LIST_PATH) as f:
        for line in f:
            path = line.strip()
            if not path:
                continue

            n_total += 1
            fname = os.path.basename(path)

            acc = extract_accession_from_fname(fname)
            if acc is None:
                print(f"[WARN] Could not extract accession from filename: {fname}")
                n_missing_acc += 1
                continue

            taxid = acc_to_taxid.get(acc)
            if taxid is None:
                print(f"[WARN] Accession not found in assembly_summary: {acc} (file: {fname})")
                n_missing_taxid += 1
                continue

            try:
                taxid_int = int(taxid)
            except ValueError:
                print(f"[WARN] Non-integer taxid '{taxid}' for accession {acc}")
                n_taxonomy_errors += 1
                continue

            # Build lineage string using ete3
            try:
                lineage = ncbi.get_lineage(taxid_int)
                ranks = ncbi.get_rank(lineage)
                names = ncbi.get_taxid_translator(lineage)

                parts = []
                for tid in lineage:
                    r = ranks.get(tid)
                    if r in wanted_ranks:
                        prefix = rank_prefix[r]
                        name = names.get(tid, "").strip()
                        if name:
                            parts.append(prefix + name)

                lineage_str = ";".join(parts)

                if not lineage_str:
                    print(f"[WARN] Empty lineage for taxid {taxid_int} (accession {acc}, file {fname})")
                    n_taxonomy_errors += 1
                    continue

                rows.append((fname, lineage_str))
                n_used += 1

            except Exception as e:
                print(f"[WARN] Failed to build lineage for taxid {taxid_int} (accession {acc}, file {fname}): {e}")
                n_taxonomy_errors += 1
                continue

    # 4. Write output TSV (no header)
    if not rows:
        raise RuntimeError("No valid genome -> lineage pairs were generated; nothing to write.")

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUTPUT_TAXONOMY_PATH, sep="\t", index=False, header=False)

    print("\n=== Done ===")
    print(f"Total genomes listed:        {n_total}")
    print(f"Genomes with valid lineage:  {n_used}")
    print(f"Missing accession in name:   {n_missing_acc}")
    print(f"Missing accession in map:    {n_missing_taxid}")
    print(f"Taxonomy / lineage errors:   {n_taxonomy_errors}")
    print(f"\nOutput written to: {OUTPUT_TAXONOMY_PATH}")


if __name__ == "__main__":
    main()
