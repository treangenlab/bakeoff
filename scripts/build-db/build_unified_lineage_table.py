#!/usr/bin/env python3
"""
Build lineage table for unified database assemblies.

Inputs:
  1) --file-list: text file with one genome file path per line.
     Example path:
       /home/Users/.../GCF_045161915.1_ASM4516191v1_genomic.fa

  2) --assembly-summary: one or more NCBI assembly_summary*.txt files.
     Must contain at least:
       assembly_accession
       taxid
       species_taxid

  3) --taxdump: directory containing nodes.dmp and names.dmp.

Output:
  - TSV with columns:
        assembly_accession
        assembly_taxid
        species_taxid
        genus_taxid
        phylum
        species_name

  - Log file with details about unmapped accessions, taxids, and lineage issues.
"""

import argparse
import csv
import os
import re
import sys
from typing import Dict, Tuple, Optional, List


def parse_args():
    p = argparse.ArgumentParser(description="Generate lineage table for unified DB assemblies.")
    p.add_argument(
        "--file-list", "-f", required=True,
        help="Text file with one genome file path per line."
    )
    p.add_argument(
        "--assembly-summary", "-a", required=True, nargs="+",
        help="One or more NCBI assembly_summary*.txt files."
    )
    p.add_argument(
        "--taxdump", "-t", required=True,
        help="Directory containing NCBI taxdump (nodes.dmp, names.dmp)."
    )
    p.add_argument(
        "--output", "-o", required=True,
        help="Output TSV file."
    )
    p.add_argument(
        "--log", "-l", required=True,
        help="Log file to record unmapped files and issues."
    )
    return p.parse_args()


def log_open(log_path: str):
    fh = open(log_path, "w")
    def _log(msg: str):
        print(msg, file=fh)
        print(msg, file=sys.stderr)
    return fh, _log


def load_assembly_summary(paths: List[str], log) -> Dict[str, Dict[str, str]]:
    """
    Load one or more assembly_summary files.

    Returns:
        mapping: assembly_accession -> {
            "taxid": taxid,
            "species_taxid": species_taxid
        }
    """
    mapping: Dict[str, Dict[str, str]] = {}

    for path in paths:
        if not os.path.exists(path):
            log(f"[ERROR] assembly_summary file not found: {path}")
            continue

        log(f"[INFO] Loading assembly_summary: {path}")
        with open(path, "r") as fh:
            header = None
            rows_iter = []

            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#assembly_accession"):
                    # This is the header line; strip leading '#' and use it
                    header = line.lstrip("#")
                elif line.startswith("#"):
                    # Comment or description line; skip
                    continue
                else:
                    # Data line
                    rows_iter.append(line)

        if header is None:
            log(f"[WARN] No '#assembly_accession' header line found in {path}. Skipping.")
            continue

        reader = csv.DictReader(rows_iter, delimiter="\t", fieldnames=header.split("\t"))

        # Validate presence of required columns
        required_cols = {"assembly_accession", "taxid", "species_taxid"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            log(f"[WARN] assembly_summary {path} missing columns: {missing}. Still using what is available.")

        for row in reader:
            acc = row.get("assembly_accession", "").strip()
            if not acc:
                continue
            taxid = (row.get("taxid") or "").strip()
            species_taxid = (row.get("species_taxid") or "").strip()

            if not taxid:
                # Sometimes taxid can be missing; log and skip
                log(f"[WARN] assembly_accession {acc} has empty taxid in {path}, skipping.")
                continue

            mapping[acc] = {
                "taxid": taxid,
                "species_taxid": species_taxid if species_taxid else taxid,
            }

    log(f"[INFO] Loaded {len(mapping)} assembly_accession → (taxid, species_taxid) mappings.")
    return mapping


def load_nodes(taxdump_dir: str, log) -> Dict[str, Tuple[str, str]]:
    """
    Load nodes.dmp: taxid -> (parent_taxid, rank)
    """
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
    if not os.path.exists(nodes_path):
        log(f"[ERROR] nodes.dmp not found in {taxdump_dir}")
        sys.exit(1)

    nodes: Dict[str, Tuple[str, str]] = {}
    with open(nodes_path, "r") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            taxid = parts[0]
            parent = parts[1]
            rank = parts[2]
            nodes[taxid] = (parent, rank)
    log(f"[INFO] Loaded {len(nodes)} taxonomy nodes from nodes.dmp.")
    return nodes


def load_names(taxdump_dir: str, log) -> Dict[str, str]:
    """
    Load names.dmp: taxid -> scientific_name (scientific name only).
    """
    names_path = os.path.join(taxdump_dir, "names.dmp")
    if not os.path.exists(names_path):
        log(f"[ERROR] names.dmp not found in {taxdump_dir}")
        sys.exit(1)

    names: Dict[str, str] = {}
    with open(names_path, "r") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            taxid = parts[0]
            name_txt = parts[1]
            name_class = parts[3]
            if name_class == "scientific name":
                names[taxid] = name_txt
    log(f"[INFO] Loaded {len(names)} scientific names from names.dmp.")
    return names


def get_lineage_ranks(
    taxid: str,
    nodes: Dict[str, Tuple[str, str]],
    target_ranks=None
) -> Dict[str, str]:
    """
    Walk up the taxonomy tree from 'taxid' to root, recording the
    first taxid encountered for each rank in target_ranks.

    Returns:
        dict: rank -> taxid
    """
    if target_ranks is None:
        target_ranks = {"species", "genus", "phylum"}

    result: Dict[str, str] = {}
    visited = set()
    current = taxid

    while current in nodes and current not in visited:
        visited.add(current)
        parent, rank = nodes[current]
        if rank in target_ranks and rank not in result:
            result[rank] = current

        if current == parent:  # reached root
            break
        current = parent

    return result


def extract_accession_from_path(path: str) -> Optional[str]:
    """
    Extract assembly accession (e.g., GCF_000001215.4, GCA_045161915.1)
    from a file path using a regex.
    """
    basename = os.path.basename(path)
    m = re.search(r"(GC[AF]_\d+\.\d+)", basename)
    if m:
        return m.group(1)
    return None


def main():
    args = parse_args()
    log_fh, log = log_open(args.log)

    # 1. Load assembly summary mappings
    asm_map = load_assembly_summary(args.assembly_summary, log)

    # 2. Load taxonomy
    nodes = load_nodes(args.taxdump, log)
    names = load_names(args.taxdump, log)

    # 3. Prepare output
    fieldnames = [
        "assembly_accession",
        "assembly_taxid",
        "species_taxid",
        "genus_taxid",
        "phylum",
        "species_name",
    ]
    try:
        out_fh = open(args.output, "w", newline="")
    except OSError as e:
        log(f"[ERROR] Cannot open output file: {e}")
        sys.exit(1)

    writer = csv.DictWriter(out_fh, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    # 4. Process file list
    total_files = 0
    written_rows = 0
    no_accession = 0
    no_asm_entry = 0
    no_taxnode = 0
    missing_species = 0
    missing_genus = 0
    missing_phylum = 0

    log(f"[INFO] Reading file list from: {args.file_list}")
    with open(args.file_list, "r") as fh:
        for line in fh:
            path = line.strip()
            if not path or path.startswith("#"):
                continue
            total_files += 1

            acc = extract_accession_from_path(path)
            if acc is None:
                no_accession += 1
                log(f"[WARN] Could not extract assembly accession from path: {path}")
                continue

            entry = asm_map.get(acc)
            if entry is None:
                no_asm_entry += 1
                log(f"[WARN] assembly_accession {acc} not found in assembly_summary mapping (file: {path})")
                continue

            assembly_taxid = entry["taxid"]
            species_taxid = entry["species_taxid"]

            if species_taxid not in nodes:
                no_taxnode += 1
                log(f"[WARN] species_taxid {species_taxid} for {acc} not found in nodes.dmp")
                lineage = {}
            else:
                lineage = get_lineage_ranks(species_taxid, nodes)

            # lineage ranks
            species_tid = lineage.get("species", species_taxid if species_taxid in nodes else "")
            genus_tid = lineage.get("genus", "")
            phylum_tid = lineage.get("phylum", "")

            species_name = names.get(species_tid, "") if species_tid else ""
            phylum_name = names.get(phylum_tid, "") if phylum_tid else ""

            if not species_tid:
                missing_species += 1
                log(f"[WARN] No species rank resolved for species_taxid {species_taxid} (accession {acc})")
            if not genus_tid:
                missing_genus += 1
                log(f"[WARN] No genus rank found for species_taxid {species_taxid} (accession {acc})")
            if not phylum_tid:
                missing_phylum += 1
                log(f"[WARN] No phylum rank found for species_taxid {species_taxid} (accession {acc})")

            out_row = {
                "assembly_accession": acc,
                "assembly_taxid": assembly_taxid,
                "species_taxid": species_tid or "",
                "genus_taxid": genus_tid or "",
                "phylum": phylum_name,
                "species_name": species_name,
            }
            writer.writerow(out_row)
            written_rows += 1

    out_fh.close()

    # 5. Summary
    log("\n[SUMMARY]")
    log(f"  Total file paths processed:              {total_files}")
    log(f"  Rows written to output:                  {written_rows}")
    log(f"  Paths with no accession parsed:          {no_accession}")
    log(f"  Accessions missing in assembly_summary:  {no_asm_entry}")
    log(f"  species_taxid missing in nodes.dmp:      {no_taxnode}")
    log(f"  Entries with no species resolved:        {missing_species}")
    log(f"  Entries with no genus resolved:          {missing_genus}")
    log(f"  Entries with no phylum resolved:         {missing_phylum}")

    log_fh.close()


if __name__ == "__main__":
    main()
