#!/usr/bin/env python3
"""
Build a ganon_lineage table from:
  1) a file list of genome paths
  2) an NCBI assembly summary file
  3) a ganon taxonomy file

Outputs:
  - ganon_lineage.tsv
  - logs:
      * ganon_lineage_missing_in_assembly_summary.tsv
      * ganon_lineage_missing_taxid_in_taxonomy.tsv
      * ganon_lineage_lineage_resolution_errors.tsv
"""

import argparse
import csv
import os
import sys
import re
from collections import defaultdict

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------


def parse_args():
    p = argparse.ArgumentParser(
        description="Build ganon_lineage.tsv by linking file list, assembly summary, and ganon taxonomy."
    )
    p.add_argument(
        "--file-list",
        required=True,
        help="Text file: one genome path per line (paths containing GCF_/GCA_ accession).",
    )
    p.add_argument(
        "--assembly-summary",
        required=True,
        help="NCBI assembly summary file (tab-delimited). Comment lines start with '#'.",
    )
    p.add_argument(
        "--ganon-taxonomy",
        required=True,
        help="Ganon taxonomy file (no header). Columns: taxid, parent_taxid, rank, name, internal_id.",
    )
    p.add_argument(
        "--out-prefix",
        default="ganon_lineage",
        help="Prefix for output files (default: ganon_lineage). "
             "Main table will be <prefix>.tsv.",
    )
    return p.parse_args()


# ---------------------------------------------------------------------
# Step 1: file list -> assembly_accession
# ---------------------------------------------------------------------

def extract_assembly_accession_from_path(path):
    """Return GC[AF]_########.# from a genome file path, or None if absent."""
    m = re.search(r'(GC[AF]_\d+\.\d+)', os.path.basename(path))
    return m.group(1) if m else None

def load_file_list(file_list_path):
    """
    Read file list, return dict:
      assembly_accession -> list of file paths
    """
    mapping = defaultdict(list)
    with open(file_list_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            acc = extract_assembly_accession_from_path(line)
            mapping[acc].append(os.path.abspath(line))
    return mapping


# ---------------------------------------------------------------------
# Step 2: assembly summary
# ---------------------------------------------------------------------


def _to_int_or_none(x):
    return int(x) if x and x != "na" else None


def load_assembly_summary(summary_path):
    """
    Load assembly summary file.

    Assumes:
      - tab-delimited
      - comment lines start with '#'
      - columns:
          1: assembly_accession
          6: taxid
          7: species_taxid
          8: organism_name

    Returns:
      dict[assembly_accession] -> {
          "assembly_taxid": int or None,
          "species_taxid": int or None,
          "organism_name": str
      }
    """
    asm_dict = {}
    with open(summary_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            acc = parts[0]
            asm_dict[acc] = {
                "assembly_taxid": _to_int_or_none(parts[5].strip()),
                "species_taxid": _to_int_or_none(parts[6].strip()),
                "organism_name": parts[7].strip(),
            }
    return asm_dict


# ---------------------------------------------------------------------
# Step 3: ganon taxonomy
# ---------------------------------------------------------------------


def load_ganon_taxonomy(tax_path):
    """
    Load ganon taxonomy file.

    Columns (no header):
      0: taxid
      1: parent_taxid
      2: rank
      3: name
      4: internal_id (ignored)

    Returns:
      parent: dict[taxid] -> parent_taxid
      rank:   dict[taxid] -> rank_string
      name:   dict[taxid] -> name_string
    """
    parent = {}
    rank = {}
    name = {}
    with open(tax_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                taxid = int(parts[0])
                parent_id = int(parts[1])
            except ValueError:
                continue
            rnk = parts[2].strip()
            nm = parts[3].strip()
            parent[taxid] = parent_id
            rank[taxid] = rnk
            name[taxid] = nm
    return parent, rank, name


# ---------------------------------------------------------------------
# Step 4: lineage walker + domain classifier
# ---------------------------------------------------------------------


def get_lineage_nodes(start_taxid, parent_map, rank_map, name_map):
    """
    Climb from start_taxid to root using parent_map.
    Collect first occurrences of species, genus, phylum, and domain/superkingdom.
    Also collect full lineage (list of taxids from leaf to root).
    """
    species_taxid = None
    genus_taxid = None
    phylum_taxid = None
    domain_taxid = None

    species_name = None
    genus_name = None
    phylum_name = None
    domain_name = None

    lineage_taxids = []

    current = start_taxid
    visited = set()

    while current is not None and current not in visited and current in parent_map:
        visited.add(current)
        lineage_taxids.append(current)

        rnk = rank_map.get(current, "").lower()
        nm = name_map.get(current, "")

        # Record first matches
        if rnk == "species" and species_taxid is None:
            species_taxid = current
            species_name = nm
        if rnk == "genus" and genus_taxid is None:
            genus_taxid = current
            genus_name = nm
        if rnk == "phylum" and phylum_taxid is None:
            phylum_taxid = current
            phylum_name = nm
        if rnk in ("superkingdom", "domain") and domain_taxid is None:
            domain_taxid = current
            domain_name = nm

        # Root detection: when parent == self or no change
        parent_id = parent_map.get(current)
        if parent_id is None or parent_id == current:
            break
        current = parent_id

    return {
        "species_taxid": species_taxid,
        "genus_taxid": genus_taxid,
        "phylum_taxid": phylum_taxid,
        "domain_taxid": domain_taxid,
        "species_name": species_name,
        "genus_name": genus_name,
        "phylum_name": phylum_name,
        "domain_name": domain_name,
        "lineage_taxids": lineage_taxids,
    }


def classify_domain(lineage_info, name_map):
    """
    Map lineage to one of: Bacteria, Archaea, Viruses, Fungi, Other.
    Uses domain_name if present, otherwise inspects names along the lineage.
    """
    domain_name = (lineage_info.get("domain_name") or "").lower()
    lineage_taxids = lineage_info["lineage_taxids"]

    # Shortcut using domain/superkingdom name if available
    if domain_name == "bacteria":
        return "Bacteria"
    if domain_name == "archaea":
        return "Archaea"

    # Collect all names in the lineage for fallback checks
    lineage_names = [name_map.get(t, "").lower() for t in lineage_taxids]

    # Viruses: any node named 'viruses' or containing 'virus'
    if any("virus" in n for n in lineage_names):
        return "Viruses"

    # Fungi: any node with 'fungi' or 'fungus'
    if any("fungi" in n or "fungus" in n for n in lineage_names):
        return "Fungi"

    return "Other"


# ---------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------


def main():
    args = parse_args()

    # 1) File list
    print(f"[INFO] Loading file list from {args.file_list}", file=sys.stderr)
    asm_to_files = load_file_list(args.file_list)
    print(f"[INFO] Found {len(asm_to_files)} distinct assembly accessions in file list.", file=sys.stderr)

    # 2) Assembly summary
    print(f"[INFO] Loading assembly summary from {args.assembly_summary}", file=sys.stderr)
    asm_summary = load_assembly_summary(args.assembly_summary)
    print(f"[INFO] Loaded {len(asm_summary)} assemblies from summary.", file=sys.stderr)

    # 3) Ganon taxonomy
    print(f"[INFO] Loading ganon taxonomy from {args.ganon_taxonomy}", file=sys.stderr)
    parent_map, rank_map, name_map = load_ganon_taxonomy(args.ganon_taxonomy)
    print(f"[INFO] Taxonomy loaded: {len(parent_map)} taxids.", file=sys.stderr)

    out_main = f"{args.out_prefix}.tsv"
    out_missing_asm = f"{args.out_prefix}_missing_in_assembly_summary.tsv"
    out_missing_tax = f"{args.out_prefix}_missing_taxid_in_taxonomy.tsv"
    out_lineage_err = f"{args.out_prefix}_lineage_resolution_errors.tsv"

    n_main_rows = 0
    n_missing_asm = 0
    n_missing_tax = 0
    n_lineage_err = 0

    with open(out_main, "w", newline="") as fout, \
         open(out_missing_asm, "w", newline="") as fmiss_asm, \
         open(out_missing_tax, "w", newline="") as fmiss_tax, \
         open(out_lineage_err, "w", newline="") as ferr:

        main_writer = csv.writer(fout, delimiter="\t")
        miss_asm_writer = csv.writer(fmiss_asm, delimiter="\t")
        miss_tax_writer = csv.writer(fmiss_tax, delimiter="\t")
        err_writer = csv.writer(ferr, delimiter="\t")

        # Headers
        main_writer.writerow([
            "assembly_accession",
            "file_path",
            "assembly_taxid",
            "species_taxid_input",
            "species_taxid_lineage",
            "genus_taxid",
            "phylum_taxid",
            "domain_taxid",
            "domain_label",
            "species_name",
            "genus_name",
            "phylum_name",
            "domain_name",
            "organism_name",
        ])

        miss_asm_writer.writerow(["assembly_accession", "file_path", "reason"])
        miss_tax_writer.writerow(["assembly_accession", "file_path", "species_taxid_input", "organism_name"])
        err_writer.writerow(["assembly_accession", "file_path", "species_taxid_input", "message"])

        # Process each assembly in file list
        for asm_acc, file_paths in asm_to_files.items():
            if asm_acc not in asm_summary:
                for fp in file_paths:
                    miss_asm_writer.writerow([asm_acc, fp, "assembly_accession_not_in_summary"])
                    n_missing_asm += 1
                continue

            asm_info = asm_summary[asm_acc]
            asm_taxid = asm_info.get("assembly_taxid")
            species_taxid_input = asm_info.get("species_taxid")
            organism_name = asm_info.get("organism_name", "")

            # Decide which taxid to use as start: prefer species_taxid, else assembly_taxid
            start_taxid = species_taxid_input or asm_taxid
            if start_taxid is None:
                # Nothing to do: no taxid at all
                for fp in file_paths:
                    err_writer.writerow([asm_acc, fp, str(species_taxid_input), "no_taxid_available"])
                    n_lineage_err += 1
                continue

            if start_taxid not in parent_map:
                for fp in file_paths:
                    miss_tax_writer.writerow([asm_acc, fp, str(species_taxid_input), organism_name])
                    n_missing_tax += 1
                continue

            # Walk lineage
            lineage = get_lineage_nodes(start_taxid, parent_map, rank_map, name_map)
            domain_label = classify_domain(lineage, name_map)

            species_tax_lineage = lineage.get("species_taxid")
            genus_taxid = lineage.get("genus_taxid")
            phylum_taxid = lineage.get("phylum_taxid")
            domain_taxid = lineage.get("domain_taxid")

            species_name = lineage.get("species_name")
            genus_name = lineage.get("genus_name")
            phylum_name = lineage.get("phylum_name")
            domain_name = lineage.get("domain_name")

            if species_tax_lineage is None and genus_taxid is None and phylum_taxid is None:
                # Could not resolve any useful rank
                for fp in file_paths:
                    err_writer.writerow([asm_acc, fp, str(species_taxid_input), "no_spec_genus_phylum_resolved"])
                    n_lineage_err += 1
                continue

            # Write one row per file_path (in case multiple files per assembly)
            for fp in file_paths:
                main_writer.writerow([
                    asm_acc,
                    fp,
                    asm_taxid if asm_taxid is not None else "",
                    species_taxid_input if species_taxid_input is not None else "",
                    species_tax_lineage if species_tax_lineage is not None else "",
                    genus_taxid if genus_taxid is not None else "",
                    phylum_taxid if phylum_taxid is not None else "",
                    domain_taxid if domain_taxid is not None else "",
                    domain_label,
                    species_name or "",
                    genus_name or "",
                    phylum_name or "",
                    domain_name or "",
                    organism_name,
                ])
                n_main_rows += 1

    print(f"[INFO] Wrote {n_main_rows} rows to {out_main}", file=sys.stderr)
    print(f"[INFO] Missing assemblies in summary: {n_missing_asm} (see {out_missing_asm})", file=sys.stderr)
    print(f"[INFO] Missing taxids in taxonomy: {n_missing_tax} (see {out_missing_tax})", file=sys.stderr)
    print(f"[INFO] Lineage resolution errors: {n_lineage_err} (see {out_lineage_err})", file=sys.stderr)


if __name__ == "__main__":
    main()
