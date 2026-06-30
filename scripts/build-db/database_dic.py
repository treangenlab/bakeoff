#!/usr/bin/env python3
"""
Map genome files to taxonomy IDs using RefSeq assembly summary.
"""

import os
import re
from pathlib import Path


def extract_accession(filepath):
    """
    Extract assembly accession number from file path.
    Example: library/bacteria/GCF_045161915.1_ASM4516191v1_genomic.fa -> GCF_045161915.1
    """
    filename = os.path.basename(filepath)
    # Match pattern: GCF_XXXXXXXXX.X or GCA_XXXXXXXXX.X
    match = re.match(r'(GC[AF]_\d+\.\d+)', filename)
    if match:
        return match.group(1)
    return None


def read_assembly_summary(summary_file):
    """
    Read assembly summary file and create a mapping of accession -> taxid.
    Assumes the file has headers starting with '#' and columns separated by tabs.
    """
    accession_to_taxid = {}
    
    with open(summary_file, 'r') as f:
        for line in f:
            # Skip comment lines except the column header line
            if line.startswith('#'):
                if 'assembly_accession' in line:
                    # This is the header line, extract column indices
                    header = line.lstrip('#').strip().split('\t')
                    try:
                        acc_idx = header.index('assembly_accession')
                        taxid_idx = header.index('taxid')
                    except ValueError as e:
                        print(f"Error: Could not find required columns in header: {e}")
                        raise
                continue
            
            # Parse data lines
            fields = line.strip().split('\t')
            if len(fields) > max(acc_idx, taxid_idx):
                accession = fields[acc_idx]
                taxid = fields[taxid_idx]
                accession_to_taxid[accession] = taxid
    
    return accession_to_taxid


def process_genome_files(path_file, accession_to_taxid, output_file, missing_file):
    """
    Process genome file paths and map them to taxids.
    """
    mapped = []
    missing = []
    
    with open(path_file, 'r') as f:
        for line in f:
            filepath = line.strip()
            if not filepath:
                continue
            
            # Extract accession number
            accession = extract_accession(filepath)
            
            if accession is None:
                missing.append((filepath, "Could not extract accession number"))
                continue
            
            # Look up taxid
            if accession in accession_to_taxid:
                taxid = accession_to_taxid[accession]
                # Convert to absolute path
                abs_path = os.path.abspath(filepath)
                mapped.append((accession, taxid, abs_path))
            else:
                missing.append((filepath, f"Accession {accession} not found in summary"))
    
    # Write mapped entries
    with open(output_file, 'w') as f:
        f.write("assembly_accession\ttaxid\tfile_path\n")
        for accession, taxid, path in mapped:
            f.write(f"{accession}\t{taxid}\t{path}\n")
    
    # Write missing entries
    with open(missing_file, 'w') as f:
        if not missing:
            f.write("All genome files are mapped to the taxid\n")
        else:
            f.write("file_path\treason\n")
            for filepath, reason in missing:
                f.write(f"{filepath}\t{reason}\n")
    
    return len(mapped), len(missing)


def main():
    # Configuration
    assembly_summary_file = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/assembly_summary_refseq.txt"  # Path to RefSeq assembly summary
    path_file = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/file.list"  # Path to file containing genome file paths
    output_file = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/refseq03032025_dict.txt"
    missing_file = "/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/unmapped_genomes.txt"
    
    print("Reading assembly summary file...")
    accession_to_taxid = read_assembly_summary(assembly_summary_file)
    print(f"Loaded {len(accession_to_taxid)} accession-taxid mappings")
    
    print("\nProcessing genome files...")
    mapped_count, missing_count = process_genome_files(
        path_file, accession_to_taxid, output_file, missing_file
    )
    
    print(f"\nResults:")
    print(f"  Mapped: {mapped_count} files -> {output_file}")
    print(f"  Missing: {missing_count} files -> {missing_file}")
    print("\nDone!")


if __name__ == "__main__":
    main()
