#!/usr/bin/env python3
import csv
import sys

# Hard-coded file paths
SEQ2TAXID_PATH = '/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/seqid2taxid.map'
IDENTS_PATH    = '/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/sm_030325/build_030325.csv'
NODES_DMP      = '/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/taxonomy/nodes.dmp'
NAMES_DMP      = '/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/taxonomy/names.dmp'
OUTPUT_PATH    = '/home/Users/pacbio_bakeoff/data/ref_db/refseq03032025/sm_030325/lineage_030325.csv'

def load_seq2taxid(path):
    """Load mapping of sequence identifiers to taxids."""
    seq2tax = {}
    with open(path) as f:
        for line in f:
            ident, taxid = line.strip().split()
            seq2tax[ident] = taxid
    return seq2tax


def load_idents(path, colname='ident'):
    """Load the list of idents from a CSV file with header 'ident'."""
    idents = []
    with open(path) as f:
        reader = csv.DictReader(f)
        if colname not in reader.fieldnames:
            raise ValueError(f"Input file {path!r} has no column {colname!r}")
        for row in reader:
            idents.append(row[colname])
    return idents


def load_names(names_dmp):
    """Load scientific names for each taxid from names.dmp."""
    taxid2name = {}
    with open(names_dmp) as f:
        for line in f:
            cols = [c.strip() for c in line.split('|')]
            tid, name_txt, _, name_class = cols[0], cols[1], cols[2], cols[3]
            if name_class == 'scientific name':
                taxid2name[tid] = name_txt
    return taxid2name


def load_nodes(nodes_dmp):
    """Load parent pointers and rank information from nodes.dmp."""
    parent = {}
    rank   = {}
    with open(nodes_dmp) as f:
        for line in f:
            cols = [c.strip() for c in line.split('|')[:3]]
            tid, p_tid, rnk = cols
            parent[tid] = p_tid
            rank[tid]   = rnk
    return parent, rank


def get_lineage(tid, parent, rank, taxid2name, ranks):
    """Walk up the taxonomic tree to assemble a lineage dict for desired ranks."""
    lin = dict.fromkeys(ranks, '')
    curr = tid
    while True:
        r = rank.get(curr)
        if r in lin:
            lin[r] = taxid2name.get(curr, '')
        if curr == '1' or curr not in parent:
            break
        curr = parent[curr]
    return lin


def main():
    # Load all inputs
    seq2tax     = load_seq2taxid(SEQ2TAXID_PATH)
    idents      = load_idents(IDENTS_PATH)
    taxid2name  = load_names(NAMES_DMP)
    parent, rank = load_nodes(NODES_DMP)

    # Desired taxonomic ranks in order
    ranks = ['superkingdom','phylum','class','order','family','genus','species','strain']

    # Write out filtered lineage CSV
    with open(OUTPUT_PATH, 'w', newline='') as out:
        writer = csv.writer(out)
        # use 'ident' for header to match your CSV
        writer.writerow(['ident','taxid'] + ranks)
        for ident in idents:
            tid = seq2tax.get(ident)
            if not tid:
                print(f"WARNING: no taxid found for {ident!r}", file=sys.stderr)
                continue
            lin = get_lineage(tid, parent, rank, taxid2name, ranks)
            writer.writerow([ident, tid] + [lin[r] for r in ranks])

if __name__ == '__main__':
    main()
