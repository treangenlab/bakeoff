"""
Helper functions for taxonomic profiling experiments:
- handy utilities for data preprocessing, evaluation, and taxid resolution
Author(s): Wenyu Huang (eh58@rice.edu)
"""

from __future__ import annotations
import pandas as pd
import re

# GENERAL HELPERS
def apply_taxid_synonyms(
    df: pd.DataFrame,
    taxid_col: str = "taxid",
    abundance_col: str = "abundance",
    synonyms: dict = None,
    collapse: bool = False,
    source_label: str = ""
):
    if synonyms is None or len(synonyms) == 0:
        return df.copy()

    df = df.copy()
    df[taxid_col] = pd.to_numeric(df[taxid_col], errors="coerce").astype("Int64")

    mask = df[taxid_col].isin(synonyms.keys())
    n_changed = mask.sum()

    if n_changed > 0:
        ex_old = df.loc[mask, taxid_col].unique()[:5]
        ex_new = [synonyms[o] for o in ex_old]
        print(
            f"[Synonyms:{source_label}] {n_changed} taxids replaced "
            f"(examples: {dict(zip(ex_old, ex_new))})"
        )
        df.loc[mask, taxid_col] = df.loc[mask, taxid_col].map(synonyms)

    if collapse:
        before = len(df)
        df = (
            df.groupby(df.columns.tolist(), dropna=False)
              .sum(numeric_only=True)
              .reset_index()
        )
        after = len(df)
        if before != after:
            print(
                f"[Synonyms:{source_label}] Collapsed duplicates after synonym mapping "
                f"({before} → {after})"
            )

    return df

# CLEAN AND NORMALIZE HELPERS

def _clean_name(n: str) -> str:
    """
    Canonicalize a taxon name string, WITHOUT touching underscores
    inside alphanumeric strain identifiers (e.g. FDAARGOS_192).

    Rules:
    - Remove s__/g__ prefix.
    - Replace underscores ONLY when both neighboring tokens are alphabetic.
    - Collapse extra spaces.
    - Keep underscores in strain names like FDAARGOS_192.
    """
    if not isinstance(n, str):
        return "Unclassified"

    s = n.strip()

    # Remove s__/g__ prefixes
    s = re.sub(r'^[a-z]__', '', s)

    # Split by underscores but only replace them if BOTH sides are alphabetic
    parts = s.split("_")
    cleaned_parts = [parts[0]]

    for i in range(1, len(parts)):
        prev = cleaned_parts[-1]
        curr = parts[i]

        # If both parts are alphabetic → treat underscore as a space
        if prev.isalpha() and curr.isalpha():
            cleaned_parts[-1] = prev + " " + curr
        else:
            cleaned_parts.append(curr)

    s = "_".join(cleaned_parts)

    # Collapse multiple spaces
    s = re.sub(r"\s+", " ", s).strip()

    return s if s else "Unclassified"

# EVALUATION HELPERS

def make_eval_key(row):
    """
    Unique, stable identity for detection metrics.
    """
    taxid = row.get("taxid", None)
    name  = row.get("name", None)

    # ----- Prefer taxid -----
    try:
        tid = int(taxid)
        if tid > 0:
            return f"tid:{tid}"
    except Exception:
        pass

    # ----- Fall back to cleaned name -----
    clean = _clean_name(name) if isinstance(name, str) else None
    if clean:
        return f"name:{clean}"

    # ----- Final fallback -----
    return f"unmapped_idx:{row.name}"


def presence_from_df(df, threshold):
    """
    Return dict eval_key -> present(0/1), from a df with at least:
      - 'taxid'
      - 'value'  (standardized numeric field)
      - 'name'   (for unmapped records)

    eval_key:
      - 'tid:{taxid}'  if taxid > 0
      - 'name:{clean_name}'  if no taxid but a usable name
      - 'unmapped_idx:{row_index}' last resort

    Threshold rule:
      - if threshold == 0:  present if value > 0
      - if threshold > 0:   present if value >= threshold
    """
    if "taxid" not in df.columns:
        raise ValueError("presence_from_df requires a 'taxid' column.")
    if "value" not in df.columns:
        raise ValueError("presence_from_df requires a 'value' column.")

    df = df.copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce").fillna(0.0)

    # Build evaluation key
    df["eval_key"] = df.apply(make_eval_key, axis=1)

    # Aggregate true duplicates at eval_key level
    agg = (
        df.groupby("eval_key", as_index=True)["value"]
          .sum(min_count=1)
    )

    if threshold == 0:
        return (agg > 0).astype(int).to_dict()
    else:
        return (agg >= threshold).astype(int).to_dict()


def prf_counts(pred_present, truth_present):
    """
    Compute TP/FP/FN/TN given two dicts eval_key -> 0/1
    over the union of keys.

    eval_key is typically:
      - "tid:{taxid}" for mapped taxa
      - "name:{clean_name}" or "unmapped_idx:{row_index}" for unmapped ones
    """
    taxa = set(pred_present) | set(truth_present)
    TP = FP = FN = TN = 0
    for t in taxa:
        p = pred_present.get(t, 0)
        r = truth_present.get(t, 0)
        if   p == 1 and r == 1: TP += 1
        elif p == 1 and r == 0: FP += 1
        elif p == 0 and r == 1: FN += 1
        else:                   TN += 1
    return TP, FP, FN, TN


def prf_scores(TP, FP, FN):
    prec = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    rec  = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1   = (2 * prec * rec / (prec + rec)) if (prec + rec) > 0 else 0.0
    return prec, rec, f1

# ETE3 HELPERS

# Raw taxon label extraction

def extract_taxon_label(raw_name: str) -> str:
    """
    Extract a cleaned *taxon label* from a FASTA/header string.

    Preserves:
      - 'uncultured', 'unknown', 'unclassified'
      - species names, 'sp.' labels, strain IDs (SCO41, DSM 3005)
      - multi-word genus names, biovar, serovar, etc.

    Removes:
      - accession prefixes (e.g. 'NZ_CP027225.1', 'NC_000913.3')
      - trailing assembly descriptors:
          chromosome(s), plasmid(s), genome(s), contig(s),
          scaffold(s), assembly/assemblies, sequence(s),
          genomic, complete
      - anything after the first comma (typical assembly info)

    Examples
    --------
    "NZ_CP027225.1 uncultured Phytobacter sp. SCO41 chromosome, complete genome"
        → "uncultured Phytobacter sp. SCO41"
    """

    # Step 1 — normalize spacing & underscores
    s = _clean_name(raw_name)
    if not s:
        return "Unclassified"

    # Step 2 — remove leading accessions
    #   e.g. "NZ_CP027225.1", "NC_000913.3", "GCA_xxxxx"
    s = re.sub(r"^[A-Za-z]{1,3}_[A-Za-z0-9.]+", "", s).strip()

    # Step 3 — drop anything after first comma:
    #   ", complete genome", ", chromosome", etc.
    if "," in s:
        s = s.split(",", 1)[0].strip()

    # Step 4 — trim trailing generic assembly descriptors
    tail_markers = {
        "chromosome", "chromosomes",
        "plasmid", "plasmids",
        "genome", "genomes",
        "contig", "contigs",
        "scaffold", "scaffolds",
        "assembly", "assemblies",
        "sequence", "sequences",
        "genomic", "complete",
    }

    toks = s.split()
    while toks and toks[-1].lower().strip(",") in tail_markers:
        toks.pop()

    cleaned = " ".join(toks).strip()

    return cleaned if cleaned else "Unclassified"

# Species label extraction

def extract_species_label(raw_name: str) -> str:
    """
    Extract a species-like label from messy strain-level names.

    Examples:
        "Blautia massiliensis (ex Durand et al. 2017) strain FDAARGOS_1576"
            → "Blautia massiliensis (ex Durand et al. 2017)"

        "Staphylococcus haemolyticus strain ATCC 29970"
            → "Staphylococcus haemolyticus"

    Strategy:
        - Clean spacing minimally, DO NOT remove parentheses.
        - Trim from keywords: strain, isolate, sample, clone, etc.
        - Trim trailing uppercase codes like _A, _B, _AQ, _XYZ.
        - Preserve GTDB-style names (no trimming for uppercase suffix in GTDB context).
    """

    if not isinstance(raw_name, str):
        return ""

    s = raw_name.strip()

    # 1. Trim on keywords indicating subspecies/strain/isolate
    # Add more if needed
    STRIP_KEYWORDS = [
        r"\bstrain\b",
        r"\bisolate\b",
        r"\bsample\b",
        r"\bclone\b",
        r"\bsubsp\b",
        r"\bsubstrain\b",
        r"\bserotype\b",
        r"\bATCC\b",
        r"\bDSM\b",
        r"\bJCM\b",
        r"\bFDAARGOS\b",
    ]

    pattern = re.compile("(" + "|".join(STRIP_KEYWORDS) + ")", flags=re.IGNORECASE)
    m = pattern.search(s)

    if m:
        # Trim before keyword
        s = s[:m.start()].strip()

    # 2. Trim trailing PURE-UPPERCASE suffixes or _ABC patterns
    #    BUT ONLY when the base is NOT gtdb-style (no genus_A format)
    # Example: "Blautia massiliensis AQ12" → "Blautia massiliensis"
    #          "Enterococcus_B lactis"     (should NOT trim in GTDB context)
    #
    # We infer non-GTDB context when name contains spaces (plain species names)
    #
    if " " in s:
        # Remove trailing strain-like codes such as "AQ12" or "FDAARGOS_1576",
        # but preserve pure numeric oral-taxon identifiers like "... taxon 807".
        tokens = s.split()
        while tokens:
            last = tokens[-1]
            prev = tokens[-2].lower().strip(".") if len(tokens) >= 2 else ""

            if prev == "taxon" and re.fullmatch(r"\d+", last):
                break

            if re.fullmatch(r"[A-Z0-9_]+", last) and re.search(r"[A-Z_]", last):
                tokens.pop()
                continue

            break

        s = " ".join(tokens).strip()

    # Final cleanup
    s = re.sub(r"\s+", " ", s).strip()

    return s

# Genus label extraction

def infer_genus_label(name: str) -> str:
    """
    Best-effort genus label from a species-like name.

    Special cases:
      - Single-word 'uncultured/unknown/unclassified' -> keep as-is.
      - 'uncultured/unknown/unclassified Genus ...' -> 'uncultured Genus'
        (e.g. 'uncultured Phytobacter sp. SCO41' -> 'uncultured Phytobacter').

    Otherwise:
      - If starts with 'Candidatus/ Ca.' and has >= 2 tokens, use first two tokens.
      - Else if contains 'sp.' or 'sp', use the token before 'sp.' as genus.
      - Else, use the first token as genus.
    """
    s = _clean_name(name)
    if not s:
        return "Unclassified"

    tokens = s.split()
    if not tokens:
        return "Unclassified"

    lower_tokens = [t.lower() for t in tokens]
    uc_tokens = {"uncultured", "unknown", "unclassified"}

    # 1) Single-word uncultured/unknown/unclassified
    if len(tokens) == 1 and lower_tokens[0] in uc_tokens:
        return tokens[0]

    # 2) 'uncultured/unknown/unclassified Genus ...' -> 'uncultured Genus'
    if lower_tokens[0] in uc_tokens and len(tokens) >= 2:
        return f"{tokens[0]} {tokens[1]}"

    # 3) 'Candidatus Genus species' or 'Ca. Genus species'
    if lower_tokens[0] in {"candidatus", "ca."} and len(tokens) >= 2:
        # two-word genus label
        return " ".join(tokens[:2])

    # 4) 'Genus sp.' patterns (no uncultured/unknown at front)
    if "sp." in lower_tokens:
        idx = lower_tokens.index("sp.")
        if idx > 0:
            return tokens[idx - 1]
    if "sp" in lower_tokens:
        idx = lower_tokens.index("sp")
        if idx > 0:
            return tokens[idx - 1]

    # 5) Default: first token as genus
    return tokens[0]

# TAXID RESOLUTION

def project_taxid_to_rank(taxid: int, target_rank: str, ncbi) -> int | None:
    """
    Project any taxid up or down its NCBI lineage to a given rank.

    Parameters
    ----------
    taxid : int
        The starting taxid (species, strain, subspecies, genus, etc.)
    target_rank : str
        Desired rank: "species" or "genus" (case-insensitive)
    ncbi : ete3.NCBITaxa instance

    Returns
    -------
    int or None
        - taxid of the lineage node at the target rank
        - None if no such rank found
    """
    if taxid is None:
        return None
    try:
        tid = int(taxid)
    except (TypeError, ValueError):
        return None
    if tid <= 0:
        return None

    target_rank = target_rank.lower()

    try:
        lineage = ncbi.get_lineage(tid)
        ranks = ncbi.get_rank(lineage)
    except Exception:
        return None

    # First: if the node *itself* is the target rank
    if ranks.get(tid, "").lower() == target_rank:
        return tid

    # Otherwise: search lineage for the first parent at that rank
    for t in lineage:
        if ranks.get(t, "").lower() == target_rank:
            return t

    return None


def canonicalize_merged(taxids, ncbi) -> dict:
    """Map retired taxids to their current id via NCBI's merged.dmp.

    ete3's get_rank/get_lineage do not auto-apply merged.dmp, so a retired taxid
    resolves to nothing and gets dropped — e.g. a retired genus node loses its
    cumulative mass, making Sigma-genus < Sigma-species. Returns {old: new} for the
    retired ids among `taxids` only (current ids are absent from the map).
    """
    ids = [int(t) for t in taxids if t and int(t) > 0]
    if not ids:
        return {}
    _, merged = ncbi._translate_merged(ids)
    return merged


def get_taxid_from_name(name: str, rank_level: str, ncbi, clean=True) -> int:
    """
    Resolve taxid from name using ete3.
    Return 0 for unmapped; do NOT drop.
    """
    if clean:
        name = _clean_name(name)
    if not name:
        return 0

    try:
        trans = ncbi.get_name_translator([name])
    except Exception:
        trans = {}

    cand = trans.get(name, [])

    if cand:
        ranks = ncbi.get_rank(cand)
        # exact rank
        for tid in cand:
            if ranks.get(tid, "").lower() == rank_level:
                return int(tid)
        # lineage projection
        for tid in cand:
            try:
                lineage = ncbi.get_lineage(tid)
                lranks = ncbi.get_rank(lineage)
                target = next((t for t,r in lranks.items() if r == rank_level), None)
                if target:
                    return int(target)
            except Exception:
                pass

    return 0

# Best effort genus resolution from species
def resolve_genus_from_species(
    species_tid,
    species_name: str,
    ncbi,
) -> tuple[int, str]:
    """
    Safe: NEVER double-maps. Species-taxid lineage resolution ALWAYS stops the process.
    Returns (genus_tid, genus_name), taxid=0 allowed, name always kept.
    """

    # Prepare
    genus_tid = 0
    genus_name = "Unclassified"

    sname = _clean_name(species_name if isinstance(species_name, str) else "")
    if not sname:
        sname = "Unclassified"

    # STEP 1. TRY SPECIES TAXID → GENUS TAXID
    try:
        stid = int(species_tid)
    except Exception:
        stid = 0

    if stid > 0:
        try:
            lineage = ncbi.get_lineage(stid)
            ranks = ncbi.get_rank(lineage)
            genus_ids = [tid for tid, rk in ranks.items() if rk == "genus"]

            if genus_ids:
                gid = int(genus_ids[0])
                genus_tid = gid

                # Try canonical genus name
                try:
                    nm = ncbi.get_taxid_translator([gid]).get(gid)
                    if isinstance(nm, str) and nm.strip():
                        genus_name = _clean_name(nm)
                    else:
                        genus_name = infer_genus_label(sname)
                except Exception:
                    genus_name = infer_genus_label(sname)

                return genus_tid, genus_name   # EARLY EXIT
        except Exception:
            pass

    # STEP 2. FALLBACK: name-based genus inference
    label = infer_genus_label(sname)
    genus_name = label

    try:
        gtid = get_taxid_from_name(label, "genus", ncbi)
    except Exception:
        gtid = 0

    if gtid and gtid > 0:
        genus_tid = int(gtid)
        try:
            nm = ncbi.get_taxid_translator([genus_tid]).get(genus_tid)
            if isinstance(nm, str) and nm.strip():
                genus_name = _clean_name(nm)
        except Exception:
            pass

    return genus_tid, genus_name
