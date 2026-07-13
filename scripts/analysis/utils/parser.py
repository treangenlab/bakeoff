"""
Parser helpers for loading reports of taxonomic profilers
- return "taxid, name, rank, abundance" with both raw and cleaned values, along with report specific fileds
- map taxid for reports without taxids utilizing ete3 with database aligned taxdump
- fix broken Centrifuge reports (missing taxa and/or abundances) with ete3
- recommend using database corresponding taxdump
Author(s): Wenyu Huang (eh58@rice.edu)
"""

from __future__ import annotations
import pandas as pd
import re
from pathlib import Path

from .util import (
    _clean_name,
    project_taxid_to_rank,
    canonicalize_merged,
    get_taxid_from_name,
    resolve_genus_from_species,
    extract_species_label,
)

# Known truth-name aliases where the simulation uses an older or alternate
# species spelling than the unified taxonomy database.
TRUTH_NAME_ALIASES = {
    "Pseudocoprococcus catus": "Coprococcus catus",
}

TRUTH_GENUS_ALIASES = {
    "Pseudocoprococcus": "Coprococcus",
}


def _iter_truth_resolution_candidates(original_name: str) -> list[str]:
    """
    Build an ordered set of candidate names for GT resolution.

    We first try the raw label, then cleaned/species-like forms, then
    known alias rewrites for taxonomy renames.
    """
    candidates: list[str] = []

    def add(name: str | None) -> None:
        if not isinstance(name, str):
            return
        name = name.strip()
        if name and name not in candidates:
            candidates.append(name)

    add(original_name)
    add(_clean_name(original_name))

    species_like = extract_species_label(original_name)
    add(species_like)
    add(_clean_name(species_like))

    # Expand exact whole-name aliases first.
    for name in list(candidates):
        add(TRUTH_NAME_ALIASES.get(name))

    # Expand genus-level aliases next.
    for name in list(candidates):
        parts = name.split(maxsplit=1)
        if len(parts) != 2:
            continue
        genus, rest = parts
        alias_genus = TRUTH_GENUS_ALIASES.get(genus)
        if alias_genus:
            add(f"{alias_genus} {rest}")

    return candidates

# Truth file parser

def build_unified_ground_truth(
    dataset_path,
    ncbi,
    out_dir,
    name_col: str,
    abund_col: str | None = None,
    sep: str = "\t",
    label: str | None = None,
):
    """
    Build a unified ground-truth table (species + genus) for a dataset.

    Returns unified_gt with columns: ['rank', 'name', 'taxid', 'abundance'].

    Output filename: `<label>_gt.csv` if `label` is given, else `<source-stem>.csv`
    (preserves the older behavior for callers that don't pass a label).
    """
    from pathlib import Path
    import pandas as pd

    dataset_path = Path(dataset_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Building unified ground truth for: {dataset_path} ===")

    # 1) Read file with user-provided separator
    try:
        gt_df = pd.read_csv(dataset_path, sep=sep)
    except Exception as e:
        print(f"Error reading ground truth file: {e}")
        raise

    print(f"Loaded {len(gt_df)} rows from {dataset_path}")
    print(f"Columns found: {list(gt_df.columns)}")

    # 2) Check name column
    if name_col not in gt_df.columns:
        raise ValueError(
            f"Ground truth file must contain a '{name_col}' column.\n"
            f"Found columns: {list(gt_df.columns)}"
        )

    # Standardize name column
    gt_df = gt_df.rename(columns={name_col: "name"})
    print(f"Using '{name_col}' as taxon name column → standardized to 'name'")

    # 3) Resolve each unique name to a clade taxid, then project to species/genus
    unique_labels = gt_df["name"].dropna().unique()
    print(f"Found {len(unique_labels)} unique names to resolve")

    have_abundance = abund_col is not None and abund_col in gt_df.columns

    species_records = []
    genus_records = []
    failed_entries = []

    # optional small cache to avoid repeated taxid→name lookups
    genus_name_cache: dict[int, str] = {}

    for raw_name in unique_labels:
        original_name = str(raw_name).strip()
        if not original_name:
            continue

        display_name = original_name  # used for species-level name
        resolved_candidates = []
        candidate_taxids = []

        for candidate_name in _iter_truth_resolution_candidates(original_name):
            resolved_candidates.append(candidate_name)
            try:
                trans = ncbi.get_name_translator([candidate_name])
            except Exception:
                trans = {}

            candidate_taxids = trans.get(candidate_name, [])
            if candidate_taxids:
                display_name = candidate_name
                break

        # If still nothing, record as unresolved
        if not candidate_taxids:
            count = int((gt_df["name"] == raw_name).sum())
            reason = (
                "Name not found in NCBI taxonomy after trying: "
                + " | ".join(resolved_candidates)
            )
            failed_entries.append({
                "name": original_name,
                "sequence_count": count,
                "reason": reason,
            })
            continue

        # Use first candidate as clade root
        clade_tid = int(candidate_taxids[0])

        abundance_value = pd.NA
        if have_abundance:
            abundance_value = pd.to_numeric(
                gt_df.loc[gt_df["name"] == raw_name, abund_col],
                errors="coerce",
            ).sum(min_count=1)

        # project to species and genus taxids
        species_tid = project_taxid_to_rank(clade_tid, "species", ncbi)
        genus_tid   = project_taxid_to_rank(clade_tid, "genus",   ncbi)

        # --- Species record: keep display_name (species-like) ---
        if species_tid > 0:
            species_records.append({
                "name": display_name,
                "taxid": species_tid,
                "genus_taxid": genus_tid if genus_tid > 0 else pd.NA,
                "abundance": abundance_value,
            })

        # --- Genus record: use canonical genus name from ete3 ---
        if genus_tid > 0:
            if genus_tid not in genus_name_cache:
                try:
                    genus_name_cache[genus_tid] = ncbi.get_taxid_translator([genus_tid]).get(genus_tid, display_name)
                except Exception:
                    genus_name_cache[genus_tid] = display_name  # fallback
            genus_name = genus_name_cache[genus_tid]

            genus_records.append({
                "name": genus_name,   # <- proper genus name here
                "taxid": genus_tid,
                "abundance": abundance_value,
            })

    # 4) Build species/genus DataFrames (deduplicated)
    species_df = pd.DataFrame(species_records)
    if not species_df.empty:
        species_agg = {
            "name": "first",
            "genus_taxid": "first",
        }
        if have_abundance:
            species_agg["abundance"] = "sum"
        species_df = (
            species_df
            .groupby("taxid", as_index=False, dropna=False)
            .agg(species_agg)
            .reset_index(drop=True)
        )
    else:
        species_df = pd.DataFrame(columns=["name", "taxid", "genus_taxid"])
        if have_abundance:
            species_df["abundance"] = pd.Series(dtype="float64")

    genus_df = pd.DataFrame(genus_records)
    if not genus_df.empty:
        genus_agg = {"name": "first"}
        if have_abundance:
            genus_agg["abundance"] = "sum"
        genus_df = (
            genus_df
            .groupby("taxid", as_index=False, dropna=False)
            .agg(genus_agg)
            .reset_index(drop=True)
        )
    else:
        genus_df = pd.DataFrame(columns=["name", "taxid"])
        if have_abundance:
            genus_df["abundance"] = pd.Series(dtype="float64")

    # 5) Attach abundances if provided (species -> genus aggregation)
    if have_abundance:
        print(f"Using '{abund_col}' as abundance column.")
    else:
        if abund_col is None:
            print("No abundance column provided (abund_col=None) → leaving abundances empty.")
        else:
            print(f"Abundance column '{abund_col}' not found → leaving abundances empty.")
        species_df["abundance"] = pd.NA
        genus_df["abundance"] = pd.NA

    # 6) Build unified table (drop internal genus_taxid column from species)
    species_out = species_df[["name", "taxid", "abundance"]].copy()
    genus_out   = genus_df[["name", "taxid", "abundance"]].copy()

    species_out = species_out.assign(rank="species")
    genus_out   = genus_out.assign(rank="genus")

    unified_gt = pd.concat([species_out, genus_out], ignore_index=True)
    unified_gt = unified_gt[["rank", "name", "taxid", "abundance"]]

    # 7) Save unified GT
    # Prefer the caller-supplied registry label (gives `<label>_gt.csv`);
    # fall back to the source file's stem for callers that don't pass one.
    if label:
        out_stem = f"{label}_gt"
    else:
        out_stem = dataset_path.stem
    out_csv = out_dir / f"{out_stem}.csv"
    unified_gt.to_csv(out_csv, index=False)
    print(f"Unified ground truth written to: {out_csv}")

    # 8) Save unresolved-name report (rewrite every time to avoid stale output)
    failed_txt = out_dir / f"{out_stem}_unresolved_names.txt"
    failed_df = pd.DataFrame(failed_entries) if failed_entries else pd.DataFrame(
        columns=["name", "sequence_count", "reason"]
    )
    total_failed = int(failed_df["sequence_count"].sum()) if not failed_df.empty else 0
    with open(failed_txt, "w") as f:
        f.write("# Names that could not be resolved to NCBI taxids\n")
        f.write(f"# Total unresolved names: {len(failed_entries)}\n")
        f.write(f"# Total entries affected: {total_failed}\n")
        f.write("# " + "="*70 + "\n\n")
        for entry in failed_entries:
            f.write(f"Name: {entry['name']}\n")
            f.write(f"  Entries: {entry['sequence_count']}\n")
            f.write(f"  Reason: {entry['reason']}\n\n")
    print(f"Unresolved names written to: {failed_txt}")

    print(
        f"Summary: resolved species={len(species_df)}, "
        f"resolved genera={len(genus_df)}, "
        f"failed names={len(failed_entries)}"
    )

    return unified_gt, failed_entries

# Kraken2 report style parser

def load_kreport_by_rank(path, rank="S", isKraken=True):
    """
    Load a Kraken-style kreport file and return a rank-specific table
    with raw + canonical taxon fields and a standardized value/value_type.

    Parameters
    ----------
    path : str
        File path to the kreport file.
    rank : str
        Single-letter rank code (e.g., "S" for species, "G" for genus).
    isKraken : bool
        If True  -> abundance column is % abundance (value_type="abund")
        If False -> abundance column is fraction of reads (value_type="fraction_reads")

    Notes
    -----
    - Does NOT add tool/db/project/sample metadata.
    - No dropping except true duplicates (not aggregated here).
    - abundance_raw is kept exactly as reported by the tool.
    """

    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["abundance", "numReads", "numDirectReads", "rank", "taxid", "name"],
        usecols=[0, 1, 2, 3, 4, 5],
    )

    # Filter to rank (S/G/etc.)
    target = str(rank).upper()
    df["rank"] = df["rank"].astype(str)
    df = df[df["rank"].str.upper().eq(target)]

    # Raw vs canonical taxid
    df["taxid_raw"] = df["taxid"]
    df["taxid"] = (
        pd.to_numeric(df["taxid"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # Raw vs canonical name
    df["name_raw"] = df["name"].astype(str)
    df["name"] = df["name_raw"].map(_clean_name)

    # Parse native numeric fields
    df["abundance_raw"] = (
        pd.to_numeric(df["abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )
    df["numReads"] = pd.to_numeric(df["numReads"], errors="coerce").fillna(0).astype("Int64")
    df["numDirectReads"] = pd.to_numeric(df["numDirectReads"], errors="coerce").fillna(0).astype("Int64")

    # Standardized output fields
    df["value"] = df["abundance_raw"]
    if isKraken:
        df["value_type"] = "abund"            # % abundance
    else:
        df["value_type"] = "fraction_reads"   # fraction of reads

    return df[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",
            "value",
            "value_type",
            "abundance_raw",
            "numReads",
            "numDirectReads",
        ]
    ].reset_index(drop=True)

# Centrifuge report parser

def parse_centrifuge_report_ete3(filepath, rank="species", ncbi=None):
    """
    Clean Centrifuge report parser using a lineage-based projection.

    Raw Centrifuge columns (TSV with header):
      - name
      - taxID
      - taxRank
      - genomeSize
      - numReads
      - numUniqueReads
      - abundance   (relative abundance; scale kept as-is)

    Logic
    -----
    1) Read ALL rows from the Centrifuge report.
    2) For each row's taxid, use ete3 to project it to the requested rank:
         project_taxid_to_rank(taxid, rank, ncbi)
       where rank ∈ {"species", "genus"}.
    3) Drop rows whose taxid cannot be projected to that rank
       (e.g. truly rank-less leaf nodes with no species/genus ancestor,
        or invalid taxid <= 0).
    4) Group by the projected taxid and sum:
         - abundance_raw
         - numReads
         - numUniqueReads
       and take first genomeSize.
       This means:
         - species/genus rows in the original report get their own
           abundance plus all children (strain, subspecies, leaf).
         - multiple children under the same parent are merged.
    5) Resolve canonical names via NCBI for the projected taxids.

    Returns
    -------
    DataFrame with NO tool/db/sample metadata yet, at the requested rank:

        taxid_raw, taxid,
        name_raw, name,
        rank,
        value, value_type,
        abundance_raw,
        genomeSize, numReads, numUniqueReads
    """

    if ncbi is None:
        raise ValueError("parse_centrifuge_report_ete3 requires ncbi=NCBITaxa(...)")

    rank = str(rank).lower()
    if rank not in ("species", "genus"):
        raise ValueError("rank must be 'species' or 'genus'")

    # 1. Read & normalize columns
    df = pd.read_csv(filepath, sep="\t")

    df = df.rename(columns={
        "taxID": "taxid",
        "name": "name",
        "taxRank": "taxRank",
        "genomeSize": "genomeSize",
        "numReads": "numReads",
        "numUniqueReads": "numUniqueReads",
        "abundance": "abundance",
    })

    # Keep essential columns; assume they all exist
    df = df[["taxid", "name", "taxRank", "genomeSize", "numReads", "numUniqueReads", "abundance"]].copy()

    # 2. Numeric cleanup
    df["taxid"] = (
        pd.to_numeric(df["taxid"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["genomeSize"] = (
        pd.to_numeric(df["genomeSize"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["numReads"] = (
        pd.to_numeric(df["numReads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["numUniqueReads"] = (
        pd.to_numeric(df["numUniqueReads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["abundance_raw"] = (
        pd.to_numeric(df["abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )

    # 3. Project each taxid to requested rank
    # Use the public helper from util.py:
    #   project_taxid_to_rank(taxid, target_rank, ncbi) -> taxid or None
    df["projected_taxid"] = df["taxid"].apply(
        lambda t: project_taxid_to_rank(t, rank, ncbi)
    )

    # Drop rows where we couldn't find a parent at that rank
    df_valid = df[df["projected_taxid"].notna()].copy()
    if df_valid.empty:
        # Return empty frame with correct columns
        return pd.DataFrame(
            columns=[
                "taxid_raw",
                "taxid",
                "name_raw",
                "name",
                "rank",
                "value",
                "value_type",
                "abundance_raw",
                "genomeSize",
                "numReads",
                "numUniqueReads",
            ]
        )

    df_valid["projected_taxid"] = df_valid["projected_taxid"].astype("Int64")

    # 4. Aggregate by projected taxid at rank
    grouped = (
        df_valid
        .groupby("projected_taxid", as_index=False)
        .agg({
            "abundance_raw": "sum",
            "genomeSize": "first",
            "numReads": "sum",
            "numUniqueReads": "sum",
        })
    )

    grouped = grouped.rename(columns={"projected_taxid": "taxid"})
    grouped["taxid"] = grouped["taxid"].astype("Int64")

    # 5. Resolve canonical names via ete3
    unique_taxids = [int(t) for t in grouped["taxid"].dropna().tolist() if int(t) > 0]
    name_map = {}
    if unique_taxids:
        try:
            # dict: taxid -> scientific name
            name_map = ncbi.get_taxid_translator(unique_taxids)
        except Exception:
            name_map = {}

    def _name_for_tid(t):
        try:
            t_int = int(t)
        except (TypeError, ValueError):
            return "Unknown"
        return name_map.get(t_int, f"taxid_{t_int}")

    grouped["name_raw"] = grouped["taxid"].apply(_name_for_tid)
    grouped["name"] = grouped["name_raw"].astype(str).map(_clean_name)

    # 6. Standardize schema
    # taxid_raw = final aggregated taxid
    grouped["taxid_raw"] = grouped["taxid"]

    grouped["value"] = grouped["abundance_raw"]*100
    grouped["value_type"] = "abund"
    grouped["rank"] = rank  # "species" or "genus"

    grouped = grouped[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",
            "value",
            "value_type",
            "abundance_raw",
            "genomeSize",
            "numReads",
            "numUniqueReads",
        ]
    ].reset_index(drop=True)

    return grouped

# Centrifuger report parser

def parse_centrifuger_report(filepath, rank="species", ncbi=None):
    """
    Parse a Centrifuger TSV report into a standardized per-rank table.

    Original Centrifuger columns (TSV with header):
      - name
      - taxID
      - taxRank        (e.g. 'species', 'genus')
      - genomeSize
      - numReads
      - numUniqueReads
      - abundance      (relative abundance; scale kept as-is)

    Rank handling: Centrifuger's `taxRank` column is unreliable — it re-derives
    ranks through a stale internal enum, so NCBI "species group" nodes print as
    "species" (see doc/centrifuger_bug_note.md, mourisl/centrifuger#80). We ignore
    it and re-rank every taxid against NCBI/ete3, keeping only nodes whose canonical
    rank == the requested rank. Because Centrifuger reports *cumulative* (clade)
    abundance, a species node already contains its strains, so we take its value
    as-is and do NOT sum descendants (unlike Centrifuge, which is per-taxon and
    rolls descendants up in parse_centrifuge_report_ete3). `taxRank_raw` is kept
    for provenance.

    Aggregates true duplicates (same taxid + name); does NOT add tool/db/sample metadata.
    """
    if ncbi is None:
        raise ValueError("parse_centrifuger_report requires ncbi=NCBITaxa(...)")

    df = pd.read_csv(filepath, sep="\t")

    # Standardize names
    df = df.rename(columns={
        "taxID": "taxid",
        "taxRank": "taxRank",
        "genomeSize": "genomeSize",
        "numReads": "numReads",
        "numUniqueReads": "numUniqueReads",
        "abundance": "abundance",
        "name": "name",
    })

    # Keep Centrifuger's raw rank label for provenance, but do not trust it.
    df["taxRank_raw"] = df["taxRank"].astype(str)

    # --- Raw vs canonical taxid ---
    df["taxid_raw"] = df["taxid"]
    df["taxid"] = (
        pd.to_numeric(df["taxid"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # Canonicalize retired taxids via merged.dmp before re-ranking — ete3's get_rank
    # does not auto-apply it, so a retired node (common in the default-DB report)
    # would resolve to nothing and be dropped, e.g. a retired genus losing its
    # cumulative mass (Σgenus < Σspecies). taxid_raw keeps the original id.
    merged = canonicalize_merged(df["taxid"].dropna().tolist(), ncbi)
    if merged:
        df["taxid"] = df["taxid"].map(
            lambda t: merged.get(int(t), int(t)) if pd.notna(t) else t
        ).astype("Int64")

    # Re-rank via NCBI and keep only nodes whose canonical rank is the target
    # (species-group / above-species nodes drop out; strains re-rank below and
    # drop out — the parent species' cumulative value already includes them).
    target = str(rank).lower()
    ids = [int(t) for t in df["taxid"].dropna().unique() if int(t) > 0]
    canon = ncbi.get_rank(ids) if ids else {}
    df = df[df["taxid"].map(
        lambda t: pd.notna(t) and int(t) > 0 and canon.get(int(t)) == target
    )].copy()
    df["rank"] = target

    # --- Raw vs canonical name ---
    df["name_raw"] = df["name"].astype(str)
    df["name"] = df["name_raw"].map(_clean_name)

    # --- Numeric cleanup ---
    df["genomeSize"] = (
        pd.to_numeric(df["genomeSize"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["numReads"] = (
        pd.to_numeric(df["numReads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["numUniqueReads"] = (
        pd.to_numeric(df["numUniqueReads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # Abundance (raw, keep Centrifuger’s scale as-is)
    df["abundance_raw"] = (
        pd.to_numeric(df["abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )
    df.drop(columns=["abundance"], inplace=True)

    # --- Standardized main numeric field ---
    df["value"] = df["abundance_raw"]*100
    df["value_type"] = "abund"

    # --- Aggregate true duplicates (same taxid + name) ---
    # This respects the "no dropping unless true duplicates" rule.
    df = (
        df.groupby(["taxid", "name"], as_index=False)
          .agg({
              "taxid_raw": "first",
              "name_raw": "first",
              "rank": "first",
              "value": "sum",
              "value_type": "first",
              "abundance_raw": "sum",
              "genomeSize": "first",
              "numReads": "sum",
              "numUniqueReads": "sum",
              "taxRank_raw": "first",
          })
    )

    # Return standardized schema (no tool/db/sample yet)
    return df[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",
            "value",
            "value_type",
            "abundance_raw",
            "genomeSize",
            "numReads",
            "numUniqueReads",
            "taxRank_raw",
        ]
    ].reset_index(drop=True)

# Ganon2 report parser

def load_ganon_tre(path, rank="species", reads=False):
    """
    Parse a ganon2 .tre report file and return a standardized per-rank table.

    ganon2 .tre format (tab-separated):
        1: rank           (e.g. 'species', 'genus', 'phylum', 'unclassified', 'root', ...)
        2: target         (taxonomic id or specialization/assembly id)
        3: lineage        (pipe-separated taxids, e.g. 1|2|1224|...)
        4: name           (taxon name)
        5: num_unique     (reads uniquely assigned to this target)
        6: num_shared     (non-unique / LCA / re-assigned matches)
        7: num_children   (assignments to all children)
        8: num_cumulative (sum up to this node)
        9: percent_cumulative (percentage of assignments / abundance)

    Returns a DataFrame WITHOUT tool/db/sample metadata, with:
        taxid_raw, taxid, name_raw, name, rank,
        value, value_type,
        abundance_raw,
        target, lineage, num_unique, num_shared, num_children, num_cumulative
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "rank_label",        # col 1
            "target",            # col 2
            "lineage",           # col 3
            "name",              # col 4
            "num_unique",        # col 5
            "num_shared",        # col 6
            "num_children",      # col 7
            "num_cumulative",    # col 8
            "percent_cumulative" # col 9
        ],
        comment="#",
        dtype=str
    )

    # Normalize rank strings
    df["rank"] = df["rank_label"].astype(str).str.strip().str.lower()

    # Map requested rank to ganon rank label; support only species/genus here
    r = str(rank).lower()
    if r not in ("species", "genus"):
        raise ValueError(
            f"load_ganon_tre currently supports rank='species' or 'genus', got {rank!r}"
        )

    # Filter to requested rank
    df = df[df["rank"] == r].copy()

    # --- Raw vs canonical taxid ---
    # 'target' may be assembly IDs or taxids; we keep original in taxid_raw_target
    df["taxid_raw"] = df["target"]
    df["taxid"] = (
        pd.to_numeric(df["target"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # --- Raw vs canonical name ---
    df["name_raw"] = df["name"].astype(str)
    df["name"] = df["name_raw"].map(_clean_name)

    # --- Numeric cleanup for ganon2-specific fields ---
    for col in ["num_unique", "num_shared", "num_children", "num_cumulative"]:
        df[col] = (
            pd.to_numeric(df[col], errors="coerce")
            .fillna(0)
            .astype("Int64")
        )

    # percent_cumulative is already in percent units in ganon2
    df["abundance_raw"] = (
        pd.to_numeric(df["percent_cumulative"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )

    # Standardized value/value_type
    # For ganon2, abundance is this percent_cumulative
    df["value"] = df["abundance_raw"]
    if reads:
        df["value_type"] = "fraction_reads"
    else:
        df["value_type"] = "abund"

    # OPTIONAL: aggregate true duplicates (same taxid + same name)
    # This respects your "no dropping unless true duplicates" rule.
    df = (
        df.groupby(["taxid", "name"], as_index=False)
          .agg({
              "taxid_raw": "first",
              "name_raw": "first",
              "rank": "first",
              "value": "sum",
              "value_type": "first",
              "abundance_raw": "sum",
              "target": "first",
              "lineage": "first",
              "num_unique": "sum",
              "num_shared": "sum",
              "num_children": "sum",
              "num_cumulative": "sum",
          })
    )

    # Return standardized schema (no tool/db/sample yet)
    return df[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",            # 'species' or 'genus'; map to your final rank labels later if needed
            "value",
            "value_type",
            "abundance_raw",
            "target",
            "lineage",
            "num_unique",
            "num_shared",
            "num_children",
            "num_cumulative",
        ]
    ].reset_index(drop=True)

# Sourmash report parser
# - kreport parser with ete3 for missing taxids
# - default report parser with ete3

def parse_sourmash_kreport_ete3(path, ncbi=None, rank="S", verbose=False):
    """
    Parse Sourmash 6-column kreport and assign taxids via ete3.

    Columns (Sourmash kreport):
        0: %clade        -> 'abundance_raw' (already in percent)
        1: clade_reads   -> 'clade_reads'
        2: direct_reads  -> 'direct_reads'
        3: rank_code     -> 'rank' (single-letter, e.g. 'S', 'G')
        4: taxid_unused  -> 'taxid_raw' (Sourmash doesn't use this; often 0/NA)
        5: name          -> 'name_raw' / 'name'

    Returns
    -------
    DataFrame with standardized schema (NO tool/db/sample metadata):

        taxid_raw, taxid,
        name_raw, name,
        rank,
        value, value_type,
        abundance_raw,
        clade_reads, direct_reads
    """

    if ncbi is None:
        raise ValueError("parse_sourmash_kreport_ete3 requires ncbi=NCBITaxa(...)")

    # ---- Load kreport ----
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "abundance",     # % of clade
            "clade_reads",
            "direct_reads",
            "rank",
            "taxid_unused",
            "name",
        ],
        usecols=[0, 1, 2, 3, 4, 5],
    )

    # ---- Filter to requested rank (e.g. 'S' or 'G') ----
    target = str(rank).upper()
    df["rank"] = df["rank"].astype(str)
    df = df[df["rank"].str.upper() == target].copy()

    if df.empty:
        if verbose:
            print(f"[sourmash kreport:{target}] No rows for this rank in {path}")
        # Return empty frame with the standardized columns
        return pd.DataFrame(
            columns=[
                "taxid_raw",
                "taxid",
                "name_raw",
                "name",
                "rank",
                "value",
                "value_type",
                "abundance_raw",
                "clade_reads",
                "direct_reads",
            ]
        )

    # ---- Raw vs canonical name ----
    df["name_raw"] = df["name"].astype(str)
    df["name"] = df["name_raw"].map(_clean_name)

    # ---- Raw vs canonical taxid ----
    # Sourmash kreport taxid column is not reliable; keep it as taxid_raw.
    df["taxid_raw"] = (
        pd.to_numeric(df["taxid_unused"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # Map taxid via ete3 from name & rank
    if target == "S":
        rank_level = "species"
    elif target == "G":
        rank_level = "genus"
    else:
        if verbose:
            print(f"[sourmash kreport] Unsupported rank '{rank}'")
        return pd.DataFrame(
            columns=[
                "taxid_raw",
                "taxid",
                "name_raw",
                "name",
                "rank",
                "value",
                "value_type",
                "abundance_raw",
                "clade_reads",
                "direct_reads",
            ]
        )

    df["taxid"] = df["name"].apply(
        lambda n: get_taxid_from_name(n, rank_level, ncbi, clean=True) if isinstance(n, str) else 0
    )
    df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").fillna(0).astype("Int64")

    # ---- Numeric fields ----
    # abundance is already percent in Sourmash kreport
    df["abundance_raw"] = (
        pd.to_numeric(df["abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )
    df["clade_reads"] = (
        pd.to_numeric(df["clade_reads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )
    df["direct_reads"] = (
        pd.to_numeric(df["direct_reads"], errors="coerce")
        .fillna(0)
        .astype("Int64")
    )

    # ---- Standardized value/value_type ----
    # For detection, we treat %clade as the primary numeric channel.
    df["value"] = df["abundance_raw"]
    df["value_type"] = "abund"   # percent of clade

    # ---- Aggregate true duplicates only (same taxid + same name) ----
    out = (
        df.groupby(["taxid", "name"], as_index=False)
          .agg({
              "taxid_raw": "first",
              "name_raw": "first",
              "rank": "first",
              "value": "sum",
              "value_type": "first",
              "abundance_raw": "sum",
              "clade_reads": "sum",
              "direct_reads": "sum",
          })
    )

    if verbose:
        total = float(out["abundance_raw"].sum())
        unmapped = int((out["taxid"] == 0).sum())
        print(
            f"[sourmash kreport:{target}] rows={len(out)}  "
            f"sum={total:.3f}%  unmapped={unmapped}"
        )

    return out[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",          # still 'S' / 'G'; normalize later if you want
            "value",
            "value_type",
            "abundance_raw",
            "clade_reads",
            "direct_reads",
        ]
    ].reset_index(drop=True)

def parse_sourmash_report_ete3(filepath, rank="species", ncbi=None, use_weighted=True):
    """
    Clean, strict Sourmash summary parser for the Bakeoff pipeline.
    Assumes all expected columns are present (no defensive if/else blocks).
    """

    if ncbi is None:
        raise ValueError("parse_sourmash_report_ete3 requires ncbi=NCBITaxa(...)")

    df = pd.read_csv(filepath)

    # ----- Filter rank -----
    target = rank.lower()
    df["rank"] = df["rank"].astype(str).str.lower()
    df = df[df["rank"] == target].copy()

    # ----- Extract name from lineage -----
    df["lineage_raw"] = df["lineage"].astype(str)
    df["name_raw"] = df["lineage_raw"].apply(lambda s: s.split(";")[-1])
    df["name"] = df["name_raw"].map(_clean_name)

    # ----- Raw numeric fields -----
    df["fraction_raw"] = df["fraction"].astype(float)
    df["f_weighted_raw"] = df["f_weighted_at_rank"].astype(float)
    df["bp_match_at_rank"] = df["bp_match_at_rank"].astype("Int64")
    df["total_weighted_hashes"] = df["total_weighted_hashes"].astype("Int64")
    df["query_ani_at_rank"] = pd.to_numeric(df["query_ani_at_rank"], errors="coerce")

    # ----- Choose main abundance channel -----
    if use_weighted:
        df["abundance_raw"] = df["f_weighted_raw"]
        df["value_type"] = "abund_weighted"
    else:
        df["abundance_raw"] = df["fraction_raw"]
        df["value_type"] = "abund_unweighted"

    df["value"] = df["abundance_raw"]*100

    # ----- Map name → taxid -----
    rank_level = target
    df["taxid"] = df["name"].apply(lambda n: get_taxid_from_name(n, rank_level, ncbi))
    df["taxid"] = df["taxid"].astype("Int64")

    # Sourmash has no taxid field → always NA
    df["taxid_raw"] = pd.Series([pd.NA] * len(df), dtype="Int64")

    # ----- Keep important metadata -----
    for col in ["query_name", "query_md5", "query_filename"]:
        df[col] = df[col].astype(str)

    return df[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",
            "value",
            "value_type",
            "abundance_raw",
            "lineage_raw",
            "query_name",
            "query_md5",
            "query_filename",
            "fraction_raw",
            "f_weighted_raw",
            "bp_match_at_rank",
            "query_ani_at_rank",
            "total_weighted_hashes",
        ]
    ].reset_index(drop=True)

# Sylph report parser
# - mpa style report (currently using)
# - default report (legacy)

# Sylph MPA-style report parser

def parse_sylph_mpa_ete3(
    filepath,
    rank="species",
    *,
    ncbi=None,
    reads=True,
    verbose=False,
):
    """
    Parse Sylph MPA-style profiling TSV and assign NCBI taxids via ete3,
    returning taxon-level abundances at species or genus rank.

    Expected columns (MPA-style Sylph output):
      - 'clade_name'          (MPA lineage, e.g. "p__...|g__X|s__X Y")
      - 'relative_abundance'  (taxonomic abundance)
      - 'sequence_abundance'  (sequence/read abundance)
      - optional: other columns (ignored)

    Logic (GTDB → NCBI):
      1) From 'clade_name', extract the *last* token for the requested rank:
           - rank == "species" -> last 's__...'
           - rank == "genus"   -> last 'g__...'
         Only keep rows where that token is literally the LAST token
         in the lineage (to avoid double-counting deeper ranks).
      2) Keep that token (minus the 's__'/'g__') as name_raw (GTDB-style).
      3) Convert GTDB-style name → NCBI-ish name with a simple GTDB
         normalization:
           - Enterococcus_B lactis    -> Enterococcus lactis
           - Bifidobacterium longum_D -> Bifidobacterium longum
           - Only strip trailing "_[A-Z]{1,3}" on the genus/species tokens.
      4) Use get_taxid_from_name(name_ncbi, rank_level, ncbi) to get a
         NCBI taxid (e.g. using ncbi_042024 for Sylph default).
      5) Build final 'name' as:
           - NCBI-style name if taxid > 0
           - otherwise fall back to name_raw
         then run _clean_name on it.
      6) Aggregate duplicates at (taxid, name) and choose numeric channel:
           - reads=True  -> value = sequence_abundance, value_type="fraction_reads"
           - reads=False -> value = relative_abundance, value_type="abund"

    Returns
    -------
    DataFrame with columns:
        taxid_raw, taxid,
        name_raw, name,
        rank,
        value, value_type,
        tax_abund_raw, seq_abund_raw
    """
    if ncbi is None:
        raise ValueError("parse_sylph_mpa_ete3 requires ncbi=NCBITaxa(...)")

    # 0. Small GTDB→NCBI name helper
    def _normalize_gtdb_label_to_ncbi(label: str, rank_level: str) -> str:
        """
        Convert a GTDB-style genus/species label to something closer to NCBI.

        Examples
        --------
        Enterococcus_B lactis    -> Enterococcus lactis
        Bifidobacterium longum_D -> Bifidobacterium longum

        Rules
        -----
        - Only operate on the genus/species tokens (first 1–2 words).
        - Strip GTDB suffixes like "_A", "_BC", "_AQ" at the end of a token:
            pattern: _[A-Z]{1,3}$
        """
        if not isinstance(label, str):
            return ""

        s = label.strip()
        if not s:
            return ""

        tokens = s.split()
        if not tokens:
            return ""

        # genus-level: only trust the first token
        if rank_level == "genus":
            genus = tokens[0]
            genus = re.sub(r"_[A-Z]{1,3}$", "", genus)  # strip _A, _BC, _AQ
            return genus

        # species-level: try to keep "Genus species"
        if len(tokens) >= 2:
            genus_tok = re.sub(r"_[A-Z]{1,3}$", "", tokens[0])
            species_tok = re.sub(r"_[A-Z]{1,3}$", "", tokens[1])
            base = f"{genus_tok} {species_tok}"
            extra = " ".join(tokens[2:])
            if extra:
                return f"{base} {extra}".strip()
            return base

        # fallback: just strip GTDB suffix on the single token
        return re.sub(r"_[A-Z]{1,3}$", "", s)

    # 1. Normalize rank input
    rank = str(rank).lower()
    if rank not in ("species", "genus"):
        raise ValueError(f"Unsupported rank '{rank}', must be 'species' or 'genus'.")

    # 2. Load MPA-style table
    df = pd.read_csv(filepath, sep="\t", comment="#")

    required_cols = {"clade_name", "relative_abundance", "sequence_abundance"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in Sylph MPA file {filepath}:\n"
            f"  Required: {sorted(required_cols)}\n"
            f"  Found:    {sorted(df.columns)}"
        )

    df["clade_name"] = df["clade_name"].astype(str)

    df["tax_abund_raw"] = (
        pd.to_numeric(df["relative_abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )
    df["seq_abund_raw"] = (
        pd.to_numeric(df["sequence_abundance"], errors="coerce")
        .fillna(0.0)
        .astype(float)
    )

    # 3. Extract target taxon
    #    (only if it's the LAST token)

    def _extract_mpa_last_target(lineage: str, target_rank: str):
        """
        From an MPA lineage string, extract the *last token* matching:
          - 's__' for species
          - 'g__' for genus

        Only return a name if that token is the LAST token in the lineage.
        Otherwise, return None (to avoid double-counting deeper levels).
        """
        if not isinstance(lineage, str):
            return None

        parts = [p.strip() for p in lineage.split("|") if p.strip()]
        if not parts:
            return None

        prefix = "s__" if target_rank == "species" else "g__"

        target_idx = None
        for i, tok in enumerate(parts):
            if tok.startswith(prefix):
                target_idx = i

        if target_idx is None:
            return None

        # Only accept if this rank is the LAST token
        if target_idx != len(parts) - 1:
            return None

        name_part = parts[target_idx].split("__", 1)[1].strip()
        return name_part or None

    # raw GTDB-style names (DO NOT clean here)
    df["name_raw"] = df["clade_name"].apply(
        lambda s: _extract_mpa_last_target(s, rank)
    )

    # keep only rows where we actually have the requested rank as LAST token
    df = df[df["name_raw"].notna() & (df["name_raw"].astype(str).str.len() > 0)].copy()
    if df.empty:
        if verbose:
            print(f"[sylph_mpa:{rank}] No final {rank} entries in {filepath}")
        return pd.DataFrame(
            columns=[
                "taxid_raw",
                "taxid",
                "name_raw",
                "name",
                "rank",
                "value",
                "value_type",
                "tax_abund_raw",
                "seq_abund_raw",
            ]
        )

    # 4. GTDB → NCBI-ish name, then taxid resolution
    rank_level = rank  # "species" or "genus"

    # NCBI-style name from GTDB label
    df["name_ncbi"] = df["name_raw"].apply(
        lambda n: _normalize_gtdb_label_to_ncbi(n, rank_level) if isinstance(n, str) else ""
    )

    # Map via NCBI ETE3
    df["taxid"] = df["name_ncbi"].apply(
        lambda n: get_taxid_from_name(n, rank_level, ncbi)
        if isinstance(n, str) and n
        else 0
    )
    df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").fillna(0).astype("Int64")

    # taxid_raw = resolved NCBI taxid (no native Sylph taxid)
    df["taxid_raw"] = df["taxid"]

    # 5. Build final 'name' (name_raw stays exact; name is normalized+cleaned)

    # name_raw must stay exactly as extracted (no cleaning here)
    df["name_raw"] = df["name_raw"].astype(str)

    # fallback display name for unmapped: normalize the GTDB label
    df["name_fallback"] = df["name_raw"].apply(
        lambda n: _normalize_gtdb_label_to_ncbi(n, rank_level) if isinstance(n, str) else ""
    )

    # prefer mapped NCBI-ish name when taxid > 0; otherwise use fallback
    mapped_mask = df["taxid"].notna() & (df["taxid"].astype("Int64") > 0)
    df["name"] = df["name_fallback"]
    df.loc[mapped_mask, "name"] = df.loc[mapped_mask, "name_ncbi"]

    # final cleaning pass (formatting only)
    df["name"] = df["name"].astype(str).map(_clean_name)

    df = df.drop(columns=["name_fallback"], errors="ignore")

    # 6. Aggregate duplicates
    out = (
        df.groupby(["taxid", "name"], as_index=False)
          .agg({
              "taxid_raw": "first",
              "name_raw": "first",
              "tax_abund_raw": "sum",
              "seq_abund_raw": "sum",
          })
    )

    # 7. Standardized numeric fields
    if reads:
        out["value"] = out["seq_abund_raw"]
        out["value_type"] = "fraction_reads"
    else:
        out["value"] = out["tax_abund_raw"]
        out["value_type"] = "abund"

    out["rank"] = rank_level  # "species" or "genus"

    if verbose:
        total_tax = float(out["tax_abund_raw"].sum())
        total_seq = float(out["seq_abund_raw"].sum())
        unmapped = int((out["taxid"] == 0).sum())
        print(
            f"[sylph_mpa:{rank_level}] rows={len(out)}  "
            f"sum_tax={total_tax:.3f}  sum_seq={total_seq:.3f}  unmapped={unmapped}"
        )

    # final schema
    return out[
        [
            "taxid_raw",
            "taxid",
            "name_raw",
            "name",
            "rank",
            "value",
            "value_type",
            "tax_abund_raw",
            "seq_abund_raw",
        ]
    ].reset_index(drop=True)

