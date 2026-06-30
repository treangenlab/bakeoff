#!/usr/bin/env bash
# Fetch, verify, and extract bakeoff bundles from Zenodo into a project root.
# Restores the original `data/...` layout that the analysis scripts expect.
#
# Usage:
#   ./fetch_bakeoff_data.sh --root /path/to/pacbio_bakeoff [--bundles "name1 name2"]
#
# Options:
#   --root DIR         Extraction root. Files land at <root>/data/... (default: $PWD)
#   --bundles "a b c"  Subset to fetch (default: all bundles in the configured records).
#                      Use --list to see available bundles.
#   --list             Print available bundles from the configured Zenodo records and exit.
#   --staging DIR      Where to download parts (default: <root>/_staging)
#   --keep-staging     Don't delete downloaded parts after a successful extract
#   --skip-verify      Skip sha256 verification (faster, not recommended)
#
# Requires: curl, jq, zstd >= 1.4, tar, sha256sum.

set -euo pipefail

# ---------------------------------------------------------------------------
# Zenodo record IDs hosting the reference-database bundles.
# Records are grouped in the project community for browsing / citation:
#     https://zenodo.org/communities/treangen_lab_bakeoff
# (Read datasets are NOT fetched by this script: mock reads are pulled from SRA
# via the instructions in data/ZymoMockD6331/pacbio/mock_pacbio_fetch.md and
# the Datasets wiki page; the simulated dataset record `20496925` is cited
# from the Datasets wiki page directly.)
# ---------------------------------------------------------------------------
ZENODO_RECORDS=(
  "20481656"   # Metadata (shared_metadata + default_metadata combined)
  "20479369"   # Unified Centrifuge index
  "20479514"   # Unified Centrifuger index
  "20481847"   # Unified Ganon2 index
  "20481715"   # Unified Sourmash sketches
  "20481766"   # Unified Sylph sketches
  "20482294"   # Default Ganon2 index
  "20484934"   # Unified Kraken2 index — part 1 of 2 (both records required)
  "20496297"   # Unified Kraken2 index — part 2 of 2
)
# ---------------------------------------------------------------------------

ROOT="$PWD"
STAGING=""
KEEP_STAGING=0
SKIP_VERIFY=0
REQUESTED=""
LIST_ONLY=0

while [[ $# -gt 0 ]]; do
  case $1 in
    --root)          ROOT=$2; shift 2;;
    --bundles)       REQUESTED=$2; shift 2;;
    --list)          LIST_ONLY=1; shift;;
    --staging)       STAGING=$2; shift 2;;
    --keep-staging)  KEEP_STAGING=1; shift;;
    --skip-verify)   SKIP_VERIFY=1; shift;;
    -h|--help)       sed -n '2,20p' "$0"; exit 0;;
    *)               echo "unknown arg: $1" >&2; exit 2;;
  esac
done

for cmd in curl jq zstd tar sha256sum; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "missing dependency: $cmd" >&2; exit 1; }
done

ROOT=$(realpath -m "$ROOT")
STAGING=${STAGING:-$ROOT/_staging}

# Build the file table: for every file in every configured Zenodo record,
# record tab-separated: <key>\t<download_url>
TMPMETA=$(mktemp)
trap 'rm -f $TMPMETA' EXIT

for rid in "${ZENODO_RECORDS[@]}"; do
  if [[ $rid == PLACEHOLDER* ]]; then
    echo "ERROR: Zenodo record '$rid' is a placeholder. Edit ZENODO_RECORDS at top of script." >&2
    exit 1
  fi
  curl -fsSL "https://zenodo.org/api/records/$rid" \
    | jq -r '.files[] | "\(.key)\t\(.links.self)"' \
    >> "$TMPMETA"
done

# Discover bundles from part filenames: anything matching <bundle>.tar.zst.part-NNN
ALL_BUNDLES=$(awk -F'\t' '{print $1}' "$TMPMETA" \
  | grep -E '\.tar\.zst\.part-[0-9]+$' \
  | sed -E 's/\.tar\.zst\.part-[0-9]+$//' \
  | sort -u)

if [[ $LIST_ONLY -eq 1 ]]; then
  echo "Available bundles:"
  printf '  %s\n' $ALL_BUNDLES
  exit 0
fi

if [[ -z $REQUESTED || $REQUESTED == all ]]; then
  BUNDLES=($ALL_BUNDLES)
else
  read -ra BUNDLES <<< "$REQUESTED"
fi

mkdir -p "$ROOT" "$STAGING"

# Helper: look up download URL for a given filename
url_for() { awk -F'\t' -v k="$1" '$1==k {print $2; exit}' "$TMPMETA"; }

fetch_file() {
  local key=$1 out=$2 url
  url=$(url_for "$key") || true
  if [[ -z $url ]]; then
    return 1
  fi
  if [[ -f $out ]]; then
    echo "  cached: $key"
    return 0
  fi
  echo "  download: $key"
  curl -fL --retry 5 --retry-delay 30 --retry-all-errors -C - -o "$out" "$url"
}

for bundle in "${BUNDLES[@]}"; do
  echo "=== $bundle ==="

  # 1) Sidecars (parts.sha256 is required for verification; others are informational)
  fetch_file "${bundle}.parts.sha256"      "$STAGING/${bundle}.parts.sha256"      || true
  fetch_file "${bundle}.recombined.sha256" "$STAGING/${bundle}.recombined.sha256" || true
  fetch_file "${bundle}.manifest.txt"      "$STAGING/${bundle}.manifest.txt"      || true

  # 2) All parts for this bundle
  while IFS=$'\t' read -r key _; do
    [[ $key =~ ^${bundle}\.tar\.zst\.part-[0-9]+$ ]] || continue
    fetch_file "$key" "$STAGING/$key"
  done < "$TMPMETA"

  # 3) Verify per-part sha256. If every part matches, the concatenated stream
  #    is correct by definition, so no separate recombined-stream check needed.
  if [[ $SKIP_VERIFY -eq 0 && -f "$STAGING/${bundle}.parts.sha256" ]]; then
    echo "  verify parts.sha256"
    ( cd "$STAGING" && sha256sum -c "${bundle}.parts.sha256" )
  fi

  # 4) Stream-reassemble + decompress + extract into the project root.
  #    Archive entries are relative ("data/..."), so they land under $ROOT.
  echo "  extract -> $ROOT"
  cat "$STAGING/${bundle}.tar.zst.part-"* \
    | zstd -dc --long=27 \
    | tar -xf - -C "$ROOT"

  # 5) Cleanup parts + sidecars unless --keep-staging
  if [[ $KEEP_STAGING -eq 0 ]]; then
    rm -f "$STAGING/${bundle}.tar.zst.part-"* \
          "$STAGING/${bundle}.parts.sha256" \
          "$STAGING/${bundle}.recombined.sha256" \
          "$STAGING/${bundle}.manifest.txt"
  fi
done

# Post-extract: NCBI taxonomy snapshot.
# The shared_metadata bundle ships taxdump as a .tar.gz to save space; tools expect
# the extracted `taxdump/` directory next to it.
TAXTGZ="$ROOT/data/ref_db/refseq03032025/taxdump0303.tar.gz"
TAXDIR="$ROOT/data/ref_db/refseq03032025/taxdump"
if [[ -f $TAXTGZ && ! -d $TAXDIR ]]; then
  echo "=== Extracting NCBI taxonomy ($TAXDIR) ==="
  mkdir -p "$TAXDIR"
  tar -C "$TAXDIR" -xzf "$TAXTGZ"
fi

# Cleanup staging if empty
if [[ $KEEP_STAGING -eq 0 ]]; then
  rmdir "$STAGING" 2>/dev/null || true
fi

echo
echo "Done. Project root: $ROOT"
echo "Verify expected layout: ls $ROOT/data/"
