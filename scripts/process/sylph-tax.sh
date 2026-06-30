#!/usr/bin/env bash
# sylph-tax.sh
#
# Bulk sylph-tax taxprof runner. For every `*/sylph-results/` directory under
# --report, runs `sylph-tax taxprof *_profiling.tsv -t <DBS>` to convert
# Sylph's native profile into a taxonomy-annotated table.
#
# Usage:
#   ./sylph-tax.sh --db {unified|default} --report <root> [--db-dir <db-root>]
#
#   --db:           Which database the Sylph profiles were produced against.
#                     unified -> use a metadata TSV under --db-dir.
#                     default -> use sylph-tax's built-in named DBs
#                                (GTDB_r220, IMGVR_4.1, FungiRefSeq-2024-07-25).
#   --report:        Existing reports tree containing one or more
#                   `sylph-results/` subdirs (typically `--out` from process.sh).
#                   New `.sylphmpa` files are written next to their inputs.
#                   Default: reports/<db>-reports.
#   --db-dir:       Parent of the per-DB index trees (default: data/ref_db).
#                   For --db unified, the script appends `refseq03032025/` to
#                   reach the sylph-tax metadata TSV. Ignored for --db default.
#   --data-groups:  Comma-separated data-group subdirectories to scan (the
#                   first path component under --report). Default:
#                   ZymoMockD6331,simulated. DYN is opt-in.
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/sylph-tax.sh --db unified
#   ./scripts/process/sylph-tax.sh --db default
#   ./scripts/process/sylph-tax.sh --db unified --data-groups ZymoMockD6331,simulated,DYN

set -euo pipefail

DB=""
WORKING_DIR=""
DATABASE_DIR="data/ref_db"
DATA_GROUPS_CSV="ZymoMockD6331,simulated"
JOBS=1

# -----------------------------------------------------------------------------
# Parse args
# -----------------------------------------------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --db)           DB="$2";              shift 2 ;;
        --report)       WORKING_DIR="$2";     shift 2 ;;
        --db-dir)       DATABASE_DIR="$2";    shift 2 ;;
        --data-groups)  DATA_GROUPS_CSV="$2"; shift 2 ;;
        --jobs|-j)      JOBS="$2";            shift 2 ;;
        -h|--help)
            sed -n '2,26p' "$0"; exit 0 ;;
        *)
            echo "Unknown parameter: $1" >&2
            echo "Usage: $0 --db {unified|default} --report <root> [--db-dir <db>]" >&2
            exit 1 ;;
    esac
done

if [[ -z "$DB" ]]; then
    echo "Error: --db is required (unified|default)." >&2; exit 1
fi

# Default --report is per-DB so unified and default outputs never collide.
if [[ -z "$WORKING_DIR" ]]; then
    WORKING_DIR="reports/${DB}-reports"
fi

if [[ ! -d "$WORKING_DIR" ]]; then
    echo "ERROR: --report does not exist: $WORKING_DIR" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Per-DB taxonomy spec. SYLPH_TAX_ARGS is what gets passed after `-t`.
# Subpath under --db-dir is editable here if your layout differs.
# -----------------------------------------------------------------------------
case "$DB" in
    unified)
        # Note: the unified-DB sylph-tax metadata is conventionally stored
        # under <db-dir>/default_db/sylph_tax/ (yes, default_db/ — the file
        # was placed there by the upstream build, despite being the unified
        # DB's metadata). Move/symlink it elsewhere if your layout differs.
        SYLPH_TAX_DB="${DATABASE_DIR}/default_db/sylph_tax/refseq030325_taxonomy_metadata.tsv"
        if [[ ! -f "$SYLPH_TAX_DB" ]]; then
            echo "ERROR: Sylph taxonomy metadata not found: $SYLPH_TAX_DB" >&2
            exit 1
        fi
        SYLPH_TAX_ARGS=( "$SYLPH_TAX_DB" )
        ;;
    default)
        # Built-in named DBs that sylph-tax downloads/manages itself.
        # --db-dir is ignored in this case.
        SYLPH_TAX_ARGS=( GTDB_r220 IMGVR_4.1 FungiRefSeq-2024-07-25 )
        ;;
    *)
        echo "Error: --db must be 'unified' or 'default', got '${DB}'." >&2
        exit 1 ;;
esac

IFS=',' read -r -a DATA_GROUPS <<< "$DATA_GROUPS_CSV"

# Return 0 if the first path component of $1 under $WORKING_DIR matches a
# member of DATA_GROUPS; 1 otherwise. Skips groups not in --data-groups.
in_data_groups() {
    local path="${1#${WORKING_DIR%/}/}"
    local group="${path%%/*}"
    local g
    for g in "${DATA_GROUPS[@]}"; do
        [[ "$g" == "$group" ]] && return 0
    done
    return 1
}

echo "=== sylph-tax taxprof runner ==="
echo "DB:           ${DB}"
echo "Tax DB(s):    ${SYLPH_TAX_ARGS[*]}"
echo "Data groups:  ${DATA_GROUPS[*]}"
echo "Searching for 'sylph-results' directories under: ${WORKING_DIR}"
echo "Parallel:     ${JOBS} concurrent job(s)"
echo

# parallel_throttle / parallel_wait_all — see generate_kreport.sh for notes.
parallel_throttle() {
    local -n _pids=$1
    while [[ "${#_pids[@]}" -ge "$JOBS" ]]; do
        wait -n "${_pids[@]}" 2>/dev/null || true
        local still=()
        for p in "${_pids[@]}"; do
            if kill -0 "$p" 2>/dev/null; then still+=( "$p" ); fi
        done
        _pids=( "${still[@]}" )
    done
}
parallel_wait_all() {
    local -n _pids=$1
    if [[ "${#_pids[@]}" -gt 0 ]]; then
        wait "${_pids[@]}" 2>/dev/null || true
    fi
    _pids=()
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
mapfile -t SYLPH_DIRS < <(find "$WORKING_DIR" -type d -name "sylph-results" | sort)

if [[ ${#SYLPH_DIRS[@]} -eq 0 ]]; then
    echo "No sylph-results directories found under $WORKING_DIR."
    exit 0
fi

shopt -s nullglob

# One sylph-tax invocation per directory (sylph-tax taxprof already accepts a
# batch of profile files at once). Parallelism is therefore across directories.
pids=()
for d in "${SYLPH_DIRS[@]}"; do
    # Skip groups not in --data-groups (e.g. DYN by default).
    if ! in_data_groups "$d"; then continue; fi
    # Pre-check inside the parent shell so empty dirs are skipped without
    # spawning a background subshell.
    prof_files=( "$d"/*_profiling.tsv )
    if [[ ${#prof_files[@]} -eq 0 ]]; then
        echo "[INFO] $d  — no *_profiling.tsv files, skipping."
        continue
    fi

    echo "[INFO] $d  — ${#prof_files[@]} profiling file(s) queued"

    # Subshell so each backgrounded job has its own cwd (and we don't race on
    # the parent's). sylph-tax writes outputs next to its inputs by default.
    (
        cd "$d"
        local_files=(*_profiling.tsv)
        if sylph-tax taxprof "${local_files[@]}" -t "${SYLPH_TAX_ARGS[@]}"; then
            echo "  DONE: sylph-tax completed in $d"
        else
            echo "  FAILED: sylph-tax in $d" >&2
        fi
    ) &
    pids+=( $! )
    parallel_throttle pids
done
parallel_wait_all pids

echo
echo "=== All done ==="
