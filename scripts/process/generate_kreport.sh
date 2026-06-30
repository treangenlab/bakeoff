#!/usr/bin/env bash
# kreport.sh
#
# Post-process per-tool profiling outputs into Kraken-style kreport files.
# Covers Centrifuge, Centrifuger, and Sourmash. Scans
#     <working-dir>/<technology>/<Tool>-results/
# and emits the kreport next to each tool's native output.
#
# Usage:
#   ./kreport.sh --db {unified|default} \
#                --report <reports-root> \
#                --db-dir <db-root> \
#                [--tools centrifuge,centrifuger,sourmash] \
#                [--techs pacbio,ont]
#
#   --db:           Which database the reports were produced against.
#   --report:        Existing reports tree to scan and augment. New kreport
#                   files are written next to their inputs (no separate output
#                   directory). Default: reports/<db>-reports.
#   --db-dir:       Parent of the per-DB index trees (default: data/ref_db).
#                   The script appends `refseq03032025/` for --db unified and
#                   `default_db/` for --db default.
#   --tools:        Comma-separated subset of {centrifuge, centrifuger, sourmash},
#                   or "all". Default: all three.
#   --techs:        Comma-separated technology subdirectories to scan.
#                   Default: pacbio,ont.
#   --data-groups:  Comma-separated data-group subdirectories to scan
#                   (the first path component under --report, e.g.
#                   ZymoMockD6331, simulated, DYN). Default:
#                   ZymoMockD6331,simulated. DYN is opt-in.
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/generate_kreport.sh --db unified --report reports/unified-reports/ZymoMockD6331
#   ./scripts/process/generate_kreport.sh --db default --tools sourmash --techs ont
#   ./scripts/process/generate_kreport.sh --db unified --data-groups ZymoMockD6331,simulated,DYN

set -euo pipefail

DB=""
WORKING_DIR=""
DATABASE_DIR="data/ref_db"
TOOLS="all"
TECHS_CSV="pacbio,ont"
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
        --tools)        TOOLS="$2";           shift 2 ;;
        --techs)        TECHS_CSV="$2";       shift 2 ;;
        --data-groups)  DATA_GROUPS_CSV="$2"; shift 2 ;;
        --jobs|-j)      JOBS="$2";            shift 2 ;;
        -h|--help)
            sed -n '2,36p' "$0"; exit 0 ;;
        *)
            echo "Unknown parameter: $1" >&2
            echo "Usage: $0 --db {unified|default} --report <root> --db-dir <db> [--tools centrifuge,centrifuger,sourmash] [--techs pacbio,ont]" >&2
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
# Per-database paths. Subpaths under --db-dir mirror the process.sh
# layout. Update here if your index dir basenames differ from the manuscript's.
# -----------------------------------------------------------------------------
case "$DB" in
    unified)
        DB_ROOT="${DATABASE_DIR}/refseq03032025"
        CENTRIFUGE_DB="${DB_ROOT}/centrifuge_abv_030325/refseq_abv"
        CENTRIFUGER_DB="${DB_ROOT}/centrifuger_abv_030325/refseq_abv"
        SOURMASH_TAX_CSV="${DB_ROOT}/sourmash_abvf_030325/lineage_030325.csv"
        ;;
    default)
        DB_ROOT="${DATABASE_DIR}/default_db"
        CENTRIFUGE_DB="${DB_ROOT}/cf_default/hpvc"
        CENTRIFUGER_DB="${DB_ROOT}/cfer_default/hpv_gbsarscov2"
        # Sourmash default uses a merged lineage CSV at kreport-stage rather
        # than the per-domain CSVs used by `process.sh`.
        SOURMASH_TAX_CSV="${DB_ROOT}/sm_default/sourmash_lineage_default_merged.csv"
        ;;
    *)
        echo "Error: --db must be 'unified' or 'default', got '${DB}'." >&2
        exit 1 ;;
esac

# -----------------------------------------------------------------------------
# Tool + tech selection
# -----------------------------------------------------------------------------
ALL_TOOLS=(centrifuge centrifuger sourmash)

if [[ "$TOOLS" == "all" ]]; then
    selected_tools=( "${ALL_TOOLS[@]}" )
else
    IFS=',' read -r -a selected_tools <<< "$TOOLS"
fi

IFS=',' read -r -a TECHS <<< "$TECHS_CSV"
IFS=',' read -r -a DATA_GROUPS <<< "$DATA_GROUPS_CSV"

# Return 0 if the first path component of $1 under $WORKING_DIR matches a
# member of DATA_GROUPS; 1 otherwise. Used to skip groups (e.g. DYN by default).
in_data_groups() {
    local path="${1#${WORKING_DIR%/}/}"
    local group="${path%%/*}"
    local g
    for g in "${DATA_GROUPS[@]}"; do
        [[ "$g" == "$group" ]] && return 0
    done
    return 1
}

echo "DB:           ${DB}"
echo "Tools:        ${selected_tools[*]}"
echo "Techs:        ${TECHS[*]}"
echo "Data groups:  ${DATA_GROUPS[*]}"
echo "Working dir:  ${WORKING_DIR}"
echo "Database dir: ${DATABASE_DIR}"
echo "Parallel:     ${JOBS} concurrent job(s)"

# -----------------------------------------------------------------------------
# parallel_throttle <pid_var>
# Helper used inside per-sample loops. After launching a background command,
# call this to block while $JOBS jobs are already in flight. Implemented with
# `wait -n` so we proceed the instant any one finishes.
#   Usage:
#       cmd_for_one "$item" & active_pids+=( $! )
#       parallel_throttle active_pids
# When JOBS=1 the loop is functionally sequential (launch, wait, launch).
# -----------------------------------------------------------------------------
parallel_throttle() {
    local -n _pids=$1
    while [[ "${#_pids[@]}" -ge "$JOBS" ]]; do
        wait -n "${_pids[@]}" 2>/dev/null || true
        # Rebuild the pid list, dropping any that have exited.
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
# Per-tool drivers. Each takes a results directory and emits kreport(s) into it.
# -----------------------------------------------------------------------------

# Centrifuge: iterate *.tsv, skip *_report.tsv / *_kreport.tsv, run
# centrifuge-kreport against --db-dir's centrifuge index.
do_centrifuge() {
    local results_dir="$1"
    shopt -s nullglob
    local files=( "${results_dir}"/*.tsv )
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
        echo "  No TSV files found."
        return
    fi

    local -a pids=()
    for tsv in "${files[@]}"; do
        local base; base="$(basename "$tsv")"
        case "$base" in
            *_report.tsv|*_kreport.tsv) continue ;;
        esac
        local sample="${base%.tsv}"
        local kreport="${results_dir}/${sample}_kreport.tsv"
        echo "  Generating kreport:"
        echo "    Input : ${tsv}"
        echo "    Output: ${kreport}"
        centrifuge-kreport -x "${CENTRIFUGE_DB}" "${tsv}" > "${kreport}" &
        pids+=( $! )
        parallel_throttle pids
    done
    parallel_wait_all pids
}

# Centrifuger: iterate *_results.tsv, skip if kreport already exists.
do_centrifuger() {
    local results_dir="$1"
    shopt -s nullglob
    local files=( "${results_dir}"/*_results.tsv )
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
        echo "  No *_results.tsv files found."
        return
    fi

    local -a pids=()
    for tsv in "${files[@]}"; do
        local sample; sample="$(basename "$tsv" _results.tsv)"
        local kreport="${results_dir}/${sample}_kreport.tsv"
        if [[ -f "$kreport" ]]; then
            echo "  Skipping existing kreport: $(basename "$kreport")"
            continue
        fi
        echo "  Generating kreport:"
        echo "    Input : ${tsv}"
        echo "    Output: ${kreport}"
        centrifuger-kreport -x "${CENTRIFUGER_DB}" "${tsv}" > "${kreport}" &
        pids+=( $! )
        parallel_throttle pids
    done
    parallel_wait_all pids
}

# Sourmash: iterate *-gather.csv, run `sourmash tax metagenome -F kreport`.
do_sourmash() {
    local results_dir="$1"

    if [[ ! -f "$SOURMASH_TAX_CSV" ]]; then
        echo "  ERROR: SOURMASH_TAX_CSV does not exist: ${SOURMASH_TAX_CSV}" >&2
        echo "  (Build/merge the per-domain lineage CSVs into one before re-running.)" >&2
        return 1
    fi

    shopt -s nullglob
    local files=( "${results_dir}"/*-gather.csv )
    shopt -u nullglob

    if [[ "${#files[@]}" -eq 0 ]]; then
        echo "  No *-gather.csv files found."
        return
    fi

    local -a pids=()
    for gather in "${files[@]}"; do
        local out_prefix="${gather%-gather.csv}"
        if compgen -G "${out_prefix}"*kreport* >/dev/null; then
            echo "  Skipping (kreport exists): $(basename "$out_prefix")*kreport*"
            continue
        fi
        echo "  Generating kreport:"
        echo "    Input  : ${gather}"
        echo "    Prefix : ${out_prefix}"
        sourmash tax metagenome \
            -t "${SOURMASH_TAX_CSV}" \
            -g "${gather}" \
            -o "${out_prefix}" \
            -F kreport &
        pids+=( $! )
        parallel_throttle pids
    done
    parallel_wait_all pids
}

# -----------------------------------------------------------------------------
# Main loop: tool × technology
# -----------------------------------------------------------------------------
declare -A TOOL_DIR=(
    [centrifuge]="Centrifuge-results"
    [centrifuger]="Centrifuger-results"
    [sourmash]="sourmash-results"
)

for tool in "${selected_tools[@]}"; do
    sub="${TOOL_DIR[$tool]:-}"
    if [[ -z "$sub" ]]; then
        echo "Unknown tool: $tool (valid: centrifuge, centrifuger, sourmash)" >&2
        continue
    fi

    for tech in "${TECHS[@]}"; do
        # Find every <data_group>/<tech>/<Tool-results>/ dir under --report
        # (works regardless of how many data groups are nested).
        mapfile -t results_dirs < <(
            find "$WORKING_DIR" -type d -path "*/${tech}/${sub}" 2>/dev/null | sort
        )
        if [[ "${#results_dirs[@]}" -eq 0 ]]; then
            echo "Skipping (no ${sub} dirs under */${tech}/): $WORKING_DIR"
            continue
        fi

        for results_dir in "${results_dirs[@]}"; do
            # Skip groups not in --data-groups (e.g. DYN by default).
            if ! in_data_groups "$results_dir"; then continue; fi
            echo
            echo "Processing ${tool} results in:"
            echo "  ${results_dir}"

            case "$tool" in
                centrifuge)  do_centrifuge  "$results_dir" ;;
                centrifuger) do_centrifuger "$results_dir" ;;
                sourmash)    do_sourmash    "$results_dir" ;;
            esac
        done
    done
done

echo
echo "All kreport outputs generated successfully."
