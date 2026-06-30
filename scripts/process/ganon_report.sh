#!/usr/bin/env bash
# ganon_report.sh
#
# Post-process Ganon2 *.rep files into per-rank .tre reports.
# Scans <out>/<technology>/ganon2-results/*.rep and runs `ganon report` on
# each to emit the {sample}_reads.tre (or {sample}.tre for abundance) that
# the analysis notebooks consume.
#
# Usage:
#   ./ganon_report.sh --db {unified|default} \
#                     [--report <reports-root>] [--db-dir <db-root>] \
#                     [--techs pacbio,ont] [--type reads|abundance]
#
#   --db:      Which database the Ganon2 .rep files were produced against.
#   --report:   Existing reports tree containing one or more `ganon2-results/`
#              subdirs (typically `--out` from process.sh). New `_reads.tre`
#              files are written next to the `.rep` inputs.
#              Default: reports/<db>-reports.
#   --db-dir:  Parent of the per-DB index trees (default: data/ref_db). The
#              script appends `refseq03032025/` for --db unified and
#              `default_db/` for --db default to reach the Ganon2 db prefix.
#   --techs:   Comma-separated technology subdirectories to scan.
#              Default: pacbio,ont.
#   --data-groups: Comma-separated data-group subdirectories to scan (the
#              first path component under --report). Default:
#              ZymoMockD6331,simulated. DYN is opt-in.
#   --type:    'reads' (read-weighted, default) or 'abundance'. Reads-weighted
#              is what the manuscript's detection benchmarks consume.
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/ganon_report.sh --db unified
#   ./scripts/process/ganon_report.sh --db default --techs ont --type abundance
#   ./scripts/process/ganon_report.sh --db unified --data-groups ZymoMockD6331,simulated,DYN

set -euo pipefail

DB=""
WORKING_DIR=""
DATABASE_DIR="data/ref_db"
TECHS_CSV="pacbio,ont"
DATA_GROUPS_CSV="ZymoMockD6331,simulated"
REPORT_TYPE="reads"
JOBS=1

# -----------------------------------------------------------------------------
# Parse args
# -----------------------------------------------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --db)           DB="$2";              shift 2 ;;
        --report)       WORKING_DIR="$2";     shift 2 ;;
        --db-dir)       DATABASE_DIR="$2";    shift 2 ;;
        --techs)        TECHS_CSV="$2";       shift 2 ;;
        --data-groups)  DATA_GROUPS_CSV="$2"; shift 2 ;;
        --type)         REPORT_TYPE="$2";     shift 2 ;;
        --jobs|-j)      JOBS="$2";            shift 2 ;;
        -h|--help)
            sed -n '2,34p' "$0"; exit 0 ;;
        *)
            echo "Unknown parameter: $1" >&2
            echo "Usage: $0 --db {unified|default} [--report <root>] [--db-dir <db>] [--techs pacbio,ont] [--type reads|abundance]" >&2
            exit 1 ;;
    esac
done

if [[ -z "$DB" ]]; then
    echo "Error: --db is required (unified|default)." >&2; exit 1
fi
if [[ "$REPORT_TYPE" != "reads" && "$REPORT_TYPE" != "abundance" ]]; then
    echo "Error: --type must be 'reads' or 'abundance', got '${REPORT_TYPE}'." >&2
    exit 1
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
# Per-database Ganon2 index prefix. Subpath under --db-dir mirrors process.sh.
# -----------------------------------------------------------------------------
case "$DB" in
    unified)
        DB_ROOT="${DATABASE_DIR}/refseq03032025"
        GANON_DB="${DB_ROOT}/ganon2_abvf_030325/ganon2_abvf_030325"
        ;;
    default)
        DB_ROOT="${DATABASE_DIR}/default_db"
        GANON_DB="${DB_ROOT}/ganon2_default/ganon2_default_abfv_rs_cg"
        ;;
    *)
        echo "Error: --db must be 'unified' or 'default', got '${DB}'." >&2
        exit 1 ;;
esac

# Best-effort check for the DB prefix.
if ! ls "${GANON_DB}"* >/dev/null 2>&1; then
    echo "WARNING: no files found with Ganon2 db prefix ${GANON_DB}" >&2
    echo "         (proceeding anyway; ganon report will surface the real error)" >&2
fi

IFS=',' read -r -a TECHS <<< "$TECHS_CSV"
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

echo "DB:           ${DB}"
echo "Techs:        ${TECHS[*]}"
echo "Data groups:  ${DATA_GROUPS[*]}"
echo "Report type:  ${REPORT_TYPE}"
echo "Working dir:  ${WORKING_DIR}"
echo "Ganon2 db:    ${GANON_DB}"
echo "Parallel:     ${JOBS} concurrent job(s)"
echo

# -----------------------------------------------------------------------------
# parallel_throttle / parallel_wait_all — see generate_kreport.sh for notes.
# When JOBS=1 the loop is functionally sequential.
# -----------------------------------------------------------------------------
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
# Main loop: find every <data_group>/<tech>/ganon2-results/ dir under --report
# (works regardless of how many data groups are nested), filtered by --techs.
# -----------------------------------------------------------------------------
for tech in "${TECHS[@]}"; do
    mapfile -t results_dirs < <(
        find "$WORKING_DIR" -type d -path "*/${tech}/ganon2-results" 2>/dev/null | sort
    )
    if [[ "${#results_dirs[@]}" -eq 0 ]]; then
        echo "Skipping (no ganon2-results dirs under */${tech}/): ${WORKING_DIR}"
        continue
    fi

    for results_dir in "${results_dirs[@]}"; do
        # Skip groups not in --data-groups (e.g. DYN by default).
        if ! in_data_groups "$results_dir"; then continue; fi
        echo "Processing ganon2 .rep files in:"
        echo "  ${results_dir}"

        shopt -s nullglob
        rep_files=( "${results_dir}"/*.rep )
        shopt -u nullglob

        if [[ "${#rep_files[@]}" -eq 0 ]]; then
            echo "  No .rep files found."
            continue
        fi

        pids=()
        for rep in "${rep_files[@]}"; do
            base="$(basename "$rep" .rep)"
            if [[ "$REPORT_TYPE" == "reads" ]]; then
                out_prefix="${results_dir}/${base}_${REPORT_TYPE}"
            else
                out_prefix="${results_dir}/${base}"
            fi

            # Skip if the final .tre already exists.
            if [[ -e "${out_prefix}.tre" ]]; then
                echo "  Skipping existing: $(basename "${out_prefix}.tre")"
                continue
            fi

            echo "  Generating ganon2 report:"
            echo "    Input : ${rep}"
            echo "    Output: ${out_prefix}.tre"

            ganon report \
                -i "${rep}" \
                -o "${out_prefix}" \
                -d "${GANON_DB}" \
                --min-count 0 \
                -t "${REPORT_TYPE}" &
            pids+=( $! )
            parallel_throttle pids
        done
        parallel_wait_all pids
        echo
    done
done

echo "All ganon2 reports generated successfully."
