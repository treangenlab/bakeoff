#!/usr/bin/env bash
# postprocess.sh
#
# Convenience wrapper that runs the three post-processing scripts in sequence
# so the typical user-facing workflow after `process.sh` is a single command:
#
#     process.sh        -> raw per-tool outputs in <out>/<data_group>/<tech>/<Tool>-{results,logs}/
#     postprocess.sh    -> kreport / .tre / .sylphmpa files the notebooks consume
#     <analysis notebooks>
#
# This script does NOT implement any post-processing logic itself. It is a
# thin orchestrator that invokes:
#
#   1) generate_kreport.sh  — Centrifuge, Centrifuger, Sourmash -> kreport
#   2) ganon_report.sh      — Ganon2 .rep -> .tre (reads-weighted)
#   3) sylph-tax.sh         — Sylph profile -> taxonomy-annotated .sylphmpa
#
# Each of those three remains usable standalone for partial / per-tool runs.
# postprocess.sh forwards --db, --report, --db-dir, --techs to all three and
# CONTINUES on per-step failure, so a missing tool binary in one sub-script
# does not block the others.
#
# Usage:
#   ./postprocess.sh --db {unified|default} \
#                    [--report <reports-root>] [--db-dir <db-root>] \
#                    [--techs pacbio,ont]
#
#   --db:      Which database the reports were produced against.
#   --report:   Existing reports tree to scan + augment. Each sub-script reads
#              tool outputs already there and writes derived files NEXT TO
#              the inputs (no separate output dir).
#              Default: reports/<db>-reports (matches process.sh --out default).
#   --db-dir:  Parent of the per-DB index trees. Default: data/ref_db. The
#              sub-scripts append `refseq03032025/` (unified) or `default_db/`
#              (default) to reach per-tool indexes / lineage tables.
#   --techs:   Comma-separated technology subdirectories to scan. Default:
#              pacbio,ont. Forwarded to all three sub-scripts.
#   --data-groups:
#              Comma-separated data-group subdirectories to scan (the first
#              path component under --report). Default:
#              ZymoMockD6331,simulated. DYN is opt-in — add it explicitly to
#              postprocess the cohort. Forwarded to all three sub-scripts.
#   --jobs N (-j N):
#              Per-tool parallelism: up to N concurrent invocations of the
#              underlying tool inside each sub-script. Forwarded to all three
#              sub-scripts. Default: 1 (sequential — current behavior).
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/postprocess.sh --db unified
#   ./scripts/process/postprocess.sh --db default --techs ont
#   ./scripts/process/postprocess.sh --db unified --jobs 4   # 4-way parallel
#   ./scripts/process/postprocess.sh --db unified --data-groups ZymoMockD6331,simulated,DYN
#   ./scripts/process/postprocess.sh --db unified \
#       --report reports/unified-reports \
#       --db-dir /scratch/refs \
#       --jobs 8

set -u

DB=""
WORKING_DIR=""
DATABASE_DIR=""
TECHS_CSV=""
DATA_GROUPS_CSV=""
JOBS=""

# -----------------------------------------------------------------------------
# Parse args. All flags are passed through to the sub-scripts so the wrapper
# has no per-script translation step.
# -----------------------------------------------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --db)           DB="$2";              shift 2 ;;
        --report)       WORKING_DIR="$2";     shift 2 ;;
        --db-dir)       DATABASE_DIR="$2";    shift 2 ;;
        --techs)        TECHS_CSV="$2";       shift 2 ;;
        --data-groups)  DATA_GROUPS_CSV="$2"; shift 2 ;;
        --jobs|-j)      JOBS="$2";            shift 2 ;;
        -h|--help)
            sed -n '2,51p' "$0"; exit 0 ;;
        *)
            echo "Unknown parameter: $1" >&2
            echo "Usage: $0 --db {unified|default} [--report <root>] [--db-dir <db>] [--techs pacbio,ont] [--data-groups ZymoMockD6331,simulated]" >&2
            exit 1 ;;
    esac
done

if [[ -z "$DB" ]]; then
    echo "Error: --db is required (unified|default)." >&2; exit 1
fi

# -----------------------------------------------------------------------------
# Build the argv list once and reuse for each sub-script.
# Only forward flags the user actually set; the sub-scripts have their own
# defaults that match this script's, so an unset flag just falls through.
# -----------------------------------------------------------------------------
COMMON_ARGS=( --db "$DB" )
[[ -n "$WORKING_DIR"  ]] && COMMON_ARGS+=( --report     "$WORKING_DIR" )
[[ -n "$DATABASE_DIR" ]] && COMMON_ARGS+=( --db-dir  "$DATABASE_DIR" )

KREPORT_ARGS=( "${COMMON_ARGS[@]}" )
GANON_ARGS=(   "${COMMON_ARGS[@]}" )
SYLPH_ARGS=(   "${COMMON_ARGS[@]}" )

# --techs is forwarded to generate_kreport.sh and ganon_report.sh
# (sylph-tax.sh walks recursively and does not use --techs).
if [[ -n "$TECHS_CSV" ]]; then
    KREPORT_ARGS+=( --techs "$TECHS_CSV" )
    GANON_ARGS+=(   --techs "$TECHS_CSV" )
fi

# --data-groups forwarded to all three sub-scripts.
if [[ -n "$DATA_GROUPS_CSV" ]]; then
    KREPORT_ARGS+=( --data-groups "$DATA_GROUPS_CSV" )
    GANON_ARGS+=(   --data-groups "$DATA_GROUPS_CSV" )
    SYLPH_ARGS+=(   --data-groups "$DATA_GROUPS_CSV" )
fi

# --jobs forwarded to all three sub-scripts.
if [[ -n "$JOBS" ]]; then
    KREPORT_ARGS+=( --jobs "$JOBS" )
    GANON_ARGS+=(   --jobs "$JOBS" )
    SYLPH_ARGS+=(   --jobs "$JOBS" )
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# -----------------------------------------------------------------------------
# Run each sub-script in sequence. Capture exit status but continue so a
# failure in one stage doesn't block the others.
# -----------------------------------------------------------------------------
overall_rc=0
run_step () {
    local label="$1"; shift
    local script="$1"; shift
    echo
    echo "════════════════════════════════════════════════════════════════════"
    echo "  postprocess: ${label}"
    echo "  ${script} $*"
    echo "════════════════════════════════════════════════════════════════════"
    if bash "$script" "$@"; then
        echo "[OK] ${label}"
    else
        local rc=$?
        echo "[FAIL rc=${rc}] ${label} — continuing with remaining steps." >&2
        overall_rc=$rc
    fi
}

run_step "Centrifuge / Centrifuger / Sourmash -> kreport" \
         "${SCRIPT_DIR}/generate_kreport.sh" "${KREPORT_ARGS[@]}"

run_step "Ganon2 -> _reads.tre" \
         "${SCRIPT_DIR}/ganon_report.sh" "${GANON_ARGS[@]}"

run_step "Sylph -> .sylphmpa (sylph-tax taxprof)" \
         "${SCRIPT_DIR}/sylph-tax.sh" "${SYLPH_ARGS[@]}"

echo
echo "════════════════════════════════════════════════════════════════════"
if [[ $overall_rc -eq 0 ]]; then
    echo "  postprocess: all 3 stages completed successfully."
else
    echo "  postprocess: completed with failures (exit $overall_rc). See output above." >&2
fi
echo "════════════════════════════════════════════════════════════════════"
exit $overall_rc
