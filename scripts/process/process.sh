#!/bin/bash
# process.sh
# Run any of the 6 paper profilers on a list of fastq files against either the
# unified RefSeq v228 database or each tool's native default database.
#
# Loop order is TOOL-outer / DATASET-inner so each tool's DB index sits warm
# in the OS page cache across all datasets in its block. Don't reorder.
#
# Usage:
#   ./process.sh --db {unified|default} \
#                [--out <output-root>] [--db-dir <db-root>] \
#                (--dataset <fastq> | --dataset-list <file>) ... \
#                --tools <tool_list> [--threads N] [--dry-run]
#
#   --db:           Which database to profile against. One of: unified | default.
#   --out:          Output root. Default: reports/<db>-reports (so unified
#                   lands in reports/unified-reports/ and default in
#                   reports/default-reports/).
#                   Per-sample outputs land under
#                   <out>/<data_group>/<technology>/<Tool>-{results,logs}/.
#   --db-dir:       Parent of the per-DB index trees.
#                   Default: /home/Users/pacbio_bakeoff/data/ref_db
#                   The script appends `refseq03032025/` for --db unified and
#                   `default_db/` for --db default to reach the per-tool indexes.
#   --dataset:      A single fastq, a glob of fastqs, or a directory (expanded
#                   to its *.fastq). Repeatable.
#   --dataset-list: File with one such path (file/dir/glob) per line
#                   (#-comments + blanks ok). --dataset and --dataset-list
#                   combine; at least one entry from either is required.
#   --tools:        Comma-separated list of tools, or "all".
#                   Options: ganon2, kraken2, centrifuge, centrifuger, sylph, sourmash
#                   "all" expands to: ganon2,kraken2,centrifuge,centrifuger,sylph,sourmash
#                   (sourmash last so it's easy to exclude / run pinned separately)
#   --threads:      Threads passed to each tool. Default: 10.
#   --dry-run:      Print the planned (tool, dataset) runs and exit.
#
# A per-run state log (START/DONE per tool+dataset) is always written to a
# uniquely-named file so concurrent runs on shared storage never clobber it:
#   <out>/run_state_<host>_<YYYYMMDD_HHMMSS>.log
#
# NUMA note: sourmash gather is single-threaded; on a multi-socket host with
# kernel.numa_balancing=1 it thrashes (~2x slower). If so, run the sourmash
# step pinned to one node (taskset/numactl) — see the Reproduction wiki page,
# "Execution notes". Not handled in-script (it's a host setting, not pipeline).
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/process.sh --db unified --tools all --threads 50 \
#       --dataset-list data/datasets_mock_sim.txt
#
#   ./scripts/process/process.sh --db default --tools sylph \
#       --dataset-list data/datasets_dyn_pacbio.txt
#
#   ./scripts/process/process.sh --db unified --tools all --threads 50 \
#       --dataset-list data/datasets_mock_sim.txt --dry-run

set -u

THREADS=10
TOOLS=""
DB=""
WORKING_DIR=""
DATABASE_DIR="/home/Users/pacbio_bakeoff/data/ref_db"
DATASETS=()
DATASET_LIST=""
DRY_RUN=0

# -----------------------------------------------------------------------------
# Parse args
# -----------------------------------------------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --db)           DB="$2";              shift 2 ;;
        --out)          WORKING_DIR="$2";     shift 2 ;;
        --db-dir)       DATABASE_DIR="$2";    shift 2 ;;
        --dataset)      DATASETS+=( "$2" );   shift 2 ;;
        --dataset-list) DATASET_LIST="$2";    shift 2 ;;
        --tools)        TOOLS="$2";           shift 2 ;;
        --threads)      THREADS="$2";         shift 2 ;;
        --dry-run)      DRY_RUN=1;            shift   ;;
        -h|--help)
            sed -n '2,54p' "$0"; exit 0 ;;
        *)
            echo "Unknown parameter: $1" >&2
            echo "Run with --help for usage." >&2
            exit 1 ;;
    esac
done

if [[ -z "$DB" ]]; then
    echo "Error: --db is required (unified|default)." >&2; exit 1
fi
if [[ -z "$TOOLS" ]]; then
    echo "Error: --tools is required (comma list or 'all')." >&2
    echo "Valid: kraken2, centrifuge, centrifuger, sylph, sourmash, ganon2" >&2
    exit 1
fi

# Default --out is per-DB so unified and default outputs never collide.
if [[ -z "$WORKING_DIR" ]]; then
    WORKING_DIR="reports/${DB}-reports"
fi

# State log: always written, uniquely named (host + start timestamp) so
# concurrent runs sharing the same OUT tree don't clobber each other's log.
STATE_LOG="${WORKING_DIR}/run_state_$(hostname -s)_$(date +%Y%m%d_%H%M%S).log"

# Append paths from --dataset-list (if any) to DATASETS. Empty lines and
# lines starting with # are skipped.
if [[ -n "$DATASET_LIST" ]]; then
    if [[ ! -f "$DATASET_LIST" ]]; then
        echo "Error: --dataset-list file does not exist: $DATASET_LIST" >&2
        exit 1
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"   # strip trailing #-comment
        # trim leading + trailing whitespace
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [[ -z "$line" ]] && continue
        DATASETS+=( "$line" )
    done < "$DATASET_LIST"
fi

if [[ "${#DATASETS[@]}" -eq 0 ]]; then
    echo "Error: no datasets provided (use --dataset and/or --dataset-list)." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Per-database index layout. Subpaths under --db-dir are EDITABLE here to
# match how you built the indexes (subdirectory names, file basenames). The
# top-level --db-dir comes from the CLI.
# Multi-DB tools (sylph profile, sourmash gather + tax) use arrays so the same
# loop body works for both unified (1 DB) and default (3-4 DBs).
# -----------------------------------------------------------------------------
case "$DB" in
    unified)
        DB_ROOT="${DATABASE_DIR}/refseq03032025"

        K2_DB="${DB_ROOT}/k2_abfv_030325"
        CF_IDX="${DB_ROOT}/centrifuge_abv_030325/refseq_abv"
        CFER_IDX="${DB_ROOT}/centrifuger_abv_030325/refseq_abv"

        SYLPH_DBS=( "${DB_ROOT}/sylph_abf_030325/database.syldb" )

        SOURMASH_GATHER_DBS=( "${DB_ROOT}/sourmash_abvf_030325/sm_abvf_030325.zip" )
        SOURMASH_LINEAGES=(  "${DB_ROOT}/sourmash_abvf_030325/lineage_030325.csv" )

        GANON_DB="${DB_ROOT}/ganon2_abvf_030325/ganon2_abvf_030325"
        ;;
    default)
        DB_ROOT="${DATABASE_DIR}/default_db"

        K2_DB="${DB_ROOT}/k2_default"
        CF_IDX="${DB_ROOT}/cf_default/hpvc"
        CFER_IDX="${DB_ROOT}/cfer_default/hpv_gbsarscov2"

        SYLPH_DBS=(
            "${DB_ROOT}/sylph_default/gtdb-r220-c200-dbv1.syldb"
            "${DB_ROOT}/sylph_default/imgvr_c200_v0.3.0.syldb"
            "${DB_ROOT}/sylph_default/fungi-refseq-2024-07-25-c200-v0.3.syldb"
        )

        # Sourmash default ships per-domain zips/lineages; list each one
        # explicitly rather than relying on brace expansion.
        SM_DEFAULT="${DB_ROOT}/sm_default"
        SOURMASH_GATHER_DBS=(
            "${SM_DEFAULT}/genbank-2022.03-archaea-k31.zip"
            "${SM_DEFAULT}/genbank-2022.03-fungi-k31.zip"
            "${SM_DEFAULT}/genbank-2022.03-bacteria-k31.zip"
            "${SM_DEFAULT}/genbank-2022.03-viral-k31.zip"
        )
        SOURMASH_LINEAGES=(
            "${SM_DEFAULT}/genbank-2022.03-archaea.lineages.csv"
            "${SM_DEFAULT}/genbank-2022.03-fungi.lineages.csv"
            "${SM_DEFAULT}/genbank-2022.03-bacteria.lineages.csv"
            "${SM_DEFAULT}/genbank-2022.03-viral.lineages.csv"
        )

        GANON_DB="${DB_ROOT}/ganon2_default/ganon2_default_abfv_rs_cg"
        ;;
    *)
        echo "Error: --db must be 'unified' or 'default', got '${DB}'." >&2
        exit 1 ;;
esac

# Tool-outer loop keeps each tool's DB index warm across its datasets. sourmash
# runs last so it's easy to drop and run separately (pinned) on NUMA hosts.
if [[ "$TOOLS" == "all" ]]; then
    selected_tools=(ganon2 kraken2 centrifuge centrifuger sylph sourmash)
else
    IFS=',' read -r -a selected_tools <<< "$TOOLS"
fi

# Pre-expand DATASETS (also gives the run count). Each --dataset / --dataset-list
# entry may be a single fastq, a glob, or a directory. Every resolved path must
# look like .../<data_group>/<technology>/<sample>.fastq so the script can lift
# <data_group>/<technology> into the output layout.
expanded_datasets=()
for entry in "${DATASETS[@]}"; do
    if [[ -d "$entry" ]]; then
        shopt -s nullglob
        matches=( "$entry"/*.fastq "$entry"/*.fastq.gz "$entry"/*.fq "$entry"/*.fq.gz )
        shopt -u nullglob
    else
        shopt -s nullglob; matches=( $entry ); shopt -u nullglob
    fi
    if [[ "${#matches[@]}" -eq 0 ]]; then
        echo "Warning: no files match: $entry" >&2
        continue
    fi
    for QUERY_FQ_PATH in "${matches[@]}"; do
        if [[ ! -f "$QUERY_FQ_PATH" ]]; then
            echo "Warning: not a file: $QUERY_FQ_PATH (from: $entry)" >&2
            continue
        fi
        expanded_datasets+=( "$QUERY_FQ_PATH" )
    done
done

if [[ "${#expanded_datasets[@]}" -eq 0 ]]; then
    echo "Error: no readable fastq files matched the inputs." >&2
    exit 1
fi

n_tools=${#selected_tools[@]}
n_datasets=${#expanded_datasets[@]}
total=$(( n_tools * n_datasets ))

echo "DB:        ${DB}"
echo "Tools:     ${selected_tools[*]}"
echo "Datasets:  ${n_datasets}"
echo "Threads:   ${THREADS}"
echo "Output:    ${WORKING_DIR}"
echo "State log: ${STATE_LOG}"
echo "Total:     ${total} runs (${n_tools} tools × ${n_datasets} datasets)"

# -----------------------------------------------------------------------------
# Dry-run: print every planned run with its input, DB index, and output/log
# paths, then exit before any mkdir / state-log writes. A [missing] marker on
# a DB line flags a wrong --db-dir before you commit to a long run.
# -----------------------------------------------------------------------------
if [[ $DRY_RUN -eq 1 ]]; then
    # Mark whether a path (file, dir, or index prefix) resolves on disk.
    _mark() { if [[ -e "$1" ]]; then echo "[ok]"; else echo "[MISSING]"; fi; }

    echo ""
    echo "=== dry-run plan ($total runs) ==="
    i=0
    for tool in "${selected_tools[@]}"; do
        for QUERY_FQ_PATH in "${expanded_datasets[@]}"; do
            i=$((i+1))
            QUERY_FQ_FILE=$(basename "$QUERY_FQ_PATH")
            QUERY_FQ_NAME="${QUERY_FQ_FILE%.*}"
            PARENT_DIR=$(dirname "$QUERY_FQ_PATH")
            TECHNOLOGY=$(basename "$PARENT_DIR")
            DATA_GROUP=$(basename "$(dirname "$PARENT_DIR")")
            OUTPUT_DIR="${WORKING_DIR}/${DATA_GROUP}/${TECHNOLOGY}"

            printf "[%3d/%d] %-12s %s/%s/%s\n" "$i" "$total" "$tool" "$DATA_GROUP" "$TECHNOLOGY" "$QUERY_FQ_NAME"
            printf "          input: %s %s\n" "$QUERY_FQ_PATH" "$(_mark "$QUERY_FQ_PATH")"
            case $tool in
                kraken2)
                    printf "          db   : %s %s\n" "$K2_DB" "$(_mark "$K2_DB")"
                    printf "          out  : %s/Kraken2-results/%s_{classification,report}.tsv\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/Kraken2-logs/%s.{log,err}\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
                centrifuge)
                    printf "          db   : %s.* %s\n" "$CF_IDX" "$(_mark "${CF_IDX}.1.cf")"
                    printf "          out  : %s/Centrifuge-results/%s{,_report}.tsv\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/Centrifuge-logs/%s.{log,err}\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
                centrifuger)
                    printf "          db   : %s.* %s\n" "$CFER_IDX" "$(_mark "${CFER_IDX}.1.cfr")"
                    printf "          out  : %s/Centrifuger-results/%s_{results,report}.tsv\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/Centrifuger-logs/%s{,_quant}.err\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
                sylph)
                    for d in "${SYLPH_DBS[@]}"; do
                        printf "          db   : %s %s\n" "$d" "$(_mark "$d")"
                    done
                    printf "          out  : %s/sylph-results/%s_profiling.tsv\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/sylph-logs/%s.{log,err}\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
                sourmash)
                    for d in "${SOURMASH_GATHER_DBS[@]}"; do
                        printf "          db   : %s %s\n" "$d" "$(_mark "$d")"
                    done
                    for l in "${SOURMASH_LINEAGES[@]}"; do
                        printf "          tax  : %s %s\n" "$l" "$(_mark "$l")"
                    done
                    printf "          out  : %s/sourmash-results/%s-{gather,report}.csv\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/sourmash-logs/%s_{sketch,gather,tax}.{log,err}\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
                ganon2)
                    printf "          db   : %s.hibf %s\n" "$GANON_DB" "$(_mark "${GANON_DB}.hibf")"
                    printf "          out  : %s/ganon2-results/%s.*\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    printf "          log  : %s/ganon2-logs/%s.{log,err}\n" "$OUTPUT_DIR" "$QUERY_FQ_NAME"
                    ;;
            esac
        done
    done
    exit 0
fi

mkdir -p "$(dirname "$STATE_LOG")"
echo "[$(date +'%F %T')] process.sh start: ${total} runs, host=$(hostname), threads=${THREADS}, db=${DB}" \
    | tee -a "$STATE_LOG"

# -----------------------------------------------------------------------------
# Main loop. Outer = TOOL (keeps each tool's DB index warm in OS page cache);
# inner = DATASET (6 datasets per tool block, first cold, rest warm).
# -----------------------------------------------------------------------------
i=0
n_pass=0
n_fail=0
for tool in "${selected_tools[@]}"; do
    for QUERY_FQ_PATH in "${expanded_datasets[@]}"; do
        i=$((i+1))

        QUERY_FQ_FILE=$(basename "$QUERY_FQ_PATH")
        QUERY_FQ_NAME="${QUERY_FQ_FILE%.*}"

        # Output layout: <WORKING_DIR>/<data_group>/<technology>/
        PARENT_DIR=$(dirname "$QUERY_FQ_PATH")
        TECHNOLOGY=$(basename "$PARENT_DIR")
        DATA_GROUP=$(basename "$(dirname "$PARENT_DIR")")

        OUTPUT_DIR="${WORKING_DIR}/${DATA_GROUP}/${TECHNOLOGY}"

        label="[${i}/${total}] ${tool} ${DATA_GROUP}/${TECHNOLOGY}/${QUERY_FQ_NAME}"
        echo "[$(date +'%F %T')] START ${label}" | tee -a "$STATE_LOG"

        case $tool in
            kraken2)
                mkdir -p "$OUTPUT_DIR/Kraken2-results" "$OUTPUT_DIR/Kraken2-logs"
                /usr/bin/time -v kraken2 \
                    --db "${K2_DB}" \
                    --threads "${THREADS}" \
                    --output "${OUTPUT_DIR}/Kraken2-results/${QUERY_FQ_NAME}_classification.tsv" \
                    --report "${OUTPUT_DIR}/Kraken2-results/${QUERY_FQ_NAME}_report.tsv" \
                    "${QUERY_FQ_PATH}" \
                    > "${OUTPUT_DIR}/Kraken2-logs/${QUERY_FQ_NAME}.log" \
                    2> "${OUTPUT_DIR}/Kraken2-logs/${QUERY_FQ_NAME}.err"
                ;;
            centrifuge)
                mkdir -p "$OUTPUT_DIR/Centrifuge-results" "$OUTPUT_DIR/Centrifuge-logs"
                /usr/bin/time -v centrifuge \
                    -U "${QUERY_FQ_PATH}" \
                    -p "${THREADS}" \
                    -x "${CF_IDX}" \
                    -S "${OUTPUT_DIR}/Centrifuge-results/${QUERY_FQ_NAME}.tsv" \
                    --report-file "${OUTPUT_DIR}/Centrifuge-results/${QUERY_FQ_NAME}_report.tsv" \
                    > "${OUTPUT_DIR}/Centrifuge-logs/${QUERY_FQ_NAME}.log" \
                    2> "${OUTPUT_DIR}/Centrifuge-logs/${QUERY_FQ_NAME}.err"
                ;;
            centrifuger)
                mkdir -p "$OUTPUT_DIR/Centrifuger-results" "$OUTPUT_DIR/Centrifuger-logs"
                RESULT_FILE="${OUTPUT_DIR}/Centrifuger-results/${QUERY_FQ_NAME}_results.tsv"
                /usr/bin/time -v centrifuger \
                    -u "${QUERY_FQ_PATH}" \
                    -t "${THREADS}" \
                    -x "${CFER_IDX}" \
                    > "${RESULT_FILE}" \
                    2> "${OUTPUT_DIR}/Centrifuger-logs/${QUERY_FQ_NAME}.err"

                REPORT_FILE="${OUTPUT_DIR}/Centrifuger-results/${QUERY_FQ_NAME}_report.tsv"
                /usr/bin/time -v centrifuger-quant \
                    -c "${RESULT_FILE}" \
                    -x "${CFER_IDX}" \
                    > "${REPORT_FILE}" \
                    2> "${OUTPUT_DIR}/Centrifuger-logs/${QUERY_FQ_NAME}_quant.err"
                ;;
            sylph)
                mkdir -p "$OUTPUT_DIR/sylph-results" "$OUTPUT_DIR/sylph-logs"
                /usr/bin/time -v sylph profile \
                    "${SYLPH_DBS[@]}" \
                    "${QUERY_FQ_PATH}" \
                    -t "${THREADS}" \
                    -o "${OUTPUT_DIR}/sylph-results/${QUERY_FQ_NAME}_profiling.tsv" \
                    > "${OUTPUT_DIR}/sylph-logs/${QUERY_FQ_NAME}.log" \
                    2> "${OUTPUT_DIR}/sylph-logs/${QUERY_FQ_NAME}.err"
                ;;
            sourmash)
                mkdir -p "$OUTPUT_DIR/sourmash-sketches" "$OUTPUT_DIR/sourmash-results" "$OUTPUT_DIR/sourmash-logs"
                /usr/bin/time -v sourmash sketch dna \
                    -p k=31,abund \
                    -o "${OUTPUT_DIR}/sourmash-sketches/${QUERY_FQ_NAME}.sketch" \
                    "${QUERY_FQ_PATH}" \
                    > "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_sketch.log" \
                    2> "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_sketch.err"

                # ( ... ) so the redirect captures both `time -v` and sourmash stderr.
                ( /usr/bin/time -v sourmash gather \
                    -k 31 \
                    "${OUTPUT_DIR}/sourmash-sketches/${QUERY_FQ_NAME}.sketch" \
                    "${SOURMASH_GATHER_DBS[@]}" \
                    -o "${OUTPUT_DIR}/sourmash-results/${QUERY_FQ_NAME}-gather.csv" \
                    > "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_gather.log" \
                ) 2> "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_gather.err"

                # sourmash tax metagenome takes --taxonomy-csv once per lineage file.
                TAX_ARGS=()
                for lin in "${SOURMASH_LINEAGES[@]}"; do
                    TAX_ARGS+=( --taxonomy-csv "${lin}" )
                done
                /usr/bin/time -v sourmash tax metagenome \
                    -o "${OUTPUT_DIR}/sourmash-results/${QUERY_FQ_NAME}-report.csv" \
                    --gather-csv "${OUTPUT_DIR}/sourmash-results/${QUERY_FQ_NAME}-gather.csv" \
                    "${TAX_ARGS[@]}" \
                    -r species \
                    -F csv_summary \
                    > "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_tax.log" \
                    2> "${OUTPUT_DIR}/sourmash-logs/${QUERY_FQ_NAME}_tax.err"
                ;;
            ganon2)
                mkdir -p "$OUTPUT_DIR/ganon2-results" "$OUTPUT_DIR/ganon2-logs"
                # ganon's EM reassign step splits the .all on tabs (expects 3
                # fields); dorado ONT headers carry tab-delimited SAM tags that
                # abort EM and corrupt per-taxon counts. Strip such headers to a
                # temp file BEFORE the timed command (so the awk pass isn't timed).
                # PacBio/illumina/simulated headers have no tabs -> no temp.
                GANON_INPUT="${QUERY_FQ_PATH}"
                GANON_TMP=""
                if head -1 "${QUERY_FQ_PATH}" | grep -q $'\t'; then
                    GANON_TMP="${OUTPUT_DIR}/ganon2-results/.${QUERY_FQ_NAME}.cleanhdr.fastq"
                    awk 'NR%4==1{sub(/\t.*/,"")} 1' "${QUERY_FQ_PATH}" > "${GANON_TMP}"
                    GANON_INPUT="${GANON_TMP}"
                fi

                /usr/bin/time -v ganon classify \
                    --db-prefix "${GANON_DB}" \
                    -s "${GANON_INPUT}" \
                    -o "${OUTPUT_DIR}/ganon2-results/${QUERY_FQ_NAME}" \
                    --min-count 0 \
                    -t "${THREADS}" \
                    --verbose \
                    > "${OUTPUT_DIR}/ganon2-logs/${QUERY_FQ_NAME}.log" \
                    2> "${OUTPUT_DIR}/ganon2-logs/${QUERY_FQ_NAME}.err"
                ganon_rc=$?
                [[ -n "${GANON_TMP}" ]] && rm -f "${GANON_TMP}"
                ( exit "$ganon_rc" )   # rc reflects ganon, not the temp cleanup
                ;;
            *)
                echo "Unknown tool specified: $tool" >&2
                echo "Valid: kraken2, centrifuge, centrifuger, sylph, sourmash, ganon2" >&2
                ;;
        esac
        rc=$?

        if [[ $rc -eq 0 ]]; then
            n_pass=$((n_pass+1))
        else
            n_fail=$((n_fail+1))
        fi
        echo "[$(date +'%F %T')] DONE  ${label} rc=${rc}" | tee -a "$STATE_LOG"
    done
done

echo "[$(date +'%F %T')] process.sh complete: ${n_pass}/${total} succeeded, ${n_fail} failed" \
    | tee -a "$STATE_LOG"

[[ $n_fail -eq 0 ]]
