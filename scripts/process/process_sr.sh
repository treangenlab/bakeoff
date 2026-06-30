#!/bin/bash
# process_sr.sh  (paired-end short reads)
#
# Paired-end counterpart of process.sh. Same model: you pass a list of fastq
# paths (here the R1 of each pair) via --dataset / --dataset-list; the script
# derives the matching R2 (_R1_001 -> _R2_001), runs the selected tools, and
# writes outputs under <out>/<data_group>/<technology>/<Tool>-{results,logs}/.
#
# Loop order is TOOL-outer / SAMPLE-inner so each tool's DB index sits warm
# in the OS page cache across all samples in its block. Don't reorder.
#
# Output basename per sample = the DYN_{patient}_{day} key parsed from the R1
# filename (so results match prior runs), falling back to the R1 stem.
#
# Tool input forms (paired):
#   - kraken2:     --paired R1 R2
#   - centrifuge:  -1 R1 -2 R2
#   - centrifuger: -1 R1 -2 R2
#   - sylph:       profile <dbs...> -1 R1 -2 R2
#   - ganon2:      -p R1 R2
#   - sourmash:    sketch dna takes BOTH reads; gather/tax unchanged
#
# Usage:
#   ./process_sr.sh --db {unified|default} \
#                   [--out <output-root>] [--db-dir <db-root>] \
#                   (--dataset <R1-fastq> | --dataset-list <file>) ... \
#                   --tools <tool_list> [--threads N] [--dry-run]
#
#   --db:           Which database to profile against. One of: unified | default.
#   --out:          Output root. Default: reports/<db>-reports.
#                   Per-sample outputs land under
#                   <out>/<data_group>/<technology>/<Tool>-{results,logs}/.
#   --db-dir:       Parent of the per-DB index trees.
#                   Default: /home/Users/pacbio_bakeoff/data/ref_db
#   --dataset:      An R1 fastq (..._R1_001.fastq), a glob of R1 fastqs, or a
#                   directory (expanded to its *_R1_001.fastq). Repeatable.
#                   The matching R2 is derived as ..._R2_001.fastq.
#   --dataset-list: File with one such path (file/dir/glob) per line.
#                   --dataset and --dataset-list combine; at least one required.
#   --tools:        Comma-separated list of tools, or "all".
#                   "all" -> ganon2,kraken2,centrifuge,centrifuger,sylph,sourmash
#   --threads:      Threads passed to each tool. Default: 10.
#   --dry-run:      Print the planned (tool, sample) runs and exit.
#
# A per-run state log (START/DONE per tool+sample) is always written to a
# uniquely-named file so concurrent runs on shared storage never clobber it:
#   <out>/run_state_<host>_<YYYYMMDD_HHMMSS>.log
#
# NUMA note: sourmash gather is single-threaded; on a multi-socket host with
# kernel.numa_balancing=1 it thrashes (~2x slower). If so, run the sourmash
# step pinned to one node (taskset/numactl) — see the Reproduction wiki page,
# "Execution notes". Not handled in-script (it's a host setting, not pipeline).
#
# Examples (run from inside bakeoff/):
#   ./scripts/process/process_sr.sh --db unified --tools all --threads 50 \
#       --dataset-list data/datasets_dyn_illumina.txt
#   ./scripts/process/process_sr.sh --db default --tools ganon2 \
#       --dataset-list data/datasets_dyn_illumina.txt --dry-run

set -u   # not -e: record per-run failures and continue to the next combo

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
        -h|--help)      sed -n '2,58p' "$0"; exit 0 ;;
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

# Append paths from --dataset-list (if any) to DATASETS.
if [[ -n "$DATASET_LIST" ]]; then
    if [[ ! -f "$DATASET_LIST" ]]; then
        echo "Error: --dataset-list file does not exist: $DATASET_LIST" >&2
        exit 1
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [[ -z "$line" ]] && continue
        DATASETS+=( "$line" )
    done < "$DATASET_LIST"
fi

if [[ "${#DATASETS[@]}" -eq 0 ]]; then
    echo "Error: no datasets provided (use --dataset and/or --dataset-list)." >&2
    echo "       Each entry is an R1 fastq path (..._R1_001.fastq) or glob." >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Per-database index layout (EDITABLE subpaths under --db-dir).
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

# Same tool order as process.sh (sourmash last, easy to split off / pin).
if [[ "$TOOLS" == "all" ]]; then
    selected_tools=(ganon2 kraken2 centrifuge centrifuger sylph sourmash)
else
    IFS=',' read -r -a selected_tools <<< "$TOOLS"
fi

# Output basename: DYN_{patient}_{day} key from the R1 filename, else R1 stem.
sample_base() {
    local bn="$1"
    if [[ "$bn" =~ ^(DYN_[^_]+_[^_]+)_ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "${bn%_R1_001.*}"
    fi
}

# -----------------------------------------------------------------------------
# Pre-expand R1 paths -> parallel arrays (R1, R2, OUT_BASE, OUTPUT_DIR). Globs
# and missing files are resolved once here, before any work, so we can count
# total runs for the [N/total] progress + dry-run plan.
# Each R1 path must look like .../<data_group>/<technology>/<sample>_R1_001.fastq
# -----------------------------------------------------------------------------
SMP_R1=(); SMP_R2=(); SMP_BASE=(); SMP_OUT=()
for entry in "${DATASETS[@]}"; do
    # Each entry may be a single R1 file, a glob, or a directory.
    if [[ -d "$entry" ]]; then
        shopt -s nullglob; r1s=( "$entry"/*_R1_001.fastq "$entry"/*_R1_001.fastq.gz ); shopt -u nullglob
    else
        shopt -s nullglob; r1s=( $entry ); shopt -u nullglob
    fi
    if [[ "${#r1s[@]}" -eq 0 ]]; then
        echo "Warning: no R1 files match: $entry" >&2
        continue
    fi
    for R1 in "${r1s[@]}"; do
        if [[ ! -f "$R1" ]]; then
            echo "Warning: not a file: $R1 (from: $entry)" >&2
            continue
        fi
        R2="${R1/_R1_001/_R2_001}"
        if [[ ! -f "$R2" ]]; then
            echo "Warning: no R2 for $R1 (expected $R2) — skipping" >&2
            continue
        fi
        bn="$(basename "$R1")"
        parent="$(dirname "$R1")"
        technology="$(basename "$parent")"
        data_group="$(basename "$(dirname "$parent")")"
        SMP_R1+=( "$R1" )
        SMP_R2+=( "$R2" )
        SMP_BASE+=( "$(sample_base "$bn")" )
        SMP_OUT+=( "${WORKING_DIR}/${data_group}/${technology}" )
    done
done

if [[ "${#SMP_R1[@]}" -eq 0 ]]; then
    echo "Error: no readable R1/R2 pairs matched the inputs." >&2
    exit 1
fi

n_tools=${#selected_tools[@]}
n_samples=${#SMP_R1[@]}
total=$(( n_tools * n_samples ))

echo "DB:        ${DB}"
echo "Tools:     ${selected_tools[*]}"
echo "Samples:   ${n_samples}"
echo "Threads:   ${THREADS}"
echo "Output:    ${WORKING_DIR}"
echo "State log: ${STATE_LOG}"
echo "Total:     ${total} runs (${n_tools} tools × ${n_samples} samples)"

# -----------------------------------------------------------------------------
# Dry-run: print every planned run with R1/R2, DB, output, then exit.
# -----------------------------------------------------------------------------
if [[ $DRY_RUN -eq 1 ]]; then
    _mark() { if [[ -e "$1" ]]; then echo "[ok]"; else echo "[MISSING]"; fi; }
    echo ""
    echo "=== dry-run plan ($total runs) ==="
    i=0
    for tool in "${selected_tools[@]}"; do
        for s in "${!SMP_R1[@]}"; do
            i=$((i+1))
            printf "[%3d/%d] %-12s %s\n" "$i" "$total" "$tool" "${SMP_BASE[$s]}"
            printf "          R1   : %s %s\n" "${SMP_R1[$s]}" "$(_mark "${SMP_R1[$s]}")"
            printf "          R2   : %s %s\n" "${SMP_R2[$s]}" "$(_mark "${SMP_R2[$s]}")"
            case $tool in
                kraken2)     printf "          db   : %s %s\n" "$K2_DB" "$(_mark "$K2_DB")" ;;
                centrifuge)  printf "          db   : %s.* %s\n" "$CF_IDX" "$(_mark "${CF_IDX}.1.cf")" ;;
                centrifuger) printf "          db   : %s.* %s\n" "$CFER_IDX" "$(_mark "${CFER_IDX}.1.cfr")" ;;
                ganon2)      printf "          db   : %s.hibf %s\n" "$GANON_DB" "$(_mark "${GANON_DB}.hibf")" ;;
                sourmash)    for d in "${SOURMASH_GATHER_DBS[@]}"; do printf "          db   : %s %s\n" "$d" "$(_mark "$d")"; done ;;
                sylph)       for d in "${SYLPH_DBS[@]}"; do printf "          db   : %s %s\n" "$d" "$(_mark "$d")"; done ;;
            esac
            printf "          out  : %s/%s-results/%s*\n" "${SMP_OUT[$s]}" "$tool" "${SMP_BASE[$s]}"
        done
    done
    exit 0
fi

mkdir -p "$(dirname "$STATE_LOG")"
echo "[$(date +'%F %T')] process_sr.sh start: ${total} runs, host=$(hostname), threads=${THREADS}, db=${DB}" \
    | tee -a "$STATE_LOG"

# -----------------------------------------------------------------------------
# Main loop. Outer = TOOL (keeps each tool's DB index warm); inner = SAMPLE.
# -----------------------------------------------------------------------------
i=0
n_pass=0
n_fail=0
for tool in "${selected_tools[@]}"; do
    for s in "${!SMP_R1[@]}"; do
        i=$((i+1))
        R1="${SMP_R1[$s]}"
        R2="${SMP_R2[$s]}"
        OUT_BASE="${SMP_BASE[$s]}"
        OUTPUT_DIR="${SMP_OUT[$s]}"

        label="[${i}/${total}] ${tool} $(basename "$(dirname "$OUTPUT_DIR")")/$(basename "$OUTPUT_DIR")/${OUT_BASE}"
        echo "[$(date +'%F %T')] START ${label}" | tee -a "$STATE_LOG"

        case $tool in
            kraken2)
                mkdir -p "$OUTPUT_DIR/Kraken2-results" "$OUTPUT_DIR/Kraken2-logs"
                /usr/bin/time -v kraken2 \
                    --db "${K2_DB}" \
                    --threads "${THREADS}" \
                    --paired "${R1}" "${R2}" \
                    --output "${OUTPUT_DIR}/Kraken2-results/${OUT_BASE}_classification.tsv" \
                    --report "${OUTPUT_DIR}/Kraken2-results/${OUT_BASE}_report.tsv" \
                    > "${OUTPUT_DIR}/Kraken2-logs/${OUT_BASE}.log" \
                    2> "${OUTPUT_DIR}/Kraken2-logs/${OUT_BASE}.err"
                ;;
            centrifuge)
                mkdir -p "$OUTPUT_DIR/Centrifuge-results" "$OUTPUT_DIR/Centrifuge-logs"
                /usr/bin/time -v centrifuge \
                    -1 "${R1}" -2 "${R2}" \
                    -p "${THREADS}" \
                    -x "${CF_IDX}" \
                    -S "${OUTPUT_DIR}/Centrifuge-results/${OUT_BASE}.tsv" \
                    --report-file "${OUTPUT_DIR}/Centrifuge-results/${OUT_BASE}_report.tsv" \
                    > "${OUTPUT_DIR}/Centrifuge-logs/${OUT_BASE}.log" \
                    2> "${OUTPUT_DIR}/Centrifuge-logs/${OUT_BASE}.err"
                ;;
            centrifuger)
                mkdir -p "$OUTPUT_DIR/Centrifuger-results" "$OUTPUT_DIR/Centrifuger-logs"
                RESULT_FILE="${OUTPUT_DIR}/Centrifuger-results/${OUT_BASE}_results.tsv"
                /usr/bin/time -v centrifuger \
                    -1 "${R1}" -2 "${R2}" \
                    -t "${THREADS}" \
                    -x "${CFER_IDX}" \
                    > "${RESULT_FILE}" \
                    2> "${OUTPUT_DIR}/Centrifuger-logs/${OUT_BASE}.err"

                /usr/bin/time -v centrifuger-quant \
                    -c "${RESULT_FILE}" \
                    -x "${CFER_IDX}" \
                    > "${OUTPUT_DIR}/Centrifuger-results/${OUT_BASE}_report.tsv" \
                    2> "${OUTPUT_DIR}/Centrifuger-logs/${OUT_BASE}_quant.err"
                ;;
            ganon2)
                # Illumina headers have no tabs, so ganon's reassign step is
                # unaffected (unlike dorado ONT in process.sh). No header strip.
                mkdir -p "$OUTPUT_DIR/ganon2-results" "$OUTPUT_DIR/ganon2-logs"
                /usr/bin/time -v ganon classify \
                    --db-prefix "${GANON_DB}" \
                    -p "${R1}" "${R2}" \
                    -o "${OUTPUT_DIR}/ganon2-results/${OUT_BASE}" \
                    --min-count 0 \
                    -t "${THREADS}" \
                    --verbose \
                    > "${OUTPUT_DIR}/ganon2-logs/${OUT_BASE}.log" \
                    2> "${OUTPUT_DIR}/ganon2-logs/${OUT_BASE}.err"
                ;;
            sylph)
                mkdir -p "$OUTPUT_DIR/sylph-results" "$OUTPUT_DIR/sylph-logs"
                /usr/bin/time -v sylph profile \
                    "${SYLPH_DBS[@]}" \
                    -1 "${R1}" -2 "${R2}" \
                    -t "${THREADS}" \
                    -o "${OUTPUT_DIR}/sylph-results/${OUT_BASE}_profiling.tsv" \
                    > "${OUTPUT_DIR}/sylph-logs/${OUT_BASE}.log" \
                    2> "${OUTPUT_DIR}/sylph-logs/${OUT_BASE}.err"
                ;;
            sourmash)
                mkdir -p "$OUTPUT_DIR/sourmash-sketches" "$OUTPUT_DIR/sourmash-results" "$OUTPUT_DIR/sourmash-logs"
                /usr/bin/time -v sourmash sketch dna \
                    -p k=31,abund \
                    "${R1}" "${R2}" \
                    --name "${OUT_BASE}" \
                    -o "${OUTPUT_DIR}/sourmash-sketches/${OUT_BASE}.sig" \
                    > "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_sketch.log" \
                    2> "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_sketch.err"

                ( /usr/bin/time -v sourmash gather \
                    --dna --ksize 31 \
                    "${OUTPUT_DIR}/sourmash-sketches/${OUT_BASE}.sig" \
                    "${SOURMASH_GATHER_DBS[@]}" \
                    -o "${OUTPUT_DIR}/sourmash-results/${OUT_BASE}-gather.csv" \
                    > "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_gather.log" \
                ) 2> "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_gather.err"

                TAX_ARGS=()
                for lin in "${SOURMASH_LINEAGES[@]}"; do
                    TAX_ARGS+=( --taxonomy-csv "${lin}" )
                done
                /usr/bin/time -v sourmash tax metagenome \
                    -o "${OUTPUT_DIR}/sourmash-results/${OUT_BASE}-report.csv" \
                    --gather-csv "${OUTPUT_DIR}/sourmash-results/${OUT_BASE}-gather.csv" \
                    "${TAX_ARGS[@]}" \
                    -r species \
                    -F csv_summary \
                    > "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_tax.log" \
                    2> "${OUTPUT_DIR}/sourmash-logs/${OUT_BASE}_tax.err"
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

echo "[$(date +'%F %T')] process_sr.sh complete: ${n_pass}/${total} succeeded, ${n_fail} failed" \
    | tee -a "$STATE_LOG"

[[ $n_fail -eq 0 ]]
