#!/usr/bin/env bash
# Sweep a figure notebook across a list of values for a single config variable,
# optionally repeating the sweep across multiple datasets and/or with fixed
# overrides for other variables.
#
# Each iteration:
#   1. patches the notebook's config cells to set the fixed overrides + the swept variable
#   2. executes the patched copy with jupyter nbconvert
#   3. the notebook's own output paths are namespaced by the swept variable,
#      so per-iteration results land in distinct subdirs and never collide.
#
# When --datasets is given, the whole sweep is repeated for each dataset key
# (i.e. a cross-product of datasets × the swept variable). For each dataset,
# DATASET is pinned in the config cell before the inner sweep runs.
#
# Examples:
#   # Default: detection threshold sweep on the notebook's current DATASET
#   scripts/analysis/sweep_notebook.sh
#
#   # Threshold sweep × every detection dataset
#   scripts/analysis/sweep_notebook.sh --datasets \
#       "pacbio_low_input pacbio_standard_input ont_zymo_kit ont_qiagen_kit sim_pacbio sim_ont"
#
#   # Threshold sweep against a pinned analysis-prep snapshot
#   scripts/analysis/sweep_notebook.sh --set 'TS="20260530_225610"'
#
#   # Abundance sweep × every mock dataset
#   scripts/analysis/sweep_notebook.sh \
#       --notebook scripts/analysis/analysis_abundance.ipynb \
#       --var MIN_ABUNDANCE \
#       --values "0 0.001 0.01 0.1" \
#       --datasets "pacbio_low_input pacbio_standard_input ont_zymo_kit ont_qiagen_kit"
#
# Options:
#   --notebook PATH    Notebook to sweep (default: scripts/analysis/analysis_detection.ipynb)
#   --var NAME         Swept config variable (default: THRESHOLD_PERCENT)
#   --values "..."     Space-separated value list (default: "0.0 0.0001 0.001 0.01 0.1")
#   --datasets "..."   Space-separated DATASET keys to iterate over (outer loop).
#                      Omit to use the DATASET already set in the notebook.
#   --set VAR=VALUE    Fixed override applied every iteration. Repeatable.
#                      VALUE is inserted verbatim into Python — quote strings:
#                          --set 'TS="20260530_225610"'
#                          --set "DATASET='sim_pacbio'"
#   --keep             Keep the executed notebooks under /tmp (otherwise removed)
#   -h, --help         Show this help
#
# Notes:
#   * Run from the bakeoff/ repo root — the notebook's import path walker
#     needs `config/bakeoff_env.yaml` reachable from CWD.
#   * The notebook's existing output paths must already namespace by the swept
#     variable (e.g. results/<TS>/detection/threshold_<x>/<dataset>/), otherwise
#     iterations clobber each other.
#   * SAVE=True is forced on every iteration (it would defeat the purpose of a
#     sweep to leave it off); --set 'SAVE=False' is ignored.

set -euo pipefail

NOTEBOOK="scripts/analysis/analysis_detection.ipynb"
VAR="THRESHOLD_PERCENT"
VALUES="0.0 0.0001 0.001 0.01 0.1"
DATASETS=""             # outer-loop DATASET keys (optional)
KEEP=0
SET_OVERRIDES=()        # array of "VAR=VALUE" strings, applied every iteration

while [[ $# -gt 0 ]]; do
  case "$1" in
    --notebook) NOTEBOOK="$2"; shift 2 ;;
    --var)      VAR="$2";      shift 2 ;;
    --values)   VALUES="$2";   shift 2 ;;
    --datasets) DATASETS="$2"; shift 2 ;;
    --set)      SET_OVERRIDES+=("$2"); shift 2 ;;
    --keep)     KEEP=1;        shift ;;
    -h|--help)  sed -n '2,52p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Sanity checks
[[ -f config/bakeoff_env.yaml ]] || {
  echo "ERROR: run from the bakeoff/ repo root (config/bakeoff_env.yaml not found in $PWD)" >&2
  exit 1
}
[[ -f "$NOTEBOOK" ]] || { echo "ERROR: notebook not found: $NOTEBOOK" >&2; exit 1; }
jupyter nbconvert --version >/dev/null 2>&1 || {
  echo "ERROR: 'jupyter nbconvert' not available in this shell." >&2
  echo "       Activate the bakeoff env (e.g. 'conda activate bakeoff') or install nbconvert." >&2
  exit 1
}

# Temp scratch dir (cleaned up on exit unless --keep)
SCRATCH=$(mktemp -d -t bakeoff_sweep.XXXXXX)
[[ $KEEP -eq 1 ]] || trap 'rm -rf "$SCRATCH"' EXIT

NB_NAME=$(basename "$NOTEBOOK" .ipynb)
TMP_NB="_sweep_${NB_NAME}.ipynb"   # lives in repo root so kernel CWD picks up utils

# Echo the fixed overrides once (helps the user see what's pinned).
if [[ ${#SET_OVERRIDES[@]} -gt 0 ]]; then
  echo "Fixed overrides applied every iteration:"
  for kv in "${SET_OVERRIDES[@]}"; do echo "    $kv"; done
  echo
fi

# Outer dataset loop: empty entry = "use the notebook's existing DATASET".
if [[ -n "$DATASETS" ]]; then
  read -r -a DATASET_ARR <<< "$DATASETS"
else
  DATASET_ARR=("")
fi

for DS in "${DATASET_ARR[@]}"; do
  if [[ -n "$DS" ]]; then
    echo "============================================================"
    echo "DATASET = ${DS}"
    echo "============================================================"
    DS_OVERRIDE=("DATASET=\"${DS}\"")
  else
    DS_OVERRIDE=()
  fi

  for VALUE in $VALUES; do
    echo "==> ${VAR} = ${VALUE}"

    # Per-iteration assignment list = fixed overrides + DATASET (if any) + swept var
    # + SAVE=True (forced last so it always wins; sweeping without writing outputs
    # to disk is never what the caller wants).
    ASSIGNMENTS=("${SET_OVERRIDES[@]}" "${DS_OVERRIDE[@]}" "${VAR}=${VALUE}" "SAVE=True")

    python3 - "$NOTEBOOK" "$TMP_NB" "${ASSIGNMENTS[@]}" <<'PY'
import json, sys, re

src_path, dst_path, *assignments = sys.argv[1:]
overrides = {}
for kv in assignments:
    var, sep, val = kv.partition('=')
    if not sep or not var:
        sys.exit(f"ERROR: --set value must be VAR=VALUE, got {kv!r}")
    overrides[var] = val

nb = json.load(open(src_path))
hits = {v: 0 for v in overrides}

for c in nb['cells']:
    if c.get('cell_type') != 'code':
        continue
    for i, line in enumerate(c.get('source', [])):
        for var, val in overrides.items():
            m = re.match(rf'^(\s*){re.escape(var)}\s*=', line)
            if m:
                indent = m.group(1)
                comment = ('    ' + line[line.index('#'):].rstrip('\n')) if '#' in line else ''
                c['source'][i] = f'{indent}{var} = {val}{comment}\n'
                hits[var] += 1
                break

missing = [v for v, n in hits.items() if n == 0]
if missing:
    sys.exit(f"ERROR: no assignment found for: {', '.join(missing)}")

json.dump(nb, open(dst_path, 'w'), indent=1)
PY

    DS_TAG=${DS:-default}
    EXECUTED="$SCRATCH/${NB_NAME}_${DS_TAG}_${VAR}_${VALUE}.executed.ipynb"
    jupyter nbconvert --to notebook --execute "$TMP_NB" --output "$EXECUTED"
    echo
  done
done

rm -f "$TMP_NB"

if [[ $KEEP -eq 1 ]]; then
  echo "Executed notebooks: $SCRATCH"
else
  echo "Done. Output tables/figures live where the notebook normally writes them"
  echo "(e.g. results/<TS>/detection/threshold_<value>/<dataset>/ or results/<TS>/abundance/...)."
fi
