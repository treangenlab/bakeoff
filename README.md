# bakeoff - Long-read taxonomic profiling of the gut microbiome

Companion repository for Huang et al., **"Reference Databases Matter: Disentangling Database and Algorithm Effects in Long-Read Gut Microbiome Profiling Performance"** (manuscript in preparation).

## Introduction

This study benchmarks six long-read-capable taxonomic profilers - **Kraken2**, **Centrifuge**, **Centrifuger**, **Ganon2**, **Sourmash**, and **Sylph** - across three dataset categories (PBSIM3 + MIMIC simulated reads, ZymoBIOMICS D6331 mock community, and the DYN clinical cohort) on both PacBio HiFi and Oxford Nanopore reads. Each method is evaluated under (a) its native default database and (b) a shared **unified RefSeq-228 build** (~69k genomes spanning archaea, bacteria, fungi, and viruses), holding the reference fixed to disentangle database effects from algorithm effects.

Three findings the paper argues:

1. **Database recency matters across the board.** Database choice materially affects accuracy regardless of the tool, especially for lower-abundance taxa. Harmonizing to a unified RefSeq reference visibly reduced cross-tool disagreement on both mock and clinical samples - keeping the reference current is best practice regardless of the algorithm.
2. **Sketch-based profilers led the controlled benchmarks.** Across simulated and mock data, **Sourmash**, **Sylph**, and **Ganon2** posted the strongest detection performance, generally driven by higher precision at comparably high recall.
3. **Controlled-benchmark rankings transfer only indirectly to real data.** On the clinical cohort, cross-tool agreement among the more accurate methods held for moderate-to-high-abundance taxa but weakened on low-abundance ones - so a multi-method composite approach may still be valuable in practice.

Prior long-read benchmarks emphasized short reads, used method-specific databases (confounding algorithm vs. DB), and relied on narrow synthetic communities. This study fixes the reference and tests across mock, simulated, and clinical data to isolate where each effect lives.

## Install

### 1. Clone the repository

```bash
git clone https://github.com/treangenlab/bakeoff.git
cd bakeoff
```

### 2. Create the conda environment

The full pipeline environment - the six paper profilers (version-pinned to manuscript Methods), profiler-adjacent helpers (`sylph-tax`, `taxonkit`, `multitax`, `mosdepth`, `lemur`), and the analysis stack (Python, Jupyter, scikit-learn, seaborn, ete3, biopython, etc.) - is declared in [`config/bakeoff_env.yaml`](config/bakeoff_env.yaml):

```bash
conda env create -n bakeoff -f config/bakeoff_env.yaml
conda activate bakeoff
```

Mamba works identically (`mamba env create ...`) and is significantly faster on this env.

Pinned profiler versions:

| Tool | Version | Notes |
| --- | ---: | --- |
| Kraken2     | 2.1.5   | |
| Centrifuge  | 1.0.4.2 | |
| Centrifuger | 1.0.9   | |
| Ganon2      | 2.1.1   | |
| Sourmash    | 4.9.0   | |
| Sylph       | 0.8.1   | |
| sylph-tax   | 1.7.0   | Other versions silently produce `NO_TAXONOMY` rows against the manuscript's metadata. |

### 3. Verify the install

```bash
for t in kraken2 centrifuge centrifuger ganon sourmash sylph; do
  command -v $t >/dev/null && echo "$t: $(command -v $t)" || echo "$t: MISSING"
done
```

All six should resolve to paths inside the active conda env.

### 4. System tools for the data fetch (typically already present)

`curl`, `jq`, `zstd ≥ 1.4`, `tar`, `sha256sum`. Quick check:

```bash
for c in curl jq zstd tar sha256sum; do
  command -v $c >/dev/null && echo "$c: OK" || echo "$c: MISSING"
done
```

## Quick reproduction

**Seven of the ten manuscript figures (Figs 2-8) plus Table T1 and Table T2 reproduce from the analysis tables shipped with this repository** - no reference databases, no read files, no profiler runs needed. After installing the env, just execute a notebook:

```bash
# cd into the bakeoff project root first
conda activate bakeoff
jupyter nbconvert --to notebook --execute scripts/analysis/analysis_detection.ipynb
```

The figure (and per-rank metric tables) land under `results/<ts>/detection/threshold_<x>/<dataset>/{tables,figures}/`. To render a different dataset or figure, change the `DATASET` (and threshold where relevant) in the notebook's config cell at the top:

| Fig / Table | Notebook | What to set |
| --- | --- | --- |
| 2, 3 | `scripts/analysis/analysis_detection.ipynb` | `DATASET` (mock or simulated key) |
| 4, Table T2 | `scripts/analysis/auxiliary/sensitivity_floor.ipynb` | - (pools all four D6331 mocks; emits Fig 4 panels + the `rare_tier_detection_matrix.csv` that backs Table T2) |
| 5, 6, S2 | `scripts/analysis/analysis_abundance.ipynb` | `DATASET` (mock key) |
| 7, S4-S9 | `scripts/analysis/dyn_heatmap.ipynb` | `DATASET`, `MIN_ABUNDANCE` |
| 8, S10 | `scripts/analysis/dyn_alpha_div.ipynb` | cohorts, `MIN_ABUNDANCE` |
| Table T1 | `scripts/analysis/auxiliary/database_comparison.ipynb` | - |

Figs 7, 8, and their supplementary panels additionally need the DYN cohort cache under `results/prepared/<ts>/dyn-prep/tables/cohorts/` - shipped with the repo. Figs 1, 9, 10, and Fig S1 (read-length, resource benchmarking, recommended-workflow diagram, error rates) require the full pipeline or are hand-drawn; they are not reproducible from the shipped tree.

To sweep a threshold (e.g. the published detection sensitivity panel `0.0 / 0.0001 / 0.001 / 0.01 / 0.1`) without editing the notebook - and across every dataset in one command:

```bash
# Detection - 6 datasets × 5 thresholds = 30 runs
scripts/analysis/sweep_notebook.sh \
    --notebook scripts/analysis/analysis_detection.ipynb \
    --var THRESHOLD_PERCENT \
    --datasets "pacbio_low_input pacbio_standard_input ont_zymo_kit ont_qiagen_kit sim_pacbio sim_ont"

# Abundance - 4 datasets × 5 thresholds = 20 runs (abundance is mock-only)
scripts/analysis/sweep_notebook.sh \
    --notebook scripts/analysis/analysis_abundance.ipynb \
    --var MIN_ABUNDANCE \
    --datasets "pacbio_low_input pacbio_standard_input ont_zymo_kit ont_qiagen_kit"
```

Both use the default `--values "0.0 0.0001 0.001 0.01 0.1"` (override with `--values "..."`).

See the wiki for full reproduction details - [Reproducing the manuscript][wiki-reproduction], including the full-pipeline path (raw reads → profilers → notebooks) for the figures above.

## Full reproduction (from raw reads)

Only needed if you want to regenerate the analysis tables themselves (or reproduce Fig 1, Fig 9, Fig S1). Reference databases and the simulated dataset live on Zenodo (grouped in the [Treangen Lab Bakeoff community](https://zenodo.org/communities/treangen_lab_bakeoff)); mock-community reads are obtained from their original SRA accessions. See the wiki for the full recipe:

- [Data fetch guide][wiki-data-fetch] - Zenodo bundles, manual reassembly, post-extract steps
- [Reproducing the manuscript][wiki-reproduction] - per-figure / per-table walkthroughs (quick and full paths)
- [Database build][wiki-db-build] - rebuilding the unified RefSeq-228 DB from scratch

[wiki-data-fetch]: https://github.com/treangenlab/bakeoff/wiki/Data-Fetch
[wiki-reproduction]: https://github.com/treangenlab/bakeoff/wiki/Reproduction
[wiki-db-build]: https://github.com/treangenlab/bakeoff/wiki/Database-Build
