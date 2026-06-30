# bakeoff/scripts — manuscript code

The scripts and notebooks that generate every figure and table in the manuscript ("Critical Evaluation of Long Read Taxonomic Profiling of the Gut Microbiome"). The pipeline is structured as **`process/`** (tool drivers) → **`analysis/`** + **`read-analysis/`** (figure/table notebooks), with reusable helpers in `analysis/utils/` and the database-construction tooling in `build-db/`.

## Layout

```
bakeoff/scripts/
├── process/        tool execution drivers (shell)
├── analysis/       Jupyter notebooks + ETL scripts that produce figures + supplementary tables
│   ├── utils/      analysis-internal Python package (parsers + taxid utils)
│   └── auxiliary/  off-core-path work: secondary analyses + advanced/provenance scripts (outputs ship)
├── build-db/       unified-database construction helpers
└── read-analysis/  read stats, error-rate analysis, taxid→genome mapping
```

## process/  — tool execution

Shell drivers that invoke each of the 6 paper methods on a fastq input. Source:
`script/process/`. Magnet and HYMET drivers are omitted (those tools are not in
the v0.1 main benchmark).

| File                      | Role                                                |
| ------------------------- | --------------------------------------------------- |
| `process.sh`              | long-read driver: `--db {unified\|default} [--out <out>] [--db-dir <db>] (--dataset <fq> \| --dataset-list <file>) --tools <list> [--threads N]`. Datasets are mandatory; at least one `--dataset` or one valid `--dataset-list` entry must be supplied. Defaults: `--out reports/<db>-reports`, `--db-dir data/ref_db`. Supersedes the original `process_uni.sh` / `process_default.sh`. |
| `process_sr.sh`           | short-read (paired Illumina) driver: same `--db`/`--out`/`--db-dir`/`--dataset`/`--dataset-list` flags (here `--dataset` is a directory or glob containing paired-end fastqs), plus `--samples KEY1,KEY2` (or `all`). Supersedes `process_uni_sr.sh` / `process_default_sr.sh`. |
| `postprocess.sh`          | **convenience wrapper** that runs the three post-processing scripts below (`generate_kreport.sh`, `ganon_report.sh`, `sylph-tax.sh`) in sequence with shared flags (`--db`, `--report`, `--db-dir`, `--techs`, `--jobs`). The typical user-facing workflow is `process.sh` then `postprocess.sh`. Continues on per-step failure so a missing tool binary in one sub-script doesn't block the others. `--jobs N` is forwarded to all three for per-tool parallelism. Each sub-script remains usable standalone. |
| `generate_kreport.sh`     | combined kreport post-processor for Centrifuge, Centrifuger, and Sourmash: `--db {unified\|default} [--report <tree>] [--db-dir <db>] [--tools centrifuge,centrifuger,sourmash] [--techs pacbio,ont]`. Reads + writes inside `--report`; new kreports land next to their inputs. Defaults: `--report reports/<db>-reports`, `--db-dir data/ref_db`. Supersedes `cf_kreport.sh` / `cfer_kreport.sh` / `sm_kreport.sh`. |
| `ganon_report.sh`         | post-process Ganon2 → `.tre`                        |
| `sylph-tax.sh`            | combined Sylph taxonomy mapping: `--db {unified\|default} [--report <tree>] [--db-dir <db>]`. Defaults: `--report reports/<db>-reports`, `--db-dir data/ref_db`. `unified` reads a metadata TSV under `--db-dir/default_db/sylph_tax/`; `default` uses sylph-tax's built-in named DBs (GTDB_r220, IMGVR_4.1, FungiRefSeq-2024-07-25). Supersedes `sylph-tax_custom.sh` / `sylph-tax_default.sh`. |

## analysis/utils/  — analysis-internal helpers

Python package imported only by the analysis notebooks + ETL scripts under
`analysis/`. Returns a standardized
`(taxid_raw, taxid, name_raw, name, rank, value, value_type, abundance_raw, ...)`
schema per tool. Lives inside `analysis/` because no other directory consumes it.
Imported as `from utils.parser import ...` / `from utils.util import ...`.

| File                                      | Role                                              |
| ----------------------------------------- | ------------------------------------------------- |
| `analysis/utils/parser.py`                | per-tool report parsers (Kraken2/Centrifuge/Centrifuger/Ganon2/Sourmash/Sylph) and unified ground-truth builder |
| `analysis/utils/util.py`                  | name cleaning, taxid lookup/projection helpers, eval-key + PRF metrics |
| `analysis/utils/results.py`               | timestamped-run discovery (`find_latest_analysis_prep_ts`, etc.) shared by the figure notebooks |
| `analysis/utils/alpha_div.py`             | Shannon / Simpson helpers consumed by `dyn_alpha_div.ipynb` |
| `analysis/utils/__init__.py`              | package init                                      |
| `analysis/utils/output_format_notes.txt`  | quick reference for each tool's native output columns |

## analysis/  — manuscript figure + table notebooks

Each notebook produces specific manuscript artifacts; their outputs land under `results/<TS>/...`.

| Notebook                            | Manuscript artifact                          |
| ----------------------------------- | -------------------------------------------- |
| `analysis_prep.py`                  | CLI: preprocesses raw per-tool reports into a timestamped run dir `results/metadata/<YYYYMMDD_HHMMSS>/analysis-prep/preprocessed/{detection,abundance}/<Tool>_<db>/<file>.csv`, builds ground-truth tables under the same run's `ground_truth/`, and emits `preprocessed/{detection,abundance}/totals.csv` (root + unclassified reads, sourmash k-mer totals). For kreport-derived tools (Kraken2 in both modes; Centrifuge / Centrifuger / Sourmash in detection), the `value` column is recomputed upstream as `(numReads / total) * 100` so the figure notebooks consume already-precise fractions — `abundance_raw` carries the original kreport-rounded value as audit trail. Upstream of every figure notebook below. `--mode both` runs both detection and abundance side-by-side; `--jobs N` for parallelism; `--out` overrides the output root; `--data-groups ZymoMockD6331,simulated[,DYN]` selects which dataset families to process. |
| `analysis_detection.ipynb`          | Reads preprocessed `detection/` CSVs from the latest `analysis_prep.py` run, applies the detection threshold to the already-precise `value` column, computes precision/recall/F1 per (tool, rank). **Figs 2, 3** (and **Fig S2** for the other two mock libraries). Select `DATASET` in the config cell. |
| `analysis_abundance.ipynb`          | Reads preprocessed `abundance/` CSVs from the latest `analysis_prep.py` run, aligns to ground truth, computes Spearman rho / 1-Bray-Curtis / 1-JSD + stacked rel-abundance bars. **Figs 5, 6** (and **Figs S3, S4** for the other two mock libraries). Select `DATASET` in the config cell. |
| `sweep_notebook.sh`                 | Parameter-sweep driver for the figure notebooks. Patches a config variable (default `THRESHOLD_PERCENT` for detection; `MIN_ABUNDANCE` for abundance) across `--values "..."` and runs the notebook once per value via `jupyter nbconvert`; per-iteration outputs land in the notebook's threshold-namespaced result dirs. Add `--datasets "key1 key2 ..."` to repeat the sweep across multiple DATASET keys (outer loop), and `--set VAR=VALUE` (repeatable) to pin extra config overrides every iteration (`SAVE=True` is forced). `scripts/analysis/sweep_notebook.sh --help` for the full option set. |
| `dyn_prep.py`                       | DYN cohort cache builder — walks per-tool reports under `--reports-{default,unified}` and writes one long-form TSV per (cohort, tool, db_mode, rank) into `results/metadata/<TS>/dyn-prep/tables/cohorts/`. Upstream of `dyn_heatmap.ipynb` and `dyn_alpha_div.ipynb` (both auto-discover the latest `<TS>` dir). Supersedes the former `dyn_load_all.py` + `dyn_alpha_div_compute.py` pair. CLI: `python dyn_prep.py --help`. |
| `dyn_heatmap.ipynb`                 | **Fig 7** (DYN species/genus heatmap; **Figs S6, S7** for ONT cohorts) — reads cohort cache directly. |
| `dyn_alpha_div.ipynb`               | **Fig 8** (DYN alpha-diversity; **Fig S8** for ONT Qiagen) — computes Shannon/Simpson from the cohort cache via `utils/alpha_div.py`. |
| `resource_benchmark.ipynb`          | **Fig 9** (runtime + memory + accuracy-cost trade-offs). |

Fig 1 (read-length distributions) is produced by `read-analysis/read_analysis.ipynb`.
Fig 10 (workflow diagram) is hand-drawn, not script-generated.

## analysis/auxiliary/  — off-core-path analyses + provenance scripts

Self-contained notebooks/scripts that are **not** part of the core three-command
reproduction path (`process.sh` → `postprocess.sh` → `analysis_prep.py`): secondary
analyses that support the Discussion but aren't tied to a main figure, plus
advanced/provenance scripts whose outputs are shipped (so a normal reproducer never
runs them). Each reuses the `analysis/utils/` package; notebook outputs land under
`results/<TS>/auxiliary/<name>/`.

| File                                   | Role                                                |
| -------------------------------------- | --------------------------------------------------- |
| `auxiliary/database_comparison.ipynb`  | **Supplementary Table S1** (taxonomic breadth across DBs) — DB-only analysis; doesn't depend on the profiler pipeline. |
| `auxiliary/sensitivity_floor.ipynb`    | **Fig 4** (detection operating characteristic: PR plane + F1 vs threshold across the sweep, pooled over all four D6331 mocks under the unified DB) and **Supplementary Table S2** (low-abundance D6331 detection across mocks). Also emits FN counts (per-tool, by abundance-tier × tool-class, per-taxon), a detection-floor table, a rare-tier detection matrix (CSV), the sensitivity-precision tradeoff at no threshold, and a PacBio-vs-ONT split of the rare-tier matrix — showing missed detections concentrate at the rarest taxa. Mock-only (the uniform-depth simulations contain no rare taxa to test a floor). Consumes preprocessed `detection/` CSVs from `analysis_prep.py`. |
| `auxiliary/abundance_sweep.ipynb`      | Abundance-accuracy counterpart to the detection sweep: scores the three concordance metrics (Spearman rho, 1-Bray-Curtis, 1-JSD) across the same threshold sweep on the four D6331 mocks. Shows that magnitude-based metrics (1-BC, 1-JSD) are threshold-invariant (removed false positives carry negligible abundance mass) while only rank-based rho responds, with Centrifuge the outlier. Three-panel figure + detail/range CSVs. Consumes preprocessed `abundance/` CSVs from `analysis_prep.py`. |
| `auxiliary/extract_simulated_truth.py` | **Provenance / advanced users.** Builds the simulated source ground-truth CSVs (`simulated_<tech>_gt.csv`) by counting pbsim reads per contig and aggregating to species-level relative abundance. The truth tables ship with the dataset, so this is only needed to *regenerate* them from the raw MIMIC/pbsim inputs (joint FASTA, ete3 sqlite, the per-contig fastqs). `python auxiliary/extract_simulated_truth.py --help` for flags. |

## build-db/  — unified-database construction + Zenodo fetch

Used once to build `data/ref_db/refseq03032025/` (the unified RefSeq v228 reference); the fetch script provides the bundled-DB shortcut for reproducers who don't want to rebuild.

| File                                  | Role                                          |
| ------------------------------------- | --------------------------------------------- |
| `fetch_bakeoff_data.sh`               | Download the published Zenodo bundles for the unified + Ganon2-default DBs and extract them into the project root. Edit `ZENODO_RECORDS` at the top of the script if record IDs change; the records currently in use are grouped at [zenodo.org/communities/treangen_lab_bakeoff](https://zenodo.org/communities/treangen_lab_bakeoff). |
| `database_dic.py`                     | accession → taxid map for selected RefSeq assemblies (writes `seqid2taxid.map`) |
| `build_unified_lineage_table.py`      | top-level lineage table for the unified DB    |
| `build_ganon_lineage_table.py`        | Ganon2-format lineage table                   |
| `build_sourmash_lineage_table.py`     | Sourmash-format lineage CSV (ident → taxid + NCBI lineage from `nodes.dmp` / `names.dmp`) |
| `build_sylph-tax_metadata_table.py`   | sylph-tax metadata for custom Sylph DB        |

## read-analysis/  — read stats, error-rate analysis, taxid→genome mapping

Read-stats outputs default to `bakeoff/data/read_stats/`; error-rate outputs to `bakeoff/results/error-rate/`.

| File                       | Role                                                      |
| -------------------------- | --------------------------------------------------------- |
| `read_stats.py`            | per-dataset read stats from a FASTQ list via NanoStat (num_reads, total_bases, mean/median/N50 length, mean Q); optional grouped CSV → `data/read_stats/`. Feeds **Table 1** and `analysis/resource_benchmark.ipynb` (**Fig 9**). |
| `read_info.py`             | per-read length tables from a FASTQ list via pysam → `data/read_stats/read-lengths/{group}__{dataset}.csv`. Parallel (`--workers`). Upstream of `read_analysis.ipynb`. |
| `read_analysis.ipynb`      | grouped read-length boxplot from `read_info.py` output. **Fig 1.** |
| `taxid2genome.py`          | NCBI taxid list → genome FASTA mapping (per-file / grouped / best). **Optional** alternative reference source for `error_rates.py` (`--ref-csv`); builds a RefSeq-by-taxid panel. *Not* used for the manuscript's D6331 Fig S1, which aligned against per-strain ZymoBIOMICS FASTAs via `--ref-dir`. |
| `error_rates.py`           | combined ref → minimap2 → pysam → `results/error-rate/error-rate-results/per_read_{dataset}.csv`. Single CLI consolidating the former `aligner.py` + `error_rate_per_read.py`. |
| `error_analysis.ipynb`     | aggregates `results/error-rate/error-rate-results/per_read_*.csv` into mismatch / insertion / deletion boxplots. **Fig S1.** |

