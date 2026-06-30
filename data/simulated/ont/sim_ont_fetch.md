# Simulated ONT — fetch

Synthetic ONT reads simulated by PBSIM3 from the same 127-species gut community as `sim_pacbio`, with read-length and quality profiles fit from a real PromethION human-gut sample (SRA accession ERR15285694, sequenced with SQK-LSK114 / R10.4.1 chemistry). Feeds the ONT panels of Fig 2 (simulated detection).

Hosted on Zenodo in the [Treangen Lab Bakeoff community](https://zenodo.org/communities/treangen_lab_bakeoff):

| File | Description | Zenodo |
| --- | --- | --- |
| `sim_ont.fastq` | PBSIM3 reads (ONT R10.4.1 profile fit from ERR15285694) | [`20496925`](https://zenodo.org/records/20496925) |
| `sim_ont_truth_abundance.tsv` | per-species realized truth (Species, Reads, Theoretical_Abundance(%)) | ships in this directory |
| `sim_ont_gt.tsv` | strain-level manifest (the 233 source assemblies) | ships in this directory |

## Fetch

From this directory:

```bash
curl -L -o sim_ont.fastq \
  "https://zenodo.org/records/20496925/files/sim_ont.fastq.gz?download=1"
```

Once present here, the path matches the entry in `data/datasets_simulated_ont.txt` and `data/datasets_all.txt`, and the truth table is already aligned for use by `analysis_prep.py`.
