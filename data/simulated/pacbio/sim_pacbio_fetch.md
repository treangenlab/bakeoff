# Simulated PacBio HiFi — fetch

Synthetic PacBio HiFi reads simulated by PBSIM3 from a 127-species gut community fitted with [MIMIC](https://github.com/treangenlab/Mimic). Feeds the PacBio panels of Fig 2 (simulated detection).

Hosted on Zenodo in the [Treangen Lab Bakeoff community](https://zenodo.org/communities/treangen_lab_bakeoff):

| File | Description | Zenodo |
| --- | --- | --- |
| `sim_pacbio.fastq` | PBSIM3 reads (PacBio HiFi profile fit from a real human-gut HiFi sample) | [`20496925`](https://zenodo.org/records/20496925) |
| `sim_pacbio_truth_abundance.tsv` | per-species realized truth (Species, Reads, Theoretical_Abundance(%)) | ships in this directory |
| `sim_pacbio_gt.tsv` | strain-level manifest (the 233 source assemblies) | ships in this directory |

## Fetch

From this directory:

```bash
curl -L -o sim_pacbio.fastq \
  "https://zenodo.org/records/20496925/files/sim_pacbio.fastq.gz?download=1"
```

Once present here, the path matches the entry in `data/datasets_simulated_pacbio.txt` and `data/datasets_all.txt`, and the truth table is already aligned for use by `analysis_prep.py`.
