# ZymoBIOMICS D6331 — PacBio HiFi fetch

ZymoBIOMICS Gut Microbiome Standard: https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard

Two SRA accessions, both prepared with the Shoreline Breaker DNA extraction:

| Accession | Description | SRA |
| --- | --- | --- |
| SRR13128013 | low input | <https://www.ncbi.nlm.nih.gov/sra/?term=SRR13128013> |
| SRR13128014 | standard input | <https://www.ncbi.nlm.nih.gov/sra/?term=SRR13128014> |

## Fetch

`sra-tools` provides `prefetch` + `fasterq-dump`. Run from this directory:

```bash
for acc in SRR13128013 SRR13128014; do
    prefetch "$acc"
    fasterq-dump --split-files --threads 8 "$acc"
    rm -rf "$acc"     # clean the .sra cache
done
```

This produces `SRR13128013.fastq` and `SRR13128014.fastq` here, matching the paths in `data/datasets_mock_pacbio.txt` / `data/datasets_all.txt`.
