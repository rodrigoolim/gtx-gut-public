# Example input

Both files are for [SRR37635745](https://www.ncbi.nlm.nih.gov/sra/?term=SRR37635745),
a publicly deposited Crohn's disease patient sample, processed through
16S rRNA amplicon sequencing (QIIME 2 + DADA2 + GreenGenes 13.8, same pipeline
described in the dissertation). Used as one of the eight comparative cases in
the dissertation's Experimentos chapter.

- `SRR37635745_relative_by_levels.csv` -- relative abundance table (all
  taxonomic levels), the `--patient-abundance` input.
- `SRR37635745_shannon_diversity2.csv` -- this sample's Shannon diversity
  index, as produced by `src/shannon.R` on the DADA2 feature table, the
  `--patient-shannon-csv` input.

Provided here purely as ready-to-use input for `generate_report_json.py` --
see the Quickstart section in the main README.

`fastq/` bundles already quality-filtered FASTQ (`*_gcfix.fastq.gz`) for all
eight public samples discussed in the dissertation's Experimentos chapter --
the raw input for the Snakefile, one file (single-end) or file pair
(paired-end) per sample. See the main README's Requirements section for why
quality control (BBDuk) is not part of this repository's workflow.

| Sample | Condition | Layout |
|---|---|---|
| SRR37635745 | Crohn's Disease | single-end |
| SRR33578315 | Type 2 Diabetes Mellitus | paired-end |
| ERR13610605 | Irritable Bowel Syndrome | single-end |
| SRR39497617 | Parkinson's Disease | paired-end |
| SRR39651412 | Depression | paired-end |
| SRR39514725 | Diarrhea | paired-end |
| DRR519923 | Obesity | paired-end |
| SRR37635673 | Ulcerative Colitis | single-end |
