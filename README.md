# GTX-GUT (reproducibility companion)

This repository accompanies the dissertation *"GTX-GUT: um pipeline automatizado
para análise do microbioma intestinal"*. It contains the statistical/methodological
core described and evaluated in the dissertation: taxonomic abundance processing,
the reference population used for population-level comparisons, and the Tukey/IQR
classifier used to flag Shannon diversity, Firmicutes/Bacteroidetes (F/B) ratio, and
individual-taxon abundances as `Very Low` / `Normal` / `Very High` relative to a
healthy reference population.

Running it on a processed sample produces a single JSON file with those results.
**It does not produce a clinical report.**

## What this is NOT

This is a stripped-down, reproducibility-focused extract of a larger commercial
pipeline (GTX-GUT / GTX-BOSTA). The following pieces belong to that commercial
product and are **not included here**:

- The curated disease/condition–taxon association knowledge base.
- AI-assisted (GPT) bacteria description and clinical-text synthesis.
- `.docx` clinical report templates and generation.
- Any real patient/client data.

If you're looking for a way to generate the same statistical outputs described in
the dissertation (Shannon, F/B ratio, enterotype, per-taxon classification) for
your own 16S rRNA amplicon samples, this repository is for that. If you're looking
for the full commercial clinical-report product, this is not it.

## Reference population

`data/` contains the reference population used throughout the dissertation:
350 healthy adults' fecal samples, derived from BioProject
[PRJEB53463](https://www.ebi.ac.uk/ena/browser/view/PRJEB53463) (ENA/SRA,
*USDA-ARS Cross-Sectional Nutritional Phenotyping Study*), after deduplication by
individual and quality control (see the dissertation, Metodologia, for the full
filtering funnel: 530 raw runs → 501 → 401 → 372 → 350 final samples).

These files are the exact ones used to generate the dissertation's results and are
what `generate_report_json.py` reads by default. You do not need to regenerate
them to use the classifier.

`src/process_population.py` is included for transparency about *how* that
population was built from the raw BioProject data, not as a required step. If you
re-run it from the raw abundance tables, note that it will **not** reproduce
exactly 350 samples: the original pipeline set a minimum-read-count QC threshold
based on real client samples processed by the commercial product (private data,
not reproducible here). This public version instead derives that threshold from
the reference population's own read-depth distribution (10th percentile), which
is methodologically defensible but numerically different. Use the pre-built files
in `data/` if you want to match the dissertation exactly.

## Requirements

- Python 3.11, `pip install -r requirements.txt`
- [QIIME 2](https://qiime2.org/) (tested with 2023.9) with the `q2-dada2` and
  `q2-feature-classifier` plugins
- [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (BBMap suite)
- R with `vegan` and `qiime2R` (for Shannon diversity calculation, `src/shannon.R`)

`resources/` bundles the two reference files needed by the steps above:

- `illumina_adapters.fa` -- standard Illumina adapter sequences, for BBDuk.
- `gg-13-8-99-nb-classifier.qza` -- pre-trained GreenGenes 13.8 (99%) Naive
  Bayes classifier for `q2-feature-classifier`, the same one used throughout
  the dissertation.

`config.example.json` already points to both; copy it to `config.json` and
adjust `singularity_builds` (and any other path) for your own environment. No
API keys are required by anything in this repository.

## Quickstart

`examples/` includes an already-processed public sample
([SRR37635745](https://www.ncbi.nlm.nih.gov/sra/?term=SRR37635745), a Crohn's
disease patient, one of the eight cases discussed in the dissertation) so you
can try the classifier immediately, without running QIIME 2/DADA2 yourself
first:

```bash
cd src
python generate_report_json.py \
    --patient-abundance ../examples/SRR37635745_relative_by_levels.csv \
    --patient-id SRR37635745 \
    --patient-shannon-csv ../examples/SRR37635745_shannon_diversity2.csv \
    --reference-dir ../data \
    --output results.json
```

This should classify 101 taxa and report Shannon `Very Low`, F/B ratio
`Normal`, and `Enterotype 1 (Bacteroides)` -- matching the dissertation's
discussion of this sample.

## Usage

1. Process your raw FASTQ through the standard 16S amplicon steps (this
   repository does not orchestrate these for you -- see the dissertation,
   Metodologia, for the exact QIIME 2/DADA2/GreenGenes commands used):
   - BBDuk quality control
   - QIIME 2 import + DADA2 denoising
   - GreenGenes taxonomic classification
   - Export to a `taxon,level,<sample_id>` relative-abundance CSV
     (`src/process_patient.py`, `src/make_pop_and_patient_df.py`)
   - Shannon diversity via `src/shannon.R` on the DADA2 feature table
     (**not** on the relative-abundance CSV above -- see note below)

2. Generate the results JSON, pointing `--patient-shannon-csv` straight at
   shannon.R's output (no need to open the file and copy the number by hand):

   ```bash
   Rscript src/shannon.R path/to/raw_derep_table.qza path/to/shannon_diversity2.csv

   python src/generate_report_json.py \
       --patient-abundance path/to/Specie_full_abundance.relative_by_levels.csv \
       --patient-id YOUR_SAMPLE_ID \
       --patient-shannon-csv path/to/shannon_diversity2.csv \
       --reference-dir data \
       --output results.json
   ```

   If you already have the Shannon value from elsewhere, `--patient-shannon 3.05`
   works too (the two flags are mutually exclusive).

   **Why Shannon can't be computed from `--patient-abundance` directly:** that
   file has already been collapsed into named taxonomic levels (phylum,
   genus, species, ...), which understates diversity relative to the raw
   ASV table DADA2 produces. We measured this on the bundled example: Shannon
   from species-level abundances gives 0.867, vs. 2.088 from the actual ASV
   table used throughout the dissertation -- a large enough gap to flip the
   `Very Low` / `Normal` classification. Always compute Shannon from the raw
   feature table via `shannon.R`, not by re-deriving it from the collapsed CSV.

   Output schema:

   ```json
   {
     "sample_id": "...",
     "shannon_diversity": {"value": 0.0, "classification": "Normal", "percentile": 0.0, "modified_zscore": 0.0},
     "firmicutes_bacteroidetes_ratio": {"value": 0.0, "classification": "Normal", "percentile": 0.0, "modified_zscore": 0.0},
     "enterotype": "Enterotype 1 (Bacteroides)",
     "reference_stats": {
       "shannon": {"n": 350, "median": 0.0, "q1": 0.0, "q3": 0.0, "iqr": 0.0},
       "fb_ratio": {"n": 350, "median": 0.0, "q1": 0.0, "q3": 0.0, "iqr": 0.0}
     },
     "taxa": [
       {"name": "...", "level": "p", "relative_abundance": 0.0, "classification": "Normal",
        "percentile": 0.0, "modified_zscore": 0.0, "reference": {"q1": 0.0, "q3": 0.0, "median": 0.0}}
     ],
     "reference_population": {"bioproject": "PRJEB53463", "n_samples": 350}
   }
   ```

   With `reference_stats` and each taxon's `reference` block, the Very Low /
   Normal / Very High classification (`< Q1 - 1.5×IQR` / between / `> Q3 +
   1.5×IQR`) can be independently verified without re-running anything.

## Known limitations (documented and discussed in the dissertation)

- `src/plotting.py`'s `plot_shannon` renders its population-comparison plot from
  three fixed anchor points (`0`, `3`, `6`) rather than the real reference
  population, even though the real population is available. This is
  intentionally kept as-is to match the current behavior of the commercial
  product being described; see the dissertation's Limitações chapter for the
  full discussion (including why it matters and what fixing it would involve).
- The Tukey/IQR criterion does not clip its fences to the valid range of a
  relative abundance (`[0%, 100%]`), which can make `Very Low` or `Very High`
  statistically unreachable for some taxa. Also discussed in Limitações.
- The F/B ratio is an unbounded ratio and can produce extreme values when
  Bacteroidetes abundance approaches zero. Also discussed in Limitações.

## License

MIT (see `LICENSE`). The reference population data (`data/`) is derived from
publicly available ENA/SRA records (BioProject PRJEB53463); see that
BioProject's own terms for the underlying sequencing data.
