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

- [Singularity/Apptainer](https://apptainer.org/) (tested with singularity-ce 3.11)
  and [Snakemake](https://snakemake.readthedocs.io/) (tested with 8.28), if you
  want the automated FASTQ→JSON pipeline (recommended)
- Python 3.11, `pip install -r requirements.txt`, if you only want to run
  `generate_report_json.py` directly on already-processed input
- R with `vegan` and `qiime2R` (for Shannon diversity, `src/shannon.R`) --
  bundled in the `qiime2r.sif` container built from `qiime2r/qiime2r.def`
- QIIME 2 (tested with 2023.9) with `q2-dada2` and `q2-feature-classifier` --
  bundled in the official `quay.io/qiime2/amplicon:2023.9` image, pulled
  automatically by the Snakefile

`resources/` bundles the GreenGenes classifier needed by taxonomic
classification: `gg-13-8-99-nb-classifier.qza`, pre-trained GreenGenes 13.8
(99%) Naive Bayes, the same one used throughout the dissertation.

`config.example.json` already points to it; copy it to `config.json` and set
`singularity_builds` to a writable directory (where the pulled/built `.sif`
containers will be cached). No API keys are required by anything in this
repository.

**Note on quality control:** this pipeline intentionally starts from
*already filtered* FASTQ (`*_gcfix.fastq.gz`). The private product's exact
BBDuk parameters are not disclosed in the dissertation (proprietary), so
there is no point pretending to reproduce that step -- `examples/fastq/`
bundles the already-filtered FASTQ for all eight samples discussed in the
dissertation's Experimentos chapter. If you're starting from your own raw
FASTQ, run your own quality control (BBDuk or otherwise) first.

## Quickstart

The fastest way to see the classifier work, using one of the eight samples
discussed in the dissertation ([SRR37635745](https://www.ncbi.nlm.nih.gov/sra/?term=SRR37635745),
a Crohn's disease patient), from FASTQ to results, fully automated:

```bash
cp config.example.json config.json   # then set singularity_builds

snakemake all --cores 1 --configfile config.json --config \
    R1=examples/fastq/SRR37635745_S1_L001_R1_001_gcfix.fastq.gz \
    output_dir=results/SRR37635745

cat results/SRR37635745/results.json
```

This runs QIIME 2 import, DADA2 denoising, GreenGenes classification, taxa
collapse, abundance export, Shannon diversity, and the final classification
against the 350-sample reference population -- the same chain of steps
described in the dissertation's Metodologia. It should report Shannon
`Very Low`, F/B ratio `Normal`, and `Enterotype 1 (Bacteroides)`, matching the
dissertation's discussion of this sample. For a paired-end example, add
`R2=examples/fastq/SRR39514725_S1_L001_R2_001_gcfix.fastq.gz` (and swap R1 to
match) to reproduce the Diarrhea sample instead.

**Note on the `taxa` count:** this script classifies every taxon present in
both the patient's abundance table and the reference population (101 for
this sample). The dissertation's comparative table reports a different count
(96) for the same sample, because that number comes from the private
product's curated disease-association database -- a name-matched subset of
taxa, not every taxon in the abundance table. That curated database is
proprietary and not included in this repository (see "What this is NOT"
above), so this script intentionally classifies the full taxon set instead.
Shannon, F/B ratio, and enterotype are unaffected and match exactly.

If you already have a processed sample (a relative-abundance CSV and a
Shannon value) and just want to run the classifier itself, skip straight to
`generate_report_json.py`:

```bash
cd src
python generate_report_json.py \
    --patient-abundance ../examples/SRR37635745_relative_by_levels.csv \
    --patient-id SRR37635745 \
    --patient-shannon-csv ../examples/SRR37635745_shannon_diversity2.csv \
    --reference-dir ../data \
    --output results.json
```

## Usage

For your own samples, the Snakefile expects already quality-filtered FASTQ
(see the note above) -- everything else is automated:

```bash
snakemake all --cores 1 --configfile config.json --config \
    R1=path/to/your_sample_R1_filtered.fastq.gz \
    R2=path/to/your_sample_R2_filtered.fastq.gz \
    output_dir=results/your_sample
```

(Drop `R2=...` for single-end data.) `results/your_sample/results.json` is the
final output. Intermediate QIIME 2 artifacts (`.qza`) and the
relative-abundance/Shannon CSVs are also left in that directory if you want
to inspect them.

If you'd rather run each step by hand (e.g. to swap in your own quality
control), the Snakefile's rules mirror the stages described in the
dissertation's Metodologia: QIIME 2 import → DADA2 denoising → GreenGenes
classification (`q2-feature-classifier`) → taxa collapse → CSV export
(`qiime2r/export.R`) → relative abundance (`src/process_patient.py`) →
Shannon diversity (`src/shannon.R`) → `src/generate_report_json.py`. Run
`snakemake --dag | dot -Tpng > dag.png` (requires graphviz) to see the exact
dependency graph.

**Why Shannon is computed from the raw ASV table, not from
`--patient-abundance`:** that CSV has already been collapsed into named
taxonomic levels (phylum, genus, species, ...), which understates diversity
relative to the raw table DADA2 produces. We measured this on the bundled
example: Shannon from species-level abundances gives 0.867, vs. 2.088 from
the actual ASV table used throughout the dissertation -- a large enough gap
to flip the `Very Low` / `Normal` classification. `shannon.R` always runs on
`raw_derep_table.qza`, before any taxonomic collapsing.

Output schema (`results.json`):

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

- `src/plotting.py`'s `plot_shannon` renders a fixed low/medium/high diversity
  gauge (anchored at `0`, `3`, and `6`) reflecting literature-consensus Shannon
  categories, rather than the empirical distribution of the 350-sample
  reference population. This is a deliberate design choice (a fixed clinical
  reference scale, not a population comparison) and is kept as-is here; the
  dissertation's Metodologia and Limitações discuss it as a methodological
  difference from the population-derived approach used for F/B ratio and
  individual taxa, not as a defect.
- The Tukey/IQR criterion does not clip its fences to the valid range of a
  relative abundance (`[0%, 100%]`), which can make `Very Low` or `Very High`
  statistically unreachable for some taxa. Also discussed in Limitações.
- The F/B ratio is an unbounded ratio and can produce extreme values when
  Bacteroidetes abundance approaches zero. Also discussed in Limitações.

## License

MIT (see `LICENSE`). The reference population data (`data/`) is derived from
publicly available ENA/SRA records (BioProject PRJEB53463); see that
BioProject's own terms for the underlying sequencing data.
