# Reference population data

Derived from BioProject [PRJEB53463](https://www.ebi.ac.uk/ena/browser/view/PRJEB53463)
(ENA/SRA), *USDA-ARS Cross-Sectional Nutritional Phenotyping Study Gut Microbiota*.
See the dissertation (Metodologia, "Construção e Uso da População de Referência")
for the full filtering funnel and methodology.

| Runs deposited | After near-empty exclusion | After dedup by individual | After QC | Final |
|---|---|---|---|---|
| 530 | 501 | 401 | 372 | **350** |

## Files

- `feces_Specie_full_abundance_v3.csv` / `.norm1k.csv` -- raw and normalized abundance
  tables, one column per individual (401 samples, post-deduplication, pre-QC).
  Provided for transparency of `process_population.py`; not required for
  `generate_report_json.py`.
- `feces_Specie_full_abundance_v3.norm1k.all_levels.csv` -- all taxonomic levels
  (kingdom through species), absolute counts, 350 final samples.
- `feces_Specie_full_abundance_v3.norm1k.relative_by_levels.csv` /
  `...taxon_indexed.csv` -- relative abundance per taxon per sample, 350 final
  samples. **This is what `generate_report_json.py` reads** to classify each
  patient taxon.
- `feces_Specie_full_abundance_v3.norm1k.fb_ratios.csv` -- Firmicutes/Bacteroidetes
  ratio for each of the 350 reference samples.
- `feces_Specie_full_abundance_v3.norm1k.relative_by_levels.shannon.csv` -- Shannon
  diversity index for each of the 350 reference samples.
- `feces_Specie_full_abundance_v3.norm1k.pop_stats.json` -- per-sample QC stats
  (read counts, outlier flags) computed by `process_population.py`.

All 350-sample files share the exact same set of samples (verified in the
dissertation, Metodologia).
