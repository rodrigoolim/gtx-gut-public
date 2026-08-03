#!/usr/bin/env python3
"""
Generates a JSON with the statistical results (Shannon diversity, Firmicutes/
Bacteroidetes ratio, enterotype, and per-taxon Tukey/IQR classification) for a
single processed sample, compared against the public reference population
(BioProject PRJEB53463, N=350 healthy adults; see data/README.md).

This script intentionally stops at the statistics layer: it does not perform
disease-association matching, AI-based text synthesis, or .docx report
generation -- those belong to the commercial product this repository was
extracted from, and are not included here.

Usage:
    python generate_report_json.py \\
        --patient-abundance path/to/Specie_full_abundance.relative_by_levels.csv \\
        --patient-id SAMPLE_ID \\
        --patient-shannon-csv path/to/shannon_diversity2.csv \\
        --output results.json

Shannon diversity must be computed beforehand with shannon.R on the patient's
DADA2 feature table (the same way the reference population's Shannon values
were computed -- see data/README.md); this script does not recompute it from
--patient-abundance, since that file has already been collapsed to named
taxonomic levels and would understate diversity relative to the raw ASV table.
Pass shannon.R's output CSV directly via --patient-shannon-csv, or a bare
number via --patient-shannon if you already have the value.
"""
import argparse
import json
import re

import numpy as np
import pandas as pd

from util import freq_stats

ENTEROTYPE_GENERA = ["bacteroides", "prevotella", "ruminococcus"]
ENTEROTYPE_NAMES = {
    "bacteroides": "Enterotype 1 (Bacteroides)",
    "prevotella": "Enterotype 2 (Prevotella)",
    "ruminococcus": "Enterotype 3 (Ruminococcus)",
}


def get_phylum_value(df, patient_col, phylum):
    phylum_rows = df[(df["level"] == "p") & (df["taxon"].str.endswith(phylum))]
    if len(phylum_rows) == 1:
        return float(phylum_rows.iloc[0][patient_col])
    return 0.0


def compute_enterotype(patient_df, patient_col):
    """Dominant genus among Bacteroides/Prevotella/Ruminococcus. Sums all
    matching lineages per genus (see dissertation Sec. 4.7 for the rationale
    -- this differs slightly from the original product, which kept only the
    last matching lineage)."""
    abundances = {g: 0.0 for g in ENTEROTYPE_GENERA}
    for _, row in patient_df.iterrows():
        taxon = str(row["taxon"])
        for genus in ENTEROTYPE_GENERA:
            if re.search(r"__" + genus + r"$", taxon, re.IGNORECASE):
                abundances[genus] += float(row[patient_col])

    dominant_genus = max(abundances, key=abundances.get)
    if abundances[dominant_genus] > 0:
        return ENTEROTYPE_NAMES[dominant_genus], abundances
    return "Other", abundances


def reference_stats_block(values):
    values = np.asarray(values, dtype=float)
    q1, q3 = np.percentile(values, [25, 75])
    return {
        "n": len(values),
        "median": round(float(np.median(values)), 3),
        "q1": round(float(q1), 3),
        "q3": round(float(q3), 3),
        "iqr": round(float(q3 - q1), 3),
    }


def read_shannon_csv(path):
    """Reads the single-value CSV produced by shannon.R (header + one
    "<row_label>,<value>" line) and returns the float. Deliberately not
    imported from plotting.py, which pulls in matplotlib -- this script has
    no other need for it."""
    import csv
    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader)  # header
        for row in reader:
            if len(row) == 2:
                return float(row[1])
    raise ValueError(f"No Shannon value found in {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--patient-abundance", required=True,
        help="Patient's relative abundance CSV (taxon, level, <patient_id> columns) -- "
             "output of process_patient.py / make_pop_and_patient_df.py.")
    parser.add_argument("--patient-id", required=True,
        help="Column name of the patient in --patient-abundance.")
    shannon_group = parser.add_mutually_exclusive_group(required=True)
    shannon_group.add_argument("--patient-shannon", type=float,
        help="Patient's Shannon diversity index, already computed.")
    shannon_group.add_argument("--patient-shannon-csv",
        help="Path to the CSV produced by shannon.R for this patient -- read directly "
             "instead of having to copy the value by hand.")
    parser.add_argument("--reference-dir", default="../data",
        help="Directory with the reference population files (default: ../data).")
    parser.add_argument("--output", required=True, help="Output JSON path.")
    args = parser.parse_args()

    patient_shannon = args.patient_shannon
    if patient_shannon is None:
        patient_shannon = read_shannon_csv(args.patient_shannon_csv)

    ref_dir = args.reference_dir
    patient_df = pd.read_csv(args.patient_abundance)
    patient_col = args.patient_id

    # --- Reference population ---
    ref_abundance = pd.read_csv(f"{ref_dir}/feces_Specie_full_abundance_v3.norm1k.relative_by_levels.csv")
    ref_sample_cols = [c for c in ref_abundance.columns if c not in ("taxon", "level")]

    ref_shannon_df = pd.read_csv(f"{ref_dir}/feces_Specie_full_abundance_v3.norm1k.relative_by_levels.shannon.csv")
    ref_shannon = ref_shannon_df[ref_shannon_df["Unnamed: 0"] != "level"]["shannon_diversity"].astype(float).tolist()

    ref_fb_df = pd.read_csv(f"{ref_dir}/feces_Specie_full_abundance_v3.norm1k.fb_ratios.csv")
    ref_fb = ref_fb_df["firm_ratio"].astype(float).tolist()

    # --- Patient metrics ---
    patient_firm = get_phylum_value(patient_df, patient_col, "Firmicutes")
    patient_bact = get_phylum_value(patient_df, patient_col, "Bacteroidetes")
    patient_fb_ratio = (patient_firm / patient_bact) if patient_bact > 0 else float("inf")

    enterotype, enterotype_abundances = compute_enterotype(patient_df, patient_col)

    # --- Classification (Tukey/IQR via freq_stats, same logic as the product's classifier) ---
    shannon_stats = freq_stats(ref_shannon, patient_shannon)
    fb_stats = freq_stats(ref_fb, patient_fb_ratio)

    taxa_results = []
    for _, row in patient_df.iterrows():
        taxon = row["taxon"]
        patient_value = float(row[patient_col])
        ref_row = ref_abundance[ref_abundance["taxon"] == taxon]
        if ref_row.empty:
            continue
        ref_values = ref_row.iloc[0][ref_sample_cols].astype(float).tolist()
        stats = freq_stats(ref_values, patient_value)
        taxa_results.append({
            "name": taxon,
            "level": ref_row.iloc[0]["level"],
            "relative_abundance": patient_value,
            "classification": stats["classification"],
            "percentile": round(stats["patient_percentile"], 2),
            "modified_zscore": round(stats["mod_zscore"], 3),
            "reference": {
                "q1": round(stats["q1"], 4),
                "q3": round(stats["q3"], 4),
                "median": round(stats["median"], 4),
            },
        })

    result = {
        "sample_id": patient_col,
        "shannon_diversity": {
            "value": round(patient_shannon, 3),
            "classification": shannon_stats["classification"],
            "percentile": round(shannon_stats["patient_percentile"], 2),
            "modified_zscore": round(shannon_stats["mod_zscore"], 3),
        },
        "firmicutes_bacteroidetes_ratio": {
            "value": None if patient_fb_ratio == float("inf") else round(patient_fb_ratio, 3),
            "classification": fb_stats["classification"],
            "percentile": round(fb_stats["patient_percentile"], 2),
            "modified_zscore": round(fb_stats["mod_zscore"], 3),
        },
        "enterotype": enterotype,
        "reference_stats": {
            "shannon": reference_stats_block(ref_shannon),
            "fb_ratio": reference_stats_block(ref_fb),
        },
        "taxa": taxa_results,
        "reference_population": {
            "bioproject": "PRJEB53463",
            "n_samples": len(ref_sample_cols),
        },
    }

    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Wrote {args.output} ({len(taxa_results)} taxa classified)")


if __name__ == "__main__":
    main()
