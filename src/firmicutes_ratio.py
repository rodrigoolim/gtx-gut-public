import sys
import pandas as pd
import numpy as np

def get_phylum_counts(df, phylum):
    phylum_level = df[df['level'] == 'p']
    phylum_rows = phylum_level[phylum_level['taxon'].str.endswith(phylum)]

    if len(phylum_rows) == 1:
        counts = {}
        for _, row in phylum_rows.iterrows():
            for key, val in row.items():
                if not key in ['taxon', 'level']:
                    counts[key] = float(val)
    
        return counts
    else:
        return {}

def remove_outliers(vals: dict, mult = 3):
    np_vec = np.array(list(vals.values()))
    mean = np.mean(np_vec)
    std = np.std(np_vec)

    uplim = mean + mult*std
    lolim = mean - mult*std

    return {k: v for k, v in vals.items() if v >= lolim and v <= uplim}

def ratio_nonewise(vals1: dict, vals2: dict):

    ratios = {}
    common_keys = set(vals1.keys()).intersection(vals2.keys())
    for k in common_keys:
        v1 = vals1[k]
        v2 = vals2[k]
        if v2 > 0:
            ratios[k] = (v1 / v2, v2, v1)
    
    return ratios

if __name__ == "__main__":
    pop_abundances_path = sys.argv[1]
    patient_abundances_path = sys.argv[2]
    patient_n = int(sys.argv[3])
    output_path = sys.argv[4]

    pop_abundances_df = pd.read_csv(pop_abundances_path)
    patient_abundances_df = pd.read_csv(patient_abundances_path)
    patient_id = patient_abundances_df.columns[patient_n+1]
    patient_df = patient_abundances_df[['taxon', 'level', patient_id]]
    
    print(pop_abundances_df.head())
    print(patient_df.head())

    pop_firm_counts = remove_outliers(get_phylum_counts(pop_abundances_df, 'Firmicutes'))
    pop_bact_counts = remove_outliers(get_phylum_counts(pop_abundances_df, 'Bacteroidetes'))

    pop_firm_bac_ratios = ratio_nonewise(pop_firm_counts, pop_bact_counts)

    patient_firm = get_phylum_counts(patient_df, 'Firmicutes')[patient_id]
    patient_bact = get_phylum_counts(patient_df, 'Bacteroidetes')
    if patient_id in patient_bact:
        patient_bact = patient_bact[patient_id]
        patient_firm_bac_ratio = patient_firm / patient_bact
    else:
        patient_bact = 0.0
        patient_firm_bac_ratio = np.inf
        
    output = open(output_path, 'w')
    output.write('"","firm_ratio","bact","firm"\n')
    for key, val in pop_firm_bac_ratios.items():
        ratio, bact, firm = val
        output.write('"'+key+'",'+str(ratio)
            +','+str(bact)+','+str(firm)+'\n')
    output.write('"patient.'+patient_id+'",'+str(patient_firm_bac_ratio)
        +','+str(patient_bact)+','+str(patient_firm)+'\n')