import sys
import pandas as pd
import numpy as np

if __name__ == "__main__":
    pop_abundances_path = sys.argv[1]
    patient_abundances_path = sys.argv[2]
    patient_n = int(sys.argv[3])
    output_path = sys.argv[4]

    pop_abundances_df = pd.read_csv(pop_abundances_path)
    pop_abundances_df = pop_abundances_df.fillna(0.0)
    patient_abundances_df = pd.read_csv(patient_abundances_path)
    patient_id = patient_abundances_df.columns[patient_n+1]
    patient_df = patient_abundances_df[['taxon', 'level', patient_id]]

    patient_percs_dict = {
        row['taxon']: row[patient_id] for index, row in patient_df.iterrows()
    }

    pop_abundances_df[patient_id] = pop_abundances_df['taxon'].apply(
        lambda t: patient_percs_dict[t] if t in patient_percs_dict else np.nan)
    
    taxa_in_patient = set(patient_df['taxon'].tolist())
    taxa_in_pop = set(pop_abundances_df['taxon'].tolist())

    patient_only_taxa = [t for t in taxa_in_patient if not t in taxa_in_pop]
    patient_only_rows = [{'taxon': row['taxon'], 'level': row['level'], patient_id: row[patient_id]} 
        for _, row in patient_df.iterrows() if row['taxon'] in patient_only_taxa]
    
    patient_only_df = pd.DataFrame(patient_only_rows)
    pop_and_patient_df = pd.concat([pop_abundances_df, patient_only_df])

    #is_in_patient = pop_and_patient_df[~pop_and_patient_df[patient_id].isna()]
    pop_and_patient_df = pop_and_patient_df.fillna(0.0)
    pop_and_patient_df.to_csv(output_path, sep=',', index=False)
    del pop_and_patient_df['level']
    pop_and_patient_df.set_index('taxon', inplace=True)
    pop_and_patient_df.to_csv(output_path.replace('.csv', '.taxon_indexed.csv'), sep=',')