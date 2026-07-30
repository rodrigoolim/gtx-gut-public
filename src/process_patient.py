
import sys

import pandas as pd

from abundance_processing import generic_counts_processing

if __name__ == "__main__":
    original_abundances_path = sys.argv[1]

    updated_df, all_levels_path, sample_names, columns_order, filtering_log_path = generic_counts_processing(
        original_abundances_path)

    updated_df.to_csv(all_levels_path, index=False)

    level_dfs = []
    for level, level_rows in updated_df.groupby('level'):
        local_total = {s: level_rows[s].sum() for s in sample_names}
        for s in sample_names:
            level_rows[s] = level_rows[s].apply(
                lambda x: float(x)/local_total[s])
        print(level_rows)
        level_dfs.append(level_rows)
    updated_df = pd.concat(level_dfs)

    #sample_totals = {s: updated_df[s].sum() for s in sample_names}
    #for sample_name in sample_names:
    #    updated_df[sample_name] = updated_df[sample_name].apply(
    #        lambda x: float(x)/sample_totals[sample_name])
    columns_order = ['taxon', 'level'] + sample_names
    updated_df = updated_df[columns_order]
    updated_df.to_csv(original_abundances_path.replace('.csv', '.relative_by_levels.csv'), index=False)