import argparse
import json

import pandas as pd

from abundance_processing import calc_pop_stats, find_pop_fb_ratios, generic_counts_processing
from taxonomy import archaea_kingdom, lower_phyla

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build the reference population statistics (Shannon, F/B ratio, "
                    "per-taxon quartiles) used by the Tukey/IQR classifier.")
    parser.add_argument("abs_abundances_path",
        help="Normalized raw abundance CSV (all reference samples as columns).")
    parser.add_argument("original_abundances_path",
        help="Raw (non-normalized) abundance CSV, same sample set as above.")
    parser.add_argument("fb_ratios_output",
        help="Output path for the population F/B ratio CSV.")
    parser.add_argument("--min-read-count", type=float, default=None,
        help="Minimum total read count for a reference sample to be kept. "
             "Defaults to the 10th percentile of the reference population's own "
             "read-depth distribution (see calc_pop_stats).")
    args = parser.parse_args()

    abs_df = pd.read_csv(args.abs_abundances_path, sep=',')
    sample_totals = {c: int(abs_df[c].sum()) for c in abs_df.columns if not c in ['taxon', 'level', 'Unnamed: 0', '']}

    original_abundances_path = args.original_abundances_path
    fb_ratios_output = args.fb_ratios_output

    updated_df, all_levels_path, sample_names, columns_order, filtering_log_path = generic_counts_processing(
        original_abundances_path)

    updated_df['is_archea'] = updated_df['taxon'].apply(
        lambda t: archaea_kingdom in str(t)
    )
    updated_df = updated_df[~updated_df['is_archea']]
    del updated_df['is_archea']

    all_levels_raw_path = original_abundances_path.replace('.csv', '.all_levels.raw.csv')
    updated_df.to_csv(all_levels_raw_path, index=False)

    pop_stats = calc_pop_stats(all_levels_raw_path, sample_totals, filtering_log_path,
                                min_read_count=args.min_read_count)
    good_samples = []
    for sample_name, stats in pop_stats.items():
        if not any([stats[bad_indicator] for bad_indicator in ['low_read_count', 'firmicute_outlier', 'bacteroidetes_outlier']]):
            good_samples.append(sample_name)
    json.dump(pop_stats, open(original_abundances_path.replace('.csv', '.pop_stats.json'), 'w'), indent=4)

    good_samples_all_levels = updated_df[['taxon', 'level'] + good_samples]
    good_samples_all_levels.to_csv(all_levels_path, index=False)

    level_dfs = []

    lower_phyla_counts = []
    for level, level_rows in good_samples_all_levels.groupby('level'):
        local_total = {s: level_rows[s].sum() for s in good_samples}
        for s in good_samples:
            level_rows[s] = level_rows[s].apply(
                lambda x: float(x)/local_total[s])
        print(level_rows)
        level_dfs.append(level_rows)

        if level == 'p':
            lower_phyla = level_rows[level_rows['taxon'].isin(lower_phyla)]
            lower_phyla_counts = {s: float(lower_phyla[s].sum())
                                  for s in good_samples}
    good_samples_all_levels = pd.concat(level_dfs)
    columns_order = ['taxon', 'level'] + good_samples
    good_samples_all_levels = good_samples_all_levels[columns_order]
    good_samples_all_levels.to_csv(original_abundances_path.replace('.csv', '.relative_by_levels.csv'), index=False)
    open(original_abundances_path.replace('.csv', '.lower_phyla_counts.txt'), 'w').write(
        '\n'.join([str(x) for x in lower_phyla_counts.values()])
    )

    for_shannon = good_samples_all_levels.set_index('taxon', inplace=False)
    for_shannon.to_csv(original_abundances_path.replace('.csv', '.relative_by_levels.taxon_indexed.csv'), sep=',')

    good_sample_stats = {s: stats for s, stats in pop_stats.items() if s in good_samples}

    find_pop_fb_ratios(good_sample_stats, fb_ratios_output, filtering_log_path)
