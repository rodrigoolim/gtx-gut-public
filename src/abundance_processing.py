import json
from os import path
import sys

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from tqdm import tqdm

from plotting import plot_val_distribution

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

def get_phylum_count(df, patient_id, phylum):
    phylum_level = df[df['level'] == 'p']
    for _, row in phylum_level.iterrows():
        if phylum in row['taxon']:
            return row[patient_id]
    
    return 0.0

def get_all_phylum_counts(df):
    phylum_level = df[df['level'] == 'p']
    counts_by_phylum = {}
    counts_by_phylum_relative = {}
    count_totals = {}
    for phylum, phylum_rows in phylum_level.groupby('taxon'):
        counts = {}
        phylum_rows = phylum_level[phylum_level['taxon'].str.endswith(phylum)]

        assert len(phylum_rows) == 1, 'Duplicate rows for phylum:' + phylum_level['taxon'].tolist() 

        for _, row in phylum_rows.iterrows():
            for key, val in row.items():
                if not key in ['taxon', 'level']:
                    if not key in count_totals:
                        count_totals[key] = 0
                    counts[key] = float(val)
                    count_totals[key] += counts[key]
        counts_by_phylum[phylum] = counts
    
    zero_total = [k for k, v in count_totals.items() if not v > 0]

    for phylum, counts in counts_by_phylum.items():
        for z in zero_total:
            del counts[z]
        
        counts_relative = {k: counts[k] / count_totals[k] for k in counts.keys()}
        counts_by_phylum_relative[phylum] = counts_relative
    
    return counts_by_phylum, counts_by_phylum_relative, count_totals

def remove_outliers(vals: dict, mult = 3, abs_mi=100, remove_low_abs=False, logger=sys.stdout):
    np_vec = np.array(list(vals.values()))
    if remove_low_abs == None:
        remove_low_abs = (len([x for x in np_vec if x > 2]) / len(np_vec)) > 0.1

    if remove_low_abs:
        prev_len = len(np_vec)
        np_vec = np.array([x for x in np_vec if x > abs_mi])
        new_len = len(np_vec)
        if new_len < prev_len:
            print('\tFound samples with very low abs counts', file=logger)
            print('\tReducing number of samples from', file=logger)
            print('\tto', len(np_vec), file=logger)

    median = np.median(np_vec)
    std = np.std(np_vec)

    uplim = median + mult*std
    lolim = median - mult*std

    if remove_low_abs:
        final = {k: v for k, v in vals.items() 
                if v >= lolim and v <= uplim and v > abs_mi}
    else:
        final = {k: v for k, v in vals.items() 
                if v >= lolim and v <= uplim}
    print('\tFinal number of samples:', len(final), file=logger)
    return final

def remove_outliers_list(vec: list, mult = 3, abs_mi=100, remove_low_abs=False, logger=sys.stdout):
    np_vec = np.array(vec)
    if remove_low_abs == None:
        remove_low_abs = (len([x for x in np_vec if x > 2]) / len(np_vec)) > 0.1

    if remove_low_abs:
        prev_len = len(np_vec)
        np_vec = np.array([x for x in np_vec if x > abs_mi])
        new_len = len(np_vec)
        if new_len < prev_len:
            print('\tFound samples with very low abs counts', file=logger)
            print('\tReducing number of samples from', file=logger)
            print('\tto', len(np_vec), file=logger)

    median = np.median(np_vec)
    std = np.std(np_vec)

    uplim = median + mult*std
    lolim = median - mult*std

    if remove_low_abs:
        final = [v for v in np_vec 
                if v >= lolim and v <= uplim and v > abs_mi]
    else:
        final = [v for v in np_vec 
                if v >= lolim and v <= uplim]
    print('\tFinal number of samples:', len(final), file=logger)
    return final

def ratio_nonewise(vals1: dict, vals2: dict):

    ratios = {}
    has_val1 = set(vals1.keys())
    has_val2 = set(vals2.keys())
    '''only_val1 = [k for k in has_val1 if not k in has_val2]
    only_val2 = [k for k in has_val2 if not k in has_val1]
    for k in only_val1:
        vals2[k] = 0.0
    for k in only_val2:
        vals1[k] = 0.0'''
    common_keys = has_val1.intersection(has_val2)
    for k in common_keys:
        v1 = vals1[k]
        v2 = vals2[k]

        if v1 > 0:
            if v2 > 0:
                ratios[k] = (v1 / v2, v2, v1)
        '''    else:
                ratios[k] = (np.inf, v2, v1)
        else:
            if v2 > 0:
                ratios[k] = (np.nextafter(0, 1), v2, v1)
            else:
                ratios[k] = (np.nan, v2, v1)'''
    
    return ratios

def boxplot_fb(pop_firm_counts, pop_bact_counts, sample_totals, output_path, min_total=None,
               min_relative=None, min_rel_firm = None, min_rel_bac = None,
               max_zscore=None):

    pop_firm_bac_ratios_tps = ratio_nonewise(pop_firm_counts, pop_bact_counts)
    pop_firm_bac_ratios = [a for a, _, _ in pop_firm_bac_ratios_tps.values()]
    pop_firm_bac_ratios = remove_outliers_list(pop_firm_bac_ratios, remove_low_abs=False)
    fig, axes = plt.subplots(2, 1, figsize=(12,6))
    ratios_ax = axes[0]
    counts_ax = axes[1]

    ratios_ax.boxplot([pop_firm_bac_ratios],
                       tick_labels=['F/B ratios'],
                       sym='gD', vert=False)
    counts_ax.boxplot(
        [list(sample_totals.values()),
        list(pop_bact_counts.values()),
        list(pop_firm_counts.values())],
        tick_labels=['Total Counts', 'Bact. Counts', 'Firm. Counts'],
        sym='gD', vert=False
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=120)

def calc_pop_stats(pop_abundances_path, sample_totals, filtering_log_path, min_read_count=None):
    """
    min_read_count: minimum total read count for a reference sample to be kept.
    If not provided, defaults to the 10th percentile of the reference population's
    own read-depth distribution (sample_totals) times 0.9, as a safety margin.

    NOTE: in the original (private) pipeline this threshold was derived from the
    minimum read count observed across real patient samples processed by the
    product. That comparison depended on private client data and is not
    reproducible here, so this public version derives the threshold from the
    reference population itself instead.
    """
    log_file = open(filtering_log_path, 'a')
    if min_read_count is None:
        min_read_count = np.percentile(list(sample_totals.values()), 10) * 0.9
    dnagtx_min_counts = min_read_count
    print('dnagtx_min_counts:', dnagtx_min_counts, file=log_file)

    output_path = pop_abundances_path + '.fb_ratios.csv'

    pop_abundances_df = pd.read_csv(pop_abundances_path)
    
    print(pop_abundances_df.head(), file=log_file)

    counts, _, local_totals = get_all_phylum_counts(pop_abundances_df)
    pop_firm_counts = counts['k__Bacteria;p__Firmicutes']
    pop_bact_counts = counts['k__Bacteria;p__Bacteroidetes']
    
    output_path = pop_abundances_path + '.fb_raw.png'
    boxplot_fb(pop_firm_counts, pop_bact_counts, local_totals, output_path)

    print('Removing samples with total less than', dnagtx_min_counts, file=log_file)
    high_count_totals = {k: v for k, v in local_totals.items() if sample_totals[k] >= dnagtx_min_counts}
    print('\tFrom', len(local_totals), 'to', len(high_count_totals), file=log_file)
    print('Removing low count samples from firm values', file=log_file)
    pop_firm_counts2 = {k: v for k, v in pop_firm_counts.items() if k in high_count_totals}
    print('\tFrom', len(pop_firm_counts), 'to', len(pop_firm_counts2), file=log_file)
    print('Removing low count samples from bact values', file=log_file)
    pop_bact_counts2 = {k: v for k, v in pop_bact_counts.items() if k in high_count_totals}
    print('\tFrom', len(pop_bact_counts), 'to', len(pop_bact_counts2), file=log_file)

    output_path2 = pop_abundances_path + '.fb_high_count.png'
    boxplot_fb(pop_firm_counts2, pop_bact_counts2, high_count_totals, output_path2)

    print('Removing outliers from firm values', file=log_file)
    pop_firm_counts_nooutlier = remove_outliers(pop_firm_counts2, mult=2.5, remove_low_abs = False, logger=log_file)
    print('Removing outliers from bact values', file=log_file)
    pop_bact_counts_nooutlier = remove_outliers(pop_bact_counts2, mult=2.5, remove_low_abs = False, logger=log_file)

    print('Samples without outliers:', file=log_file)
    normal_firm_and_bact = {k for k in pop_firm_counts_nooutlier.keys() if k in pop_bact_counts_nooutlier}
    print('\t', len(normal_firm_and_bact), file=log_file)
    pop_firm_counts_nooutlier = {k: v for k, v in pop_firm_counts_nooutlier.items() if k in normal_firm_and_bact}
    pop_bact_counts_nooutlier = {k: v for k, v in pop_bact_counts_nooutlier.items() if k in normal_firm_and_bact}
    total_counts_nooutlier = {k: v for k, v in high_count_totals.items() if k in normal_firm_and_bact}
    
    output_path3 = pop_abundances_path + '.fb_no_outliers.png'
    boxplot_fb(pop_firm_counts_nooutlier, pop_bact_counts_nooutlier, total_counts_nooutlier, output_path3)

    pop_stats = {}

    print('Collecting population stats', file=log_file)
    for sample_name in local_totals.keys():
        #if sample_name in pop_firm_counts and sample_name in pop_bact_counts:
        stats = {
            'total_count': sample_totals[sample_name],
            'total_count_norm1k': local_totals[sample_name],
            'low_read_count': not sample_name in high_count_totals,
            'firmicutes_count': pop_firm_counts[sample_name],
            'bacteroidetes_count': pop_bact_counts[sample_name],
            'firmicute_outlier': not sample_name in pop_firm_counts_nooutlier,
            'bacteroidetes_outlier': not sample_name in pop_bact_counts_nooutlier
        }
        pop_stats[sample_name] = stats
    
    return pop_stats

def find_pop_fb_ratios(pop_stats, output_path, filtering_log_path):
    print('Writing population F/B ratios file')
    log_file = open(filtering_log_path, 'a')
    pop_firm_counts = {}
    pop_bact_counts = {}
    for sample_name, stats in pop_stats.items():
        pop_firm_counts[sample_name] = stats['firmicutes_count']
        pop_bact_counts[sample_name] = stats['bacteroidetes_count']
    
    pop_firm_bac_ratios = ratio_nonewise(pop_firm_counts, pop_bact_counts)
    pop_firm_bac_ratios_only = {sample: v[0] for sample, v in pop_firm_bac_ratios.items()}
    pop_firm_bac_ratios_only = remove_outliers(pop_firm_bac_ratios_only, mult=3, remove_low_abs = False, logger=log_file)
    pop_firm_bac_ratios_only = {sample: v for sample, v in pop_firm_bac_ratios_only.items()
        if v == v and v < np.inf}

    pop_firm_bac_ratios = {sample: (pop_firm_bac_ratios_only[sample], v[1], v[2])
                           for sample, v in pop_firm_bac_ratios.items() 
                           if sample in pop_firm_bac_ratios_only}
        
    output = open(output_path, 'w')
    output.write('"","firm_ratio","bact","firm"\n')
    
    ratios = [v[0] for sample, v in pop_firm_bac_ratios.items()]
    for key, val in pop_firm_bac_ratios.items():
        ratio, bact, firm = val
        output.write('"'+key+'",'+str(ratio)
                +','+str(bact)+','+str(firm)+'\n')
           
    plot_val_distribution(ratios, output_path + '.png', 
        outlier_filter=20, x_label = "F/B Ratios", plot_title = "Population F/B Ratios - No Outliers")
    #ratios = [r for r in ratios if r == r and r < np.inf]

def read_pop_fb_ratio(ratio_csv_path):

    pop_ratios = []
    for rawline in open(ratio_csv_path, 'r'):
        cells = [a.strip('"') for a in rawline.rstrip('\n').split(',')]
        if len(cells[0]) > 0:
            if len(cells) == 4:
                sample_id, ratio, bact, firm = cells
                if ratio == 'inf':
                    pass
                else:
                    ratio = float(ratio)
                    pop_ratios.append(ratio)
            else:
                print('Wrong number of cells:')
                print(len(cells), cells)
        
    return pop_ratios

def generic_counts_processing(original_abundances_path):
    filtering_log_path = original_abundances_path.replace('.csv', '.filtering_log.txt')

    pop_level = 's'

    original_abundances = pd.read_csv(original_abundances_path, sep=',')
    original_abundances['taxon'] = [t.replace('[', '').replace(']', '') for t in original_abundances['Unnamed: 0'].tolist()]
    del original_abundances['Unnamed: 0']
    filtering_log = open(filtering_log_path, 'w')
    print(original_abundances.columns.tolist(), file = filtering_log)
    print(original_abundances.head(), file = filtering_log)

    all_species = original_abundances['taxon'].tolist()
    sample_names = [s for s in original_abundances.columns.tolist() if not s in ['taxon']]
    if pop_level == 's':
        taxa_level_prefixes = ['p__', 'c__', 'o__', 'f__', 'g__']
        new_abundance_lines = []
        for prefix in taxa_level_prefixes:
            taxa_tuples = set([(s.split(';'+prefix)[0], prefix + s.split(';'+prefix)[-1].split(';')[0] )
                            for s in all_species])
            taxa_set = set([parent_taxon + ';' + taxon 
                            for parent_taxon, taxon in taxa_tuples if len(taxon) > 3])
            print(prefix, file = filtering_log)
            print(str(len(taxa_set)), file = filtering_log)

            level = prefix.rstrip('_')
            for t in tqdm(taxa_set):
                string_to_search = t
                species_contained_in_taxon = {s for s in all_species if string_to_search in s}
                taxon_df = original_abundances[original_abundances['taxon'].isin(species_contained_in_taxon)]
                perc_totals = {s: taxon_df[s].sum() for s in sample_names}
                perc_totals['level'] = level
                perc_totals['taxon'] = string_to_search
            
                new_abundance_lines.append(perc_totals)

        original_abundances['level'] = 's'
        original_abundances = original_abundances[~original_abundances['taxon'].str.endswith('s__')]
        new_df = pd.DataFrame(new_abundance_lines)
        updated_df = pd.concat([original_abundances, new_df])
    elif pop_level == 'p':
        original_abundances['level'] = 'p'
        updated_df = original_abundances
    
    columns_order = ['taxon', 'level'] + sample_names
    updated_df = updated_df[columns_order]

    all_levels_path = original_abundances_path.replace('.csv', '.all_levels.csv')

    filtering_log.close()

    return updated_df, all_levels_path, sample_names, columns_order, filtering_log_path