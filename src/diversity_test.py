import sys
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Rectangle

from plotting import plot_fb_ratios, plot_shannon, plot_colorbar_with_indicator
    
if __name__ == "__main__":
    #fb_ratio_csv_path = sys.argv[1]
    fb_ratio_csv_path = "~/Documents/feces_Specie_full_abundance_v3.norm1k.fb_ratios.csv"
    #shannon_index_csv_path = sys.argv[2]
    shannon_index_csv_path = "~/Documents/feces_Specie_full_abundance_v3.norm1k.relative_by_levels.shannon.csv"
    shannon_df = pd.read_csv(shannon_index_csv_path)
    shannon_df.rename(columns = {"Unnamed: 0": "Sample"}, inplace=True)
    fb_df = pd.read_csv(fb_ratio_csv_path)
    fb_df.rename(columns = {"Unnamed: 0": "Sample"}, inplace=True)
    fb_df = fb_df[fb_df['firm_ratio'] <= 16]
    fb_vec = fb_df['firm_ratio'].tolist()
    plot_fb_ratios(1.6, fb_vec)
    plot_shannon(3.2, shannon_df['shannon_diversity'].tolist())
    plot_colorbar_with_indicator(2.5, fb_vec, print_patient=True)