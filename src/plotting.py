from os import path
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import csv
from scipy.stats import gaussian_kde
import math

figsize = (8,3.6)

def read_taxa_percentages(csv_path):
    #"501",50.6067890015829
    patient_ratio = None

    pop_ratios = []
    for rawline in open(csv_path, 'r'):
        cells = [a.strip('"') for a in rawline.rstrip('\n').split('",')]
        if len(cells[0]) > 0 and len(cells) == 2:
            _, ratio = cells
            ratio = float(ratio)
            pop_ratios.append(ratio)
    patient_ratio = pop_ratios[-1]
    pop_ratios = pop_ratios[:-1]
    return pop_ratios, patient_ratio

def read_patient_firm_ratio(ratio_csv_path):
    # expected row format: "patient.<sample_id>",<ratio>
    patient_id = None
    patient_ratio = None
    patient_bact = None
    patient_firm = None

    pop_ratios = []
    for rawline in open(ratio_csv_path, 'r'):
        cells = [a.strip('"') for a in rawline.rstrip('\n').split(',')]
        if len(cells[0]) > 0:
            if len(cells) == 4:
                sample_id, ratio, bact, firm = cells
                if ratio == 'inf':
                    ratio = None
                else:
                    ratio = float(ratio)
                if sample_id.startswith('patient.'):
                    patient_id = sample_id.replace('patient.', '')
                    patient_ratio = ratio
                    patient_bact = bact
                    patient_firm = firm
                else:
                    pop_ratios.append(ratio)
            else:
                print('Wrong number of cells:')
                print(len(cells), cells)
        
    return patient_id, pop_ratios, patient_ratio, patient_bact, patient_firm

def read_patient_shannon_diversity(ratio_csv_path):
    """
    Reads the patient's own Shannon diversity value, from the single-sample
    CSV produced by the shannon_vegan_r rule (shannon.R run on the patient's
    own raw_derep_table.qza). Returns a single float, or None if unavailable.
    """
    patient_ratio = None

    with open(ratio_csv_path, mode='r', newline='') as csvfile:
        reader = csv.reader(csvfile)

        try:
            next(reader)
        except StopIteration:
            return patient_ratio

        for row in reader:
            if len(row) == 2:
                sample_id, ratio_str = row

                try:
                    patient_ratio = float(ratio_str)
                except ValueError:
                    print(f"Aviso: Não foi possível converter '{ratio_str}' para float na linha: {row}")

    return patient_ratio

def read_pop_shannon_diversity(shannon_csv_path):
    """
    Reads the Shannon diversity of every sample in the reference population,
    from the pop_shannon_vegan_r output (data/*.relative_by_levels.shannon.csv).
    Skips the "level" placeholder row written alongside the real per-sample values.
    """
    pop_ratios = []

    with open(shannon_csv_path, mode='r', newline='') as csvfile:
        reader = csv.reader(csvfile)

        try:
            next(reader)
        except StopIteration:
            return pop_ratios

        for row in reader:
            if len(row) == 2:
                sample_id, ratio_str = row
                if sample_id == 'level':
                    continue

                try:
                    pop_ratios.append(float(ratio_str))
                except ValueError:
                    print(f"Aviso: Não foi possível converter '{ratio_str}' para float na linha: {row}")

    return pop_ratios

def plot_val_distribution(vals, 
                          plot_path, 
                          outlier_filter=None,
                          x_label = "",
                          plot_title = ""):
    vals.sort()
    if outlier_filter:
        vals = [r if r < outlier_filter else outlier_filter for r in vals]

    fig, ax = plt.subplots(1, 1, figsize=figsize) # Make it 14x7 inch
    #plt.style.use('Solarize_Light2') # nice and clean grid
    n, bins, patches = ax.hist(vals, bins=60, 
        facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
    #print(n)
    #print(patient_ratio)
    #print(bins)
    bin_width = (bins[1]-bins[0])
    #print('bin width', bin_width)
    n = n.astype('int') # it MUST be integer
    # Good old loop. Choose colormap of your taste
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.hsv((n[i]/max(n))*0.5))

    ax.set_title(plot_title)
    ax.set_xlabel(x_label)
    ax.get_yaxis().set_visible(False)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    #s_cmap = plt.cm.ScalarMappable(cmap = plt.get_cmap('viridis'), norm = Normalize(min(n), max(n)))
    #fig.colorbar(s_cmap, ax=ax, fraction=0.04, label='Frequency')
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)

    return plot_path

def plot_firmi_ratios(pop_ratios, patient_ratio, 
                      patient_str, output_dir, outlier_filter=1.3):
    pop_ratios.sort()
    #mean = np.median(pop_ratios)
    q1 = np.quantile(pop_ratios, 0.25)
    q3 = np.quantile(pop_ratios, 0.75)
    iqr = q3-q1
    fence_high = q3 + (outlier_filter*iqr)
    #std = np.std(pop_ratios)
    #print(min(pop_ratios), q1, mean, q3, fence_high, max(pop_ratios))
    #print(iqr, std)
    #max_ratio = mean*4
    pop_ratios = [r for r in pop_ratios if r < fence_high]

    fig, ax = plt.subplots(1, 1, figsize=figsize) # Make it 14x7 inch
    #plt.style.use('Solarize_Light2') # nice and clean grid
    n, bins, patches = ax.hist(pop_ratios, bins=50, 
        facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
    #print(n)
    #print(patient_ratio)
    #print(bins)
    bin_width = (bins[1]-bins[0])
    #print('bin width', bin_width)
    patient_bin = len(bins)-1
    for i in range(len(bins)-1):
        if patient_ratio >= bins[i] and patient_ratio < bins[i+1]:
            patient_bin = i
            break
    print(patient_bin)
    patient_bin_x = (bin_width)*patient_bin + bins[0]
    patient_bin_x_center = patient_bin_x + bin_width*0.5
    #patient_bin_x2 = patient_bin_x + bin_width
    #print(patient_bin, bin_width, patient_bin_x)
    patient_bin_y = n[patient_bin] + 1
    #patient_bin_y2 = 0
    n = n.astype('int') # it MUST be integer
    # Good old loop. Choose colormap of your taste
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.viridis(n[i]/max(n)))
    #patches[patient_bin].set_fc('red') # Set color
    #patches[patient_bin].set_alpha(1) # Set opacity
    
    #patches[patient_bin].set_linestyle('--')
    #patches[patient_bin].set_linewidth(2.5)
    #patches[patient_bin].set_edgecolor('black')

    rect = Rectangle((patient_bin_x, 0.1), bin_width, patient_bin_y-0.9)
    rect.set_linestyle('--')
    rect.set_linewidth(3)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)

    bbox = dict(boxstyle="round", fc="0.8")
    arrowprops = dict(
        arrowstyle="->",
        connectionstyle="angle,angleA=0,angleB=90,rad=10",
        linewidth=3)

    ax.annotate(patient_str, (patient_bin_x_center, patient_bin_y), 
                 xytext=((bins[-1]-bins[0])/2 + bins[0], max(n)*(0.95)),
                 bbox=bbox, arrowprops=arrowprops, ha='center',
                 fontsize=12)

    ax.set_title('Patient Firm Ratio Compared to Population') 
    ax.set_xlabel('Firm Ratios')
    ax.get_yaxis().set_visible(False)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    plot_path = output_dir + '/firm_ratio.png'
    s_cmap = plt.cm.ScalarMappable(cmap = plt.get_cmap('viridis'), norm = Normalize(min(n), max(n)))
    fig.colorbar(s_cmap, ax=ax, fraction=0.04, label='Frequency')
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)


    return plot_path

def plot_shannon_diversity(pop_ratios, patient_ratio, 
                      patient_str, output_dir, outlier_filter=1.5):
    pop_ratios.sort()
    q1 = np.quantile(pop_ratios, 0.25)
    q3 = np.quantile(pop_ratios, 0.75)
    iqr = q3-q1
    fence_high = q3 + (outlier_filter*iqr)
    pop_ratios = [r for r in pop_ratios if r < fence_high]

    fig, ax = plt.subplots(1, 1, figsize=figsize) # Make it 14x7 inch
    #plt.style.use('Solarize_Light2') # nice and clean grid
    n, bins, patches = ax.hist(pop_ratios, bins=50, 
        facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
    
    bin_width = (bins[1]-bins[0])
    patient_bin = len(bins)-1
    for i in range(len(bins)-1):
        if patient_ratio >= bins[i] and patient_ratio < bins[i+1]:
            patient_bin = i
            break
    patient_bin_x = (bin_width)*patient_bin + bins[0]
    patient_bin_x_center = patient_bin_x + bin_width*0.5
    patient_bin_y = n[patient_bin] + 1
    n = n.astype('int') # it MUST be integer
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.viridis(n[i]/max(n)))
    
    rect = Rectangle((patient_bin_x, 0.1), bin_width, patient_bin_y-0.9)
    rect.set_linestyle('--')
    rect.set_linewidth(3)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)

    bbox = dict(boxstyle="round", fc="0.8")
    arrowprops = dict(
        arrowstyle="->",
        connectionstyle="angle,angleA=0,angleB=90,rad=10",
        linewidth=3)

    #print(bins)
    #print(bins[-1], bins[0])
    #print((bins[-1]-bins[0])/2 + bins[0])
    ax.annotate(patient_str, (patient_bin_x_center, patient_bin_y), 
                 xytext=((bins[-1]-bins[0])/2  + bins[0], max(n)*(0.95)),
                 bbox=bbox, arrowprops=arrowprops, ha='center',
                 fontsize=12)

    ax.set_title('Patient Shannon Diversity Compared to Population') 
    ax.set_xlabel('Shannon Diversities')
    ax.get_yaxis().set_visible(False)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    
    plot_path = output_dir + '/shannon_diversity.png'
    s_cmap = plt.cm.ScalarMappable(cmap = plt.get_cmap('viridis'), 
                                   norm = Normalize(min(n), max(n)))
    fig.colorbar(s_cmap, ax=ax, fraction=0.04, label='Frequency')
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)


    return plot_path

def plot_single_bact(pop_percs, patient_perc, bact_pop_graph, stats):
    low_plot_lim = 0.0
    high_plot_lim = 100.0
    pop_percs.append(patient_perc)
    pop_percs.sort()

    if max(pop_percs) < 100:
        high_plot_lim = max(pop_percs)
    
    '''if high_plot_lim > stats['min_high_outlier']:
        high_plot_lim = stats['min_high_outlier']

    if patient_perc > high_plot_lim:
        high_plot_lim = patient_perc + stats['mad']*0.01'''
    
    '''q1 = np.quantile(pop_percs, 0.25)
    q3 = np.quantile(pop_percs, 0.75)
    iqr = q3-q1
    fence_high = min(100, q3 + (1.7*iqr))
    fence_low = max(0, q1 - (1.7*iqr))
    if fence_high < patient_perc:
        fence_high = patient_perc + iqr*0.5
    if fence_low > patient_perc:
        fence_low = patient_perc - iqr*0.5
    if pop_percs[0] == 0.0:
        fence_low = 0.0'''
    total_width = high_plot_lim - low_plot_lim

    #if patient_perc >= 0.0049832913173 and patient_perc <= 0.0049832913174:
    #    print(fence_low, q1, np.median(pop_percs), q3, fence_high)
    #    print(pop_percs)

    plot_figsize = (5, 2.1)
    fig, ax = plt.subplots(1, 1, figsize=plot_figsize)
    n, bins, patches = ax.hist(pop_percs, bins=20, range=(low_plot_lim, high_plot_lim),
        facecolor = '#2ab0ff', edgecolor='#169acf', linewidth=0.5)
    #print('histogram made')

    n = n.astype('int') # it MUST be integer
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.viridis(n[i]/max(n)))
        patches[i].set_linewidth(0)
    triangle_top = max(n)*1.15
    ax.vlines([patient_perc], -1, max(n)*0.93, linewidth=4, linestyles=['--'])
    indicator_triangle = Polygon([(patient_perc,  max(n)*0.95), 
                                  (patient_perc + (total_width*0.05), triangle_top),
                                  (patient_perc - (total_width*0.05), triangle_top)])
    ax.add_artist(indicator_triangle)
    ax.get_yaxis().set_visible(False)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    if (high_plot_lim - low_plot_lim) <= 0.01:
        #print(fence_high - fence_low)
        ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    fig.tight_layout()
    #print('saved')
    fig.savefig(bact_pop_graph, dpi=120)

    return plot_figsize[1] / plot_figsize[0], high_plot_lim

def custom_boxplot(ax, stats, box_linewidth, median_linewidth, box_height = 1.0):
    rect = Rectangle((stats['q1'], 0.5), stats['iqr'], 1.0, zorder=10)
    rect.set_linestyle('-')
    rect.set_linewidth(box_linewidth)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)

    ax.hlines([1, 1], [stats['l_outlier'], stats['q3']], [stats['q1'], stats['r_outlier']],
        linewidth=box_linewidth/2, linestyles=['-', '-'], colors='black', zorder=8)
    ax.vlines([stats['l_outlier'], stats['r_outlier']], 
        0.75, 1.25, linewidth=box_linewidth/2, linestyles=['-', '-'], colors='black', zorder=7)
    ax.vlines([stats['median']], 0.5, 1.5,
        linewidth=median_linewidth, linestyles=['-'], colors='red', zorder=11)

def boxplot_single_bact(pop_percs, patient_perc, bact_pop_graph, stats):
    all_percs = [patient_perc] + pop_percs
    first_line = stats['low_ths'][0]
    last_line = stats['high_ths'][1]
    if patient_perc < first_line:
        first_line = patient_perc
    if first_line < 0:
        first_line = 0.0    
    if patient_perc > last_line:
        last_line = patient_perc
    visual_content_range = last_line - first_line

    bottom_value = first_line - visual_content_range*0.1
    top_value = last_line + visual_content_range*0.1
    '''iqr_higher = stats['q3'] + stats['iqr']*6
    if iqr_higher < top_value:
        top_value = iqr_higher
    if patient_perc > top_value:
        top_value = patient_perc + (patient_perc - bottom_value)*0.05
    bottom_value -= top_value*0.025'''
    value_range = top_value - bottom_value
    plot_figsize = (6, 1.2)
    fig, ax = plt.subplots(1, 1, figsize=plot_figsize)
    y_center = 1.0
    box_height = 1.0
    box_top_y = y_center + box_height*0.5
    box_bottom_y = y_center - box_height*0.5

    scatter_x_low = [bottom_value for x in pop_percs if x <= bottom_value]
    scatter_x_high = [top_value for x in pop_percs if x >= top_value]
    scatter_x = [x for x in pop_percs if x > bottom_value and x < top_value]

    for x_list, marker in [(scatter_x, '.'), (scatter_x_low, '<'), (scatter_x_high, '>')]:
        scatter_y = [y if (y >= 0.5 and y <= 1.5) else 1.0
            for y in np.random.normal(1.0, 0.15, size=len(x_list))]
        ax.plot(x_list, scatter_y, marker, alpha=0.2, color='#21a960', zorder=1)
    
    custom_boxplot(ax, stats, 3.0, 4.0, box_height=box_height)
    #ax.boxplot([pop_percs], orientation='horizontal', notch=False, 
    #    sym='', widths=[box_height],
    #    boxprops={'linewidth': 2.0}, medianprops={'linewidth': 3.0})
    
    triangle_top = box_top_y + 0.06
    triangle_bottom = box_top_y - 0.15
    triangle_width = value_range*0.05
    indicator_color = 'orange'
    ax.vlines([patient_perc], box_bottom_y-0.1, triangle_bottom-0.05, linewidth=3, linestyles=[':'], 
              colors=indicator_color, zorder=9)
    indicator_triangle = Polygon([(patient_perc,  triangle_bottom), 
                                  (patient_perc + (triangle_width/2), triangle_top),
                                  (patient_perc - (triangle_width/2), triangle_top)],
                                  color=indicator_color, zorder=9)
    ax.add_artist(indicator_triangle)
    ax.get_yaxis().set_visible(False)
    final_x_start = bottom_value - value_range*0.005
    final_x_end = top_value + value_range*0.005
    ax.set_xlim(final_x_start, final_x_end)
    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    '''if value_range <= 0.01:
        #print(fence_high - fence_low)
        ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))'''
    fig.tight_layout()
    #print('saved')
    fig.savefig(bact_pop_graph, dpi=120)
    fig.clear()
    plt.close()

    plot_details = {
        'first_line': first_line,
        'last_line': last_line,
        'top_value': top_value,
        'bottom_value': bottom_value,
        'final_x_start': final_x_start,
        'final_x_end': final_x_end
    }

    return plot_figsize[1] / plot_figsize[0], plot_details


cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["#d73027","limegreen","palegreen"])

def plot_short_colorbar(to_indicate, vec, points = 1000, 
        output_path=None, print_start_end=True,
        start_from_zero = True):
    vec2 = sorted([to_indicate] + vec)
    q1 = np.percentile(vec2, 25)
    q3 = np.percentile(vec2, 75)
    iqr = q3-q1
    low_limit = q1 - 1.4*iqr
    high_limit = q3 + 1.4*iqr
    
    if start_from_zero:
        mi = 0.0
    else:
        mi = max(0.0, low_limit)
        if to_indicate < mi:
            mi = to_indicate
    
    ma = high_limit
    if max(vec2) < high_limit:
        ma = max(vec2)
    if ma < to_indicate:
        ma = to_indicate + ((to_indicate - mi)*0.08)
    if ma == 0:
        ma = 0.0025
    range_size = ma - mi
    
    indicator_at_beggining = round(to_indicate, 4) == round(mi, 4)

    values = np.linspace(mi, ma, num=points)
    try:
        kde = gaussian_kde(vec2)(values)
    except Exception as err:
        print(err)
        print('Error on KDE gaussian')
        print(min(vec2))
        print(q1)
        print(np.median(vec2))
        print(q3)
        print(max(vec2))
        quit(1)

    kde = kde / max(kde)
    gradient = np.vstack((kde, kde))
    fig_h = 0.5
    fig_w = 5
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.imshow(gradient, cmap=cmap, aspect='auto', )
    
    try:
        patient_offset = to_indicate - mi
        pos_float = patient_offset / range_size
        pos_real = pos_float*points
        indicator_width = points / 20
    except Exception as e:
        print("ERRROOOOOOOOOOOOOOOR::: " + e)

    #if not indicator_at_beggining:
    if indicator_at_beggining:
        pos_real = points*0.03
    ax.vlines([pos_real], 1.85, -0.85, colors='#00000001')
    rect = Rectangle((pos_real-(indicator_width*0.5), -0.9), indicator_width, 2.65)
    rect.set_linestyle('-')
    rect.set_linewidth(2.5)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)
    #else:
    #    ax.vlines([pos_real], 1.85, -0.85, colors='#000000ff')
    
    #ax.set_xticks([min_fb, len(fb_values)], [str(round(min_fb, 1)), str(round(max_fb, 1))])
    ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
    ax.get_yaxis().set_visible(False)
    #plt.yticks(fontsize=14)
    if print_start_end:
        #print(range_size)
        #print(mi)
        if mi == 0.0:
            min_str = "0"
        else:
            min_str = str(round(mi, 3)).rstrip('0').rstrip(".")
        
        ax.text(-(points/100), 0.5, min_str, 
                va='center', ha='right', fontsize=12)
        max_str = str(round(ma*100, 3))
        if not '.' in max_str:
            max_str += '.0'
        while len(max_str) < 5:
            max_str = max_str + '0'
        ax.text(points+(points/100), 0.5, max_str, 
                va='center', ha='left', fontsize=12)
        
    ax.get_xaxis().set_visible(False)
    #fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
    else:
        fig.savefig('./output.png', dpi=200, bbox_inches='tight')
    
    fig_details = {
        'min': min(vec2),
        'mi': mi,
        'q1': q1,
        'median': np.median(vec2),
        'q3': q3,
        'ma': ma,
        'max': max(vec2),
        'to_indicate': to_indicate,
        'range_size': range_size,
        'patient_offset': patient_offset,
        'pos_float': pos_float,
        'pos_real': pos_real,
        'indicator_width': indicator_width
    }

    return fig_h / fig_w, fig_details

def plot_colorbar_with_indicator(to_indicate, vec, points = 1000, 
                                 output_path=None, print_start_end=True, 
                                 print_patient=True, start_from_zero = True):
    
    debug_param_filename = f"{output_path}_debug.txt"

    with open(debug_param_filename, 'w') as f:
        f.write(f"{to_indicate}\n")
        f.write(f"{vec}\n")
    
    to_indicate = to_indicate * 100
    #vec_np = np.array([v * 100 for v in vec]) # Aplicando a mesma escala de to_indicate

    vec2 = sorted([to_indicate] + [v*100 for v in vec])
    #vec2 = sorted([to_indicate] + list(vec)) # Use 'list(vec)' caso 'vec' seja numpy array

    q1 = np.percentile(vec2, 25)
    q3 = np.percentile(vec2, 75)
    iqr = q3-q1
    low_limit = q1 - 1.4*iqr
    high_limit = q3 + 1.4*iqr
    
    if start_from_zero:
        mi = 0.0
    else:
        mi = max(0.0, low_limit)
        if to_indicate < mi:
            mi = to_indicate
    
    ma = high_limit
    if max(vec2) < high_limit:
        ma = max(vec2)
    if ma < to_indicate:
        ma = to_indicate + ((to_indicate - mi)*0.08)
    range_size = ma - mi
    
    indicator_at_beggining = round(to_indicate, 4) == round(mi, 4)

    values = np.linspace(mi, ma, num=points)
    try:
        kde = gaussian_kde(vec2)(values)
    except Exception as err:
        print(err)
        print('Error on KDE gaussian')
        print(min(vec2))
        print(q1)
        print(np.median(vec2))
        print(q3)
        print(max(vec2))
        quit(1)

    kde = kde / max(kde)
    gradient = np.vstack((kde, kde))
    fig_h = 0.6
    fig_w = 5
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.imshow(gradient, cmap=cmap, aspect='auto', )
    
    patient_offset = to_indicate - mi
    pos_float = patient_offset / range_size
    pos_real = pos_float*points
    indicator_width = points / 20
    
    if not indicator_at_beggining:
        ax.vlines([pos_real], 1.85, -0.85, colors='#00000001')
        rect = Rectangle((pos_real-(indicator_width*0.5), -0.9), indicator_width, 2.65)
        rect.set_linestyle('-')
        rect.set_linewidth(2.5)
        rect.set_edgecolor('black')
        rect.set_facecolor('#00000000')
        ax.add_artist(rect)
    else:
        ax.vlines([pos_real], 1.85, -0.85, colors='#000000ff')
    
    #ax.set_xticks([min_fb, len(fb_values)], [str(round(min_fb, 1)), str(round(max_fb, 1))])
    ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
    ax.get_yaxis().set_visible(False)
    #plt.yticks(fontsize=14)
    if print_start_end:
        #print(range_size)
        #print(mi)
        if mi == 0.0:
            min_str = "0"
        else:
            min_str = str(round(mi, 3)).rstrip('0').rstrip(".")
        
        ax.text(-(points/100), 0.5, min_str, 
                va='center', ha='right', fontsize=14)
        max_str = str( round(ma, 3) ).rstrip('0').rstrip(".")
        ax.text(points+(points/100), 0.5, max_str, 
                va='center', ha='left', fontsize=14)
        #ax.set_yticks([0.5], [str(mi)])
    if print_patient:
        p_str = str(round(to_indicate, 3)).rstrip('0').rstrip(".")
        ax.text(pos_real, -0.9, p_str, 
                va='bottom', ha='center', fontsize=14)
    ax.get_xaxis().set_visible(False)
    #fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
    else:
        fig.savefig('./output.png', dpi=200, bbox_inches='tight')
    
    fig_details = {
        'min': min(vec2),
        'mi': mi,
        'q1': q1,
        'median': np.median(vec2),
        'q3': q3,
        'ma': ma,
        'max': max(vec2),
        'to_indicate': to_indicate,
        'range_size': range_size,
        'patient_offset': patient_offset,
        'pos_float': pos_float,
        'pos_real': pos_real,
        'indicator_width': indicator_width
    }

    return fig_h / fig_w, fig_details

def plot_fb_ratios(patient_fb, fb_vec, points = 1000, output_path=None, print_values=True):
    max_fb = max(fb_vec)
    min_fb = 0.0
    fb_range_size = max_fb - min_fb
    fb_vec2 = sorted([patient_fb] + fb_vec)
    fb_values = np.linspace(min_fb, max_fb, num=points)
    fb_kde = gaussian_kde(fb_vec2)(fb_values)
    
    fb_kde = fb_kde / max(fb_kde)
    max_index = np.argmax(fb_kde, axis=0)
    kde1 = 1.0-fb_kde[:max_index]
    kde1 += 1.0
    kde2 = fb_kde[max_index:]
    kde_alt = np.concatenate((kde1, kde2), axis=0)
    kde_alt = kde_alt / max(kde_alt)
    
    gradient = np.vstack((kde_alt, kde_alt))
    fig_h = 1
    fig_w = 5
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.imshow(gradient, cmap=cmap, aspect='auto', )
    
    patient_offset = patient_fb - min_fb
    pos_float = patient_offset / fb_range_size
    pos_real = pos_float*points
    indicator_width = points / 20
    
    ax.vlines([pos_real], 1.85, -0.85, colors='#00000001')
    rect = Rectangle((pos_real-(indicator_width*0.5), -0.9), indicator_width, 2.65)
    rect.set_linestyle('-')
    rect.set_linewidth(2.5)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)
    
    #ax.set_xticks([min_fb, len(fb_values)], [str(round(min_fb, 1)), str(round(max_fb, 1))])
    ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    if print_values:
        ax.text(-(points/100), 0.5, str(round(min_fb, 2)).rstrip('0').rstrip("."), 
                va='center', ha='right', fontsize=14)
        ax.text(points+(points/100), 0.5, str(round(max_fb, 1)).rstrip('0').rstrip("."), 
                va='center', ha='left', fontsize=14)
        ax.text(pos_real, -0.9, str(round(patient_fb, 3)).rstrip('0').rstrip("."), 
                va='bottom', ha='center', fontsize=14)
    
    ax.get_xaxis().set_visible(False)
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200,bbox_inches='tight')
    else:
        fig.savefig('./output.png', dpi=200, bbox_inches='tight')
        
    return fig_h / fig_w
        
def plot_shannon(patient_sh, shannon_vec, points = 1000, output_path=None,
                 print_values=True):
    

    if not patient_sh or not shannon_vec:
        print( "ERRRRRROR BBBBBB", patient_sh, shannon_vec )
        return 1/5

    shannon_vec2 = sorted([patient_sh] + [ 0, 3, 6 ])
    max_shannon = max(shannon_vec2)
    min_shannon = min(shannon_vec2)
    shannon_range_size = max_shannon - min_shannon
    shannon_values = np.linspace(min_shannon, max_shannon, num=points)

    print( "DEBUGGGGGIN CCCCC:::: ", shannon_vec2, shannon_values )
    shannon_kde = gaussian_kde(shannon_vec2)(shannon_values)

    shannon_kde = shannon_kde / max(shannon_kde)
    max_index = np.argmax(shannon_kde, axis=0)
    shannon_kde1 = shannon_kde[:max_index]
    shannon_kde2 = 1.0-shannon_kde[max_index:]
    shannon_kde2 += 1.0
    shannon_kde_alt = np.concatenate((shannon_kde1, shannon_kde2), axis=0)
    shannon_kde_alt = shannon_kde_alt / max(shannon_kde_alt)

    gradient = np.vstack((shannon_kde_alt, shannon_kde_alt))
    fig_h = 1
    fig_w = 5
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.imshow(gradient, cmap=cmap, aspect='auto', )
    patient_offset = patient_sh - min_shannon
    pos_float = patient_offset / shannon_range_size
    pos_real = pos_float*points
    indicator_width = points / 20
    ax.vlines([pos_real], 1.85, -0.85, colors='#00000001')
    rect = Rectangle((pos_real-(indicator_width*0.5), -0.9), indicator_width, 2.65)
    rect.set_linestyle('-')
    rect.set_linewidth(2.5)
    rect.set_edgecolor('black')
    rect.set_facecolor('#00000000')
    ax.add_artist(rect)
    #ax.set_xticks([min_shannon, len(shannon_values)], [str(round(min_shannon, 1)), str(round(max_shannon, 1))])
    ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    if print_values:
        ax.text(-(points/100), 0.5, str(round(min_shannon, 1)).rstrip('0').rstrip("."), 
                va='center', ha='right', fontsize=14)
        ax.text(points+(points/100), 0.5, str(round(max_shannon, 1)).rstrip('0').rstrip("."), 
                va='center', ha='left', fontsize=14)
        ax.text(pos_real, -0.9, str(round(patient_sh, 3)).rstrip('0').rstrip("."), 
                va='bottom', ha='center', fontsize=14)
    ax.get_xaxis().set_visible(False)
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches='tight')
    else:
        fig.savefig('./output.png', dpi=200, bbox_inches='tight')
        
    return fig_h / fig_w

if __name__ == "__main__":
    ratio_csv_path = sys.argv[1] + '/firm_ratio2.csv'
    shannon_csv_path = sys.argv[1] + '/shannon_diversity2.csv'
    patient_id, pop_ratios, patient_ratios = read_patient_firm_ratio(ratio_csv_path)
    #firm_ratio_plot_path = plot_firmi_ratios(pop_ratios, patient_ratios, patient_id, path.dirname(ratio_csv_path))
    #pop_div, patient_div = read_patient_shannon_diversity(shannon_csv_path, patient_id)
    #firm_ratio_plot_path = plot_shannon_diversity(pop_div, patient_div, patient_id, path.dirname(ratio_csv_path))

    low_freq_bact = 'alactolyticus'
    patient_and_pop_freqs = pd.read_csv(sys.argv[1] + '/pop_and_patient.taxon_indexed.csv')
    patient_and_pop_freqs = patient_and_pop_freqs[patient_and_pop_freqs['taxon'].str.contains(low_freq_bact)]

    first_row = [r for _, r in patient_and_pop_freqs.iterrows()][0]
    patient_freq = first_row[patient_id]
    pop_freqs = [f for k, f in first_row.items() if type(f) != str and k != patient_id]
    output_fig = sys.argv[1] + '/test_fig.png'

    plot_single_bact(pop_freqs, patient_freq, output_fig)
