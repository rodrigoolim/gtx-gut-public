from glob import glob
import gzip
from typing import List
import subprocess
import os
from os import path, mkdir, remove
import json
from time import sleep
import numpy as np
from scipy import stats


def run_command(cmd_vec: List[str], stdin="", no_output=True):
    '''Executa um comando no shell e retorna a saída (stdout) dele.'''
    cmd_vec = " ".join(cmd_vec)
    #logging.info(cmd_vec)
    if no_output:
        #print(cmd_vec)
        result = subprocess.run(cmd_vec, shell=True)
        return ""
    else:
        result = subprocess.run(cmd_vec, capture_output=True, 
            text=True, input=stdin, shell=True)
        return result.stdout

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def file_length(filepath):
    return os.path.getsize(filepath)

def open_writer(filepath):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'wt')
    else:
        return open(filepath, 'w')

def open_reader(filepath):
    if filepath.endswith('.gz'):
        stream = None
        try:
            stream = gzip.open(filepath, 'rt')
            testline = stream.readline()
        except OSError as err:
            print(filepath, 'is not a gzip file')
            print(err)
            if stream:
                stream.close()
            return open(filepath, 'r')
        return gzip.open(filepath, 'rt')
    else:
        stream = None
        try:
            stream = open(filepath, 'r')
            testline = stream.readline()
        except UnicodeDecodeError as err:
            print(filepath, 'is not a text file')
            print(err)
            if stream:
                stream.close()
            return gzip.open(filepath, 'rt')

        return open(filepath, 'r')

def get_sample_label(bam_path):
    return os.path.basename(bam_path).replace(".bam", "")

def chr_translation_dict():
    chr_names_path = "chr_names.tsv"
    download_cmd = ['wget', '-O '+chr_names_path, 
        'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/'
        +'GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report+ucsc_names.txt']
    run_command(download_cmd, no_output=True)
    chr_name_aliases = {}
    invalid = 0
    for rawline in open(chr_names_path,'r'):
        if not rawline.startswith('#'):
            cells = rawline.rstrip('\n').split('\t')
            realname = cells[-1]
            for index in [0, 4, 6, -1]:
                chr_name_aliases[cells[index]] = realname
    chr_name_aliases["chrMT"] = "chrM"
    for key in list(chr_name_aliases.keys()):
        chr_name_aliases['CHR_'+key] = chr_name_aliases[key]
    #logging.info("Invalid contig names: "+str(invalid))

    assert 'chr1' in chr_name_aliases
    os.remove(chr_names_path)
    return chr_name_aliases

def list_files(dir: str, extension: str, not_include=[]):
    res = []
    # Iterate directory
    for file in os.listdir(dir):
        # check only text files
        print(file)
        if file.endswith(extension) and not file in not_include:
            res.append(file)
            print('included')
        else:
            print('not included')
    return [os.path.join(dir, x) for x in res]

def find_overlapping (start1: int, end1: int, start2: int, end2: int):
    return max(0, min(end1, end2) - max(start1, start2) + 1)

def all_exist(path_list):
    for x in path_list:
        if not os.path.exists(x):
            return False
    return True

def one_exists(path_list):
    for x in path_list:
        if os.path.exists(x):
            return True
    return False

def gc_content(seq: str):
    return (seq.count('G')+seq.count('C')) / len(seq)


def find_closest(x, search_list):
    dists = [abs(x-y) for y in search_list]
    min_i = 0
    for i in range(len(dists)):
        if dists[min_i] > dists[i]:
            min_i = i
    return min_i

def read_tsv(file_path: str):
    input_stream = open_reader(file_path)
    columns = input_stream.readline().rstrip('\n').split('\t')
    lines = []
    for rawline in input_stream:
        cells = rawline.rstrip('\n').split('\t')
        newline = {columns[i]: cells[i] for i in range(len(columns))}
        lines.append(newline)

    return lines, columns


def load_config():
    proj_dir = path.dirname(path.dirname(__file__))
    templates = glob(proj_dir+'/templates/*.docx')
    data_files = glob(proj_dir+'/data/*')
    config = {
        'proj_dir': proj_dir,
        'templates': {path.basename(p): p for p in templates},
        'data_files': {path.basename(p): p for p in data_files}
    }

    try:
        other_configs = json.load(open(proj_dir+'/config.json'))
        for key, val in other_configs.items():
            config[key] = val
    except FileNotFoundError:
        # config.json is optional in this public repository (see config.example.json);
        # none of the functions kept here (e.g. freq_stats) require it.
        pass
    return config

def value_at_modzscore(mean, mad, modzscore):
    return (modzscore*mad + 0.6745*mean) / 0.6745

def mod_zscore_calc(mean, mad, x):
    if mad > 0:
        return (0.6745*(x-mean))/mad
    else:
        return 0.0

# def freq_stats(pop, patient):
    #increased_names = ['Increased', 'Very High']
    #decreased_names = ['Very Low', 'Decreased']
    percs = np.array(sorted(pop))
    all_percs = np.array(sorted(pop + [patient]))
    patient_i = all_percs.tolist().index(patient)
    patient_percentile = (patient_i / len(all_percs))*100

    q1 = np.percentile(percs, 25)
    mean = percs.mean()
    median = np.median(percs)
    q3 = np.percentile(percs, 75)
    iqr = q3-q1

    low_limit = q1 - 1.4*iqr
    high_limit = q3 + 1.4*iqr
    
    low_ths  = [q1 - 1.5*iqr, q1]
    high_ths = [q3, q3 + 1.5*iqr]

    if patient < low_limit:
        classification = 'Very Low'
    elif patient  > high_limit:
        classification = 'Very High'
    else:
        classification = 'Normal'
   
    mad = np.absolute(percs - mean).mean()

    mod_zscore = mod_zscore_calc(mean, mad, patient)
    increased = patient > median

    labels = []

    if patient > 0:
        labels.append('present')
    else:
        labels.append('absent')
    
    if classification == 'Very Low':
        labels.append('decrease')
    elif classification == 'Very High':
        labels.append('increase')

    try:
        stats = {
            'patient_perc': patient,
            'std': np.std(percs),
            'median': np.median(percs),
            'mean': mean,
            'mad': mad,
            'l_outlier': q1 - iqr*1.5,
            'q1': q1,
            'iqr': iqr,
            'iqr_1.5': iqr*1.5,
            'q3': q3,
            'r_outlier': q3 + iqr*1.5,
            'mod_zscore': mod_zscore,
            'increased': bool(increased),
            'patient_percentile': patient_percentile,
            'classification': classification,
            #'classification_bounds': bounds,
            #'class_index': class_index,
            'is_present': bool(patient > 0),
            'low_ths': low_ths,
            'high_ths': high_ths,
            #'classifications': classifications,
            'labels': labels
        }

        return stats
    except UnboundLocalError as err:
        print(err)
        print(mod_zscore)
        print(mean, mad, patient)
        quit(1)


def read_json_locked(json_path: str, miliseconds_timeout: int = 500):
    lock_path = json_path + '.lock'
    while path.exists(lock_path):
        sleep(miliseconds_timeout/1000)
    open(lock_path, 'w').write('locked')
    json_content = json.load(open(json_path, 'r'))
    remove(lock_path)
    return json_content

def write_json_locked(json_content: dict, json_path: str, miliseconds_timeout: int = 500, max_timeout: int = 4000):
    lock_path = json_path + '.lock'
    total_timeout = 0
    while path.exists(lock_path) and total_timeout <= max_timeout:
        sleep(miliseconds_timeout/1000)
        total_timeout += miliseconds_timeout
    open(lock_path, 'w').write('locked')
    json.dump(json_content, open(json_path, 'w'), indent=4)
    remove(lock_path)

config = load_config()




##########################

def detect_distribution_type(pop):
    """Detect the type of distribution based on population data"""
    pop = np.array(pop)
    zero_proportion = np.sum(pop == 0) / len(pop)
    non_zero = pop[pop > 0]
    
    if len(non_zero) == 0:
        return 'zero_only'
    
    # If no zeros, check for normality using Shapiro-Wilk
    if zero_proportion == 0:
        if len(pop) >= 3 and np.std(pop) > 1e-10:  # Use full population since no zeros
            try:
                _, p_value = stats.shapiro(pop)
                if p_value <= 0.05:
                    return 'normal'
            except ValueError:
                pass
        return 'zero_inflated'  # Default for non-normal or insufficient data
    
    # For populations with zeros, check for bimodality or zero-inflated
    if len(non_zero) > 10:
        kernel = stats.gaussian_kde(non_zero)
        x = np.linspace(min(non_zero), max(non_zero), 100)
        density = kernel(x)
        peaks = np.where((density[1:-1] > density[:-2]) & (density[1:-1] > density[2:]))[0]
        
        if len(peaks) > 1 and zero_proportion > 0.1:
            return 'bimodal'
        elif zero_proportion > 0.3:
            return 'zero_inflated'
    
    # If Shapiro-Wilk is feasible for non-zero values, check normality
    if len(non_zero) >= 3 and np.std(non_zero) > 1e-10:
        try:
            _, p_value = stats.shapiro(non_zero)
            if p_value <= 0.05:
                return 'normal'
        except ValueError:
            pass
    
    return 'zero_inflated'

def mod_zscore_calc(mean, mad, value):
    if mad == 0:
        return 0
    return 0.6745 * (value - mean) / mad

def freq_stats(pop, patient):
    pop = np.array(pop)
    distribution_type = detect_distribution_type(pop)
    
    # Initialize labels
    labels = ['present' if patient > 0 else 'absent']
    
    # Default stats for all cases
    percs = np.array(sorted(pop))
    #all_percs = np.array(sorted(pop + [patient]))

    # with open( "/tmp/debug.txt", "a" ) as f:
    #     f.write( f"{all_percs}\n{ patient }\n\n" )
    
    # Find patient index using tolerance-based comparison
    # patient_i = np.where(np.isclose(all_percs, patient, rtol=1e-8, atol=1e-8))[0]
    # if len(patient_i) == 0:
    #     patient_i = np.searchsorted(all_percs, patient)
    #     #patient_i = min(patient_i, len(all_percs) - 1)
    # else:
    #     patient_i = patient_i[0]
    # patient_percentile = (patient_i / len(all_percs)) * 100

    patient_i = np.searchsorted(percs, patient, side='left')
    patient_percentile = (patient_i / (len(percs) + 1)) * 100

    
    if distribution_type == 'zero_only':
        return {
            'patient_perc': patient,
            'std': 0.0,
            'median': 0.0,
            'mean': 0.0,
            'mad': 0.0,
            'l_outlier': 0.0,
            'q1': 0.0,
            'iqr': 0.0,
            'iqr_1.5': 0.0,
            'q3': 0.0,
            'r_outlier': 0.0,
            'mod_zscore': 0.0,
            'increased': False,
            'patient_percentile': patient_percentile,
            'classification': 'Normal',
            'is_present': bool(patient > 0),
            'low_ths': [0.0, 0.0],
            'high_ths': [0.0, 0.0],
            'labels': labels
        }
    
    if distribution_type == 'normal':
        # Normal distribution handling
        q1 = np.percentile(percs, 25)
        mean = percs.mean()
        median = np.median(percs)
        q3 = np.percentile(percs, 75)
        iqr = q3 - q1
        
        # Use standard 1.5*IQR for outlier detection
        low_limit = q1 - 1.5 * iqr
        high_limit = q3 + 1.5 * iqr
        low_ths = [q1 - 1.5 * iqr, q1]
        high_ths = [q3, q3 + 1.5 * iqr]
        
        if patient < low_limit:
            classification = 'Very Low'
            labels.append('decrease')
        elif patient > high_limit:
            classification = 'Very High'
            labels.append('increase')
        else:
            classification = 'Normal'
        
        mad = np.abs(percs - mean).mean()
        mod_zscore = mod_zscore_calc(mean, mad, patient)
        increased = patient > median
    else:
        # Handle zero-inflated or bimodal distributions
        non_zero = percs[percs > 0]
        if len(non_zero) == 0:
            return {
                'patient_perc': patient,
                'std': 0.0,
                'median': 0.0,
                'mean': 0.0,
                'mad': 0.0,
                'l_outlier': 0.0,
                'q1': 0.0,
                'iqr': 0.0,
                'iqr_1.5': 0.0,
                'q3': 0.0,
                'r_outlier': 0.0,
                'mod_zscore': 0.0,
                'increased': False,
                'patient_percentile': patient_percentile,
                'classification': 'Normal',
                'is_present': bool(patient > 0),
                'low_ths': [0.0, 0.0],
                'high_ths': [0.0, 0.0],
                'labels': labels
            }
        
        # For bimodal, focus on secondary peak
        if distribution_type == 'bimodal':
            kernel = stats.gaussian_kde(non_zero)
            x = np.linspace(min(non_zero), max(non_zero), 100)
            density = kernel(x)
            peaks = x[np.where((density[1:-1] > density[:-2]) & (density[1:-1] > density[2:]))[0]]
            reference_point = peaks[peaks > 0].min() if len(peaks) > 1 else np.median(non_zero)
        else:
            reference_point = np.median(non_zero)
        
        # Calculate stats on non-zero values
        q1 = np.percentile(non_zero, 25)
        mean = np.mean(non_zero)
        median = np.median(non_zero)
        q3 = np.percentile(non_zero, 75)
        iqr = q3 - q1
        mad = np.abs(non_zero - mean).mean()
        
        # Modified classification: only "Very High" or "Normal"
        high_limit = reference_point + 1.4 * iqr
        low_ths = [reference_point - 1.5 * iqr, reference_point]
        high_ths = [reference_point, reference_point + 1.5 * iqr]
        
        if patient == 0 or patient <= reference_point:
            classification = 'Normal'
        else:
            classification = 'Very High' if patient > high_limit else 'Normal'
            if classification == 'Very High':
                labels.append('increase')
        
        mod_zscore = mod_zscore_calc(mean, mad, patient)
        increased = patient > median
    
    return {
        'patient_perc': patient,
        'std': np.std(percs),
        'median': median,
        'mean': mean,
        'mad': mad,
        'l_outlier': q1 - iqr * 1.5,
        'q1': q1,
        'iqr': iqr,
        'iqr_1.5': iqr * 1.5,
        'q3': q3,
        'r_outlier': q3 + iqr * 1.5,
        'mod_zscore': mod_zscore,
        'increased': bool(increased),
        'patient_percentile': patient_percentile,
        'classification': classification,
        'is_present': bool(patient > 0),
        'low_ths': low_ths,
        'high_ths': high_ths,
        'labels': labels
    }