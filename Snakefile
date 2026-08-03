"""
Snakefile for the GTX-GUT reproducibility companion.

Takes already quality-filtered FASTQ (see note below) through QIIME 2 import,
DADA2 denoising, GreenGenes taxonomic classification, and abundance export,
then calls generate_report_json.py for the final Tukey/IQR classification
against the public reference population.

Quality filtering (BBDuk) is deliberately NOT part of this workflow: the
private product's exact BBDuk parameters are not disclosed (see the
dissertation, Metodologia), so this Snakefile instead expects the *already
filtered* FASTQ as input -- e.g. the `*_gcfix.fastq.gz` files bundled in
examples/fastq/ for the eight samples discussed in the dissertation's
Experimentos chapter. If you are starting from raw FASTQ of your own, run
your own quality control first.

Usage:
    snakemake --cores 1 --configfile config.json --config \\
        R1=examples/fastq/SAMPLE_R1_gcfix.fastq.gz \\
        R2=examples/fastq/SAMPLE_R2_gcfix.fastq.gz \\
        output_dir=results/SAMPLE
"""
from os import path

configfile: "config.json"
builds_path = config['singularity_builds']

SINGLE_END = not config.get('R2')


def singularity_exec_cmd(img_name, volumes, cmd):
    volumes_str = ""
    volumes = sorted(volumes, key=lambda tp: (len(tp), tp[-1]))
    for v in volumes:
        v = list(v)
        if v[0][0] != '/':
            v[0] = path.abspath(v[0])
        if len(v) == 2:
            volumes_str += "," + v[0].rstrip('/') + ':' + v[1]
        elif len(v) == 1:
            volumes_str += "," + v[0].rstrip('/')
    return ("singularity exec --bind " + volumes_str.strip(',')
            + " " + builds_path + "/" + img_name + " " + cmd)


# --------------------------------------------------------------------------
# Container builds (all from public sources: the official QIIME 2 image, and
# the two .def files bundled in this repo)
# --------------------------------------------------------------------------

rule build_qiime2_amplicon:
    output:
        builds_path + "/amplicon_2023.9.sif"
    shell:
        "singularity pull " + builds_path + "/amplicon_2023.9.sif docker://quay.io/qiime2/amplicon:2023.9"
        + " && chmod a+rwx " + builds_path + "/amplicon_2023.9.sif"

rule build_qiime2r:
    input:
        "qiime2r/qiime2r.def"
    output:
        builds_path + "/qiime2r.sif"
    shell:
        "singularity build --fakeroot " + builds_path + "/qiime2r.sif qiime2r/qiime2r.def"
        + " && chmod a+rwx " + builds_path + "/qiime2r.sif"

rule build_microbio_py:
    input:
        "microbio_py/microbio_py.def"
    output:
        builds_path + "/microbio_py.sif"
    shell:
        "singularity build --fakeroot " + builds_path + "/microbio_py.sif microbio_py/microbio_py.def"
        + " && chmod a+rwx " + builds_path + "/microbio_py.sif"

rule build:
    input:
        builds_path + "/amplicon_2023.9.sif",
        builds_path + "/qiime2r.sif",
        builds_path + "/microbio_py.sif"


# --------------------------------------------------------------------------
# Per-sample pipeline
# --------------------------------------------------------------------------

if all(x in config for x in ['R1', 'output_dir']):
    # QIIME 2 requires absolute paths in the manifest; accept relative paths
    # on the command line for convenience, but resolve them here.
    config['R1'] = path.abspath(config['R1'])
    config['output_dir'] = path.abspath(config['output_dir'])
    if config.get('R2'):
        config['R2'] = path.abspath(config['R2'])

    manifest_path = config['output_dir'] + '/manifest.csv'
    metadata_path = config['output_dir'] + '/metadados.csv'

    qza_raw = config['output_dir'] + "/rawdata.qza"
    raw_derep_table = config['output_dir'] + "/raw_derep_table.qza"
    raw_derep_seqs = config['output_dir'] + "/raw_derep_seqs.qza"
    raw_derep_stats = config['output_dir'] + "/raw_derep_stats.qza"
    raw_taxonomy = config['output_dir'] + "/raw_taxonomy.qza"
    raw_derep_table_collapsed7 = config['output_dir'] + "/raw_derep_table_collapsed7.qza"

    specie_full_abundance = config['output_dir'] + "/Specie_full_abundance.csv"
    relative_by_levels = specie_full_abundance.replace('.csv', '.relative_by_levels.csv')
    shannon_csv = config['output_dir'] + "/shannon_diversity2.csv"
    results_json = config['output_dir'] + "/results.json"

    rule make_manifest:
        input:
            [config['R1']] + ([config['R2']] if not SINGLE_END else [])
        output:
            manifest_path,
            metadata_path
        shell:
            "mkdir -p " + config['output_dir']
            + " && python3 pre_processing/manifester.py"
            + " " + manifest_path
            + " " + metadata_path
            + " " + config['R1']
            + ("" if SINGLE_END else " " + config['R2'])

    rule import_qza:
        input:
            builds_path + "/amplicon_2023.9.sif",
            manifest_path,
            config['R1']
        output: qza_raw
        shell:
            singularity_exec_cmd(
                img_name="amplicon_2023.9.sif",
                volumes=[[config['output_dir'], "/results"], [path.dirname(config['R1'])]],
                cmd="qiime tools import"
                    + " --input-format " + ("SingleEndFastqManifestPhred33" if SINGLE_END else "PairedEndFastqManifestPhred33")
                    + " --input-path /results/" + path.basename(manifest_path)
                    + " --output-path /results/" + path.basename(qza_raw)
                    + " --type 'SampleData[" + ("SequencesWithQuality" if SINGLE_END else "PairedEndSequencesWithQuality") + "]'"
            )

    rule denoise:
        input: qza_raw
        output:
            raw_derep_table,
            raw_derep_seqs,
            raw_derep_stats
        threads: 1
        shell:
            singularity_exec_cmd(
                img_name="amplicon_2023.9.sif",
                volumes=[[config['output_dir'], "/results"]],
                cmd=("qiime dada2 denoise-single --verbose --p-n-threads {threads}"
                     + " --i-demultiplexed-seqs /results/" + path.basename(qza_raw)
                     + " --p-trunc-len 0"
                     if SINGLE_END else
                     "qiime dada2 denoise-paired --verbose --p-n-threads {threads}"
                     + " --i-demultiplexed-seqs /results/" + path.basename(qza_raw)
                     + " --p-trunc-len-f 0 --p-trunc-len-r 0")
                    + " --o-table /results/" + path.basename(raw_derep_table)
                    + " --o-representative-sequences /results/" + path.basename(raw_derep_seqs)
                    + " --o-denoising-stats /results/" + path.basename(raw_derep_stats)
            )

    rule classifysklearn:
        input:
            raw_derep_seqs,
            config['greengenes_classifier_path']
        output:
            raw_taxonomy
        threads: 1
        shell:
            singularity_exec_cmd(
                img_name="amplicon_2023.9.sif",
                volumes=[[config['output_dir'], "/results"], [config['greengenes_classifier_path']]],
                cmd="qiime feature-classifier classify-sklearn"
                    + " --i-classifier " + path.abspath(config['greengenes_classifier_path'])
                    + " --i-reads /results/" + path.basename(raw_derep_seqs)
                    + " --p-confidence 0"
                    + " --p-read-orientation same"
                    + " --o-classification /results/" + path.basename(raw_taxonomy)
            )

    rule taxa_collapse_species:
        input:
            raw_derep_table,
            raw_taxonomy
        output:
            raw_derep_table_collapsed7
        shell:
            singularity_exec_cmd(
                img_name="amplicon_2023.9.sif",
                volumes=[[config['output_dir'], "/results"]],
                cmd="qiime taxa collapse"
                    + " --i-table /results/" + path.basename(raw_derep_table)
                    + " --i-taxonomy /results/" + path.basename(raw_taxonomy)
                    + " --p-level 7"
                    + " --o-collapsed-table /results/" + path.basename(raw_derep_table_collapsed7)
            )

    rule export_species_csv:
        input:
            builds_path + "/qiime2r.sif",
            raw_derep_table_collapsed7
        output:
            specie_full_abundance
        shell:
            singularity_exec_cmd(
                img_name="qiime2r.sif",
                volumes=[[config['output_dir'], "/results"], ["qiime2r", "/qiime2r"]],
                cmd="Rscript /qiime2r/export.R"
                    + " /results/" + path.basename(raw_derep_table_collapsed7)
                    + " /results/" + path.basename(specie_full_abundance)
            )

    rule preprocess_patient_abundances:
        input:
            builds_path + "/microbio_py.sif",
            specie_full_abundance,
            "src/abundance_processing.py",
            "src/process_patient.py"
        output:
            relative_by_levels
        shell:
            singularity_exec_cmd(
                img_name="microbio_py.sif",
                volumes=[["src", "/src"], [config['output_dir'], "/output"]],
                cmd="python /src/process_patient.py /output/" + path.basename(specie_full_abundance)
            )

    rule shannon:
        input:
            builds_path + "/qiime2r.sif",
            raw_derep_table
        output:
            shannon_csv
        shell:
            singularity_exec_cmd(
                img_name="qiime2r.sif",
                volumes=[["src", "/src"], [config['output_dir'], "/results"]],
                cmd="Rscript /src/shannon.R"
                    + " /results/" + path.basename(raw_derep_table)
                    + " /results/" + path.basename(shannon_csv)
            )

    rule generate_json:
        input:
            builds_path + "/microbio_py.sif",
            relative_by_levels,
            shannon_csv,
            "src/generate_report_json.py",
            "src/util.py"
        output:
            results_json
        params:
            patient_id=lambda wc: path.basename(config['R1']).split('_')[0]
        shell:
            singularity_exec_cmd(
                img_name="microbio_py.sif",
                volumes=[["src", "/src"], ["data", "/data"], [config['output_dir'], "/output"]],
                cmd="python /src/generate_report_json.py"
                    + " --patient-abundance /output/" + path.basename(relative_by_levels)
                    + " --patient-id {params.patient_id}"
                    + " --patient-shannon-csv /output/" + path.basename(shannon_csv)
                    + " --reference-dir /data"
                    + " --output /output/" + path.basename(results_json)
            )

    rule all:
        input: results_json
