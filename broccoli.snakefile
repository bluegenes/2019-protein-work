"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  broccoli.snakefile --use-conda --dryrun  # remove `--dryrun` to actually run
"""

import os
import pandas as pd

def read_samples(samples_file):
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator).set_index(["sample", "unit"], drop=False)
            samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str).set_index(["sample", "unit"], drop=False)
            samples['name'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples

# read in elvers csv, grab sample names
#samples_csv = "/home/ntpierce/2019-protein-work/haptophyta_samples_wref.csv"

configfile: "haptophyta_broccoli_config.yml"
samples_csv = config["samples"]
samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist() 
pep_dir = config["pep_dir"]
out_dir = config.get("out_dir", "haptophyta_broccoli")
pep_fasta  = os.path.join(out_dir, "pep_fasta")

#input files = pep_dir/{sample}.pep
#MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
# pass pep_dir to broccoli

#pep_files = expand("{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep", sample=SAMPLES)
#MAYBE A PROBLEM: peptide names (before first space) need to be unique. Is this just within the file, or across files? If across, add MMETSP sample number to each?
#import pdb;pdb.set_trace()

rule all:
    input: os.path.join(out_dir, "dir_step3/orthologous_groups.txt")


# if the peptide files don't end in ".fasta", we need to add this. Also, since the input to broccoli is a directory,
# we need to copy *just* the relevant pep files to a new location
rule prep_pepdir:
    output: expand(os.path.join(pep_fasta, "{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta"), sample=SAMPLES)
    params:
        dest=pep_fasta
    shell:
        """
        python copy_samplefiles.py haptophyta.txt --source {input} --destination {params.dest}
        """

rule run_broccoli: 
    input: expand(os.path.join(pep_fasta, "{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta"), sample=SAMPLES)
    output: 
        step3="dir_step3/orthologous_groups.txt",
        step4="dir_step4/orthologous_pairs.txt",
        #"dir_step3/chimeric_proteins.txt",
    log: os.path.join(out_dir, "logs/broccoli.log") 
    params:
        pep_dir=pep_fasta
    conda: "broccoli.yml"
    threads: 10
    shell: 
        """
        python Broccoli/broccoli.py -dir {params.pep_dir} -threads {threads} 2> {log}
        """
        #python Broccoli/broccoli.py -dir Broccoli/example_dataset

rule move_broccoli_results:
    input: rules.run_broccoli.output 
    output: os.path.join(out_dir, "dir_step3/orthologous_groups.txt")
    params:
        dest= out_dir
    shell:
        """
        mv dir_step1 dir_step2 dir_step3 dir_step4 {params.dest}
        """

