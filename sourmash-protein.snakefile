"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import re
import pandas as pd
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# testing with some hardcoded vals
genomic_ext = "_genomic.fna.gz"
protein_ext = "_protein.faa.gz"
rna_ext = "_rna_from_genomic.fna.gz"
cds_ext = "_cds_from_genomic.fna.gz"

outbase = "smash_testing"
    
rule all:
    input: sig_files

sig_files = []

csv_list = ["../2018-test_datasets/gingivalis.csv", "../2018-test_datasets/bacteroides.csv", "../2018-test_datasets/denticola.csv"]
### TBD: # build sig outputs
#genome_exts = [".sig"]
#ksizes = ["21", "31", "51"]

#genome_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]
#protein_sigs = [x.split('.faa.gz')[0] + ".sig" for x in genome_files]

#rna_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]
#cds_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]

# compute nucleotide sigs
rule compute_genomic:
    input: os.path.join(outbase, "{csv_name}", genomic_dir, "{sample}" + genomic_ext)
    output: os.path.join(outbase, "{csv_name}", genomic_dir, "sigs", "{sample}_k{k}_scaled{scaled}_genomic.sig" )
    params:
        scaled=2000,
        k=[21,31,51],
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute rna sigs
rule compute_rna:
    input: os.path.join(outbase, "{csv_name}", rna_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{csv_name}", rna_dir, "sigs", "{sample}_k{k}_scaled{scaled}_rna_from_genomic.sig" )
    params:
        scaled=2000,
        k=[21,31,51],
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

rule compute_cds:
    input: os.path.join(outbase, "{csv_name}", cds_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{csv_name}", cds_dir, "sigs", "{sample}_k{k}_scaled{scaled}_cds_from_genomic.sig" )
    params:
        scaled=2000,
        k=[21,31,51],
        extra=" --protein "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# 6 frame translation --> signatures
### why are we using rna_from_genomic over cds .. check, bc cds = only coding = cleaner!?
rule compute_translated_rna:
    input: os.path.join(outbase, "{csv_name}", rna_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{csv_name}", rna_dir, "sigs", "{sample}_k{k}_scaled{scaled}_rna_from_genomic_translated.sig" )
    params:
        scaled=2000,
        k=[21,33,51],
        extra=" --protein "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

rule compute_translated_cds:
    input: os.path.join(outbase, "{csv_name}", cds_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{csv_name}", cds_dir, "sigs", "{sample}_k{k}_scaled{scaled}_cds_from_genomic_translated.sig" )
    params:
        scaled=2000,
        k=[21,33,51],
        extra=" --protein "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute protein sigs
rule compute_prot:
    input: os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{csv_name}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_protein.sig" )
    params:
        scaled=2000,
        k=[7,11,17],
        extra=" --input-is-protein --protein "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"


# compute dayhoff sigs
rule compute_dayhoff:
    input: os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{csv_name}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_dayhoff.sig" )
    params:
        scaled=2000,
        k=[7_11_17],
        extra=" --input-is-protein --dayhoff "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute hp sigs
rule compute_hp:
    input: os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{csv_name}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_hp.sig" )
    params:
        scaled=2000,
        k=[7,11,17],
        extra=" --input-is-protein --hp "
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# build all compare matrices: csv output
# build all compare matrices: numpy output

# sourmash plot each compare matrix numpy output


# run regex on this pandas series instead:
#sample_names = genomeInfo.iloc[:,2].str.split('/').str[-1]
#sample_info = sample_names.str.extract(r"([A-Z]+)_(\d{3})(\d{3})(\d{3})")
#sample_info.columns = ["alpha", "first", "second", "third"]
#sample_info['name'] = sample_names 

