"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -k -n

Note: some protein, rna, and cds downloads will fail because those datasets
are unavailable for that sample. Use the `-k` flag to tell snakemake to keep
going despite these failures. 
"""

import os
import re
import pandas as pd
#from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
#HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

def build_genbank_url(row, base_url):
    sample_name = row[2].split('/')[-1]
    alpha,first,second,third = re.match("([A-Z]+)_(\d{3})(\d{3})(\d{3})", sample_name).groups()
    name = sample_name.split('_genomic.fna.gz')[0]
    # build folder path
    genbank_path =os.path.join(base_url,alpha,first,second,third,name)
    # add to df
    row["name"] = name
    row["base_url"] = os.path.join(genbank_path, name)
    return row 

# testing with some hardcoded vals
outbase = "smash_testing"
genome, protein, rna, cds = True, True, True, True, #False, False, False 
csv_list = ["../2018-test_datasets/gingivalis.csv", "../2018-test_datasets/bacteroides.csv", "../2018-test_datasets/denticola.csv"]

# some hardcoded vars that can probably stay
genbank_url = 'https://ftp.ncbi.nih.gov/genomes/all'
genomic_dir, protein_dir, rna_dir, cds_dir = "genomic", "protein", "rna", "cds"
genomic_ext = "_genomic.fna.gz"
protein_ext = "_protein.faa.gz"
rna_ext = "_rna_from_genomic.fna.gz"
cds_ext = "_cds_from_genomic.fna.gz"

# build output files and sample:genbank_url dict
data_files=[]
linkdb={}
for csv in csv_list:
    genomeInfo = pd.read_csv(csv, header=None)
    csv_basename = os.path.basename(csv.rsplit('.', 1)[0])#assuming good csv naming
    outdir = os.path.join(outbase, csv_basename)
    genomeInfo = genomeInfo.apply(build_genbank_url, axis=1, base_url=genbank_url)
    csv_linkdb = genomeInfo.set_index('name').to_dict()["base_url"]
    samples = csv_linkdb.keys()
    csv_outfiles=[]
    if genome:
        csv_outfiles += [os.path.join(outdir, genomic_dir, x + genomic_ext) for x in samples]
    if protein:
        csv_outfiles += [os.path.join(outdir, protein_dir, x + protein_ext) for x in samples]
    if rna:
        csv_outfiles += [os.path.join(outdir, rna_dir, x + rna_ext) for x in samples]
    if cds:
        csv_outfiles += [os.path.join(outdir, cds_dir, x + cds_ext) for x in samples]
    # add this csv's info to outfile list, linkdb
    data_files+=csv_outfiles
    linkdb[csv_basename] = csv_linkdb
    
rule all:
    input: data_files 

# download datasets
rule get_genomic_datasets:
    input: lambda wildcards: FTP.remote(linkdb[wildcards.csv_name][wildcards.sample] + genomic_ext, static=True, keep_local=True, immediate_close=True)
    #output: f"{outbase}/{{csv_name}}/{genomic_dir}/{{sample}}_genomic.fna.gz"
    output: os.path.join(outbase, "{csv_name}", genomic_dir, "{sample}" + genomic_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_protein_datasets:
    input: lambda wildcards: FTP.remote(linkdb[wildcards.csv_name][wildcards.sample] + protein_ext, static=True, keep_local=True, immediate_close=True)
    output: os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_rna_datasets:
    input: lambda wildcards: FTP.remote(linkdb[wildcards.csv_name][wildcards.sample] + rna_ext, static=True, keep_local=True, immediate_close=True)
    output: os.path.join(outbase, "{csv_name}", rna_dir, "{sample}" + rna_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_cds_datasets:
    input: lambda wildcards: FTP.remote(linkdb[wildcards.csv_name][wildcards.sample], static=True, keep_local=True, immediate_close=True)
    output: os.path.join(outbase, "{csv_name}", cds_dir, "{sample}" + cds_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"
