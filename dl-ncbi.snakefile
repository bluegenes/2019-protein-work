"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -k -n
Note: some downloads will fail bc the rna/cds/protein/etc files do not
actually exist.
"""

import os
import re
import pandas as pd

def grab_name(row):
    sample_name = row[2].split('/')[-1]
    alpha,first,second,third = re.match("([A-Z]+)_(\d{3})(\d{3})(\d{3})", sample_name).groups()
    name = sample_name.split('_genomic.fna.gz')[0]
    accession = re.match("^([A-Z]+_[^_]*)", sample_name).groups()[0]
    # add to df
    row["name"] = name
    row["accession"] = accession
    return row 

# testing with some hardcoded vals
outbase = "smash-testing"
genome, protein, rna, cds = True, True, True, True #False, False, False 
csv_list = ["../2018-test_datasets/gingivalis.csv", "../2018-test_datasets/bacteroides.csv", "../2018-test_datasets/denticola.csv"]

# some hardcoded vars that can probably stay
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
    genomeInfo = genomeInfo.apply(grab_name, axis=1)
    csv_linkdb = genomeInfo.set_index('name').to_dict()["accession"]
    samples = genomeInfo["name"].tolist()
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

rule ncbi_genome_download:
    output: 
        genomic = os.path.join(outbase, "{csv_name}", genomic_dir, "{sample}" + genomic_ext),
        protein = os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext),
        rna = os.path.join(outbase, "{csv_name}", rna_dir, "{sample}" + rna_ext),
        cds = os.path.join(outbase, "{csv_name}", cds_dir, "{sample}" + cds_ext),
        failed = os.path.join(outbase, "ncbi_failed_downloads.txt"),
    params: 
        accession= lambda w: linkdb[w.csv_name][w.sample],
        ghost_files=True # make ghost files for failed downloads so snakemake is happy
    conda: "dl-test-datasets.yml"
    script: "ncbi-genome-download.wrapper.py"
