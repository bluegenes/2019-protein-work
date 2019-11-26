"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import re
import pandas as pd
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

def build_genbank_urls(row, base_url, outdir):
    sample_name = row[2].split('/')[-1]
    alpha,first,second,third = re.match("([A-Z]+)_(\d{3})(\d{3})(\d{3})", sample_name).groups()
    name = sample_name.split('_genomic.fna.gz')[0]
    # add to df
    genbank_path =os.path.join(base_url,alpha,first,second,third,name)
    row["name"] = name
    row["genome_file"] = os.path.join(outdir, genomic_dir, name + genomic_ext)
    row["genome_url"] = os.path.join(genbank_path, name + genomic_ext)
    
    row["protein_file"] = os.path.join(outdir, protein_dir, name + protein_ext)
    row["protein_url"] = os.path.join(genbank_path, name + protein_ext)
    
    row["rna_file"] = os.path.join(outdir, rna_dir, name + rna_ext)
    row["rna_url"] = os.path.join(genbank_path, name + rna_ext)
    
    row["cds_file"] = os.path.join(outdir, cds_dir, name + cds_ext)
    row["cds_url"] = os.path.join(genbank_path, name + cds_ext)
    return row 

# testing with some hardcoded vals
csv_list = ["../2018-test_datasets/gingivalis.csv", "../2018-test_datasets/bacteroides.csv", "../2018-test_datasets/denticola.csv"]

genbank_url = 'https://ftp.ncbi.nih.gov/genomes/all'
genome, protein, rna, cds = True, True, True, True
genomic_dir, protein_dir, rna_dir, cds_dir = "genomic", "protein", "rna", "cds"
genomic_ext = "_genomic.fna.gz"
protein_ext = "_protein.faa.gz"
rna_ext = "_rna_from_genomic.fna.gz"
cds_ext = "_cds_from_genomic.fna.gz"

outbase = "smash_testing"
genome_files, protein_files, rna_files, cds_files = [],[],[],[]
glinkdb, plinkdb, rlinkdb, clinkdb = {}, {}, {}, {}
for csv in csv_list:
    genomeInfo = pd.read_csv(csv, header=None)
    csv_basename = os.path.basename(csv.rsplit('.', 1)[0])#assuming good csv naming
    out_dir = os.path.join(outbase, csv_basename)
    genomeInfo = genomeInfo.apply(build_genbank_urls, axis=1, base_url=genbank_url, outdir=out_dir)
    # build dicts
    if genome:
        genome_linkdb = genomeInfo.set_index('name').to_dict()["genome_url"]
        genome_files += genomeInfo["genome_file"].tolist()
        glinkdb.update(genome_linkdb) 
    if protein:
        protein_linkdb = genomeInfo.set_index("name").to_dict()["protein_url"]
        protein_files += genomeInfo["protein_file"].tolist() 
        plinkdb.update(protein_linkdb) 
    if rna:
        rna_linkdb = genomeInfo.set_index("name").to_dict()["rna_url"]
        rna_files += genomeInfo["rna_file"].tolist()
        rlinkdb.update(rna_linkdb) 
    if cds:
        cds_linkdb = genomeInfo.set_index("name").to_dict()["cds_url"]
        cds_files += genomeInfo["cds_file"].tolist()
        clinkdb.update(cds_linkdb) 
    
rule all:
    input: genome_files + protein_files + rna_files + cds_files 


### TBD: # build sig outputs
#genome_exts = [".sig"]
#ksizes = ["21", "31", "51"]

#genome_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]
#protein_sigs = [x.split('.faa.gz')[0] + ".sig" for x in genome_files]

#rna_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]
#cds_sigs = [x.split('.fna.gz')[0] + ".sig" for x in genome_files]



# download rna, protein, genomic sequences
rule get_genomic_datasets:
    input: lambda wildcards: HTTP.remote(glinkdb[wildcards.sample])
    #output: f"{outbase}/{{csv_name}}/{genomic_dir}/{{sample}}_genomic.fna.gz"
    output: os.path.join(outbase, "{csv_name}", genomic_dir, "{sample}" + genomic_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_protein_datasets:
    input: lambda wildcards: HTTP.remote(plinkdb[wildcards.sample])
    output: os.path.join(outbase, "{csv_name}", protein_dir, "{sample}" + protein_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_rna_datasets:
    input: lambda wildcards: HTTP.remote(rlinkdb[wildcards.sample])
    output: os.path.join(out_dir, "{csv_name}", rna_dir, "{sample}" + rna_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"

rule get_cds_datasets:
    input: lambda wildcards: HTTP.remote(clinkdb[wildcards.sample])
    output: os.path.join(out_dir, "{csv_name}", rna_dir, "{sample}" + cds_ext)
    conda: "dl-test-datasets.yml"
    shell: "mv {input} {output} 2> {log}"


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

