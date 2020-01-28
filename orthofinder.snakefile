"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s orthofinder.snakefile --use-conda #--cores 26 --cluster "sbatch -t 10:00:00 -N 11 -n 26 -p bmm --mem=60gb" --jobs 5 -k --rerun-incomplete 
"""

import os
import pandas as pd

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

# 10 subset
samplelist = "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/ten_haptophytes.txt"
samplelist_just_samples = "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/ten_haptophytes_samplelist.txt"
SAMPLES = [x.strip().split('\t')[0] for x in open(samplelist, "r")]

pep_src = config.get("pep_src", "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep")
out_dir = config.get("out_dir", "mmetsp_orthofinder")
pep_fasta = os.path.join(out_dir, "pep_fasta")
name = "ten_hapto"

pepfiles = []
if "MMETSP0251" in SAMPLES:
    SAMPLES.remove("MMETSP0251")
    pepfiles = [os.path.join(pep_fasta, "MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep.fasta")] 
if "MMETSP0754" in SAMPLES:
    SAMPLES.remove("MMETSP0754")
pepfiles += expand(os.path.join(pep_fasta, "{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta"), sample=SAMPLES)

rule all:
    input: os.path.join(out_dir, name, "Orthogroups.csv")


# orthofinder runs off of entire directory
rule prep_pepdir:
    input: samplelist_just_samples
    output: pepfiles 
    params:
        pep_src=pep_src,
        pep_dest=pep_fasta
    shell:
        """
        python copy_samplefiles.py {input} --source {params.pep_src} --destination {params.pep_dest}
        """

rule run_orthofinder:
    input:
        pep=pepfiles
    output:
        os.path.join(out_dir, name, "Orthogroups.csv")
    params:
        search="diamond",
        tmp_dir=os.path.join(out_dir, "tmp_output"),
        out_dir=out_dir,
        pep_dir=pep_fasta,
        run_name=name
    log: os.path.join(out_dir, "logs", name + "_orthofinder.log")
    benchmark: os.path.join(out_dir, "logs", name +"_orthofinder.benchmark")
    threads: 28
    shadow: "shallow" ## this means tmpdir will be on local scratch (using --shadow-prefix) --> faster? 
    conda: "orthofinder-env.yml"
    shell:
        #orthofinder -b {params.tmp_dir} -t {threads} -a {threads} -n {params.run_name} > {log} 2>&1  
        """
        orthofinder -f {params.pep_dir} -S {params.search} -t {threads} -a {threads} -o {params.tmp_dir} -n {params.run_name} > {log} 2>&1  
        mv {params.tmp_dir}/* {params.out_dir} 
        """

#  this is not used for LM - # nodes depends on memory asked for SBATCH --ntasks-per-node 16. each job = 1 core/48GB requested

## folder "haptophyta_pep" has all the peptide sequences for the MMETSP dataset, haptophyte samples
#SIMPLE USAGE:
#Run full OrthoFinder analysis on FASTA format proteomes in <dir>
#  orthofinder [options] -f <dir>
# -t = num parallel search threads
# -t <int>          Number of parallel sequence search threads [Default = 28]
# -a <int>          Number of parallel analysis threads [Default = 1]
# -M <txt>          Method for gene tree inference. Options 'dendroblast' & 'msa'
#                   [Default = dendroblast]
# -S <txt>          Sequence search program [Default = diamond]
#                   Options: blast, mmseqs, blast_gz, diamond
# -p <dir>          Write the temporary pickle files to <dir>
# -o <txt>          Non-default results directory
# -n <txt>          Name to append to the results directory

#orthofinder -f haptophyta_pep -S diamond -t 16 -a 16 -p $PROJECT -o /pylon5/mc5phkp/ntpierce/haptophyta_orthofinder_results -n _june_2019
