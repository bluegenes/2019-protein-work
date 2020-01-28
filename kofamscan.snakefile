"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
Function: Download and run kofamscan 
"""


import os
import sys
import itertools

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()


samplelist = "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/ten_haptophytes.txt"
SAMPLES = [x.strip().split('\t')[0] for x in open(samplelist, "r")]

pep_dir = "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep" #config["pep_dir"]
out_dir = config.get("out_dir", "mmetsp_kofamscan")
kofam_dir = os.path.join(out_dir, "kofam_dbs")
logs_dir = os.path.join(out_dir, "logs")

rule all:
    input: 
        expand(os.path.join(out_dir, "{sample}.txt"), sample = SAMPLES)

def get_pep(w):
#most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
#odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")


localrules: download_ko_list, download_ko_profiles, download_kofam_readme, download_kofamscan_program
# download ko list
rule download_ko_list:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(kofam_dir,"ko_list.gz"),
        unzipped = os.path.join(kofam_dir, "ko_list"),
    log: os.path.join(logs_dir,"download_ko_list.log")
    shell: 
        """
        mv {input} {output.gz_file}
        gunzip -c {output.gz_file} > {output.unzipped} 
        """

# download hmm profiles
rule download_ko_profiles:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(kofam_dir, "profiles.tar.gz"),
        unzipped = directory(os.path.join(kofam_dir,"profiles"))
    log: os.path.join(logs_dir,"download_profiles.log")
    params:
        kofam_dir = kofam_dir,
    shell: 
        """
        mv {input} {output.gz_file} 2> {log}
        tar xzf {output.gz_file} --directory {params.kofam_dir}
        """

rule download_kofam_readme:
    input: FTP.remote("ftp://ftp.genome.jp/pub/tools/kofamscan/README.md", static=True, keep_local=True, immediate_close=True)
    output: os.path.join(kofam_dir,"README.md")
    log: os.path.join(logs_dir, "download_kofam_readme.log")
    shell: "mv {input} {output} 2> {log}"

# download program (replace with conda install in the future!?)
rule download_kofamscan_program:
    input: FTP.remote("ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan-1.1.0.tar.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        unzipped = directory(os.path.join(kofam_dir,"kofamscan-1.1.0")),
        gz_file = os.path.join(kofam_dir, "kofamscan-1.1.0.tar.gz"),
        executable = os.path.join(kofam_dir,"kofamscan-1.1.0/exec_annotation")
    log: os.path.join(kofam_dir,"logs/download_kofamscan_program.log")
    params:
        kofam_dir = kofam_dir,
    shell: 
        """
        mv {input} {output.gz_file} 
        tar xzf {output.gz_file} --directory {params.kofam_dir}
        """

rule run_kofamscan:
    input:
        fasta = get_pep,
        executable = os.path.join(kofam_dir, "kofamscan-1.1.0/exec_annotation"),
        profile_dir = os.path.join(kofam_dir,"profiles"),
        ko_list = os.path.join(kofam_dir, "ko_list") 
    output:
        os.path.join(out_dir, "{sample}.txt")
    log:
        os.path.join(logs_dir, "kofamscan_results_{sample}.log")
    shadow: "shallow"
    threads: 28
    conda:
        "kofamscan-env.yml"
    shell:
        """
        {input.executable} --ko-list {input.ko_list} --profile {input.profile_dir} --cpu {threads} -f mapper -o {output} {input.fasta}
        """
        #{input.executable} --ko-list {input.ko_list} --profile {input.profile_dir} --cpu {threads} -f detail -o {output} {input.fasta}
