"""
Author: N. Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s eggnog.snakefile --use-conda #--cores 26 --cluster "sbatch -t 10:00:00 -N 11 -n 26 -p bmm --mem=60gb" --jobs 5 -k --rerun-incomplete 
"""

import os
import pandas as pd

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

db_version="5.0.0"

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
samples_csv = config.get("samples", None)
if not samples_csv:
    sys.stderr.write(f"\n\tError: Please provide samples file via yml config\n\n")
    sys.exit(-1)

samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist()
pep_dir = config["pep_dir"]
out_dir = config.get("out_dir", "mmetsp_eggnog")
database_dir = config.get("db_dir", "eggnog_data")

rule anNOGtate:
    input: 
         #expand(os.path.join(out_dir, "{sample}.emapper.annotations"), sample=SAMPLES), 
         expand(os.path.join(out_dir, "{sample}.emapper.seed_orthologs"), sample=SAMPLES)

rule get_eggnog_dbs:
    output:
        os.path.join(database_dir, "eggnog.db"),
        os.path.join(database_dir, "eggnog_proteins.dmnd"),
    conda: "eggnog.yml"
    params:
        data_dir = database_dir
    shell:
        """
        mkdir -p {params.data_dir}
        download_eggnog_data.py -y -f --data_dir {params.data_dir}
        """
# grab a gz version so I can more easily read it to local for I/O intensive step
#def download_eggnog_db_gz:
#    input: lambda wildcards: HTTP.remote(f"http://eggnogdb.embl.de/download/emapperdb-{db_version}s/eggnog.db.gz", static=True, keep_local=True, allow_redirects=True)
#    output: os.path.join(database_dir, "eggnog.db.gz")
#    shell:
#        """
#        mv {input} {output}
#        """

def get_pep(w):
#most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
#odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

## to optimize on cluster, run eggnog mapper steps separately
# https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
# https://github.com/eggnogdb/eggnog-mapper/issues/80

# this part is cpu intensive
rule run_eggnog_mapper_dmnd:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        nog_dmnd=os.path.join(database_dir, "eggnog_proteins.dmnd"),
        pep=get_pep
    output:
        os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    params:
        mode="diamond",
        data_dir=directory("eggnog_data"),
        out_dir=out_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_dmnd.benchmark")
    threads: 20
    shadow: "shallow" ## this means tmpdir will be on local scratch (using --shadow-prefix) --> faster? 
    conda: "eggnog.yml"
    shell:
        """
        emapper.py --no_annot --no_file_comments -i {input.pep} --output {params.out_dir}/{wildcards.sample} -m {params.mode} --cpu {threads} --data_dir {params.data_dir} > {log} 2>&1 
        """

# this part is I/O intensive --> USE LOCAL tmpdir! (use shadow, with --shadow-prefix /$SCRATCH/ntpierce)
rule run_eggnog_mapper_annotate:
    input:
        nog_db=os.path.join(database_dir, "eggnog.db"),
        #gz_nog_db=os.path.join(database_dir, "eggnog.db.gz"),
        #nog_dmnd="eggnog_data/eggnog_proteins.dmnd",
        seed_orthologs=os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
    output:
        os.path.join(out_dir, "{sample}.emapper.annotations"),
    params:
        mode="diamond",
        #data_dir=directory("eggnog_data")
        data_dir=database_dir
    log: os.path.join(out_dir, "logs", "{sample}_emapper_annot.log")
    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper_annot.benchmark")
    threads: 1
    conda: "eggnog.yml"
    shadow: "shallow" ## here, it helps to have local eggnog db, write locally, then copy out to final location
    shell:
        """
        emapper.py --annotate_hits_table {input.seed_orthologs} --no_file_comments -o {output} --data_dir {params.data_dir} --cpu {threads} > {log} 2>&1
        """

## if running on single node, can run all at once
#rule run_eggnog_mapper:
#    input: 
#        nog_db="eggnog_data/eggnog.db",
#        nog_dmnd="eggnog_data/eggnog_proteins.dmnd",
#        pep=get_pep
#    output:
#        os.path.join(out_dir, "{sample}.emapper.annotations"),
#        os.path.join(out_dir, "{sample}.emapper.seed_orthologs")
#    params:
#        mode="diamond",
#        data_dir=directory("eggnog_data")
#    log: os.path.join(out_dir, "logs", "{sample}_emapper.log")
#    benchmark: os.path.join(out_dir, "logs", "{sample}_emapper.benchmark")
#    threads: 10
#    conda: "eggnog.yml"
#    shell:
#        """
#        emapper.py -i {input.pep} --output {wildcards.sample} -m {params.mode} --data_dir {params.data_dir} --cpu {threads} > {log} 2>&1
#        """
#    #script: "eggnog.wrapper.py"


## don't need to do this for txomes - anything < 1million contigs/file = small enough (apparently)!
#rule split_fasta_into_chunks:
#    input: "input_file.faa""
#    output:
#    params:
#        num_lines_per_chunk=2000000 
#    conda: "eggnog.yml"
#    shell:
#        """
#        split -l {params.num_lines_per_chunk} -a 3 -d {input} input_file.chunk_
#        """
#
