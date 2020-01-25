"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  swiftortho.snakefile --use-conda --dryrun  # remove `--dryrun` to actually run
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

samples_csv = config["samples"]
samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist() 
pep_dir = config["pep_dir"]
out_dir = config.get("out_dir", ".")
pep_fasta  = os.path.join(out_dir, "pep_fasta")

if "MMETSP0251" in SAMPLES:
    SAMPLES.remove("MMETSP0251")

localrules: prep_pepdir

rule all:
    #input: os.path.join(out_dir, "dir_step3/orthologous_groups.txt")
    #input: os.path.join("SwiftOrtho", "example/run.sh")
    input: "SwiftOrtho/bin/find_hit.py"

# download program (replace with conda install in the future!?)
rule download_swiftortho_program:
    output: 
        "SwiftOrtho/install.sh"
    log: os.path.join(out_dir, "logs/download_swiftortho_program.log")
    conda: "swiftortho.yml"
    shell:
        """
        git clone https://github.com/Rinoahu/SwiftOrtho.git
        """

rule install_swiftortho:
    input: "SwiftOrtho/install.sh"
    output: "SwiftOrtho/bin/find_hit.py" 
    log: os.path.join(out_dir, "logs/swiftortho_install.log")
    conda: "swiftortho.yml"
    shell:
       """
       # 1. download portable pypy
       wget -c https://bitbucket.org/squeaky/portable-pypy/downloads/pypy-7.1.0-linux_x86_64-portable.tar.bz2

       # 2. unzip the compressed file
       tar xvf pypy-7.1.0-linux_x86_64-portable.tar.bz2

       # 3. remove the compressed file
       rm pypy-7.1.0-linux_x86_64-portable.tar.bz2

       # 4. install pypy

       mv ./pypy-7.1.0-linux_x86_64-portable/ SwiftOrtho/pypy
       
       mkdir -p SwiftOrtho/pypy/install_dir
       SwiftOrtho/pypy/bin/virtualenv-pypy -p SwiftOrtho/pypy/bin/pypy SwiftOrtho/pypy/install_dir
       SwiftOrtho/pypy/install_dir/bin/pypy -mpip install -U pip rpython
#       SwiftOrtho/pypy/install_dir/bin/pypy SwiftOrtho/bin/find_hit.py
       """

#./install_dir/bin/pypy -mpip install -U pip rpython


rule prep_pepdir:
    input: config["samplelist"]
    output: expand(os.path.join(pep_fasta, "{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta"), sample=SAMPLES), 
            os.path.join(pep_fasta, "MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep.fasta")
    params:
        pep_src=pep_dir,
        pep_dest=pep_fasta
    shell:
        """
        python copy_samplefiles.py {input} --source {params.pep_src} --destination {params.pep_dest}
        """

rule run_swiftortho: 
    input: 
        expand(os.path.join(pep_fasta, "{sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep.fasta"), sample=SAMPLES),
        os.path.join(pep_fasta, "MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep.fasta")
    output: 
        step3="dir_step3/orthologous_groups.txt",
        step4="dir_step4/orthologous_pairs.txt",
        #"dir_step3/chimeric_proteins.txt",
    log: os.path.join(out_dir, "logs/broccoli.log") 
    params:
        pep_dir=pep_fasta
    conda: "broccoli.yml"
    threads: 60
    shell: 
        """
        python Broccoli/broccoli.py -dir {params.pep_dir} -threads {threads} > {log} 2>&1
        """
        #python Broccoli/broccoli.py -dir Broccoli/example_dataset
