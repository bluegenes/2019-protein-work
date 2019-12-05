"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  rnasamba.snakefile --use-conda --dryrun  # remove `--dryrun` to actually run
"""

import os
import re
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

MODELS=["full_length_weights", "partial_length_weights"]
SAMPLES=["test-transcript"]

outbase = "rnasamba_work"
logs_dir = os.path.join(outbase, "logs")
#data_dir = os.path.join(outbase, "data")
data_dir = "/Users/tessa/Dropbox/dib-lab/2019-protein-work/test-data" 
models_dir = os.path.join(outbase, "models")
results_dir = os.path.join(outbase, "results")

downloaded_models = expand(os.path.join(models_dir, "{model}.hdf5"), model=MODELS)

classified_samples = expand(os.path.join(results_dir, "{sample}", "{sample}_x_{model}_classification.tsv"), sample=SAMPLES, model=MODELS)

rule all:
    #input: downloaded_models
    input: classified_samples 

rule download_pretrained_models:
    input: lambda w: HTTP.remote(f"https://raw.githubusercontent.com/apcamargo/RNAsamba/master/data/{w.model}.hdf5", static=True, keep_local=True, allow_redirects=True) 
    output: os.path.join(models_dir, "{model}.hdf5")
    log: os.path.join(logs_dir, "dl_trained_models_{model}.log")
    shell: "mv {input} {output} 2> {log}"
    

rule rnasamba_classify:
    input: 
        model=os.path.join(models_dir, "{model}.hdf5"),
        seqs=os.path.join(data_dir, "{sample}.fa"),
    output:
        predicted_proteins= os.path.join(results_dir, "{sample}", "{sample}_x_{model}_predicted_proteins.fa"),
        classification= os.path.join(results_dir, "{sample}", "{sample}_x_{model}_classification.tsv")
    params:
        verbose=1, # 0=silent, 1= current step
    log: os.path.join(logs_dir, "classify_{sample}_x_{model}.log")
    conda: "rnasamba-env.yml"
    shell:
        """
        rnasamba classify -p {output.predicted_proteins} {output.classification} {input.seqs} {input.model} > {log} 2>&1
        """

#rule rnasamba_train:
#    input:
#    output:
#    params:
#    conda: "rnasamba-env.yml"
#    shell:
#        """
#        rnasamba train -v 2 mouse_model.hdf5 gencode.vM21.pc_transcripts.fa gencode.vM21.lncRNA_transcripts.fa
#        """
#
