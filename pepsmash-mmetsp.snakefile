"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  mmetsp-pepsmash.snakefile --use-conda -n
"""

import os
import re
from read_samples import *

# read in elvers csv, grab sample names
samples_csv = config.get("samples", None)
if not samples_csv:
    sys.stderr.write(f"\n\tError: Please provide samples file via yml config file\n\n")
    sys.exit(-1)

samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist()
SAMPLES.remove("MMETSP0754") # pep file is empty!!!
pep_dir = config["pep_dir"]
scaled_val=config.get("scaled", "2000")
ksizes=config.get("ksizes", ["7", "11", "17"])
encodings=config.get("encodings", ["protein", "hp", "dayhoff"])

#build some output dirs
out_dir = config.get("out_dir", "mmetsp_pepsmash")
logs_dir = os.path.join(out_dir, "logs")
compute_dir = os.path.join(out_dir, "sigs")
compare_dir = os.path.join(out_dir, "compare")
plots_dir = os.path.join(out_dir, "plots")

compare_exts= [".np", ".csv"]

sigfiles = expand(os.path.join(compute_dir, "{sample}.sig"), sample= SAMPLES)
comparefiles = expand(os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare{end}"), k=ksizes, encoding=encodings, end=compare_exts)
plotfiles = expand(os.path.join(plots_dir, "mmetsp_k{k}_{encoding}_compare.np.matrix.pdf"), k=ksizes, encoding=encodings, end=compare_exts)

rule all:
    input: plotfiles
#    input: sigfiles, comparefiles, plotfiles


def get_pep(w):
#most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
#odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

# compute protein sigs
rule compute_protein:
    input: get_pep 
    output: os.path.join(compute_dir, "{sample}.sig")
    params:
        k=ksizes,
        scaled=scaled_val,
        compute_moltypes=["protein", "dayhoff", "hp"],
        input_is_protein=True,
        track_abundance=True,
    log: os.path.join(logs_dir, "{sample}_compute.log")
    benchmark: os.path.join(logs_dir, "{sample}_compute.benchmark")
    conda: "sourmash-3.1.0.yml"
    script: "sourmash-compute.wrapper.py"

# build all compare matrices: np and csv output
rule sourmash_compare:
    input: sigs=expand(os.path.join(compute_dir, "{sample}.sig"), sample= SAMPLES)
    output: 
        np=os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare.np"),
        csv=os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare.csv")
    params:
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: f"{w.k}",
    log: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_compare.log")
    benchmark: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_compare.benchmark")
    conda: "sourmash-3.1.0.yml"
    script: "sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
rule sourmash_plot:
    input: os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare.np")
    output: os.path.join(plots_dir, "mmetsp_k{k}_{encoding}_compare.np.matrix.pdf")
    params:
        plot_dir=plots_dir 
    log: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_plot.log")
    benchmark: os.path.join(logs_dir, "mmetsp_k{k}_{encoding}_plot.benchmark")
    conda: "sourmash-3.1.0.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
