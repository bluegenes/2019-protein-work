"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  mmetsp-pepsmash.snakefile --use-conda -n
"""

import os
import re
from read_samples import *


configfile: "config_mmetsp_pepsmash_singleton.yml"

# read in elvers csv, grab sample names
samples_csv = config.get("samples", None)
if not samples_csv:
    sys.stderr.write(f"\n\tError: Please provide samples file via yml config file\n\n")
    sys.exit(-1)

samplesDF = read_samples(samples_csv)
SAMPLES = samplesDF["sample"].tolist()
SAMPLES.remove("MMETSP0754") # pep file is empty!!!
pep_dir = config["pep_dir"]
scaled_vals=config.get("scaled", [1])
ksizes=config.get("ksizes", ["7", "11", "17"])
encodings=config.get("encodings", ["protein", "hp", "dayhoff"])

#build some output dirs
out_dir = config.get("out_dir", "mmetsp_pepsmash")
logs_dir = os.path.join(out_dir, "logs")
compute_dir = os.path.join(out_dir, "johnson_pep", "singleton_sigs")
compare_dir = os.path.join(out_dir, "johnson_pep", "singleton_compare")
plots_dir = os.path.join(out_dir, "johnson_pep", "singleton_plots")
wrappers_dir = config.get("wrappers_dir", "wrappers") 

compare_exts= [".np", ".csv"]

sigfiles = expand(os.path.join(compute_dir, "{sample}.sig"), sample= SAMPLES)
comparefiles = expand(os.path.join(compare_dir, "mmetsp_k{k}_{encoding}_compare{end}"), k=ksizes, encoding=encodings, end=compare_exts)
plotfiles = expand(os.path.join(plots_dir, "mmetsp_k{k}_{encoding}_compare_{type}.np.matrix.pdf"), k=ksizes, encoding=encodings, type=["cosine", "jaccard"], end=compare_exts)

rule all:
    #input: plotfiles
    #input: sigfiles, comparefiles, plotfiles
    #input: expand(os.path.join(out_dir, "johnson_pep", "sigs", "{sample}_jpep_scaled{scaled}_{encoding}_k{k}.sig"), sample= SAMPLES, k=[5,7,11,13,15,17,19,21,25,31,35,41], scaled= ["50", "200", "500", "1000","2000"],encoding=["protein", "dayhoff", "hp"])
    input: 
        #expand(os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_compare_jaccard.np.matrix.pdf"), sample=SAMPLES, k=ksizes, scaled=scaled_vals,encoding=encodings),
        #expand(os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"), sample=SAMPLES, k=ksizes, scaled=scaled_vals,encoding=encodings, min_count=2)
        expand(os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.hashes.txt"), sample=SAMPLES, k=ksizes, scaled=scaled_vals,encoding=encodings, min_count=2)

def get_pep(w):
#most: MMETSP0224.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep
#odd one: MMETSP0251.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep
    if w.sample == "MMETSP0251":
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.3.2.Trinity.fasta.transdecoder.pep")
    else:
        return os.path.join(pep_dir, f"{w.sample}.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep")

# compute protein sigs
#rule compute_protein:
#    input: get_pep 
#    output: os.path.join(compute_dir, "{sample}.sig")
#    params:
#        k=ksizes,
#        scaled=scaled_val,
#        compute_moltypes=["protein", "dayhoff", "hp"],
#        input_is_protein=True,
#        track_abundance=True,
#    log: os.path.join(logs_dir, "compute", "{sample}_compute.log")
#    benchmark: os.path.join(logs_dir, "{sample}_compute.benchmark")
#    conda: "sourmash-3.1.0.yml"
#    script: "sourmash-compute.wrapper.py"

rule johnson_pep_sourmash_compute:
    input: get_pep 
    output: os.path.join(compute_dir, "{sample}_jpep_scaled{scaled}_{encoding}_k{k}.sig")
    params:
        k= lambda w: w.k,
        scaled= lambda w: w.scaled,
        compute_moltypes= lambda w: w.encoding,
        input_is_protein=True,
        track_abundance=True,
        singleton=True,
    log: os.path.join(logs_dir, "sourmash", "{sample}_jpep_scaled{scaled}_{encoding}_k{k}_compute.log")
    benchmark: os.path.join(logs_dir, "sourmash", "{sample}_jpep_scaled{scaled}_{encoding}_k{k}_compute.benchmark")
    conda: "sourmash-3.2.2.yml"
    script: "sourmash-compute.wrapper.py"


rule drop_unique_hashes_jpep:
# rule modified from @taylorreiter get_greater_than_1_filt_sigs at https://github.com/taylorreiter/ibd/blob/master/Snakefile
    input: expand(os.path.join(compute_dir, "{sample}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}.sig"), sample=SAMPLES)
    output:
        hashes=os.path.join(compute_dir, "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.hashes.txt"),
        sig=os.path.join(compute_dir, "all_hashes_above_{min_count}_jpepscaled{scaled}_{encoding}_k{k}.sig"),
    log: os.path.join(logs_dir, "sourmash", "allhashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.log")
    benchmark: os.path.join(logs_dir, "sourmash", "allhashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.benchmark")
    params:
        scaled= lambda w: w.scaled,
        ksize= lambda w: w.k,
        molecule= lambda w: w.encoding,
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    script: os.path.join(wrappers_dir, "drop_unique_hashes.py")


# how does intersect work with --singleton??
rule intersect_to_drop_unique:
    input:
        keep_hashes=os.path.join(out_dir, "hashclust", "all_hashes_above_{min_count}_jpep_scaled{scaled}_{encoding}_k{k}.sig"),
        sig=os.path.join(compute_dir, "{sample}_jpep_scaled{scaled}_{encoding}_k{k}.sig")
    output:
        filt=os.path.join(compute_dir, "{sample}_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}.sig"),
        filt_renamed=os.path.join(compute_dir, "{sample}_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_renamed.sig")
    conda: os.path.join(wrappers_dir, "sourmash-3.2.2.yml")
    shell:
        """
        sourmash signature intersect -o {output.filt} --abundances-from {input.sig} -k {wildcards.k} {input.sig} {input.keep_hashes}
        sourmash signature rename -o {output.filt_renamed} -k {wildcards.k} {input.sig} {wildcards.sample}_above{wildcards.min_count}_renamed
        """

rule sourmash_compare_cosine_nounique:
    input: sigs= expand(os.path.join(compute_dir, "{sample}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output: 
        np=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"),
        csv=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.csv"),
    params:
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.log")
    benchmark: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.benchmark")
    conda: "sourmash-3.2.2.yml"
    script: "sourmash-compare.wrapper.py"

rule sourmash_compare_jaccard_nounique:
    input: sigs= expand(os.path.join(compute_dir, "{sample}_jpep_scaled{{scaled}}_{{encoding}}_k{{k}}_above{{min_count}}_renamed.sig"), sample=SAMPLES)
    output: 
        np=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"),
        csv=os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.csv"),
    params:
        ignore_abundance=True,
        include_encodings = lambda w: w.encoding,
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will exclude everything except for included encoding
        k = lambda w: w.k,
    log: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.log")
    benchmark: os.path.join(logs_dir, "compare", "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.benchmark")
    conda: "sourmash-3.2.2.yml"
    script: "sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
rule sourmash_plot_cosine_nounique:
    input: os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np"), 
    output: os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_cosine.np.matrix.pdf"),
    params:
        plot_dir=plots_dir 
    log: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.log")
    benchmark: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_cosine_plot.benchmark")
    conda: "sourmash-3.2.2.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """

rule sourmash_plot_jaccard_nounique:
    input: os.path.join(compare_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np"), 
    output: os.path.join(plots_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_compare_jaccard.np.matrix.pdf"),
    params:
        plot_dir=plots_dir 
    log: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.log")
    benchmark: os.path.join(logs_dir, "mmetsp_jpep_scaled{scaled}_{encoding}_k{k}_above{min_count}_jaccard_plot.benchmark")
    conda: "sourmash-3.1.0.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
