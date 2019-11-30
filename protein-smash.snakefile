"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  protein-smash.snakefile --use-conda -n
"""

import os
import re

# set some useful vars 
genomic_dir, protein_dir, rna_dir, cds_dir = "genomic", "protein", "rna", "cds"
outbase = "smash-testing"
compare_dir = os.path.join(outbase, "compare")
plots_dir = os.path.join(outbase, "plots")

exts = {"genomic": "_genomic.fna.gz", "protein": "_protein.faa.gz", "rna": "_rna_from_genomic.fna.gz", "cds": "_cds_from_genomic.fna.gz"}
mol_dir = {"genomic": "genomic", "protein": "protein", "rna": "rna", "cds": "cds"}

compare_exts= [".np", ".csv"]

SUBSETS=["bacteroides", "denticola", "gingivalis"]
prot_moltypes=["protein"]
prot_encodings=["protein", "dayhoff", "hp"]
prot_ksizes = ["7", "11", "17"]

nucl_moltypes=["genomic", "rna", "cds"]
nucl_encodings=["nucl"]
nucl_ksizes = ["21", "31", "51"]

translate_moltypes=["rna", "cds"]
translate_ksizes = ["21","33","51"] # need to be divisible by 3. these divide to the prot_ksizes, above. need to also do 31, but might run into /3 warning/error

prot_compare_files = expand(os.path.join(compare_dir,"{subset}_k{k}_{mol}_{enc}_compare{ext}"), subset=SUBSETS, k=prot_ksizes, mol=prot_moltypes, enc=prot_encodings, ext=compare_exts)
nucl_compare_files = expand(os.path.join(compare_dir,"{subset}_k{k}_{mol}_{enc}_compare{ext}"), subset=SUBSETS, k=nucl_ksizes, mol=nucl_moltypes, enc=nucl_encodings, ext=compare_exts)
translate_compare_files = expand(os.path.join(compare_dir,"{subset}_k{k}_{mol}_{enc}_compare{ext}"), subset=SUBSETS, k=prot_ksizes, mol=translate_moltypes, enc=prot_encodings, ext=compare_exts)

prot_compare_plots = expand(os.path.join(plots_dir, "{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=prot_ksizes, mol=prot_moltypes, enc=prot_encodings, ext=compare_exts)
nucl_compare_plots = expand(os.path.join(plots_dir, "{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=nucl_ksizes, mol=nucl_moltypes, enc=nucl_encodings, ext=compare_exts)
translate_compare_plots = expand(os.path.join(plots_dir,"{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=prot_ksizes, mol=translate_moltypes, enc=prot_encodings, ext=compare_exts)

rule all:
    #input: prot_compare_files , nucl_compare_files, translate_compare_files
    input: prot_compare_plots, nucl_compare_plots, translate_compare_plots

# compute nucleotide sigs
rule compute_genomic:
    input: os.path.join(outbase, "{subset}", "genomic", "{sample}.fna.gz") 
    output: os.path.join(outbase, "{subset}", "genomic", "sigs", "{sample}.sig")
    params: 
        k=nucl_ksizes,
        scaled=2000,
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute protein sigs
rule compute_protein:
    input: os.path.join(outbase, "{subset}", "protein", "{sample}.faa.gz")
    output: os.path.join(outbase, "{subset}", "protein", "sigs", "{sample}.sig" )
    params:
        k=prot_ksizes,
        scaled=2000,
        extra=" --input-is-protein --protein --dayhoff --hp",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute nucl and protein sigs from 6-frame translation of rna, cds data
rule compute_rna:
    input: os.path.join(outbase, "{subset}", "rna", "{sample}.fna.gz") 
    output: os.path.join(outbase, "{subset}", "rna", "sigs", "{sample}.sig" )
    params:
        scaled=2000,
        k=translate_ksizes,
        extra=" --protein --dayhoff --hp",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

rule compute_cds:
    input: os.path.join(outbase, "{subset}", "cds", "{sample}.fna.gz") 
    output: os.path.join(outbase, "{subset}", "cds", "sigs", "{sample}.sig" )
    params:
        scaled=2000,
        k=translate_ksizes,
        extra=" --protein --dayhoff --hp",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

def aggregate_sigs(w):
    if (w.molecule == "protein"):
        path = os.path.join(outbase, w.subset, mol_dir[w.molecule], "{sample}.faa.gz")
    else:
        path = os.path.join(outbase, w.subset, mol_dir[w.molecule], "{sample}.fna.gz")
    sigbase= os.path.join(outbase, w.subset, mol_dir[w.molecule], "sigs","{sample}.sig")
    sigfiles = expand(sigbase, sample=glob_wildcards(path).sample)
    return sigfiles

# build all compare matrices: np and csv output
all_encodings = ["nucl", "protein", "dayhoff", "hp"]
rule sourmash_compare:
    input: sigs=aggregate_sigs
    output: 
        np=os.path.join(compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.np"),
        csv=os.path.join(compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.csv")
    params:
        include_encodings = lambda w: "{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will excude everything except for included encoding
        k = lambda w: f"{w.k}",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
rule sourmash_plot:
    input: os.path.join(compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.np")
    output: os.path.join(plots_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.pdf")
    params:
        plot_dir= os.path.join(plots_dir),
    conda: "sourmash-2.3.0.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
