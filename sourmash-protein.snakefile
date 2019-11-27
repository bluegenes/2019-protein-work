"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import re

# set some useful vars 
genomic_dir, protein_dir, rna_dir, cds_dir = "genomic", "protein", "rna", "cds"
outbase = "smash_testing"
compare_dir = "compare"
plots_dir = "plots"

genomic_ext = "_genomic.fna.gz"
protein_ext = "_protein.faa.gz"
rna_ext = "_rna_from_genomic.fna.gz"
cds_ext = "_cds_from_genomic.fna.gz"
    

# Infer sample IDs from previously downloaded files (dl-datasets.snakefile), rather than getting them from the csvs
#GENOMIC = glob_wildcards(os.path.join(outbase, "*", genomic_dir, "{sample}.fna.gz" )
#PROTEIN = glob_wildcards(os.path.join(outbase, "*", protein_dir, "{sample}.faa.gz")
#RNA = glob_wildcards(os.path.join(outbase, "*", rna_dir, "{sample}.fna.gz")
#CDS = glob_wildcards(os.path.join(outbase, "*", cds_dir, "{sample}.fna.gz")

compare_extensions= [".np", ".csv"]

SUBSETS=["bacteroides", "denticola", "gingivalis"]
prot_moltypes=["protein", "rna", "cds"]
prot_encodings=["protein", "dayhoff", "hp"]
prot_ksizes = ["7", "11", "17"]

nucl_moltypes=["genomic", "rna", "cds"]
nucl_encodings=["nucl"]
nucl_ksizes = ["21", "31", "51"]

translate_moltypes=["rna", "cds"]
translate_ksizes = ["21", "33", "51"] # need to be divisible by 3

prot_compare_files = expand(os.path.join(outbase,compare_dir,"{subset}_k{k}_{mol}_{enc}_compare{ext}"), subset=SUBSETS, k=prot_ksizes, mol=prot_moltypes, enc=prot_encodings, ext=compare_extensions)
nucl_compare_files = expand(os.path.join(outbase,compare_dir,"{subset}_k{k}_{mol}_{enc}_compare.np"), subset=SUBSETS, k=prot_ksizes, mol=nucl_moltypes, enc=nucl_encodings, ext=compare_extensions)
translate_compare_files = expand(os.path.join(outbase,compare_dir,"{subset}_k{k}_{mol}_{enc}_compare.np"), subset=SUBSETS, k=translate_ksizes, mol=translate_moltypes, enc=prot_encodings, ext=compare_extensions)

prot_compare_plots = expand(os.path.join(outbase,compare_dir, plots_dir, "{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=prot_ksizes, mol=prot_moltypes, enc=prot_encodings, ext=compare_extensions)
nucl_compare_plots = expand(os.path.join(outbase,compare_dir,plots_dir, "{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=prot_ksizes, mol=nucl_moltypes, enc=nucl_encodings, ext=compare_extensions)
translate_compare_plots = expand(os.path.join(outbase,compare_dir,plots_dir,"{subset}_k{k}_{mol}_{enc}_compare.pdf"), subset=SUBSETS, k=translate_ksizes, mol=translate_moltypes, enc=prot_encodings, ext=compare_extensions)

rule all:
    #input: prot_compare_files, nucl_compare_files, translate_compare_files
    input: prot_compare_plots, nucl_compare_plots, translate_compare_plots
    #expand(os.path.join(outbase, compare_dir, "{subset}_k{k}_genomic_compare.csv"), subset=SUBSETS, k=nucl_ksizes)

#subworkflow download_datasets:
#    snakefile:
#        "dl-ncbi.snakefile"

# compute nucleotide sigs
rule compute_genomic:
    input: os.path.join(outbase, "{subset}", genomic_dir, "{sample}" + genomic_ext)
    output: os.path.join(outbase, "{subset}", genomic_dir, "sigs", "{sample}_k{k}_scaled{scaled}.sig" )
    params:
        scaled=2000,
        k=nucl_ksizes,
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute rna sigs
rule compute_rna:
    input: os.path.join(outbase, "{subset}", rna_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{subset}", rna_dir, "sigs", "{sample}_k{k}_scaled{scaled}.sig" )
    params:
        scaled=2000,
        k=nucl_ksizes,
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

rule compute_cds:
    input: os.path.join(outbase, "{subset}", cds_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{subset}", cds_dir, "sigs", "{sample}_k{k}_scaled{scaled}.sig" )
    params:
        scaled=2000,
        k=nucl_ksizes,
        extra=" --protein ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# 6 frame translation --> signatures
### why are we using rna_from_genomic over cds .. check, bc cds = only coding = cleaner!?
rule compute_translated_rna:
    input: os.path.join(outbase, "{subset}", rna_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{subset}", rna_dir, "sigs", "{sample}_k{k}_scaled{scaled}_translated.sig" )
    params:
        scaled=2000,
        k=translate_ksizes,
        extra=" --protein ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

rule compute_translated_cds:
    input: os.path.join(outbase, "{subset}", cds_dir, "{sample}" + rna_ext)
    output: os.path.join(outbase, "{subset}", cds_dir, "sigs", "{sample}_k{k}_scaled{scaled}_translated.sig" )
    params:
        scaled=2000,
        k=translate_ksizes,
        extra=" --protein ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute protein sigs
rule compute_prot:
    input: os.path.join(outbase, "{subset}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{subset}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_protein.sig" )
    params:
        scaled=2000,
        k=prot_ksizes,
        extra=" --input-is-protein --protein ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute dayhoff sigs
rule compute_dayhoff:
    input: os.path.join(outbase, "{subset}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{subset}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_dayhoff.sig" )
    params:
        scaled=2000,
        k=prot_ksizes,
        extra=" --input-is-protein --dayhoff ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

# compute hp sigs
rule compute_hp:
    input: os.path.join(outbase, "{subset}", protein_dir, "{sample}" + protein_ext)
    output: os.path.join(outbase, "{subset}", protein_dir, "sigs", "{sample}_k{k}_scaled{scaled}_hp.sig" )
    params:
        scaled=2000,
        k=prot_ksizes,
        extra=" --input-is-protein --hp ",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"

def aggregate_sigs(w):
    if (w.molecule == "protein"):
        samples = glob_wildcards(os.path.join(outbase, w.subset, f"{w.molecule}_dir", "*.faa.gz"))
    else:
        samples = glob_wildcards(os.path.join(outbase, w.subset, f"{w.molecule}_dir", "*.fna.gz"))
    sigfiles = expand(os.path.join(outbase, w.subset, f"{w.molecule}_dir", "sigs", f"{{sample}}_k{w.k}_scaled2000_{w.encoding}.sig"), sample=samples)
    return sigfiles

# build all compare matrices: np and csv output
rule sourmash_compare:
    input: aggregate_sigs
    output: 
        np=os.path.join(outbase, compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.np"),
        csv=os.path.join(outbase, compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.csv")
    params:
        include_encodings = lambda w: f"{w.encoding}",
        k = lambda w: f"{w.k}",
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
rule sourmash_plot:
    input: os.path.join(outbase, compare_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.np")
    output: os.path.join(outbase, compare_dir, plots_dir, "{subset}_k{k}_{molecule}_{encoding}_compare.pdf")
    params:
        plot_dir= os.path.join(outbase, compare_dir, plots_dir),
    conda: "sourmash-2.3.0.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """

#################

# run regex on this pandas series instead:
#sample_names = genomeInfo.iloc[:,2].str.split('/').str[-1]
#sample_info = sample_names.str.extract(r"([A-Z]+)_(\d{3})(\d{3})(\d{3})")
#sample_info.columns = ["alpha", "first", "second", "third"]
#sample_info['name'] = sample_names 
        

        #k=prot_ksizes,
        #subsets=SUBSETS,
        #molecule=["protein", "rna", "cds"], 
        #encoding=["protein", "dayhoff", "hp"], # do this in the desired outputs, rather than here!

#shell:
#        """
#        #sourmash compare -k {wildcards.k} -o {output.np} {input} 
#        #sourmash compare -k {wildcards.k} ---csv {output.csv} {input} 
#        """

