
import os
import pandas as pd

samplelist = "/Users/tessa/Dropbox/dib-lab/2019-burgers-shrooms/mmetsp_info/ten_haptophytes.txt"
SAMPLES = [x.strip().split('\t')[0] for x in open(samplelist, "r")]

ksizes=["17"]
encodings=["protein"]
compare_exts=[".np", ".csv", ".np.matrix.pdf"]

out_dir = "/Users/tessa/Dropbox/dib-lab/2019-burgers-shrooms/mmetsp_info/sourmash_sigs/mmetsp_pepsmash/ten_hapto"
compare_dir = os.path.join(out_dir, "compare")
plots_dir = os.path.join(out_dir, "plots")


sigfiles = expand("/Users/tessa/Dropbox/dib-lab/2019-burgers-shrooms/mmetsp_info/sourmash_sigs/mmetsp_pepsmash/sigs/{sample}.sig", sample=SAMPLES)
comparefiles = expand(os.path.join(compare_dir, "ten_hapto_k{k}_{encoding}_compare{end}"), k=ksizes, encoding=encodings, end = compare_exts)
plotfiles = expand(os.path.join(plots_dir, "ten_hapto_k{k}_{encoding}_compare.np.matrix.pdf"), k=ksizes, encoding=encodings)

rule all:
    input: plotfiles

# build compare matrices: np and csv output
rule sourmash_compare:
    input: sigs=sigfiles
    output:
        np=os.path.join(compare_dir, "ten_hapto_k{k}_{encoding}_compare.np"),
        csv=os.path.join(compare_dir, "ten_hapto_k{k}_{encoding}_compare.csv")
    params:
        include_encodings = lambda w: f"{w.encoding}",
        exclude_encodings = ["nucl", "protein", "dayhoff", "hp"], # this will excude everything except for included encoding
        k = lambda w: f"{w.k}",
    conda: "sourmash-3.1.0.yml"
    script: "sourmash-compare.wrapper.py"

# sourmash plot each compare matrix numpy output
rule sourmash_plot:
    input: os.path.join(compare_dir, "ten_hapto_k{k}_{encoding}_compare.np")
    output: os.path.join(plots_dir, "ten_hapto_k{k}_{encoding}_compare.np.matrix.pdf")
    params:
        plot_dir= os.path.join(plots_dir),
    conda: "sourmash-3.1.0.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
