"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s cleanup-dl.snakefile -n 

Note: some protein, rna, and cds downloads will fail because those datasets
are unavailable for that sample. This snakefile will delete the empty files.
"""

import os
outbase = "smash-testing"

rule rm_failed_dl:
    input: os.path.join(outbase, "ncbi_failed_downloads.txt")
    output: os.path.join(outbase, "ncbi_failed_downloads.txt.rm")
    shell:
        """
        xargs rm < {input}
        mv {input} {output}
        """
