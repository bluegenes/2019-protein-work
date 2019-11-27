"""Snakemake wrapper for ncbi-genome-download"""

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

import os

# handle parameters
extra = snakemake.params.get("extra", "")
k = snakemake.params.get("k", "31")

accession = snakemake.params.get("accession")
acc_cmd = ""
if accession:
    acc_cmd = f" -A {accession} "

# handle outputs
genomic_out = snakemake.output.get("genomic")
protein_out = snakemake.output.get("protein")
rna_out = snakemake.output.get("rna")
cds_out = snakemake.output.get("cds")
failed_out= snakemake.params.get("failed", "ncbi_failed_downloads.txt")
ghost= snakemake.params.get("ghost_files", False)

outfiles = [genomic_out, protein_out, rna_out, cds_out]
assert (list(map(bool, outfiles)).count(True) >= 1), "please specify at least one output format by using the 'genomic', 'protein', 'rna', 'cds' keywords in the output field)"

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

dl_dir=os.path.join("genbank","*",accession)

# run ncbi-genome-download and if it works,  move output to the right spot
if genomic_out:
    tmp_genomic=os.path.join(dl_dir,"*_genomic.fna.gz")
    try:
        shell("ncbi-genome-download all -A {accession} --format 'fasta' -s genbank -p {snakemake.threads} {log}")
        shell("mv {tmp_genomic} {genomic_out}")
    except:
        shell("echo {genomic_out} >> {failed_out}")
        if ghost:
            shell("touch {genomic_out}")

if protein_out:
    tmp_protein=os.path.join(dl_dir,"*_protein.faa.gz")
    try:
        shell("ncbi-genome-download all -A {accession} --format 'protein-fasta' -s genbank -p {snakemake.threads} {log}")
        shell("mv {tmp_protein} {protein_out}")
    except:
        shell("echo {protein_out} >> {failed_out}")
        if ghost:
            shell("touch {protein_out}")

if rna_out:
    tmp_rna=os.path.join(dl_dir,"*_rna_from_genomic.fna.gz")
    try:
        shell("ncbi-genome-download all -A {accession} --format 'rna-fna' -s genbank -p {snakemake.threads} {log}")
        shell("mv {tmp_rna} {rna_out}")
    except:
        shell("echo {rna_out} >> {failed_out}")
        if ghost:
            shell("touch {rna_out}")

if cds_out:
    tmp_cds= os.path.join(dl_dir,"*_cds_from_genomic.fna.gz")
    try:
        shell("ncbi-genome-download all -A {accession} --format 'cds-fasta' -s genbank -p {snakemake.threads} {log}")
        shell("mv {tmp_cds} {cds_out}")
    except:
        shell("echo {cds_out} >> {failed_out}")
        if ghost:
            shell("touch {cds_out}")
