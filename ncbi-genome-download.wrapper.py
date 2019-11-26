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
prot_out = snakemake.output.get("protein")
rna_out = snakemake.output.get("rna")
cds_out = snakemake.output.get("cds")
outfiles = [genomic_out, prot_out, rna_out, cds_out]
assert (list(map(bool, outfiles)).count(True) >= 1), "please specify at least one output format by using the 'genomic', 'protein', 'rna', 'cds' keywords in the output field)"

format_list = []
if genomic_out:
    format_list.append("fasta")
if prot_out:
    format_list.append("protein-fasta")
if rna_out:
    format_list.append("rna-fna")
if cds_out:
    format_list.append("cds-fasta")

formats = ",".join(format_list)
format_cmd= f" --format {formats} "

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

dl_dir=os.path.join("genbank","*",accession)

# run ncbi-genome-download
shell("ncbi-genome-download all -A {accession} {format_cmd} -s genbank -p {snakemake.threads} {log}")

# move output to the right spot
if genomic_out:
    tmp_genomic=os.path.join(dl_dir,"*_genomic.fna.gz")
    shell("mv {tmp_genomic} {genomic_out}")
else:
    shell("touch {genomic_out}")
    shell("echo {genomic_out} >> ncbi_failed_downloads.txt")

if protein_out:
    tmp_protein=os.path.join(dl_dir,"*_protein.faa.gz")
    shell("mv {tmp_protein} {protein_out}")
else:
    shell("touch {protein_out}")
    shell("echo {protein_out} >> ncbi_failed_downloads.txt")

if rna_out:
    tmp_rna=os.path.join(dl_dir,"*_rna_from_genomic.fna.gz")
    shell("mv {tmp_rna} {rna_out}")
else:
    shell("touch {rna_out}")
    shell("echo {rna_out} >> ncbi_failed_downloads.txt")

if cds:
    tmp_cds= os.path.join(dl_dir,"*_cds_from_genomic.fna.gz")
    shell("mv {tmp_cds} {cds_out}")
else:
    shell("touch {cds_out}")
    shell("echo {cds_out} >> ncbi_failed_downloads.txt")
