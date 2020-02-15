"""Snakemake wrapper for sourmash compare."""

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

# handle inputs
sigs = snakemake.input.get("sigs")
traverse_dir = snakemake.input.get("traverse_sig_dir")

assert (sigs is not None or traverse_dir is not None), "please specify either a directory to traverse (traverse_sig_dir) OR a list of signatures (sigs)"

if sigs:
    input_cmd = sigs
elif traverse_dir:
    input_cmd = f" --traverse-directory {traverse_dir} "

# handle outputs
np_out = snakemake.output.get("np", "")
csv_out = snakemake.output.get("csv", "")

# handle parameters
extra = snakemake.params.get("extra", "")
k = snakemake.params.get("k", "31")

ignore_abund = snakemake.params.get("ignore_abundance", False)
abund_cmd = ""
if ignore_abund:
    abund_cmd = " --ignore-abundance "

# encodings
include = snakemake.params.get("include_encodings")
exclude = snakemake.params.get("exclude_encodings")
encoding_cmd = ""

include = [include] if isinstance(include, str) else include
exclude = [exclude] if isinstance(exclude, str) else exclude

for encoding in include:
    if encoding == "nucl":
        encoding_cmd +=" --dna "
        if 'nucl' in exclude:
            exclude.remove('nucl')
    if encoding == "protein":
        encoding_cmd += " --protein "
        if 'protein' in exclude:
            exclude.remove('protein')
    if encoding == "dayhoff":
        encoding_cmd += " --dayhoff "
        if 'dayhoff' in exclude:
            exclude.remove('dayhoff')
    if encoding == "hp":
        if 'hp' in exclude:
            exclude.remove('hp')
        encoding_cmd += " --hp "

for encoding in exclude:
    if encoding == "nucl":
        encoding_cmd +=" --no-dna "
    if encoding == "protein":
        encoding_cmd += " --no-protein "
    if encoding == "dayhoff":
        encoding_cmd += " --no-dayhoff "
    if encoding == "hp":
        encoding_cmd += " --no-hp "


log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if np_out:
    shell("sourmash compare {encoding_cmd} {abund_cmd} -k {k} -p {snakemake.threads} -o {np_out} {extra} {input_cmd} {log}")
if csv_out:
    shell("sourmash compare {encoding_cmd} {abund_cmd} -k {k} -p {snakemake.threads} --csv {csv_out} {extra} {input_cmd} {log}")
