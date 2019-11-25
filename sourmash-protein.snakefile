"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import re
import pandas as pd

# testing with hardcoded csvs
csv= "../2018-test_datasets/gingivalis.csv"
#csv= "../2018-test_datasets/"
#csv= "../2018-test_datasets/gingivalis.csv"

genomeInfo = pd.read_csv(csv)
basedir = csv.rsplit('.', 1)[0]#assuming good csv naming
sample_names = genomeInfo.iloc[:,2].str.split('/').str[-1]

# run regex on this pandas series instead:
sample_info = sample_names.str.extract(r"([A-Z]+)_(\d{3})(\d{3})(\d{3})")
sample_info.columns = ["alpha", "first", "second", "third"]

genbank_url = 'https://ftp.ncbi.nih.gov/genomes/all'

#alpha,first,second,third = re.match("([A-Z]+)_(\d{3})(\d{3})(\d{3})", name).groups()

# use sample names for this bit
folder_name = name.split('_genomic.fna.gz')[0]
genbank_path =os.path.join(genbank_url,alpha,first,second,third,folder_name)
genome_url = os.path.join(genbank_path,name)


    if protein:
        protein_name = folder_name + '_protein.faa.gz'
        protein_url = os.path.join(genbank_path,protein_name)
        if not quiet:
            print('protein: ' + protein_url)
        prot_name = name.split('_genomic.fna.gz')[0] + '_protein.faa.gz'
        protein_outF = os.path.join(protein_outD, prot_name)
        if not os.path.exists(protein_outD):
            os.makedirs(protein_outD, exist_ok=True)
        get_genbank_file(protein_url, protein_outF)

    if rna:
        rna_name = folder_name + '_rna_from_genomic.fna.gz'
        #rna_name = folder_name + '_cds_from_genomic.fna.gz'
        rna_url = os.path.join(genbank_path,rna_name)
        if not quiet:
            print('rna: ' + rna_url)
        rna_name = name.split('_genomic.fna.gz')[0] + '_rna_from_genomic.fna.gz'
        rna_outF = os.path.join(rna_outD, rna_name)
        if not os.path.exists(rna_outD):
            os.makedirs(rna_outD, exist_ok=True)
        #outP = outF.split('_genomic.fna.gz')[0] + '_cds_from_genomic.fna.gz'
        get_genbank_file(rna_url, rna_outF)


# download rna, protein, genomic sequences
rule get_test_datasets:
    output:
        data_dir = directory("eggnog_data")
    conda: "sourmash-env.yml"
    run:
        """
        """

# compute nucleotide sigs
rule compute_nucl:
    input:
        fasta="test/p53.fa"
        data_dir=directory("eggnog_data")
    output: "p53_maNOG"
    params:
        scaled=2000,
        k= [21,31,51],
    conda: "sourmash-2.3.0.yml"
    script: "sourmash-compute.wrapper.py"



# compute protein sigs

# compute dayhoff sigs

# compute hp sigs

# build all compare matrices: csv output
# build all compare matrices: numpy output

# sourmash plot each compare matrix numpy output



