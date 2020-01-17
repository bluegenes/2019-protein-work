#!/bin/bash -login
#SBATCH -D /home/ntpierce/2019-protein-work
#SBATCH -J brocc 
#SBATCH -p bmh
#SBATCH -t 3-0:00:00
#SBATCH -N 1
#SBATCH --output /home/ntpierce/2019-protein-work/haptophyta_broccoli/logs/sbatch-brocc-%j.out
#SBATCH --error /home/ntpierce/2019-protein-work/haptophyta_broccoli/logs/sbatch-brocc-%j.err

# activate conda in general
source /home/ntpierce/.bashrc # if you have the conda init setting

# activate a specific conda environment, if you so choose
conda activate snakemake 

# go to a particular directory
cd /home/ntpierce/2019-protein-work

# make things fail on errors
set -o nounset
set -o errexit
set -x

### run your commands here!
snakemake -s broccoli.snakefile --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -n {cluster.ntasks} -o {cluster.output} -e {cluster.error} -D {cluster.chdir} " --cluster-config farm_config.yml --jobs 2 --latency-wait=15 --use-conda

