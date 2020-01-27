
#snakemake -s eggnog.snakefile --configfile config_eggnog.yml  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.output} -e {cluster.error} -D {cluster.chdir} " --cluster-config farm_config.yml --jobs 20 --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k -n

snakemake -s eggnog.snakefile --configfile config_eggnog.yml  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -D {cluster.chdir} " --cluster-config farm_config.yml --jobs 2 --latency-wait=25 --use-conda -p --shadow-prefix /scratch/ntpierce --rerun-incomplete -k

