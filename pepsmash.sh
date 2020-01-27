

snakemake -s pepsmash-mmetsp.snakefile --configfile config_mmetsp_pepsmash.yml  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.output} -e {cluster.error} -D {cluster.chdir} " --cluster-config farm_config.yml --jobs 20 --latency-wait=25 --use-conda -p
