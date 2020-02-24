
#snakemake -s pepsmash-mmetsp.snakefile --configfile config_mmetsp_pepsmash.yml  --cluster "sbatch -t {cluster.time} -p {cluster.partition} -N {cluster.nodes} --mem {cluster.mem} -c {cluster.cpus_per_task} -o {cluster.stdout} -e {cluster.stderr} -J {cluster.jobname} -D {cluster.chdir} " --cluster-config farm_config.yml --jobs 5 --latency-wait=25 --use-conda -p --rerun-incomplete

#snakemake -s pepsmash-mmetsp.snakefile --configfile config_mmetsp_pepsmash.yml --cores 1 --latency-wait=25 --use-conda -p -n

snakemake -s pepsmash-mmetsp.snakefile --configfile config_mmetsp_pepsmash.yml --cores 1 -p -n

