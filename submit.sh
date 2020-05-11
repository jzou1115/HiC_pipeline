#main process run interactively
#this works
snakemake -j 100 --cluster-config config/cluster.json --cluster "qsub -N {cluster.name} -cwd -V -o {cluster.output} -e {cluster.error} -l h_data={cluster.memory},highp,h_rt={cluster.time} -m {cluster.mail}"



#submit snakemake to cluster and allow snakemake to submit additional jobs
#these do not work, unclear why
#qsub -N PIPE -cwd -j yes -b y python snakemake --cluster "ssh jzou1115@hoffman2.idre.ucla.edu 'qsub -N pipe_task -j yes -cwd -S /bin/sh ' " -j
#qsub -N PIPE -cwd -V -j yes -b y snakemake --cluster "ssh jzou1115@hoffman2.idre.ucla.edu 'qsub -N pipe_task -j yes -cwd -S /bin/sh ' " -j
