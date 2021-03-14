TODO: BITP snakemake pipeline howto

0. Get reqs and connect to BITP cluster:

```
ssh {username}@clusterui.bitp.kiev.ua
```

1. Downloading hg38 and correct annotation is done by:
```
chmod +x download_prereq.sh
download_prereq.sh
```

2. Install snakemake, set up conda and environment
```
chmod +x install.sh
install.sh
```

3. Generate hg38 index
```
qsub pbs_index_hg_38_star
```

4. Activate conda env and run pipeline
```
conda activate snakemake-rnaseq

# you will need to run this serveral times until all files will be done correctly
nohup snakemake --use-conda --cluster "qsub -l walltime=5:00:00 -l nodes=1:ppn=24 -k oe -N snakemake" -k -j 4 --latency-wait 60000 --rerun-incomplete
```
Some notes for command above: 

 - nohup - prevents from stopping in case of disconnect
 - --useconda - says to snakemake to use conda in steps where environment is specified
 - --cluster - what command to run on cluster
	- -l nodes=1:ppn=24 asks for 24 core machines
 - -k - continue jobs even if one fails (it will)
 - -j number of jobs
 - --latency-wait - because steps with downloading data takes too long, wait 100 minutes before marking job as incomplete

