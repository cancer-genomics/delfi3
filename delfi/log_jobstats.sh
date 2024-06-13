#!/bin/bash
#SBATCH --job-name=log_jobstats
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=cancergen,shared
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

outdir="log_jobstats"
mkdir -p $outdir

job_id=$1
job_name=$(sacct -j $job_id -o JobName%50 | awk 'NR%3==0{gsub(" ","",$0);print $0}' | uniq)

sacct -j ${job_id}.batch -o JobID%25,NodeList,ExitCode,MaxRSS,AveRSS,Elapsed,AveDiskRead,AveDiskWrite --delimiter "\t" --units="G" > $outdir/${job_name}.${job_id}.txt
