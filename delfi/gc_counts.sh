#!/bin/bash
#SBATCH --job-name gc_counts
#SBATCH --mem=15G
#SBATCH --time=24:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

PATH="$3:$PATH"

module load conda_R/4.3

samples_file=$1
gnom_build=$2
fragdir="../../starch"
outdir="../../gc-counts"
mkdir -p $outdir
bedops_path=$3

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample

in_file=${fragdir}/${sample}.${gnom_build}.starch
out_file=${outdir}/${sample}_gc.${gnom_build}.csv

if [ ! -f $out_file ]; then
    Rscript 02-gc_counts.r --id $in_file --outfile $out_file --bedops_path $bedops_path
	if [ $? -eq 0 ];then
		echo "Finished creating fragments!"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
	else
		echo "Error"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
		exit 100
	fi
fi ;
