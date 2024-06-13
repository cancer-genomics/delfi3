#!/bin/bash
#SBATCH --job-name qc_stats 
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

module load conda_R/4.3

samples_file=$1
gnom_build=$2
bedops_path=$3
delfi_dir="../.."
qcdir="${delfi_dir}/qc"
mkdir -p $qcdir

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample

# Check index file exists and is nonzero
if [ -f $delfi_dir/starch/${sample}.${gnom_build}.starch ]; then
	Rscript 01-qcmeasures.r --id $sample --gnom_build $gnom_build --src_dir $delfi_dir --outdir $qcdir --bedops_path $bedops_path
	if [ $? -eq 0 ];then
		echo "Done!"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
	else
		echo "Error"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
		exit 100
	fi
else 
	echo -e "Error: File does not exist $delfi_dir/starch/${sample}.${gnom_build}.starch\nExiting"
	exit 0
fi


