#!/bin/bash
#SBATCH --job-name bin_corrected
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

PATH="$5:$PATH"

module load conda_R/4.3

samples_file=$1
gnom_build=$2
bin_size=$3
sequencer=$4
bedops_path=$5
CWD=$PWD
fragdir="../../starch"
outdir="../../bins/${bin_size}_bins"
ref_bins="../references/${gnom_build}/bins_${bin_size}.csv"
#target="../references/${gnom_build}/target-gc-marginal.csv"
mkdir -p $outdir

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample

frag_file=${fragdir}/${sample}.${gnom_build}.starch
bin_file=${outdir}/${sample}_${bin_size}.${gnom_build}.csv

# Check sequencer
sequencer_lower=$(echo $sequencer | awk '{print tolower($0)}')
if [[ $sequencer_lower =~ ^nova ]];then
	# Produced in R package PlasmaToolsNovaseq.hg19 (object name: target20)
	target="../references/${gnom_build}/target-gc-novaseq.csv"
elif [[ $sequencer_lower =~ ^hiseq ]];then
	# Produced in R package PlasmaToolsHiseq.hg19 (object name: target55)
	target="../references/${gnom_build}/target-gc-hiseq.csv"
else
    echo "Unrecognized sequencer"
    exit 100
fi

if [ ! -f $bin_file ]; then
	Rscript 03-bin_corrected.r --id $frag_file --outfile $bin_file --bins $ref_bins --target $target --bedops_path $bedops_path
	if [ $? -eq 0 ];then
		echo "Done!"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
	else
 		echo "Error"
		#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
		exit 100
	fi
fi ;

