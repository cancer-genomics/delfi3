#!/bin/bash
#SBATCH --job-name fastp 
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=730:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

FASTP=$2

umask g+w

# Inputs
samples_file=$1
fqdir="../../fastq"
outdir="../../cram"
scratchdir="../../cram/fastq_trimmed"
outdirfastp=$outdir/fastp

# Main
mkdir -p $outdir $outdir/fastp $scratchdir

sample_in=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file)
sample_out=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample_out

read1fq=$fqdir/${sample_in}_R1.fastq.gz
read2fq=$fqdir/${sample_in}_R2.fastq.gz

# Trimming with fastp
read1fqpaired=$scratchdir/${sample_out}.paired.R1.fastq.gz
read2fqpaired=$scratchdir/${sample_out}.paired.R2.fastq.gz

json=$outdirfastp/${sample_out}.json
html=$outdirfastp/${sample_out}.html

if [ -s $json ] && [ -s $html ]; then
   echo Sample $sample_out has been processed. Exiting.
   exit 0
fi

echo "Trimming $sample_out with fastp"
$FASTP -i $read1fq -I $read2fq -o $read1fqpaired -O $read2fqpaired \
    -Q -z 1 --detect_adapter_for_pe -j $json -h $html -w $SLURM_CPUS_PER_TASK

#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
