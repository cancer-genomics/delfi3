#!/bin/bash
#SBATCH --job-name align
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=730:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

module load bowtie/2.5.1
module load samtools/1.18 

umask g+w
set -euo pipefail

# Inputs
samples_file=$1
gnom_build=$2
sequencer=$3
fqdir="../../fastq"
outdir="../../cram"
scratchdir="tmp"
fqtrimmeddir="../../cram/fastq_trimmed"
## align with alternates?
refgenome_basename="../references/${gnom_build}/${gnom_build}"

# Check sequencer
sequencer_lower=$(echo $sequencer | awk '{print tolower($0)}')
if [[ $sequencer_lower =~ ^nova ]];then
    dist=2500
elif [[ $sequencer_lower =~ ^hiseq ]];then
    dist=100
else
    echo "Unrecognized sequencer"
    exit 100
fi

# Main
mkdir -p $outdir $outdir/dup_metrics $outdir/flagstat \
    $outdir/idxstat $outdir/bt2stats $scratchdir

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample
echo $sequencer

# Check index file exists and is nonzero
if [ -f $outdir/${sample}.${gnom_build}.cram.crai ] && [ -s $outdir/${sample}.${gnom_build}.cram.crai ]; then
   # Quickcheck makes sure cram file is not corrupted
   samtools quickcheck $outdir/${sample}.${gnom_build}.cram
   if [ $? == "0" ];then
   	echo Sample ${sample} has been processed. Exiting.
   	exit 0
   fi
fi

read1fqpaired=$fqtrimmeddir/${sample}.paired.R1.fastq.gz
read2fqpaired=$fqtrimmeddir/${sample}.paired.R2.fastq.gz

# Check if file is empty
if [ ! -s $read1fqpaired ] || [ ! -s $read2fqpaired ];then
	echo "Error: Fastq file is empty. Exiting."
	exit 100
fi

echo "Aligning ${sample} with bowtie2"
let nthreads=$SLURM_CPUS_PER_TASK-1 # samtools -@ argument equals 1 - # cores requested
bowtie2 --seed 42 --very-fast --end-to-end --threads $nthreads -x $refgenome_basename \
        -1 $read1fqpaired -2 $read2fqpaired \
         2> $outdir/bt2stats/${sample}_bt2.${gnom_build}.txt |
    samtools fixmate -u -@$nthreads -m - - | \
    samtools sort -u -@$nthreads -T $scratchdir - | \
    samtools markdup -@$nthreads -d $dist \
        -f $outdir/dup_metrics/${sample}_dupstat.txt \
        --reference ${refgenome_basename}.fa \
        - $outdir/${sample}.${gnom_build}.cram

samtools index -@$SLURM_CPUS_PER_TASK $outdir/${sample}.${gnom_build}.cram
samtools flagstat $outdir/${sample}.${gnom_build}.cram > $outdir/flagstat/${sample}_flagstat.${gnom_build}.txt
samtools idxstats $outdir/${sample}.${gnom_build}.cram > $outdir/idxstat/${sample}_idxstat.${gnom_build}.txt

rm $read1fqpaired
rm $read2fqpaired

#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
