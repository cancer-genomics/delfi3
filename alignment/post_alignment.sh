#!/bin/bash
#SBATCH --job-name post_alignment 
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
#SBATCH --time=730:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

## the following line is needed for recent bedtools version
PATH="$3:$PATH"
module load bedtools/2.31.0 
module load samtools/1.18

umask g+w

# Inputs
samples_file=$1
gnom_build=$2
fqdir="../../fastq"
cramdir="../../cram"
beddir="../../bed"
scratchdir="tmp"
refgenome_basename="../references/${gnom_build}/${gnom_build}"

mkdir -p $beddir $scratchdir

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')
echo $sample

outfile=$beddir/${sample}.${gnom_build}.starch

## Also check of bedfile nonempty?
if [ -f $outfile ]; then
   echo Sample $sample has been processed. Exiting.
   exit 0
fi

let nthreads=$SLURM_CPUS_PER_TASK-1 # samtools -@ argument equals 1 - # cores requested
samtools view -@ $nthreads -h -T ${refgenome_basename}.fa $cramdir/${sample}.${gnom_build}.cram | \
    samtools markdup -@$nthreads -r - - | \
    samtools sort -@$nthreads -T $scratchdir -n - | \
    samtools fixmate -u -@$nthreads - - | \
    bamToBed -bedpe -i - | \
    awk '($2 != -1) && ($5 != -1) && ($1 == $4)' | \
    cut -f 1,2,6,8,9 | \
    sort-bed --tmpdir $scratchdir --max-mem 5G - | \
    starch - > $outfile

## Can pipe the above into bedtools nuc to get nucleotide frequencies for each
## base fragment. Takes a few hours and will be much faster if we use
## bedtools intersect -c compared to a sorted bed file of every G/C position.

#bedtools nuc -fi ${refgenome_basename}.fa -bed - > $beddir/${sample}.${gnom_build}.bed

echo Done processing $sample!
#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
