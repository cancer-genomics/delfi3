#!/bin/bash
#SBATCH --job-name FeatureMatrix
#SBATCH --mem=5G
#SBATCH --time=96:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH -o logs/%x.o%j.txt
#SBATCH -e logs/%x.e%j.txt

gnom_build=$1
bin_size=$2

outdir="../../features"
mkdir -p $outdir

module load conda_R/4.3
Rscript 06-FeatureMatrix.r $gnom_build $bin_size $outdir
