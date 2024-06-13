#!/bin/bash
#SBATCH --job-name get_coverage
#SBATCH --mem=5G
#SBATCH --time=24:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%j.txt
#SBATCH -e logs/%x.e%j.txt

gnom_build=$1
bin_size=$2

bin_dir="../../bins/${bin_size}_bins"
outdir="../../features"
mkdir -p $outdir

outfile="$outdir/coverage.${bin_size}.${gnom_build}.csv"


module load conda_R/4.3
Rscript 05-getCoverage.r $gnom_build $bin_dir $outfile

