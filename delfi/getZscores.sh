#!/bin/bash
#SBATCH --job-name zscores 
#SBATCH --mem=5G
#SBATCH --time=24:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH -o logs/%x.o%j.txt
#SBATCH -e logs/%x.e%j.txt

module load conda_R/4.3

gnom_build=$1
bin_size=$2
sequencer=$3

# Check sequencer
sequencer_lower=$(echo $sequencer | awk '{print tolower($0)}')
if [[ $sequencer_lower =~ ^nova ]];then
    sequencer_standardized="novaseq"
elif [[ $sequencer_lower =~ ^hiseq ]];then
    sequencer_standardized="hiseq"
else
    echo "Unrecognized sequencer"
    exit 100
fi

bindir="../../bins/${bin_size}_bins"
armref="../references/${gnom_build}/armsummaries_${sequencer_standardized}.csv"
outdir="../../features"
mkdir -p $outdir

Rscript 04-getZscores.r $bindir $armref $outdir $bin_size $gnom_build
