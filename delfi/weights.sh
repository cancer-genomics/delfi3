#!/bin/bash
#SBATCH --mem=5G
#SBATCH --time=96:00:00
#SBATCH --mail-type=TIME_LIMIT_90

module load conda_R/4.0.x
if [ ! -f "target.csv" ]; then
    Rscript --vanilla weights.r ../gc-counts target.csv
fi
