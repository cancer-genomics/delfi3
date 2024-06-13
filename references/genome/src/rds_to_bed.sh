#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_fsize=100G
#$ -N rds_to_bed
#$ -o $JOB_NAME.o$JOB_ID.txt

module load conda_R/4.0.x
Rscript rds_to_bed.R
