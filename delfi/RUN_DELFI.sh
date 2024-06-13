#!/bin/bash
#SBATCH --mem=1G
#SBATCH -o %x.o%j
#SBATCH -e %x.o%j

set -eu

bedops_path="/dcs04/scharpf/data/nvulpesc/tools/bedops-v2.4.41"
fastp_path=/users/scristia/fastp
gnom_build="hg19"
sequencer="NovaSeq" #NovaSeq or HiSeq
bin_size="5mb"
samples_file="samples.txt"
n_samples=$(wc -l $samples_file | sed 's/ .*//')

mkdir -p logs
fastp_id=$(sbatch --parsable --array=1-${n_samples}%10 ../alignment/fastp.sh $samples_file $fastp_path)
align_id=$(sbatch --parsable --array=1-${n_samples}%5 -d aftercorr:$fastp_id ../alignment/align.sh $samples_file $gnom_build $sequencer)
post_align_id=$(sbatch --parsable --array=1-${n_samples}%20 -d aftercorr:$align_id ../alignment/post_alignment.sh $samples_file $gnom_build $bedops_path)
map_gc_id=$(sbatch --parsable  --array=1-${n_samples}%20 -d aftercorr:$post_align_id ../alignment/map_gc.sh $samples_file $gnom_build $bedops_path)
qc=$(sbatch --parsable --array=1-${n_samples}%5 -d aftercorr:$map_gc_id qcmeasures.sh $samples_file $gnom_build $bedops_path)
gc_counts_id=$(sbatch --parsable --array=1-${n_samples}%5 -d aftercorr:$map_gc_id gc_counts.sh $samples_file $gnom_build $bedops_path)
bin_corrected_id=$(sbatch --parsable --array=1-${n_samples}%5 -d aftercorr:$gc_counts_id bin_corrected.sh $samples_file $gnom_build $bin_size $sequencer $bedops_path)
zscore_id=$(sbatch --parsable -d afterok:$bin_corrected_id getZscores.sh $gnom_build $bin_size $sequencer)
coverage_id=$(sbatch --parsable -d afterok:$bin_corrected_id getCoverage.sh $gnom_build $bin_size)
feature_matrix_id=$(sbatch --parsable -d afterok:$coverage_id,$zscore_id FeatureMatrix.sh $gnom_build $bin_size)

## Log memory usage for all jobs
job_id_array=($fastp_id $align_id $post_align_id $map_gc_id $qc $gc_counts_id $bin_corrected_id $zscore_id $coverage_id $feature_matrix_id)
for job_id in ${job_id_array[@]};do
	sbatch --quiet -d afterany:$job_id log_jobstats.sh $job_id
done

echo "All scripts submitted successfully!!"
