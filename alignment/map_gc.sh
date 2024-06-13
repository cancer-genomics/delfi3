#!/bin/bash
#SBATCH --job-name map_gc 
#SBATCH --mem=1G
#SBATCH --time=730:00:00
#SBATCH --mail-type=TIME_LIMIT_90
#SBATCH --partition=cancergen
#SBATCH -o logs/%x.o%A.%a.txt
#SBATCH -e logs/%x.e%A.%a.txt

PATH="$3:$PATH"

umask g+w

# Inputs
samples_file=$1
gnom_build=$2
beddir="../../bed"
starchdir="../../starch"
ref_dir="../references/genome"

mkdir -p $starchdir 

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples_file | awk '{ gsub("_WGS.*", ""); print $0 }')

outfile="$starchdir/${sample}.${gnom_build}.starch"
if [ -f $outfile ] && [ -s $outfile ];then
	# --is-starch makes sure the file is not corrupted
	is_starch=$(unstarch --is-starch $outfile)
	if [ $is_starch == "1" ];then
		echo "Starch file already exists: $outfile"
		echo "Exiting!"
		exit 0
	fi
	echo "Starch file is corrupted: $outfile"
	echo "Regenerating..."
fi

unstarch ${beddir}/${sample}.${gnom_build}.starch | \
awk -F"\t" '{frag_len=$3-$2;mapq=$4;if(frag_len<=1000 && mapq>=30){print $0}}' | \
bedmap --echo --count --delim "\t" - ${ref_dir}/cytosine_ref.${gnom_build}.starch | \
bedmap --echo --count --delim "\t" - ${ref_dir}/all_CpG_sites.${gnom_build}.starch | \
bedmap --echo --echo-map-id --delim "\t" - ${ref_dir}/cpg_domains.${gnom_build}.starch | \
awk -F"\t" 'BEGIN{OFS="\t";print "#chrom","chromStart","chromEnd","id","mapq","strand","gc_count","cpg_count","cpg_type"}{print $1,$2,$3,".",$4,$5,$6,$7,$8}' | \
starch --header --note="$gnom_build" - > $outfile

echo Done processing $sample!

#sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxRSS,NodeList
