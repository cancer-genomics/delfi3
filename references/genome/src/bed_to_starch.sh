#!/bin/bash
#$ -cwd
#$ -j y
#$ -R y
#$ -l mem_free=30G
#$ -l h_vmem=30G
#$ -l h_fsize=100G
#$ -N beds_to_starch
#$ -hold_jid rds_to_bed
#$ -o $JOB_NAME.o$JOB_ID.txt

sort_tmp="tmp.sort_files"
mkdir -p $sort_tmp

#sort-bed --max-mem 15G --tmpdir $sort_tmp cytosine_ref.hg19.bed | \
# Default is bzip2, but gzip extracts faster
# Currently a bug with extracting gzip archives
#starch --bzip2 --note="hg19" - > ../cytosine_ref.hg19.starch


sort-bed --max-mem 15G --tmpdir $sort_tmp all_CpG_sites.hg19.bed | \
# Default is bzip2, but gzip extracts faster
# Currently a bug with extracting gzip archives
starch --bzip2 --note="hg19" - > ../all_CpG_sites.hg19.starch


sort-bed --max-mem 15G --tmpdir $sort_tmp cpg_domains.hg19.bed | \
# Default is bzip2, but gzip extracts faster
# Currently a bug with extracting gzip archives
starch --bzip2 --note="hg19" - > ../cpg_domains.hg19.starch

