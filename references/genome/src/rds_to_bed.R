
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

refdir <- "/dcl01/scharpf/data/data-warehouse/references/genome"
#cytosines <- readRDS(file.path(refdir, "cytosine_ref.rds"))
#export.bed(cytosines,con='cytosine_ref.hg19.bed')

all_CpG_sites <- readRDS(file.path(refdir, "all_CpG_sites.rds"))
export.bed(all_CpG_sites,con='all_CpG_sites.hg19.bed')

annotations <- readRDS(file.path(refdir, "cpg_domains.rds"))
export.bed(annotations,con='cpg_domains.hg19.bed')