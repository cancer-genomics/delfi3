library(getopt)
source("PlasmaTools_Functions.R")

args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
  unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
  option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
sample <- options.args[1]
gnom_build <- options.args[2]
src_dir <- options.args[3]
outdir <- options.args[4]
bedops_path <- options.args[5]

dir.create(outdir, showWarnings = F, recursive = T)

library(data.table)
setDTthreads(1)
library(tidyverse)
library(jsonlite)
source("../references/temp_bedops.R")


jsonpath <- file.path(src_dir, "cram", "fastp", paste0(sample, ".json"))
fragpath <- file.path(src_dir, "starch", paste0(sample,".", gnom_build,".starch"))
duppath <- file.path(src_dir, "cram", "dup_metrics", paste0(sample, "_dupstat.txt"))
flagstatpath <- file.path(src_dir, "cram", "flagstat", paste0(sample, "_flagstat.",gnom_build,".txt"))
idxstatpath <- file.path(src_dir, "cram", "idxstat", paste0(sample, "_idxstat.",gnom_build,".txt"))

x <- read_json(jsonpath)

x1 <- x$summary %>%
  as_tibble() %>%
  mutate(fastp_measure = names(before_filtering)) %>%
  unnest(cols=c(before_filtering, after_filtering)) %>%
  filter(grepl("rate|gc_", fastp_measure))

x2 <- c(x1$before_filtering)
names(x2) <- c(x1$fastp_measure)

frags <- unstarch(fragpath)
frags[,w:=chromEnd-chromStart]
frags[,gc:=gc_count/w]
setnames(frags,  c("chrom","chromStart",  "chromEnd"), c("chr", "start","end"))
b <-frag.stats(frags)

flagstat <- readLines(flagstatpath)[c(1,5,7,12,14,15)]
flagstat <- as.integer(sapply(strsplit(flagstat, " "), "[[", 1))
names(flagstat) <- c("total", "duplicates", "mapped", "proper",
                     "singletons", "mate_diff_chr")


### Save QC file
stats <- cbind(sample, b, t(flagstat),  t(x2))
setnames(stats, "gc_content", "read_gc_perc")
stats[,unmapped:=total-mapped-singletons]

filename <- file.path(outdir, paste0(sample, "_qc.",gnom_build,".csv"))

fwrite(stats, filename)



