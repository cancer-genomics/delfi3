## Reweight each fragment by mapping to target distribution based on GC
#
library(getopt)

source("PlasmaTools_Functions.R")
load("../references/hg19/filters.hg19.rda")

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
fragfile <- options.args[1]
outfile <- options.args[2]
binfile <- options.args[3]
targetfile <- options.args[4]
bedops_path <- options.args[5]
source("../references/temp_bedops.R")

if(file.exists(outfile)){
  print(paste0("File exists:",outfile,"\nExiting!"))
  q("no")
}

sample <- gsub(".hg\\d\\d.starch|.CHM13.starch", "", basename(fragfile))
print(sample)

library(data.table)
suppressMessages(library(GenomicRanges))

setDTthreads(1) 

target <- fread(targetfile)
bins <- fread(binfile)
bins <- bins[map >= 0.90 & gc >= 0.30]

fragments <- unstarch(fragfile)
setnames(fragments,  c("chrom","chromStart",  "chromEnd"), c("chr", "start","end"))
fragments <- fragments[chr%in%paste0("chr", c(1:22, "X"))]
fragments <- fragments[,w:=end-start] ## Added since we now do not save width in starch files
fragments <- fragments[,width:=w] ## REDUNDANT but both column names are required in this script
fragments <- fragments[,gc:=gc_count/w] ## Added since we now do not save width in starch files
#
### Filter blacklisted regions

filters <- as.data.table(filters.hg19)
setnames(filters, "seqnames",  "chr")
setkey(filters, chr, start, end)
fragdt <- foverlaps(fragments, filters, type="any")
fragdt <- fragdt[is.na(start)][,`:=`(start=NULL, end=NULL, width=NULL,
                                     strand=NULL, name=NULL, score=NULL)]

setnames(fragdt, c("i.start", "i.end", "i.width",  "i.strand"),
         c("start", "end", "width", "strand"))

fragdt <- gcCorrectTarget(fragdt, target, bychr=TRUE)


bins <- bins[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]
fragdt <- fragdt[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]

setnames(fragdt, "gc", "fraggc")
setkey(fragdt, chr, start, end)
setkey(bins, chr, start, end)

bins[,bin:=1:.N]
bins2 <- binFrags(fragdt, bins)
bins2[,`:=`(ratio=short/(short+long), ratio.cor=short.cor/(short.cor+long.cor))]

bins2[,id:=gsub("_downsamp", "", sample)]
setcolorder(bins2, c("id", "chr", "start", "end", "bin"))

filename <- file.path(outfile)

fwrite(bins2, filename)
q('no')

