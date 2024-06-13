suppressMessages(library(getopt))
suppressMessages(library(data.table))


data.table::setDTthreads(1)

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
out_file <- options.args[2]
bedops_path <- options.args[3]
source("../references/temp_bedops.R")

sample <- gsub(".hg\\d\\d.starch|.CHM13.starch", "", basename(fragfile))

DT <- unstarch(fragfile)
DT[,id:=sample]
DT[,w:=chromEnd-chromStart]
DT[,gc:=gc_count/w]

DT <- DT[chrom != "chrM" & w >= 100 & w <= 300]
DT[, gc := round(gc, 2)]

DT.gc <- DT[,.(n=.N), by=.(id, gc, w)]
DT.gc <- DT.gc[gc >= .20 & gc <= .80]
DT.gc <- DT.gc[order(w, gc)]

filename <- file.path(out_file)
fwrite(DT.gc, file=filename)
q('no')
