args <- commandArgs(trailingOnly = TRUE)
bindir <- args[1]
armsreffile <- args[2]
outdir <- args[3]
bin_size <- args[4]
gnom_build <- args[5]
print(bindir)
print(armsreffile)

library(data.table)
source("PlasmaTools_Functions.R")


binslist <- list.files(bindir, full.names=TRUE)
armsref <- fread(armsreffile)
binsforzscores <- rbindlist(lapply(binslist, fread))
binsforzscores <- binsforzscores[chr!="chrX"]

binsforzscores[,cov:=short+long]
zscores <- getZscores(binsforzscores, refpanel=armsref, measure="cov")

zscores[,normal.indices:=1]

zscores[,normal.indices:=NULL]

fwrite(zscores, paste0(outdir, "/zscores.",bin_size,".",gnom_build,".csv"))
q('no')
