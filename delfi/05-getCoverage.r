
args <- commandArgs(trailingOnly = TRUE)

gnom_build <- args[1]
bindir <- args[2]
outfile <- args[3]

library(data.table)
library(readr)

print("loaded!")
filenames<-list.files(bindir,pattern=paste0(gnom_build,".csv"),full.names=TRUE)

x <- rbindlist(lapply(filenames, fread))
x_f <- x[chr!="chrX",]

setkey(x_f, id)
x_f[,cov:=short+long]
x_f[,ratio.scaled:=scale(ratio.cor), by=id]
x_f[,ratiovar:=factor(paste0("ratio_", bin), paste0("ratio_", 1:.N)), by=id]

features.ratios <- dcast(x_f, id ~ ratiovar, value.var="ratio.scaled")


## Create coverage ##
x_f[,cov.cor:=short.cor+long.cor]
x_f[,cov.scaled:=scale(cov.cor), by=id]
x_f[,covvar:=factor(paste0("cov_", bin), paste0("cov_", 1:.N)), by=id]

features.covs <- dcast(x_f, id ~ covvar, value.var="cov.scaled")
setkey(features.ratios, id)
setkey(features.covs, id)
features.full <- features.covs[features.ratios]

write_csv(features.full,outfile)
