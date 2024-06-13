
args <- commandArgs(trailingOnly = TRUE)
gnom_build <- args[1]
bin_size <- args[2]
outdir <- args[3]

zpath <- paste0(outdir,"/zscores.",bin_size,".",gnom_build,".csv")
covpath <- paste0(outdir,"/coverage.",bin_size,".",gnom_build,".csv")
outfile <- paste0(outdir,"/allfeatures.",bin_size,".",gnom_build,".csv")

#collect all the bins
library(data.table)
suppressMessages(library(devtools))
library(dplyr)
library(tidyverse)

zscores<- read_csv(zpath)
cov <- read_csv(covpath)

setDT(zscores)
setDT(cov)
setkey(cov, id)
setkey(zscores, id)
full_features <- zscores[cov]


write_csv(full_features,outfile)
