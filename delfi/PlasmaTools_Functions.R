gcCorrectTarget <- function(fragments, ref, bychr=TRUE){
  fragments[, gc := round(gc, 2)]
  if(bychr) {
    DT.gc <- fragments[,.(n=.N), by=.(gc, chr)]
    DT.gc <- DT.gc[gc >= .20 & gc <= .80]
    DT.gc <- DT.gc[order(gc, chr)]
  } else {
    DT.gc <- fragments[,.(n=.N), by=gc]
    DT.gc <- DT.gc[gc >= .20 & gc <= .80]
    DT.gc <- DT.gc[order(gc)]
  }
  # setkey(mediandt, gc, seqnames)
  if(bychr) {
    setkey(DT.gc, gc, chr)
    setkey(ref, gc, chr)
  } else {
    setkey(DT.gc, gc)
    setkey(ref, gc)
  }
  #     DT.gc <- DT.gc[ref][order(chr, gc)]
  DT.gc <- DT.gc[ref]
  DT.gc[,w:=target/n]
  if(bychr) {
    fragments[DT.gc, on= .(chr, gc), weight := i.w]
  }
  else fragments[DT.gc, on= .(gc), weight := i.w]
  fragments <- fragments[!is.na(weight)] 
  #fragments[,weight := weight * .N/sum(weight)] # Only difference from PlasmaTools version  of function
  fragments[]
}


binFrags <- function(fragments, bins, cutoff=150,
                     chromosomes=paste0("chr",c(1:22, "X"))) {
  setkey(bins, chr, start, end)
  fragbins <- foverlaps(fragments[chr %in% chromosomes],
                        bins, type="within", nomatch=NULL)
  bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                       short = sum(w >= 100 & w <= cutoff ),
                       long = sum(w > cutoff & w <= 250),
                       short.cor = sum(weight[w >= 100 & w <= cutoff]),
                       long.cor = sum(weight[w > cutoff & w <= 250]),
                       ultrashort = sum(w < 100),
                       ultrashort.cor = sum(weight[w < 100]),
                       multinucs = sum(w > 250),
                       multinucs.cor = sum(weight[w > 250]),
                       mediansize = as.integer(median(w)),
                       frag.gc = mean(fraggc)),
                    by=.(chr, start, end)]
  
  setkey(bins2, chr, start, end)
  bins2 <- bins2[bins]
  bins2 <- bins2[is.na(i.gc), which(grepl("short|long|multi", colnames(bins2))):=0]
  bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
  bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
  bins2[, bin:=1:.N]
  setcolorder(bins2, c("chr", "start", "end", "bin"))
  bins2[]
}



getZscores <- function(bins, normal.ids=NULL, loo = FALSE,
                       refpanel=NULL, measure="cov") {
  bins2 <- copy(bins)
  countAndNormalize(bins2, measure=measure)
  armmeansdt <- getArmMeans(bins2)
  
  armmeansdt <- armmeansdt[refpanel, on="arm"]
  armmeansdt[,zscore:=(armmean-mu)/sigma]
  armmeansdt[,`:=`(armmean=NULL, mu=NULL, sigma=NULL)]
  armmeansdt[]
  
  #     if(is.null(refpanel) && !is.null(normal.ids)) {
  #         armmeansdt[,normal.indices:=createNormalIndex(id, normal.ids)]
  #         #     armmeans[,normal.means:=mean(armmean[normal.indices]),by=.(arm)]
  #         armmeansdt[,zscore:=.loo.zscore(armmean, normal.indices), by=arm]
  #         armmeansdt[]
  #     }
}


countAndNormalize <- function(bins, measure) {
  ### normalize for filtered bases
  ### weight is to get count per expected width of basepairs
  ## Maybe
  bins[,weight:=(end - start + 1)/(end - start + 1 - filtered.bases)]
  bins[,counts.mult:=log2((get(measure) * weight + 1))]
  bins[,loess.pred:=gcCorrectLoess(counts.mult, gc),by=id]
  bins[,adjusted:=counts.mult - loess.pred]
  bins[]
}


gcCorrectLoess <- function(counts.mult, gc) {
  lower <- 0
  upper <- quantile(counts.mult, .99, na.rm=TRUE)
  trend.points <- counts.mult > lower & counts.mult < upper
  trend.counts <- counts.mult[trend.points]
  trend.gc <- gc[trend.points]
  num.sample.points <- min(length(trend.counts), 10E3L)
  samp <- sample(1:length(trend.counts), num.sample.points)
  #    pad counts and GC
  include <- c(which(gc==min(gc)), which(gc==max(gc)))
  trend.counts <- c(counts.mult[include], trend.counts[samp])
  trend.gc <- c(gc[include], trend.gc[samp])
  initial.trend <- loess(trend.counts ~ trend.gc)
  i <- seq(min(gc, na.rm=TRUE), max(gc, na.rm=TRUE), by = 0.001)
  final.trend <- loess(predict(initial.trend, i) ~ i)
  counts.pred <- predict(final.trend, gc)
  return(counts.pred)
}


getArmMeans <- function(bins) {
  bins[,adjusted.cent:=adjusted-median(adjusted, na.rm=TRUE), by=id]
  arms <- bins[,.(armmean=mean(adjusted.cent, na.rm=TRUE)),
               by=.(id, arm)]
  arms <- arms[,armmean:=armmean - mean(armmean, na.rm = TRUE)]
  arms[]
}



.loo.zscore <-  function(x, indices) {
  zscores <- rep(NA, length(x))
  for(i in seq_along(indices)) {
    if(indices[i]) {
      keep <- indices==1
      keep[i] <- FALSE
    } else {
      keep <- indices==1
    }
    mu <- mean(x[keep])
    sigma <- sd(x[keep])
    zscores[i] <- (x[i] - mu)/sigma
  }
  zscores
}



frag.stats <- function(fragments, cutoff=150, groups=NULL) {
  Mode <- function(x){
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  fragments[,.(mode=Mode(w), median=as.integer(median(w)),
               mean=mean(w), size_iqr=matrixStats::iqr(w),
               nfrags=.N, mononucs=sum(w>=100 & w<=250),
               multinucs = sum(w>=250), ultrashort = sum(w<100),
               slratio = sum(w<=cutoff)/sum(w > cutoff & w <= 250),
               meangc=mean(gc), mediangc=median(gc),
               gc_iqr=matrixStats::iqr(gc), chrMrep=-log(sum(chr=="chrM")/.N)),
            keyby=groups][]
}
