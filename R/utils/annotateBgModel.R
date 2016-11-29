 annotateBgModel <- function (hotspots.df, allBins, exp.col) {


     for (hi in 1:nrow(hotspots.df)) {

         bin.overlap <-(pmin(allBins$chromEnd,   hotspots.df$end.bp[hi])-pmax(allBins$chromStart,   hotspots.df$start.bp[hi]))
         bin.overlap[(as.character(allBins$chr)!=as.character(hotspots.df$chr[hi]))] <- 0
         bin.overlap[bin.overlap<0] <- 0
         
         matching.bins <- which(bin.overlap >0)
         
         hotspots.df$d.bg[hi] <- sum(allBins[matching.bins, exp.col])/sum(allBins$chromEnd[matching.bins]-allBins$chromStart[matching.bins]  )
         # hotspots.df$d.bg.acc[hi] <- sum((bin.overlap/5e5)*allBins[, exp.col])/sum(bin.overlap)
         
     }


     hotspots.df$d.obs.exp <-  hotspots.df$d.seg/hotspots.df$d.bg
     #hotspots.df$d.obs.exp.acc <-  hotspots.df$d.seg/hotspots.df$d.bg.acc
     
     hotspots.df


}
