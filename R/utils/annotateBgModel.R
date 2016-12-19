 annotateBgModel <- function (hotspots.df, allBins, exp.col) {
                                        # Calculates expected breakpoint density in each segment. Informed by negative binomial regression.
                                        # Args:
                                        #   hotspots.df: data frame. One hotspot per row.
                                        #   allBins: result of fitting the negative regression
                                        #   exp.col: which column to use from allBins (which rearrangement signature.)

     for (hi in 1:nrow(hotspots.df)) { # loop over hotspots

         bin.overlap <-(pmin(allBins$chromEnd,   hotspots.df$end.bp[hi])-pmax(allBins$chromStart,   hotspots.df$start.bp[hi]))
         bin.overlap[(as.character(allBins$chr)!=as.character(hotspots.df$chr[hi]))] <- 0
         bin.overlap[bin.overlap<0] <- 0
         
         matching.bins <- which(bin.overlap >0)
         
         hotspots.df$d.bg[hi] <- sum(allBins[matching.bins, exp.col])/sum(allBins$chromEnd[matching.bins]-allBins$chromStart[matching.bins]  )        
     }

     hotspots.df$d.obs.exp <-  hotspots.df$d.seg/hotspots.df$d.bg

     
     hotspots.df


}
