calc.gw.mut.rate <- function(bps) {

    total.size <- 0

 
    bps <- subset(bps, chr!='Y')
    chroms <-  unique(bps$chr)
    
    chr.rate <- list()
    
    for (loop.chr in chroms) {
        cat(paste('calculating rate on chromosome', loop.chr), '\n')
        bps.chr <- subset(bps,chr==loop.chr)
     
                                        # right bp

        # number of non-N bases within this range
        chr.seq <- getSeq(Hsapiens, paste0('chr',loop.chr), start=min(bps.chr$pos), end=max(bps.chr$pos))
        f <- alphabetFrequency(chr.seq, baseOnly=TRUE)
        chr.mappable.size <-  max(bps.chr$pos) - min(bps.chr$pos) - f[['other']] # subtracting the positions in the genome not mappable
        
        total.size <- total.size + chr.mappable.size

        chr.rate[[loop.chr]] <- nrow(bps.chr)/chr.mappable.size 
    }

    names(chr.rate)[names(chr.rate)=='X'] <- '23'
    chr.rate <- chr.rate[order(as.integer(names(chr.rate)))]
    result <- list()
    result$gw.rate <- nrow(bps)/total.size
    result$chr.rate <- chr.rate
    result
}
