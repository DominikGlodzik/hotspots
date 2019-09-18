downsampleRearrBins <- function(rearr.bps){
    # generate bins
    binList <- list()
    for (c in c(as.character(1:22), 'X')) {       
        chr.length <- seqlengths(Hsapiens)[[paste('chr',c,sep='')]]       
        binStarts <- seq(from=1, to=chr.length-binSize, by=binSize ) # ER: edit so that bin size is not hard-coded in "to"
        binEnds <- binStarts + binSize-1
        binList[[c]] <- data.frame(chr=c, chromStart=binStarts, chromEnd=binEnds)
    }
    bins <- do.call('rbind', binList)
    bins$bin <- row.names(bins)
    bins$chr <- as.character(bins$chr)
    
    # randomly sample one breakpoint per bin
    set.seed(42)
    out <- rearr.bps %>%
        left_join(bins, by = "chr") %>%
        filter(position > chromStart & position < chromEnd) %>%
        group_by(sample, bin) %>%
        sample_n(size = 1) %>%
        as.data.frame() %>%
        dplyr::select(-bin, -chromStart, -chromEnd)
        
    return(out)
}
    