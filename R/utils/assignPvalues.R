assignPvalues <- function(kat.regions, chrom.bps, bp.rate=NA) {

                                        # annotate the hotspots:
                                        # d.seg (density of breakpoints),
                                        # rate.factor (density of breakpoint in hospots compared to chormosome)
                                        # pvalue (probability of observint this many breakpoints, based on chromsomome expectation, binomial)

if (is.na(bp.rate)) { # estimate the chromosome rate
    left.bp <- min(chrom.bps$pos)
    right.bp <-  max(chrom.bps$pos)
    bp.rate <- nrow(chrom.bps)/ (right.bp - left.bp)
}


# assume binomial distribution
kat.regions$pvalue <- 1-pbinom(kat.regions$number.bps, kat.regions$end.bp - kat.regions$start.bp, bp.rate)

kat.regions$d.seg<- (kat.regions$number.bps/( kat.regions$end.bp - kat.regions$start.bp))

kat.regions$rate.factor <- kat.regions$d.seg/bp.rate

kat.regions
}
