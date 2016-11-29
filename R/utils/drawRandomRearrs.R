drawRandomRearrs <- function(allBins, rearrs.df, binMeanCol, binSize) {

    
    simulated.rearrs.df <- rearrs.df[,c('sample', 'Chromosome.1', 'Chromosome.2', 'pf', 'length')]

    no.rearrs <- nrow(rearrs.df)


    coord.bp1 <- vector()
    coord.bp2 <- vector()

    is.within.ref.bp1 <- vector()
    is.within.ref.bp2 <- vector()
    
    ref.gen.bp1 <- rep('N',no.rearrs) # store base at the first breakpoint
    ref.gen.bp2 <- rep('N',no.rearrs) # store base at the first breakpoint
    
    isNonMapping <- ref.gen.bp1=='N' | ref.gen.bp2=='N'
    
    while (sum(isNonMapping)>0) { # while any of the breakpoints are not mapping
        


        
        # sampling proportionally to the mean
        bin.draws <- sample(nrow(  allBins), size=no.rearrs , replace=TRUE, prob=allBins[,binMeanCol]/sum(allBins[,binMeanCol]))
        bin.draws2 <- sample(nrow(  allBins), size=no.rearrs, replace=TRUE, prob=allBins[,binMeanCol]/sum(allBins[,binMeanCol])) # in case of translocation

                                        # sampling proportional to the nb draw
        #nb.draw <- vector()
        #for (bi in 1:nrow(  allBins)) { # draw from negative-bionmial distribution at each bin
        #    nb.draw[bi] <- rnbinom(1, mu=allBins[bi, binMeanCol], size=2.2)
        #}
        #bin.draws <- sample(nrow(  allBins), size=no.rearrs , replace=TRUE, prob=nb.draw/sum(nb.draw))
                                        #bin.draws2 <- sample(nrow(  allBins), size=no.rearrs, replace=TRUE, prob=nb.draw/sum(nb.draw)) # in case of translocation
        
        bin.bp <- sample.int(binSize-1, no.rearrs, replace=TRUE) # random position within each bin
        
                                        # draw the first and second coordinate
        coord.bp1[isNonMapping] <- allBins[bin.draws[isNonMapping], 'chromStart'] + bin.bp[isNonMapping]
        simulated.rearrs.df$Chromosome.1[isNonMapping] <- allBins[bin.draws[isNonMapping], 'chr']
                                        # non-translocations
        coord.bp2[isNonMapping & simulated.rearrs.df$pf!=32] <- coord.bp1[isNonMapping & simulated.rearrs.df$pf!=32] + rearrs.df$length[isNonMapping & simulated.rearrs.df$pf!=32]
        simulated.rearrs.df$Chromosome.2[isNonMapping & simulated.rearrs.df$pf!=32] <- allBins[bin.draws[isNonMapping & simulated.rearrs.df$pf!=32], 'chr']
                                        # translocations
        coord.bp2[isNonMapping & simulated.rearrs.df$pf==32] <-  allBins[bin.draws2[isNonMapping & simulated.rearrs.df$pf==32], 'chromStart'] + bin.bp[isNonMapping & simulated.rearrs.df$pf==32]     
        simulated.rearrs.df$Chromosome.2[isNonMapping & simulated.rearrs.df$pf==32] <- allBins[bin.draws2[isNonMapping & simulated.rearrs.df$pf==32], 'chr']
        
                                        # reference genome at the drawn sites
        # update the base for the breakpoints that fall in the reference genome and are not mapping yet
        is.within.ref.bp1 <- coord.bp1<=seqlengths(Hsapiens)[paste0('chr', simulated.rearrs.df$Chromosome.1)] & isNonMapping
        if (sum(is.within.ref.bp1)>0) {
            ref.gen.bp1[is.within.ref.bp1 ] <- as.character(getSeq(Hsapiens, paste0('chr',simulated.rearrs.df$Chromosome.1[is.within.ref.bp1 ]), start=coord.bp1[is.within.ref.bp1 ], end=coord.bp1[is.within.ref.bp1 ]))
        }
                                        # update the base for the breakpoints that fall in the reference genome and are not mapping yet        

        is.within.ref.bp2 <- coord.bp2<=seqlengths(Hsapiens)[paste0('chr', simulated.rearrs.df$Chromosome.2)] & isNonMapping
        if (sum(is.within.ref.bp2)>0) {
            ref.gen.bp2[is.within.ref.bp2 ] <- as.character(getSeq(Hsapiens, paste0('chr',simulated.rearrs.df$Chromosome.2[is.within.ref.bp2 ]), start=coord.bp2[is.within.ref.bp2], end=coord.bp2[is.within.ref.bp2 ]))
        }
        
        isNonMapping <- ref.gen.bp1=='N' | ref.gen.bp2=='N'
        print (paste(sum(isNonMapping), 'rearrnagements dont map'))
        
    }

    simulated.rearrs.df$pos.1 <- coord.bp1
    simulated.rearrs.df$pos.2 <- coord.bp2
    simulated.rearrs.df
}
