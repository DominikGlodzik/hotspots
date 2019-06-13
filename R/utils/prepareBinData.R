prepareBinData <- function(maxBins=Inf,
	bedList=c(), # vector of paths to bed files, with names of variables
	bpList=NULL,
	cnDf=NULL,
	repTimingBed=NULL,
    genes.high.expr=NULL,
    genes.low.expr=NULL) {
	
	# design the bins: 500kb bins 
    binList <- list()
    for (c in c(as.character(1:22),'X')) {       
        chr.length <- seqlengths(Hsapiens)[[paste('chr',c,sep='')]]       
        binStarts <- seq(from=1, to=chr.length-5e5, by=binSize )
        binEnds <- binStarts + binSize-1
        binList[[c]] <- data.frame(chr=c, chromStart=binStarts, chromEnd=binEnds)
    }
    allBins <- do.call('rbind', binList)

    gr.allBins <- GRanges(seqnames=Rle(paste0('chr',allBins$chr)),
                          ranges=IRanges(allBins$chromStart, allBins$chromEnd),
                          strand=rep(c("*"), nrow(allBins)))
    	# down-sample the bins
	if (maxBins!=Inf) {
		si <- sample(nrow(allBins))
		allBins <- allBins[si,]
		gr.allBins <- gr.allBins[si,]
	}


	if (!is.null(genes.high.expr) & !is.null(genes.low.expr)) {
	    gr.highly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.high.expr$chr)),
	                         ranges=IRanges(genes.high.expr$chromStart+1, genes.high.expr$chromEnd),
	                         strand=rep(c("*"), nrow(genes.high.expr)) ,
	                             seqlengths=seqlengths(Hsapiens)                        
	                     )
	    gr.lowly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.low.expr$chr)),
	                                    ranges=IRanges(genes.low.expr$chromStart+1, genes.low.expr$chromEnd),
	                                    strand=rep(c("*"), nrow(genes.low.expr)) ,
	                                    seqlengths=seqlengths(Hsapiens)                        
	                                   )

	    coverage.expr.genes <- coverage(gr.highly.expr.genes)
    	coverage.low.expr.genes <- coverage( gr.lowly.expr.genes)
        allBins$highExpGenes <-  rep(NA, nrow(allBins))
    	allBins$lowExpGenes <-  rep(NA, nrow(allBins))
	}

    # replication timing: medianRepTime
    if (!is.null(repTimingBed)) {
    	rep.timing.df <- read.csv(repTimingBed, sep='\t', header=FALSE)
		colnames(rep.timing.df) <- c('chr', 'chromStart', 'chromEnd', 'timing')
		gr.timing <- GRanges(seqnames=Rle(paste0(rep.timing.df$chr)),
        	             ranges=IRanges(rep.timing.df$chromStart+1, rep.timing.df$chromEnd),
            	         strand=rep(c("*"), nrow(rep.timing.df )),
                	     timing=rep.timing.df$timing,
                    		         seqlengths=seqlengths(Hsapiens)
                     	)
		overlaps.rep.domains <- as.data.frame(findOverlaps(gr.allBins, gr.timing ))
	    allBins$medianRepTime <- rep(NA, nrow(allBins))
	}	

	# extract sequences of the bins, in preparation to count the non-mapping bases
    binSequences <-     as.character(getSeq(Hsapiens, paste0('chr',allBins$chr),
                                            start=allBins$chromStart,
                                            end=allBins$chromEnd))

    # copy number data
    if (!is.null(cnDf)) {
    	gr.ascat <-  GRanges(seqnames=Rle(paste0(cnDf$Chromosome)),
                      ranges=IRanges(cnDf$chromStart, cnDf$chromEnd),
                      strand=rep(c("*"), nrow(cnDf)),
                     seqlengths=seqlengths(Hsapiens),
                     totalCn=cnDf$total.copy.number.inTumour
                     )
	    bin.midpoints <- rowMeans(allBins[,c('chromStart', 'chromEnd')])
    	gr.bin.midpoints <- GRanges(seqnames=Rle(paste0('chr',allBins$chr)),
                          ranges=IRanges(bin.midpoints , bin.midpoints ),
                          strand=rep(c("*"), nrow(allBins)))
    	cn.overlaps <- as.data.frame(findOverlaps(gr.bin.midpoints, gr.ascat))
    	allBins$meanCn <- NA
    }


	covList <- list()
	# all the other bed files
	for (bi in 1:length(bedList)) {
		print(names(bedList)[bi])
		aBed <- read.table(bedList[bi])
		colnames(aBed) <- c('chr', 'chromStart', 'chromEnd')
		gr.bed <- GRanges(seqnames=Rle(paste0(aBed[,1])),
                    ranges=IRanges(pmax(aBed[,2],1), aBed[,3]),
                    strand=rep(c("*"), nrow(aBed)),
                    seqlengths=seqlengths(Hsapiens)
                    )
        coverage.bed <-coverage(gr.bed)
        covList[[bi]] <- coverage.bed

        allBins[,names(bedList)[bi]] <- NA

	}


	# counting the breakpoints
	if (!is.null(bpList)) {
		# for all breakpoint lists
		for (bps.i in 1:length(bpList)) {
			allBins[,names(bpList)[bps.i]] <- NA
			gr.bps <- GRanges(seqnames=Rle(paste0('chr',bpList[[bps.i ]]$chr)),
                    ranges=IRanges(bpList[[bps.i ]]$position, bpList[[bps.i ]]$position),
                        strand=rep(c("*"), nrow(bpList[[bps.i ]])))
			allBins[,names(bpList)[bps.i]] <-  countOverlaps(gr.allBins,gr.bps)
		}
	}

	# summarize all bins
    for (bi in 1:min(maxBins,nrow(allBins))) {
    	print(bi)

  		if (!is.null(genes.high.expr) & !is.null(genes.low.expr)) {
			allBins$highExpGenes[bi] <- sum(coverage.expr.genes [[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
    	    allBins$lowExpGenes[bi] <- sum(coverage.low.expr.genes[[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
       	}

        matching.domains <- subset(overlaps.rep.domains, queryHits==bi)
        if (nrow(matching.domains)>0) {
            allBins$medianRepTime[bi] <-  median(mcols(gr.timing[matching.domains$subjectHits,])$timing)
        }

        # count the N bases
        allBins$noNbases[bi] <- str_count(binSequences[bi],'N')

        # copy number
        if (!is.null(cnDf)) {
        	cn.overlaps.bin <- subset(cn.overlaps, queryHits==bi)
        	allBins$meanCn[bi] <- mean(mcols(gr.ascat[cn.overlaps.bin$subjectHits,])$totalCn)
        }

        # loop over the other bed files
       	for (bl in 1:length(bedList)) {
       		allBins[bi,names(bedList)[bl]] <-  sum(covList[[bl]][[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
       	}

    }
    return(allBins)

}