prepareBinData <- function(maxBins=Inf,
	bedList=c(), # vector of paths to bed files, with names of variables
	bpList=NULL
	) {
	
	# design the bins: 500kb bins every 250mb
    binList <- list()
    for (c in as.character(1:22)) {       
        chr.length <- seqlengths(Hsapiens)[[paste('chr',c,sep='')]]       
        binStarts <- seq(from=1, to=chr.length-5e5, by=binSize )
        binEnds <- binStarts + binSize-1
        binList[[c]] <- data.frame(chr=c, chromStart=binStarts, chromEnd=binEnds)
    }
    allBins <- do.call('rbind', binList)
    allBins$highExpGenes <-  rep(NA, nrow(allBins))
    allBins$lowExpGenes <-  rep(NA, nrow(allBins))
    allBins$medianRepTime <- rep(NA, nrow(allBins))

    gr.allBins <- GRanges(seqnames=Rle(paste0('chr',allBins$chr)),
                          ranges=IRanges(allBins$chromStart, allBins$chromEnd),
                          strand=rep(c("*"), nrow(allBins)))
    	# down-sample the bins
	if (maxBins!=Inf) {
		si <- sample(nrow(allBins))
		allBins <- allBins[si,]
		gr.allBins <- gr.allBins[si,]
	}

    # load the list of genes first
	load('../../../data/biomart.19.05.2016.RData') # incl. genes.table

    # gene expression: highExpGenes', 'lowExpGenes'


	if (!exists('expr.data')) {
	    expr.matrix <- read.table('/nfs/cancer_archive04/dg17/BASIS/expression/FPKM/QuantNorm_log2_FPKM_n342.txt', header=TRUE, sep='\t')
	    expr.matrix.raw <- expr.matrix[, 6:ncol(expr.matrix )]
	    expr.matrix.raw.t <- t(expr.matrix[, 6:ncol(expr.matrix )])
	    colnames(expr.matrix.raw.t) <-  expr.matrix$Ensembl

	    medianGeneExpression <- data.frame(EnsemblID=expr.matrix$Ensembl, medianExpression=apply(expr.matrix.raw.t, 2, median, na.rm=TRUE))
	    medianGeneExpression <- subset(medianGeneExpression, !is.na(medianExpression))
	    medianGeneExpression$exprQuantile <- with(medianGeneExpression, factor(findInterval(medianExpression, c(-Inf, quantile(medianExpression, probs=c(0.25, .5, .75)), Inf) ), labels=c("gene.expr.Q1","gene.expr.Q2","gene.expr.Q3","gene.expr.Q4")))
	    
	    
	    genes.high.expr <- genes.table
	    genes.high.expr <- subset(genes.high.expr,ensembl_gene_id %in% subset(medianGeneExpression, exprQuantile=="gene.expr.Q4")$EnsemblID)

	    gr.highly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.high.expr$chr)),
	                         ranges=IRanges(genes.high.expr$chromStart+1, genes.high.expr$chromEnd),
	                         strand=rep(c("*"), nrow(genes.high.expr)) ,
	                             seqlengths=seqlengths(Hsapiens)                        
	                     )

	    genes.low.expr <- subset(genes.table,ensembl_gene_id %in% subset(medianGeneExpression, exprQuantile!="gene.expr.Q4")$EnsemblID)
	    gr.lowly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.low.expr$chr)),
	                                    ranges=IRanges(genes.low.expr$chromStart+1, genes.low.expr$chromEnd),
	                                    strand=rep(c("*"), nrow(genes.low.expr)) ,
	                                    seqlengths=seqlengths(Hsapiens)                        
	                                   )

	    coverage.expr.genes <- coverage(gr.highly.expr.genes)
    	coverage.low.expr.genes <- coverage( gr.lowly.expr.genes)
	}

    # replication timing: medianRepTime
    rep.timing.df <- read.csv('/lustre/scratch116/casm/cgp/users/sm22/nfs/cancer_archive04/sm22/BASIS/Repliseq_data/RepliTime/MCF7_data/MCF7_RepliSeq.bedGraph', sep='\t', header=FALSE)
	colnames(rep.timing.df) <- c('chr', 'chromStart', 'chromEnd', 'timing')
	gr.timing <- GRanges(seqnames=Rle(paste0(rep.timing.df$chr)),
                     ranges=IRanges(rep.timing.df$chromStart+1, rep.timing.df$chromEnd),
                     strand=rep(c("*"), nrow(rep.timing.df )),
                     timing=rep.timing.df$timing,
                             seqlengths=seqlengths(Hsapiens)
                     )
	overlaps.rep.domains <- as.data.frame(findOverlaps(gr.allBins, gr.timing ))


	covList <- list()
	# all the other bed files
	for (bi in 1:length(bedList)) {
		aBed <- read.table(bedList[bi])
		colnames(aBed) <- c('chr', 'chromStart', 'chromEnd')
		gr.bed <- GRanges(seqnames=Rle(paste0(aBed[,1])),
                    ranges=IRanges(aBed[,2], aBed[,3]),
                    strand=rep(c("*"), nrow(aBed)),
                    seqlengths=seqlengths(Hsapiens)
                    )
        coverage.bed <-coverage(gr.bed)
        covList[[bi]] <- coverage.bed

        allBins[,names(bedList)[bi]] <- NA

	}





	# breakpoints
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

	# Ssummarize all bins
    for (bi in 1:min(maxBins,nrow(allBins))) {
    	print(bi)
		allBins$highExpGenes[bi] <- sum(coverage.expr.genes [[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
        allBins$lowExpGenes[bi] <- sum(coverage.low.expr.genes[[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
       
        matching.domains <- subset(overlaps.rep.domains, queryHits==bi)
        if (nrow(matching.domains)>0) {
            allBins$medianRepTime[bi] <-  median(mcols(gr.timing[matching.domains$subjectHits,])$timing)
        }

        # loop over the other bed files
       	for (bl in 1:length(bedList)) {
       		allBins[bi,names(bedList)[bl]] <-  sum(covList[[bl]][[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
       	}

    }
    return(allBins)

}