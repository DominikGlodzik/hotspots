prepareBinData <- function() {
	
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

	# Ssummarize all bins
    for (bi in 1:nrow(allBins)) {
		allBins$highExpGenes[bi] <- sum(coverage.expr.genes [[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
        allBins$lowExpGenes[bi] <- sum(coverage.low.expr.genes[[paste0('chr', allBins$chr[bi])]][allBins$chromStart[bi]: allBins$chromEnd[bi]]>0)
       
        matching.domains <- subset(overlaps.rep.domains, queryHits==bi)
        if (nrow(matching.domains)>0) {
            allBins$medianRepTime[bi] <-  median(mcols(gr.timing[matching.domains$subjectHits,])$timing)
        }

    }
    return(allBins)

}