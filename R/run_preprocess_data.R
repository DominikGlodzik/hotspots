    # process the gene expression data
genes.table <- read.csv('../data/genes.table.csv')
expr.matrix <- read.table('/nfs/cancer_archive04/dg17/BASIS/expression/FPKM/QuantNorm_log2_FPKM_n342.txt', header=TRUE, sep='\t')
expr.matrix.raw <- expr.matrix[, 6:ncol(expr.matrix )]
expr.matrix.raw.t <- t(expr.matrix[, 6:ncol(expr.matrix )])
colnames(expr.matrix.raw.t) <-  expr.matrix$Ensembl
medianGeneExpression <- data.frame(EnsemblID=expr.matrix$Ensembl, medianExpression=apply(expr.matrix.raw.t, 2, median, na.rm=TRUE))
medianGeneExpression <- subset(medianGeneExpression, !is.na(medianExpression))
medianGeneExpression$exprQuantile <- with(medianGeneExpression, factor(findInterval(medianExpression, c(-Inf, quantile(medianExpression, probs=c(0.25, .5, .75)), Inf) ), labels=c("gene.expr.Q1","gene.expr.Q2","gene.expr.Q3","gene.expr.Q4")))      
genes.high.expr <- genes.table
genes.high.expr <- subset(genes.high.expr,ensembl_gene_id %in% subset(medianGeneExpression, exprQuantile=="gene.expr.Q4")$EnsemblID)
genes.low.expr <- subset(genes.table,ensembl_gene_id %in% subset(medianGeneExpression, exprQuantile!="gene.expr.Q4")$EnsemblID)
        
write.csv(genes.high.expr, file='../data/genes.high.expr.csv')
write.csv(genes.low.expr , file='../data/genes.low.expr.csv')


   # chromosome staining
    #   
recreateBands <- FALSE
if (recreateBands) {
    FILE.bands <- '../data/chromBands'
    bands <- read.table(FILE.bands, header=TRUE)
    bands$chrom <- bands$chrom,4,100000L
    bands$chromStart <- bands$chromStart +1
    bands <- bands [order(bands$chrom, bands$chromStart ),] # order the bands
    names(bands) <- c('chr', 'chromStart',  'chromEnd', 'name', 'gieStain')
    bands$isStaining <- bands$gieStain %in% c('gpos100', 'gpos25', 'gpos50', 'gpos75')
    bands <- subset(bands, isStaining==TRUE)
    write.table(bands, file= '/nfs/users/nfs_d/dg17/breast_rearr/data/chromBands/chromBandsPos', sep='\t', quote=FALSE)
}

        # load the breakpoints
DATA.VERSION <- '29.11'
REARRANGEMENTS.PATH=paste0('/nfs/cancer_archive04/dg17/BASIS/versions/',DATA.VERSION,'/rearrangements')
PRELOAD.REARR.PATH<-paste(REARRANGEMENTS.PATH, '/allRearrangements.RData',sep='')
load(PRELOAD.REARR.PATH)
# prepare the breakpoint data
# data frame with columns 'chr' (values without chr prefix) and 'posiiton'
bps.sv1 <-  subset(r$sample.bps.all, Signature.SV1>0.5)
write.csv(bps.sv1, file='../data/bps.sv1.csv')