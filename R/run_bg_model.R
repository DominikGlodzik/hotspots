library(lmtest)
library(MASS)
library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
source('utils/plotCoefficients.R')
source('utils/prepareBinData.R')

normaliseVector <- function(v, refv=NULL) {
    if (is.null(refv)) {
        r <- (v - mean(v, na.rm=TRUE))/sd(v, na.rm=TRUE)
    } else {
        r <- (v - mean(refv, na.rm=TRUE))/sd(refv, na.rm=TRUE)
    }
    return(r)
}


# define a data frame with relevant features
binSize = 5e5
prepareData <- TRUE
if (prepareData) {
    print('Loading the data ...')
    bedList <- c('/nfs/users/nfs_d/dg17/breast_rearr/data/metadata/dnase/MCF7DukeDNaseSeq.pk',
                '/nfs/users/nfs_d/dg17/breast_rearr/data/chromBands/chromBandsPos',
                '/nfs/users/nfs_d/dg17/breast_rearr/data/fragileSites/broadFragile37Chrom.txt',
                '/nfs/users/nfs_d/dg17/breast_rearr/data/metadata/repeats/alu.repeatMasker',
                '/nfs/users/nfs_d/dg17/breast_rearr/data/metadata/repeats/repeatMasker',
                '/nfs/users/nfs_d/dg17/breast_rearr/data/metadata/segDup/GRCh37GenomicSuperDup.nohead.tab'
        )
    names(bedList) <- c('dnase', 'staining', 'fragile', 'counts.alu', 'counts.rep','segDup')

    repTimingBed <- '/lustre/scratch116/casm/cgp/users/sm22/nfs/cancer_archive04/sm22/BASIS/Repliseq_data/RepliTime/MCF7_data/MCF7_RepliSeq.bedGraph'

    # load the breakpoints
    DATA.VERSION <- '29.11'
    REARRANGEMENTS.PATH=paste0('/nfs/cancer_archive04/dg17/BASIS/versions/',DATA.VERSION,'/rearrangements')
    PRELOAD.REARR.PATH<-paste(REARRANGEMENTS.PATH, '/allRearrangements.RData',sep='')
    load(PRELOAD.REARR.PATH)
    # prepare the breakpoint data
    # data frame with columns 'chr' (values without chr prefix) and 'posiiton'
    bps.sv1 <-  subset(r$sample.bps.all, Signature.SV1>0.5)
    bpList=list()
    bpList[[1]] <- bps.sv1
    names(bpList)[1] <- 'counts.sv1'

    # prepare the copy number data
    # data frame with columns Chromosome chromStart  chromEnd  total.copy.number.inTumour
    ascat.df <- read.csv('../data/ascat.df.csv')

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

	allBins <- prepareBinData(
        maxBins=100, 
        bedList=bedList, 
        bpList=bpList, 
        cnDf=ascat.df , 
        repTimingBed=repTimingBed,
        genes.high.expr=genes.high.expr,
        genes.low.expr=genes.low.expr
        )

} else {
    print('loading the example data')
    load(paste0('../data/regressionData08.06-',binSize,'.RData'))
}

                                        # regression part
                                        # prepare bins for regression
#bins.regression <- subset(allBins, noNbases<0.1*binSize &
#                              noCensusGenes==0 & !is.nan(meanCn))
bins.regression <- allBins

### run the regression model
# choose the variables
variables.included <- 1:11
model.variables <- c('noNbases', 'medianRepTime', 'meanCn', 'dnase', 'highExpGenes', 'lowExpGenes', 
    'staining', 'fragile', 'counts.alu', 'counts.rep','segDup')[variables.included ]
model.variable.names <- c('N bases', 'early replication', 'copy number', 'DNAse', 'highly expressed genes', 'lowly expressed genes', 'chromatin staining', 
'fragile sites', 'ALU elements', 'repeats', 'segmental duplications', 'copy number ovarian')[variables.included ]#, 'histone mod H3K4me3 (MCF7)', 'histone mod H3K9me3 (MCF7)', 'histone mod H3K27ac (MCF7)', 'histone mod me1K4 (K562)')
names(model.variable.names) <- paste0(model.variables , '.norm')


# normalize the data
names(model.variable.names) <- paste0(model.variables, '.norm')
for (vi in 1:length(model.variables)) {
    normalisedV <-  normaliseVector(bins.regression[,as.character(model.variables[vi])])
    vn <- paste0(model.variables[vi], '.norm')
    bins.regression[[vn]] <- normalisedV
}

# fit the regression model
print('Fitting regression model...')
formula.sv1 <- as.formula(paste('counts.sv1 ~ ' , paste(paste0(model.variables,'.norm'), collapse=' + ')))
fm_nbin.sv1 <- glm.nb(formula.sv1 , data = bins.regression) # negative binomial fit

# make plots
pdf(paste0('../data/regressionCoefficients-bin-',binSize,'.pdf'), width=12, height=15)
layout(matrix(c(1,2,3,4, 5, 6), 3, 2, byrow = TRUE))
par(mar=c(3,13,3,3))

c <- plotCoefficients(fm_nbin.sv1, nbfit.ref=NULL, 'RS1 breast', model.variable.names)
dev.off()
