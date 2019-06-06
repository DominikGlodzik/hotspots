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


print('Loading the data ...')
bedList <- c(            '../data/breastData/genes.high.expr.tsv',
            '../data/breastData/genes.low.expr.tsv',
            '../data/breastData/MCF7DukeDNaseSeq.pk',
            '../data/breastData/chromBandsPos',
            '../data/breastData/broadFragile37Chrom.txt',
            '../data/breastData/alu.repeatMasker',
            '../data/breastData/repeatMasker1000',
            '../data/breastData/GRCh37GenomicSuperDup.nohead.tab'
    )
names(bedList) <- c('highExpGenes', 'lowExpGenes', 'dnase', 'staining', 'fragile', 'counts.alu', 'counts.rep','segDup' )

# the following columns are required: c('chr', 'chromStart', 'chromEnd', 'timing')
repTimingBed <- '../data/breastData/MCF7_RepliSeq.bedGraph'

# load the breakpoints
bps.sv1 <- read.csv('../data/breastData/bps.sv1.csv')
bpList=list()
bpList[[1]] <- bps.sv1
names(bpList)[1] <- 'counts.sv1'

# prepare the copy number data
# data frame with columns Chromosome chromStart  chromEnd  total.copy.number.inTumour
ascat.df <- read.csv('../data/breastData/ascat.df.csv')

# process the gene expression data


allBins <- prepareBinData(
    maxBins=100, 
    bedList=bedList, 
    bpList=bpList, 
    cnDf=ascat.df , 
    repTimingBed=repTimingBed
    )



# regression part
# prepare bins for regression
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
