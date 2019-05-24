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
    bedList <- c('/nfs/users/nfs_d/dg17/breast_rearr/data/metadata/dnase/MCF7DukeDNaseSeq.pk')
    names(bedList) <- c('dnase')

    # load the breakpoints
    DATA.VERSION <- '29.11'
    REARRANGEMENTS.PATH=paste0('/nfs/cancer_archive04/dg17/BASIS/versions/',DATA.VERSION,'/rearrangements')
    PRELOAD.REARR.PATH<-paste(REARRANGEMENTS.PATH, '/allRearrangements.RData',sep='')
    load(PRELOAD.REARR.PATH)
    bps.sv1 <-  subset(r$sample.bps.all, Signature.SV1>0.5)
    bpList=list()
    bpList[[1]] <- bps.sv1
    names(bpList)[1] <- 'counts.sv1'

	allBins <- prepareBinData(maxBins=100, bedList=bedList, bpList=bpList)
} else {
    load(paste0('../data/regressionData08.06-',binSize,'.RData'))
}

                                        # regression part
                                        # prepare bins for regression
#bins.regression <- subset(allBins, noNbases<0.1*binSize &
#                              noCensusGenes==0 & !is.nan(meanCn))
bins.regression <- allBins

### run the regression model
# choose the variables
variables.included <- c(2,4,5,6)
model.variables <- c('noNbases', 'medianRepTime', 'meanCn', 'dnase', 'highExpGenes', 'lowExpGenes', 'staining', 'fragile', 'counts.alu', 'counts.rep',
                     'segDup')[variables.included ]
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


