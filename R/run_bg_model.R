library(lmtest)
library(MASS)
library(stringr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
source('utils/plotCoefficients.R')
source('utils/prepareBinData.R')
#library(randomForest)
#library(reprtree)
#library(randomForestExplainer)


normaliseVector <- function(v, refv=NULL) {
    if (is.null(refv)) {
        r <- (v - mean(v, na.rm=TRUE))/sd(v, na.rm=TRUE)
    } else {
        r <- (v - mean(refv, na.rm=TRUE))/sd(refv, na.rm=TRUE)
    }
    return(r)
}


# define a data frame with relevant features
binSize = 1e4

load('../data/bgModel-13.06-5e+05.RData')
load('../data/allRearrangements.29.11.RData')


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
repTimingBed <- '../data/breastData/MCF7_RepliSeq.bedGraph.n1000'

# load the breakpoints
bps.sv1 <- read.csv('../data/breastData/bps.sv1.csv')
bpList=list()
# random SV1
random.sv1.bed <- subset( r$random.bps.all, Signature.SV1>0.5)
random.sv1.bed <- random.sv1.bed[,c('chr', 'position',  'position')]
bpList[[1]] <- random.sv1.bed
names(bpList)[1] <- 'counts.random.sv1'
# random SV3
random.sv3.bed <- subset( r$random.bps.all, Signature.SV3>0.5)
random.sv3.bed <- random.sv3.bed[,c('chr', 'position',  'position')]
bpList[[2]] <- random.sv3.bed
names(bpList)[2] <- 'counts.random.sv3'



allBins2 <- prepareBinData(
    bedList=c(), 
    bpList=bpList, 
    cnDf=NULL , 
    repTimingBed=NULL,
    addSummary=FALSE
    )



# prepare the copy number data
# data frame with columns Chromosome chromStart  chromEnd  total.copy.number.inTumour
ascat.df <- read.csv('../data/breastData/ascat.df.csv')

# process the gene expression data


allBins <- prepareBinData(
    #maxBins=100, 
    bedList=bedList, 
    bpList=bpList, 
    cnDf=ascat.df , 
    repTimingBed=repTimingBed
    )
save(paste0('~/repos/hotspots/data/allBins_',binSize, '.RData' ))






# regression part
# prepare bins for regression
bins.regression <- allBins

### run the regression model
# choose the variables
variables.included <- 1:11
model.variables <- c('noNbases', 'medianRepTime', 'meanCn', 'dnase', 'highExpGenes', 'lowExpGenes', 
    'staining', 'fragile', 'counts.alu', 'counts.rep','segDup', 'early.repl.x.high.exp')[variables.included ]
model.variable.names <- c('N bases', 'early replication', 'copy number', 'DNAse', 'highly expressed genes', 'lowly expressed genes', 'chromatin staining', 
'fragile sites', 'ALU elements', 'repeats', 'segmental duplications',  'repl expr interct')[variables.included ]#, 'histone mod H3K4me3 (MCF7)', 'histone mod H3K9me3 (MCF7)', 'histone mod H3K27ac (MCF7)', 'histone mod me1K4 (K562)')
names(model.variable.names) <- paste0(model.variables , '.norm')

bins.regression$early.repl.x.high.exp <- bins.regression$medianRepTime * bins.regression$highExpGenes


# normalize the data
names(model.variable.names) <- paste0(model.variables, '.norm')
for (vi in 1:length(model.variables)) {
    normalisedV <-  normaliseVector(bins.regression[,as.character(model.variables[vi])])
    vn <- paste0(model.variables[vi], '.norm')
    bins.regression[[vn]] <- normalisedV
}

bins.regression.m <- merge(bins.regression, allBins2, by=c('chr', 'chromStart', 'chromEnd'))

# fit the regression model
print('Fitting regression model...')
formula.sv1 <- as.formula(paste('counts.sv1 ~ ' , paste(paste0(model.variables,'.norm'), collapse=' + ')))
fm_nbin.sv1 <- glm.nb(formula.sv1 , data = bins.regression.m ) # negative binomial fit

formula.sv3 <- as.formula(paste('counts.sv3 ~ ' , paste(paste0(model.variables,'.norm'), collapse=' + ')))
fm_nbin.sv3 <- glm.nb(formula.sv3 , data = bins.regression.m ) # negative binomial fit


formula.random.sv1 <- as.formula(paste('counts.random.sv1 ~ ' , paste(paste0(model.variables,'.norm'), collapse=' + ')))
fm_nbin.random.sv1 <- glm.nb(formula.random.sv1 , data = bins.regression.m) # negative binomial fit

formula.random.sv3 <- as.formula(paste('counts.random.sv3 ~ ' , paste(paste0(model.variables,'.norm'), collapse=' + ')))
fm_nbin.random.sv3 <- glm.nb(formula.random.sv3 , data = bins.regression.m) # negative binomial fit

# make plots
pdf(paste0('../data/regressionCoefficients-bin-',binSize,'.pdf'), width=12, height=15)
layout(matrix(c(1,2,3,4, 5, 6), 3, 2, byrow = TRUE))
par(mar=c(3,13,3,3))

c <- plotCoefficients(fm_nbin.sv1, nbfit.ref=NULL, 'RS1 breast', model.variable.names)
c <- plotCoefficients(fm_nbin.random.sv1, nbfit.ref=NULL, 'Scattered RS1 breast', model.variable.names)
c <- plotCoefficients(fm_nbin.sv3, nbfit.ref=NULL, 'RS3 breast', model.variable.names)
c <- plotCoefficients(fm_nbin.random.sv3, nbfit.ref=NULL, 'Scattered RS3 breast', model.variable.names)
dev.off()

bins.regression.m.sel <- subset(bins.regression.m, rowSums(is.na(bins.regression.m))==0)
save(bins.regression.m.sel, file='../data/bins.regression.m.sel.RData')


doInteractions <- FALSE
model <- randomForest(formula.sv1, 
                          data = bins.regression.m.sel, importance=TRUE, localImp = TRUE) 

# https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html

if (doInteractions) {
    pdf('interact.pdf')
    plot_predict_interaction(model, bins.regression.m.sel, "medianRepTime.norm", "highExpGenes.norm", grid=20)
    dev.off()

    qcut = function(x, n) {
      quantiles = seq(0, 1, length.out = n+1)
      cutpoints = unname(quantile(x, quantiles, na.rm = TRUE))
      as.character(cut(x, cutpoints, include.lowest = TRUE, labels=FALSE))
    }

    bins.regression.m.sel$medianRepTime.quantile<- qcut(bins.regression.m.sel$medianRepTime, 10)
    bins.regression.m.sel$highExpGenes.quantile<- qcut(bins.regression.m.sel$highExpGenes, 10)    
    }

