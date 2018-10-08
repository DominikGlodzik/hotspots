# annotate rearrangements as to whether they are clustered or not
.libPaths(c("/software/R-3.0.0/lib/R/library", '/nfs/users/nfs_d/dg17/R/circos-library/3.0'))

source('utils/calcIntermutDist.R')
source('utils/plotScatterCirco.R')
source('utils/rearrangement.clustering.demo.R')
source('utils/fastPCF.R')
source('utils/extract.kat.regions.R')


library(VariantAnnotation) # for checking if genome regions are mappable
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)

sh <- read.table('../../sign-tracker/data/signatures.txt', header=TRUE, sep='\t')
mut.order <- (sh[,'Somatic.Mutation.Type'])

args <- commandArgs(trailingOnly = TRUE)
if (length(args)>0) {
    PROJECT.ID <- as.integer(args[[1]])
    SAMPLE.ID <- args[[2]]
    output.path <- args[[3]]
    output.path.subs <- args[[4]]
    sample.path <- args[[5]]
    onlyAssembled <- as.logical(as.character(args[6]))
    LOAD.DATA <- TRUE
    print(paste('arguments: ', output.path, output.path.subs, sample.path))
    print((sample.path!='NA' || is.na(sample.path)))
    if (length(args)>6) {
        SAMPLE.PANCAN.ID <- args[[7]]
    } else {
        SAMPLE.PANCAN.ID <- NA
    }      
} 

cat(paste('preparing rearrangements ... \n'))


rearrs <- read.table(sample.path, header=TRUE)
                                        # prepare a data frame with rearrangements
rearrs$Chromosome.1 <- rearrs$chromosome.1
rearrs$Chromosome.2 <- rearrs$chromosome.2
rearrs$pos.1 <- rearrs$pos.1.min
rearrs$pos.2 <- rearrs$pos.2.min
rearrs$id <- rearrs$ID
rearrs$pf <- rearrs$svclass
sample.rearrs <- rearrs

                                        # prepare a data frame with breakpoints
rearrs.left <- rearrs[,c('chromosome.1','pos.1.min', 'svclass','pf', 'sample', 'ID','nts', 'mh')]; names(rearrs.left ) <- NA
rearrs.right <- rearrs[,c('chromosome.2','pos.2.min', 'svclass', 'pf', 'sample', 'ID','nts', 'mh')];names(rearrs.right ) <- NA
rearrs.cncd <- rbind(rearrs.left , rearrs.right  );
rearrs.cncd$isLeft <- c(rep(TRUE, nrow(rearrs.left)), rep(FALSE, nrow(rearrs.left)))
colnames(rearrs.cncd) <- c('chr', 'position', 'svclass', 'pf', 'sample', 'id', 'nts', 'mh', 'isLeft')
sample.bps <- rearrs.cncd


   
                                        # annotate each rearrangement/bp whether it is clustered
cat('pcf: rearrangements \n')
plot.path <- paste(output.path,'/', SAMPLE.ID, sep='') #  set tp NA of plot not needed
if (nrow(sample.bps)>0) { # if there are any rearrangements
    
                                        # annotate each rearrangement whether it falls into repeats
    #sample.rearrs <- annotateRearrsMh(sample.rearrs, metadata)
    sample.bps <- sample.bps[order(sample.bps$chr, sample.bps$position),]
    
    clustering.result <- rearrangement.clustering(sample.bps,
                                                  plot.path = paste0(output.path, SAMPLE.ID,'.png'),
                                                  kmin=10,                                                  
                                                  kmin.samples=1,
                                                  gamma.sdev=25,
                                                  PEAK.FACTOR=10,
                                                  thresh.dist=NA)
    
    sample.bps <- clustering.result$sample.bps
    # mark both breakpoints of a rearrangement as clustered if any is
    sample.rearrs$is.clustered <- sample.rearrs$id %in% sample.bps$id[sample.bps$is.clustered.single]
        
}

result.path <-  paste(output.path, '/', SAMPLE.ID, '.RData', sep='')
cat(paste('Writing prepared rearrangements to ', result.path))
save(sample.bps, sample.rearrs,  file=result.path)

write.table( sample.rearrs, file=paste0(sample.path,'.prepared'), sep='\t', row.names=FALSE, quote = FALSE,
                 col.names = TRUE)

