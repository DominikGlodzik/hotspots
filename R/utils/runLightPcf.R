runLightPcf <- function (rearr.df,
                         kmin,
                         gamma,
                         rate.factor.thresh=NA,
                         bg.rate=NA,
                         logScale,
                         doMerging,
                         obs.exp.thresh=NA,
                         allBins = NULL,
                         exp.col = NA, # which column of allBins to use as expectation,
                         plot.path='../plots/simulations/'
                         ) {

    pvalue.thresh <- 1 #1e-10
    imd.global <- 1e7

    rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$pf),
                           data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$pf)
                       )
    rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]
    
    
    intermut.dist <- calcIntermutDist(rearr.bps)

    
    kat.regions.list <- list()
    for (loop.chr in unique(intermut.dist$chr)) {
    #    print(loop.chr)
    #for (loop.chr in c('8')) {
        intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)
        pdf(paste0(plot.path, 'RS3.imd.pdf'))
        hist(log10(intermut.dist.chr$distPrev), main='imd RS3 chromosome 1')
        dev.off()
        
        if (logScale) {
            interDist = log10(intermut.dist.chr$distPrev);
           
        } else {
            interDist = intermut.dist.chr$distPrev;
           
        }

        
        sdev <- getMad(interDist,k=25);
        res= exactPcf(interDist,kmin,gamma*sdev,T)# do the segmentation
                                        #

        
        kat.regions.list[[loop.chr]] <- extract.kat.regions2(res,
                                                             subs=subset(rearr.bps, chr==loop.chr),
                                                             kmin.samples=kmin,
                                                             rate.factor.thresh=rate.factor.thresh,
                                                             doMerging=TRUE,
                                                             bp.rate=bg.rate,
                                                             obs.exp.thresh= obs.exp.thresh,
                                                             allBins=allBins,
                                                             exp.col=exp.col
                                                             )

        
    }
    
    kat.regions <- do.call('rbind', kat.regions.list)
    kat.regions
}
