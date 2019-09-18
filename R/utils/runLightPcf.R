runLightPcf <- function (rearr.df,
                         kmin,
                         gamma,
                         rate.factor.thresh=NA,
                         bg.rate=NA,
                         logScale,
                         doMerging,
                         obs.exp.thresh=NA,
                         allBins = NULL,
                         exp.col = NA, 
                         binSize = 5e5,
                         downsampleBins = FALSE,
                         plot.path = '../plots/simulations/'
                         ) {

                                        # Identifies regions of chromosomes with increased frequencies of breakpoints
                                        #
                                        # Args:
                                        #   rearr.df: data frame with one rearrangement per row (column names see below)
                                        #             required columns: sample, Chromosome.1, pos.1, pf, Chromosome.2, pos.2, pf
                                        #   kmin: parameter of PCF algorithm
                                        #   gamma: parameter of PCF
                                        #   rate.factor.thresh: hotspot threshold compare to genome-wide average
                                        #   bg.rate: genome-wide average
                                        #   logScale: TRUE/FALSE whether inter-breakpoint distances are log-transformed
                                        #   doMerging: TRUE/FALSE whether to merge adjacent hotspots
                                        # obs.exp.thresh: hotspots threshold: ratio of observed to expected breakpoint frequency
                                        # allBins: expected breakpoint frequency according to negative binomial regression
                                        # exp.col: which column of allBins to use as expectation
                                        # binSize: bin size for downsampling
                                        # downsampleBins: if true, randomly sample one breakpoint per bin
                                        # plot.path: where to store plots

                                        # make a dataframe with one row per breakpoint (as opposed to one row per rearrangement)
    rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$pf),
                       data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$pf)
                       )
    
    print(paste0("Total number of breakpoints: ", nrow(rearr.bps)))
    
    if(downsampleBins){
        rearr.bps <- downsampleRearrBins(rearr.bps)
        print(paste0("Number of breakpoints after downsampling: ", nrow(rearr.bps)))
    }

    
                                        # order by chromosome and position
    rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]
                                        # calculate distances between neighbouring breakpoints
    intermut.dist <- calcIntermutDist(rearr.bps)
    
                                        # loop over chromosomes
    kat.regions.list <- list()
    for (loop.chr in unique(intermut.dist$chr)) {

        intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)
          
        if (logScale) {
            interDist = log10(intermut.dist.chr$distPrev);
           
        } else {
            interDist = intermut.dist.chr$distPrev;
           
        }
        # call the PCF functions
        sdev <- getMad(interDist,k=25); # 
        res= exactPcf(interDist,kmin,gamma*sdev,T)# PCF call
                                        #
        
                                        # extract hotspots from the PCF segmentation
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
            
                                        # make a dataframe with hotspots and return it
    kat.regions <- do.call('rbind', kat.regions.list)
    kat.regions
}
