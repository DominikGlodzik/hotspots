source('utils/assignPvalues.R')

extract.kat.regions2 <- function (res, # segmentation result
                                  subs, # breakpoints (chr, position)
                                  kmin.samples=10, # minimum number of samples to dorm a hotspot
                                  rate.factor.thresh=1,
                                  doMerging=FALSE,
                                  bp.rate=NA, # baseline rate
                                  obs.exp.thresh=NA,
                                  allBins = NULL,
                                  exp.col = NA # which column of allBins to use as expectation
                                  ) {

                                        # Extracts hotspots from a PCF segmentation.
                                        # Args:
                                        #   res: PCF segmentation
                                        #   kmin.samples: minimum number of samples for a hotspots
                                        #   doMerging: if merining neighbouring hotspots
                                        #   bp.rate: baseline rate
                                        #   obs.exp.thresh: threshold - observed to expected breakpoints
                                        #   allBins: expected number of breakpoints per bin according to negative binomial
                                        #   exp.col: relevant column of the regression data frame
    
	segInterDist <-  res$yhat
	kat.regions.all = data.frame()	
	chr <- as.character(subs$chr[1])
	positions <- subs$pos

        # positions where segmentation changes: possible starts of hotspots
        start.regions = which(segInterDist[-1] != segInterDist[-length(segInterDist)] ) +1 # endpoints between peaks
        start.regions <- c(1, start.regions)

         # positions where segmentation changes: possible ends of hotspots
        end.regions = which( segInterDist[-1] != segInterDist[-length(segInterDist)] )
        end.regions <- c( end.regions, length(segInterDist))

        
        start.regions.init <- start.regions
        end.regions.init <- end.regions 
        
        if (length(start.regions)!=length(end.regions)) {
            print('starts dont have the same length as ends')
            browser()
        }
        
                                # prepare a data structure that will be later filled up
        kat.regions.all <- data.frame(
            chr=subs$chr[1],
            start.bp=rep(NA,length(start.regions)), # start coordinate [bp]
            end.bp=rep(NA,length(start.regions)), # end coordinate [bp]
            length.bp=rep(NA,length(start.regions)), # length [bp]
            number.bps=rep(NA,length(start.regions)),
            number.bps.clustered=rep(NA,length(start.regions)),
            avgDist.bp=rep(NA,length(start.regions)),
            no.samples=rep(NA,length(start.regions)),
            no.del =rep(NA,length(start.regions)),
            no.dup =rep(NA,length(start.regions)),
            no.inv= rep(NA,length(start.regions)),
            no.trn = rep(NA,length(start.regions)),
            firstBp=start.regions,
            lastBp=end.regions                                    )
        
        kat.regions.all <- hotspotInfo(kat.regions.all, subs, segInterDist)
        
        step.segInterDist.left <- rep(NA, length(segInterDist))
        step.segInterDist.left[2:length(segInterDist)] <- segInterDist[2:length(segInterDist)]- segInterDist[1:(length(segInterDist)-1)]       
        step.segInterDist.right <- rep(NA, length(segInterDist))
        step.segInterDist.right[1:(length(segInterDist)-1)] <- segInterDist[1:(length(segInterDist)-1)]- segInterDist[2:(length(segInterDist))]          
        kat.regions.all$step.left <-  step.segInterDist.left[start.regions]
        kat.regions.all$step.right <-  step.segInterDist.right[end.regions]
            
                                       
        if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
            kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate) # assign side info to the hotspots

            if (!is.na(rate.factor.thresh)) {
                                        # only keep the hotspots that exceed the theshold
                kat.regions.all <- subset(kat.regions.all, rate.factor>=rate.factor.thresh)
            }

            if (!is.na(obs.exp.thresh) & !is.null(allBins)) {
                kat.regions.all <- annotateBgModel(kat.regions.all, allBins, exp.col)
                kat.regions.all <- subset(kat.regions.all, d.obs.exp>obs.exp.thresh)
            }

            if (nrow(kat.regions.all)>0) { # if there are any hotspots left 

                # merging of adjacent hotspots
                if (doMerging) {
                    if(nrow(kat.regions.all)>1){
                        for(r in 2:nrow(kat.regions.all)){
                            if (kat.regions.all$lastBp[r-1] == (kat.regions.all$firstBp[r]-1)) {
                                        # merge two segments
                                kat.regions.all$firstBp[r] <- kat.regions.all$firstBp[r-1]
                                kat.regions.all$firstBp[r-1] <- NA
                                kat.regions.all$lastBp[r-1] <- NA
                                kat.regions.all$avgDist.bp[r] <- NA # this will need to be updated as segments are being merged
                            }
                        }        
                    }
                                        # remove some of the merged segments
                    kat.regions.all <- subset(kat.regions.all, !is.na(firstBp) & !is.na(lastBp))
                    
                                        # update the info on hotspots that might have changed when they were merged
                    kat.regions.all <- hotspotInfo( kat.regions.all ,  subs, segInterDist)
                    kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate)

                    if (!is.null(allBins)) {
                        kat.regions.all <- annotateBgModel(kat.regions.all, allBins, exp.col=exp.col)
                    }
                } # end merging
            }

            
        }  

        # make sure there are enough samples in each hotspot
        if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
            kat.regions.all <- subset(kat.regions.all, no.samples>=kmin.samples)
        }
        

	kat.regions.all
        
    }

# Annotates the hotspots: cooridnates in chromosome, length, how many rearrangements, what type and from which samples
hotspotInfo <- function(kat.regions.all, subs, segInterDist=c()) {
    if(nrow(kat.regions.all)>0){
        for(r in 1:nrow(kat.regions.all)){
            
                                        # indices of the breakpoints in the hotspot
            subs.hotspot <-subs[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r],]
            
            kat.regions.all[r,'start.bp'] <- min(subs.hotspot$pos)
            kat.regions.all[r,'end.bp'] <- max(subs.hotspot$pos)
            kat.regions.all[r,'length.bp'] <-  kat.regions.all[r,'end.bp'] - kat.regions.all[r,'start.bp'] 
            kat.regions.all[r,'number.bps'] <- nrow(subs.hotspot)
            kat.regions.all[r,'number.bps.clustered'] <- sum(subs.hotspot$is.clustered)
            
            if (length(segInterDist)>0 & is.na(kat.regions.all[r,'avgDist.bp'])) {
                kat.regions.all[r,'avgDist.bp'] <- mean(segInterDist[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r]])
            }
            kat.regions.all[r,'no.samples'] <- length(unique(subs.hotspot$sample))

            if ('pf' %in% colnames(subs.hotspot)){
                kat.regions.all[r,'no.del'] <- nrow(subset(subs.hotspot, pf==2))
                kat.regions.all[r,'no.dup'] <- nrow(subset(subs.hotspot, pf==4))
                kat.regions.all[r,'no.inv'] <- nrow(subset(subs.hotspot, pf==1 | pf==8))
                kat.regions.all[r,'no.trn'] <- nrow(subset(subs.hotspot, pf==32))
            }
            
        } # for all peaks
    } # if there is at least one peak
    kat.regions.all
}
