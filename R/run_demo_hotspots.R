# find hotspots of rearrangements attributed to signatures RS1 and RS3

source('utils/drawRandomRearrs.R')
source('utils/runLightPcf.R')
source('utils/drawSegmentation.R')
source('utils/fastPCF.R')
source('utils/extract.kat.regions2.R')
source('utils/annotateBgModel.R')
source('utils/calcIntermutDist.R')

# where data is stored and output written
fp <- '../data/'

# parameters of the PCF algorithm
gamma <- 8
imd.factor <- 2
kmin.my <- 8

# prepare the observed rearrangements
PRELOAD.REARR.PATH<-paste0(fp,'allRearrangements.29.11.RData')
load(PRELOAD.REARR.PATH)

# identify RS1 and RS3 rearrangements
rs1.obs.df <- subset(r$sample.rearrs.all, Signature.SV1>0.5 & is.clustered==FALSE) # select rearrangements confidently attributed to signature RS1
rs3.obs.df <- subset(r$sample.rearrs.all, Signature.SV3>0.5 & is.clustered==FALSE) # select rearrangements confidently attributed to signature RS3

# load background model: negative binomial regression
binSize <- 5e5 # bin size: 0.5 mb
load(paste0(fp,'bgModel-13.06-',binSize,'.RData')) # allBins, with columns fm_nbin.sv1, fm_nbin.sv3 

# background breakpoint rates per signature, bp.rate.clust.sv
load(paste0(fp, 'rates.RData'))

# Hotspots - RS1 signature
hotspots.obs.rs1 <- runLightPcf(rearr.df=rs1.obs.df, # RS1 rearrangements 
                                kmin=kmin.my,# kmin parameter for PCF
                                gamma=gamma, # gamma parameter for PCF
                                bg.rate=bp.rate.clust.sv[[1]]$gw.rate, # genome-wide rate of RS1 breakpoints
                                logScale=TRUE, # log-scale transformation of inter-mutation disances
                                doMerging=TRUE, # merge adjacent hotspots
                                obs.exp.thresh=imd.factor, # paramter for PCF
                                allBins = allBins, # expected breakpoints/bin
                                exp.col='exp.nb.sv1', # use the RS1 column from allBins
                                plot.path='../data/' # where to save a plot
                                )

hotspots.obs.rs1$hotspot.id <- paste0('RS1_chr',hotspots.obs.rs1$chr,'_', round(hotspots.obs.rs1$start.bp/1e6,1), 'Mb' ) # give IDs to each hotspots
drawSegmentation(rs1.obs.df,  hotspots.obs.rs1 , fn=paste0(fp, 'obs.RS1.segmentation.i=',imd.factor , '-g=',gamma,'.pdf'), 'observed RS1') # draw the rainfall plots 
write.csv(hotspots.obs.rs1, file=paste0(fp, 'obs.RS1.hotspots.i=',imd.factor, '-g=',gamma,'.csv'), row.names = FALSE) # save table with hotspots

# Hotspots - RS3 signature
hotspots.obs.rs3 <- runLightPcf(rearr.df=rs3.obs.df , # RS3 rearrangements 
                                kmin=kmin.my,
                                gamma=gamma,
                                bg.rate=bp.rate.clust.sv[[3]]$gw.rate,
                                logScale=TRUE,
                                doMerging=TRUE,
                                obs.exp.thresh=imd.factor,
                                allBins = allBins,
                                exp.col='exp.nb.sv3',
                                plot.path='../data/')

hotspots.obs.rs3$hotspot.id <- paste0('RS3_chr',hotspots.obs.rs3$chr,'_', round(hotspots.obs.rs3$start.bp/1e6,1), 'Mb' ) # give IDs to each hotspots
drawSegmentation(rs3.obs.df,  hotspots.obs.rs3 , fn=paste0(fp, 'obs.R3.segmentation.i=',imd.factor, '-g=',gamma,'.pdf'), 'observed RS3') # draw the rainfall plots 
write.csv(hotspots.obs.rs3, file=paste0(fp, 'obs.RS3.hotspots.i=',imd.factor, '-g=',gamma,'.csv'), row.names = FALSE) # save table with hotspots    

