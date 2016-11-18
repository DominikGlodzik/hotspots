source('utils/drawRandomRearrs.R')
source('utils/runLightPcf.R')
source('utils/drawSegmentation.R')
#source('load_dependencies.R')
source('pcf/fastPCF.R')
source('pcf/extract.kat.regions2.R')
source('utils/annotateBgModel.R')
source('../../breast_circo/code/utils/calcIntermutDist.R')

fp <- '/nfs/cancer_archive04/dg17/BASIS/versions/10.06.2016/hotspots/'

gamma <- 8
imd.factor <- 2
kmin.my <- 8

# prepare the observed rearrangements
DATA.VERSION <- '29.11'
RESULT.VERSION <- '29.11'
REARRANGEMENTS.PATH=paste0('/nfs/cancer_archive04/dg17/BASIS/versions/',DATA.VERSION,'/rearrangements')
PRELOAD.REARR.PATH<-paste(REARRANGEMENTS.PATH, '/allRearrangements.RData',sep='')
if (!exists('r')) {
    load(PRELOAD.REARR.PATH)
}
rs1.obs.df <- subset(r$sample.rearrs.all, Signature.SV1>0.5 & is.clustered==FALSE)
rs3.obs.df <- subset(r$sample.rearrs.all, Signature.SV3>0.5 & is.clustered==FALSE)

if (!exists('allBins')) {
    binSize <- 5e5
    load(paste0('../data/regression/bgModel-13.06-',binSize,'.RData')) # allBins, fm_nbin.sv1, fm_nbin.sv3
   
}
# background rates
load(paste0('/nfs/cancer_archive04/dg17/BASIS/versions/','21.01.2016','/rearr.clustering.global', "/rates.RData"))

hotspots.obs.rs1 <- runLightPcf(rearr.df=rs1.obs.df ,
                                kmin=kmin.my,
                                gamma=gamma,
                                rate.factor.thresh=NA,
                                bg.rate=bp.rate.clust.sv[[1]]$gw.rate,
                                logScale=TRUE,
                                doMerging=TRUE,
                                obs.exp.thresh=imd.factor,
                                allBins = allBins,
                                exp.col='exp.nb.sv1')

hotspots.obs.rs1$hotspot.id <- paste0('RS1_chr',hotspots.obs.rs1$chr,'_', round(hotspots.obs.rs1$start.bp/1e6,1), 'Mb' )
drawSegmentation(rs1.obs.df,  hotspots.obs.rs1 , fn=paste0(fp, 'obs.RS1.segmentation.i=',imd.factor , '-g=',gamma,'.pdf'), 'observed RS1')
#hotspots.obs.rs1.annotated <- annotatePcfPeaks(hotspots.obs.rs1, type='RS1',  bps=subset(r$sample.bps.all, Signature.SV1>0.5))
write.csv(hotspots.obs.rs1.annotated, file=paste0(fp, 'obs.RS1.hotspots.i=',imd.factor, '-g=',gamma,'.csv'), row.names = FALSE)

hotspots.obs.rs3 <- runLightPcf(rearr.df=rs3.obs.df ,
                                kmin=kmin.my,
                                gamma=gamma,
                                rate.factor.thresh=NA,
                                bg.rate=bp.rate.clust.sv[[3]]$gw.rate,
                                logScale=TRUE,
                                doMerging=TRUE,
                                obs.exp.thresh=imd.factor,
                                allBins = allBins,
                                exp.col='exp.nb.sv3')
hotspots.obs.rs3$hotspot.id <- paste0('RS3_chr',hotspots.obs.rs3$chr,'_', round(hotspots.obs.rs3$start.bp/1e6,1), 'Mb' )
drawSegmentation(rs3.obs.df,  hotspots.obs.rs3 , fn=paste0(fp, 'obs.R3.segmentation.i=',imd.factor, '-g=',gamma,'.pdf'), 'observed RS3')
#hotspots.obs.rs3.annotated <- annotatePcfPeaks(hotspots.obs.rs3, type='RS3',  bps=subset(r$sample.bps.all, Signature.SV3>0.5))
write.csv(hotspots.obs.rs3.annotated, file=paste0(fp, 'obs.RS3.hotspots.i=',imd.factor, '-g=',gamma,'.csv'), row.names = FALSE)    

