plotScatterCirco <- function(pos, plot.path) {

    png(file=plot.path , height=3000, width=3000, res=550)


    data(UCSC.HG19.Human.CytoBandIdeogram);
    hg19.cyto <- UCSC.HG19.Human.CytoBandIdeogram;
    RCircos.Set.Core.Components(cyto.info=hg19.cyto, chr.exclude=NULL, 	tracks.inside=10, tracks.outside=1);

    params <- RCircos.Get.Plot.Parameters();
    params[names(params.my)] <- params.my
    RCircos.Reset.Plot.Parameters(params)
    
    RCircos.Set.Plot.Area();
    RCircos.Chromosome.Ideogram.Plot.my();
    params <- RCircos.Get.Plot.Parameters();
    params$sub.tracks <- 1
    RCircos.Reset.Plot.Parameters(params)


    
    rearrs.data <- data.frame(Chromosome=pos$chr, chromStart=pos$position, chromEnd=pos$position, pf=pos$pf, intermut.dist=pos$intermut.dist, mean.intermut.dist=pos$mean.intermut.dist, is.clustered=pos$is.clustered)
    
    rearrs.data$Chromosome <- as.character(rearrs.data$Chromosome)
    rearrs.data$Chromosome[rearrs.data$Chromosome=='23'] <- 'X'
    rearrs.data$Chromosome[rearrs.data$Chromosome=='24'] <- 'Y'
    
    
    if (nrow(rearrs.data )>0) {
        
        data.ceiling <- max(rearrs.data$intermut.dist)

        if (sum(!rearrs.data$is.clustered)>0) {
            RCircos.Scatter.Plot.color(rearrs.data[!rearrs.data$is.clustered,], 5, 1, 'in', scatter.colors=rep('black', sum(!rearrs.data$is.clustered)), draw.bg =TRUE, draw.scale=TRUE, no.sort=FALSE,  data.ceiling = data.ceiling )
        }
        if (sum(rearrs.data$is.clustered)>0) {
            RCircos.Scatter.Plot.color(rearrs.data[rearrs.data$is.clustered,], 5, 1, 'in', scatter.colors=rep('red',  sum(rearrs.data$is.clustered)), draw.bg =FALSE, draw.scale=FALSE, no.sort=FALSE,  data.ceiling = data.ceiling )
    }        
        RCircos.Scatter.Plot.color(rearrs.data, 6, 1, 'in', scatter.colors=rep('grey', nrow(rearrs.data)), draw.bg =FALSE, draw.scale=FALSE, no.sort=FALSE,  data.ceiling = data.ceiling) 
        
    }
        
    
       dev.off()
}

.libPaths('/nfs/users/nfs_d/dg17/R/x86_64-unknown-linux-gnu-library/3.0')
library(RCircos);
library(scales)

params.my <- list()
params.my$track.background <- 'white'
params.my$highlight.width <- 0.2
params.my$point.size <- 0.01
params.my$point.type <- 19
params.my$radius.len <- 3
params.my$chr.ideog.pos <- 3.2
params.my$highlight.pos <- 3.35
params.my$chr.name.pos <- 3.45
params.my$plot.radius <- 3.55
params.my$track.in.start <- 3.05
params.my$track.out.start <- 3.2
params.my$text.size <-  0.6
params.my$track.padding <- 0.07
params.my$grid.line.color <- 'lightgrey'
params.my$track.heights <- c(0.7, 0.4, 0.4, 0.4, 100)
params.my$track.height <- 0.1
params.my$sub.tracks <- 1
params.my$heatmap.cols <- c(alpha('lightcoral', 1), alpha('lightcoral', 0.5), alpha('lightgrey',0.10), alpha('olivedrab2', 0.3),  alpha('olivedrab2', 0.5), alpha('olivedrab2',.7), alpha('olivedrab2', 0.75),
                         alpha('olivedrab3', 0.9), alpha('olivedrab4', 0.9))
params.my$heatmap.ranges <- c(0,1,2,4,8,16, 32,64,1000)



setwd('../../breast_circo/code/')
source('utils/RCircos.Scatter.Plot.color.R')
source('utils/my.RCircos.Tile.Plot.R')
source('utils/calcIntermutDist.R')
source('utils/RCircos.Get.Plot.Data.nosort.R')
source('utils/RCircos.Chromosome.Ideogram.Plot.my.R')
source('utils/RCircos.Track.Outline.my.R')
source('utils/RCircos.Line.Plot.my.R')
source('utils/RCircos.Link.Plot.my.R')
source('utils/RCircos.Track.Positions.my.R')
source('utils/RCircos.Gene.Connector.Plot.my.R')
source('utils/RCircos.Heatmap.Plot.my.R')
source('utils/RCircos.Validate.Genomic.Data.my.R')
source('utils/toPyr.R')
source('utils/processSubs.R')
source('utils/RCircos.Gene.Name.Plot.my.R')
source('../../sign-tracker/code/utils/merge.with.order.R')
source('utils/RCircos.Get.Plot.Data.segment.R')
source('utils/RCircos.Scatter.Plot.cn.R')
source('utils/RCircos.Scatter.Plot.ra.R')
setwd('../../breast_rearr/code/')
