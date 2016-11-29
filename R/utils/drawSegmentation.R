 drawSegmentation <- function (rearr.df,  hotspots.df,  fn, main.text) {

     rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$pf),
                        data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$pf)
                        )
     rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]
     

     intermut.dist <- calcIntermutDist(rearr.bps)

     pdf(fn, width=20, height=5)
     
     for (loop.chr in as.character(1:22)) {
         intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)
         
         interDist = log10(intermut.dist.chr$distPrev);


         plot(intermut.dist.chr$position, interDist, main=paste(main.text, 'chromosome', loop.chr), pch=19, col='aquamarine3', cex=0.5)


         hotspots.chr <- subset(hotspots.df, chr==loop.chr )

         if (nrow(hotspots.chr)>0) {             
             segments(x0=hotspots.chr$start.bp, y0=hotspots.chr$avgDist.bp,x1=hotspots.chr$end.bp, y1=hotspots.chr$avgDist.bp, lwd=2 )
         }
         
     }

     dev.off()


}
