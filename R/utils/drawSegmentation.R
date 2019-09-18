drawSegmentation <- function (rearr.df,  
                              hotspots.df,  
                              fn, 
                              main.text, 
                              chroms=c(as.character(1:22), 'X'), 
                              xWindow=NULL, 
                              bp.col='aquamarine3', 
                              annot=NULL, 
                              hotspots=NULL
                              ) {

                                        # Draw the rainfall plots for rearrangement breakpoints
                                        #   rearr.df: data frame with rearrangements, one rearrangement per row
                                        #   hotspots.df: hotspots identified
                                        #   fn: file for saving PDF
                                        #   main.text: tile for PDF
    
    # make a data frame: one breakpoint per line 
    rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$pf),
                       data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$pf)
                       )
    
    rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]
    
    intermut.dist <- calcIntermutDist(rearr.bps)

    pdf(fn, width=20, height=5, useDingbats=FALSE)
     
    for (loop.chr in chroms) {
        intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)
        
        interDist = log10(intermut.dist.chr$distPrev);

        if (!is.null(bp.col)) {
            intermut.dist.chr$col<- bp.col
        } else {
            intermut.dist.chr$col <- NA
            
            
            intermut.dist.chr$col[intermut.dist.chr$pf==1] <- 'dodgerblue2'
            intermut.dist.chr$col[intermut.dist.chr$pf==2] <- 'coral2'
            intermut.dist.chr$col[intermut.dist.chr$pf==4] <- 'olivedrab3'
            intermut.dist.chr$col[intermut.dist.chr$pf==8] <- 'dodgerblue2'
            intermut.dist.chr$col[intermut.dist.chr$pf==32] <- 'gray35'
        }
        
        if (is.null(xWindow)) {
           
            plot(intermut.dist.chr$position, interDist, main=paste(main.text, 'chromosome', loop.chr), pch=19, col=intermut.dist.chr$col, cex=0.5, xlab='chromosome position', ylab='inter-breakpoint distance', ylim=c(0,7))     
        } else {
             browser()
            plot(intermut.dist.chr$position, interDist, main=paste(main.text, 'chromosome', loop.chr), pch=19, col=intermut.dist.chr$col, cex=0.5, xlim=xWindow, xlab='chromosome position', ylab='inter-breakpoint distance', ylim=c(0,7))     
        }
        if (nrow(hotspots.df)>0) {  
            hotspots.chr <- subset(hotspots.df, chr==loop.chr )
            if (nrow(hotspots.chr)>0) {             
                segments(x0=hotspots.chr$start.bp, y0=hotspots.chr$avgDist.bp,x1=hotspots.chr$end.bp, y1=hotspots.chr$avgDist.bp, lwd=2 )
            }
        }
        

         
    }


                                            # Annotations on top
     if (!is.null(annot) ) {

         arrow.start<-vector()
         arrow.end <- vector()
         arrow.start[annot$strand==1]<- annot[annot$strand==1, 2]
         arrow.end[annot$strand==1] <- annot[annot$strand==1, 3]
         arrow.start[annot$strand==-1]<- annot[annot$strand==-1, 3]
         arrow.end[annot$strand==-1] <- annot[annot$strand==-1, 2]
         
        arrows(
            x0 = arrow.start,
            y0 = annot[, 'height'] ,#+ (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            x1 =  arrow.end,
            lwd = 1,
            col=annot[, 'col'],
            length=.1
        )
        
        text(
            rowMeans(annot[, 2:3]),
            annot[, 'height'], #+ (0:(sum(annot[,1] %in% chrs)-1)) * 0.05 * yrange_size,
            labels = annot[,4],
            cex=1,
            col=annot[, 'col'],
            pos=3
        )
     }

    if (!is.null(hotspots)) {
        segments(
            x0 = hotspots$start.bp,
            y0 = 0 ,
            x1 =  hotspots$end.bp,
            lwd = 2,
            col='black'
        )
    }

    
    dev.off()

}
