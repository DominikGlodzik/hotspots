plotHotspots <- function(hotspots_toplot, sv, window = 500000, gene_path = NULL, enhancer_path = NULL, driver_path = NULL){
    "
    Utility function to visualize the SV breakpoint distribution in hotspots
    Mandatory input:
        hotspots_toplot: hotspot object from runLightPcf
        sv: all the SV breakpoints to plot. Mandatory columns: sample, chrom1 (1, 2,.. X), pos1, chrom2, pos2, sv_pcawg
            sv_pcawg should contain sv class annotations in the following categories: 
            DEL, DUP, INV, TRA, Templated_insertion, Chromothripsis, Chromoplexy, Complex
            window: number of bases to plot on either side of the hotspot margins

    Optional annotation tracks:
        gene_path: path to gene reference with the following mandatory columns: 
            hg19.knownGene.chrom (chr1, chr2, ... chrX), hg19.knownGene.txStart (start pos), 
            hg19.knownGene.txEnd (end pos), hg19.kgXref.geneSymbol (gene name to plot)
        enhancer_path: path to enhancer bed file with the following columns (without headers):
            chromosome (chr1, chr2, ... chrX), start position, end position
        driver_path: path to driver gene list. Mandatory columns: 
            gene (driver gene names corresponding to the names in gene_path)
    "
    
    if(is.null(enhancer_path)){
        enhancer <- data.frame(chr = NA, start = NA, end = NA, driver = NA)
    } else {
        ## Import enhancers
        enhancer <- read.table(enhancer_path, stringsAsFactors = F, header = F, sep = '\t')
        names(enhancer) <- c('chr', 'start', 'end')
        enhancer$driver <- 'Enhancer'
    }
    
    ## Setup genes
    if(is.null(gene_path)){
        gene <- data.frame(chr = NA, start_gene = NA, end_gene = NA, gene = NA, driver = NA)
    } else {
        gene<- read.delim(gene_path, stringsAsFactors = F)
        gene<- select(gene, hg19.knownGene.chrom, hg19.knownGene.txStart, hg19.knownGene.txEnd, hg19.kgXref.geneSymbol)
        names(gene)<-c("chr", "start_gene", "end_gene", "gene")
        
        gene <- gene %>%
            distinct() %>%
            group_by(gene) %>%
            filter((end_gene-start_gene) == max(end_gene-start_gene)) %>%
            ungroup()
    }
    
    ## Setup drivers
    if(is.null(driver_path)){
        gene$driver <- "Unknown"
    } else {
        drivers <- read.delim(driver_path, stringsAsFactors = F)
        gene <- mutate(gene, driver = ifelse(gene %in% drivers$gene, "Driver/COSMIC", "Unknown"))
    }
    
    # Palettes
    
    scale_fill_pcawg <- function(...){
        ggplot2:::manual_scale('fill', 
                               values = setNames(c('#FF7F00', 'black', '#377EB8', '#4DAF4A', '#E41A1C', '#A9A9A9', 
                                                   '#984EA3', '#add8e6', '#A65628'),
                                                 c('Complex', 'TRA', 'INV', 'DUP', 'DEL', 'NO_SV',
                                                   'Chromothripsis', 'Chromoplexy', 'Templated_insertion')), 
                               ...)
    }
    
    scale_col_driver_gray <- function(...){
        ggplot2:::manual_scale('colour', 
                               values = setNames(c('darkgray', 'red', 'brown'),
                                                 c('Unknown', 'Driver/COSMIC', 'Enhancer')), 
                               ...)
    }
    
    ### PREPARE ALL SVS FOR PLOTTING
    sv1 <- sv[,c("sample","chrom1","pos1", 'sv_pcawg')]
    sv2 <- sv[,c("sample","chrom2","pos2", 'sv_pcawg')]
    colnames(sv2) <- colnames(sv1) <- c("sample","chr","pos", 'sv_pcawg')
    sv_comb <- rbind(sv1,sv2) %>%
        mutate(chr = paste0("chr", chr))
    
    ### PREPARE HOTSPOTS FOR PLOTTING
    
    hotspots_toplot <- mutate(hotspots_toplot,
                              chr = paste0('chr', chr),
                              start = start.bp - window, 
                              end = end.bp + window)
    
    sv_hot_expl_fig <- list()
    for(i in 1:nrow(hotspots_toplot)){
        # Defining for each hotspot region the chromosome, start and end position
        c <- as.character(hotspots_toplot$chr[i])
        id <- hotspots_toplot$hotspot.id[i]
        st <- hotspots_toplot$start[i]
        en <- hotspots_toplot$end[i]
        
        # Generating subset
        sub <- sv_comb %>%
            filter(chr == c, pos > st & pos < en)
        
        gene_sub <- gene %>%
            filter(chr == c, end_gene > st & start_gene < en) %>%
            mutate(start_gene = ifelse(start_gene < st, st, start_gene),
                   end_gene = ifelse(end_gene > en, en, end_gene))
        
        enhancer_sub <- enhancer %>%
            filter(chr == c, start > st & end < en)
        
        # Generate plots
        
        p1 <- ggplot()+
            geom_histogram(data = sub, aes(pos/1000000, fill = sv_pcawg), bins = 200)+
            geom_vline(xintercept = (st+window)/1000000, size = 1, linetype = 2)+
            geom_vline(xintercept = (en-window)/1000000, size = 1, linetype = 2)+
            scale_fill_pcawg() +
            scale_col_driver_gray()+
            scale_x_continuous(limits = c(st/1000000, en/1000000))+
            expand_limits(y=-3.5) +
            labs(x = 'Position (Mb)',
                 y = 'Breakpoints (n)',
                 fill = 'SV type',
                 col = 'Annotation',
                 title = 'Distribution of SVs within hotspot region') +
            theme_bw()+
            theme(text = element_text(size = 10),
                  plot.title = element_text(hjust = 0.5))
        
        # If enhancer in region, add annotation track
        if(nrow(enhancer_sub)>0){
            p1 <- p1 + geom_segment(data = enhancer_sub, aes(x=start/1000000, xend=end/1000000, y=-1, yend=-1, col = driver), size=3)
        }
        
        # If genes in region, add gene annotation tracks
        if(nrow(gene_sub)>0){
            p1 <- p1 + geom_segment(data = gene_sub, aes(x=start_gene/1000000, xend=end_gene/1000000, y=-0.5, yend=-0.5, col = driver), size=3)
            
            if(nrow(filter(gene_sub, driver == 'Unknown'))>0){
                p1 <- p1 + geom_text(data = filter(gene_sub, driver == 'Unknown'), mapping = aes(x = (start_gene+end_gene)/2000000, y=-2.5, label = gene, col = driver), size=3, angle = 90)
            }
            if(nrow(filter(gene_sub, driver != 'Unknown'))>0){
                p1 <- p1 + geom_text(data = filter(gene_sub, driver != 'Unknown'), mapping = aes(x = (start_gene+end_gene)/2000000, y=-2.5, label = gene, col = driver), size=3, angle = 90)
            }
        }
        
        # Combine relevant layers of data into final plot for hotspot
        
        pcomb <- suppressWarnings(ggarrange(p1))
        pcomb <- annotate_figure(pcomb, top = text_grob(paste('Hotspot:', id, '|', st/1000000, '-', en/1000000, 'Mb'), 
                                                        face = "bold", size = 14, hjust = 0.5))
        sv_hot_expl_fig[[i]] <- pcomb
    }
    return(sv_hot_expl_fig)
}