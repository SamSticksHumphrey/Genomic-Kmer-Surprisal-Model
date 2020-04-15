

library(plyr)
library(data.table)
library(tidyverse)
library(purrr)
library(gplots)
library(ggpubr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(Logolas)
library(cowplot)
library(matrixStats)


n <- 16  
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]

CRUKcolours <- c(rgb(red = 200, green = 201, blue = 199, alpha = 255, names = "CRUK_Grey", maxColorValue = 255), 
                  rgb(red = 0, green = 182, blue = 237, alpha = 255, names = "CRUK_LightBlue", maxColorValue = 255),
                  rgb(red = 236, green = 0, blue = 140, alpha = 255, names = "CRUK_Magenta", maxColorValue = 255), 
                  rgb(red = 46, green = 0, blue = 139, alpha = 255, names = "CRUK_Blue", maxColorValue = 255), 
                  rgb(red = 0, green = 0, blue = 0, alpha = 255, names = "CRUK_Black", maxColorValue = 255))

mytheme <- theme(
      			line = element_line(colour = "black", size = 1, lineend = "square"), 
            axis.line = element_line(colour = "black", size = 1, lineend = "square"), 
      			text = element_text(family = "Helvetica", face = "bold", colour = "black", size = 14, hjust = 0.5, vjust = 0),
      			panel.grid.major = element_blank(), 
      			panel.grid.minor = element_blank(), 
      			panel.background = element_blank(), 
      			legend.position='none')

dinucleotides <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

aminoAcids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

trinucleotides <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                    "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                    "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

lowLimit <- function(x){floor(min(x))}
highLimit <- function(x){ceiling(max(x))}

get_plot_limits <- function(plot) {
    gb = ggplot_build(plot)
    xmin = gb$layout$panel_params[[1]]$x.range[1]
    xmax = gb$layout$panel_params[[1]]$x.range[2]
    ymin = gb$layout$panel_params[[1]]$y.range[1]
    ymax = gb$layout$panel_params[[1]]$y.range[2]
    list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

# ---------------------------------------------------------------
# plotSpectra takes in a distribution of the form {k-mer occurance}, a
#	name for the distribution and a numeric depth to evaluate the spectra to.  
#
# It returns a figure of the k-mer occurance on the x-axis and the frequency
#	that they occur on the y.  
# ---------------------------------------------------------------
plotDistributionSpectra <- function(dist, distName, spectraDepth) {
  
    data <- dist
    colnames(data) <- c("seq", "occurance")
    data <- subset(data, occurance > 0)
    occurSpectra <- as.data.frame(table(data$occurance))
    colnames(occurSpectra) <- c("occuranceNames", "freqSpectra")
    occurSpectra$occuranceNames <- as.numeric(occurSpectra$occuranceNames)
    occurSpectra <- subset(occurSpectra, occuranceNames <= spectraDepth)

    spectraFigure <- ggplot(data = occurSpectra, aes(x = occuranceNames, y = freqSpectra)) +
        geom_bar(stat = "identity") + mytheme + 
        labs( x = "",  y = distName) +
        scale_x_continuous(expand=c(0.01,0)) +
        scale_y_continuous(expand=c(0.01,0))

    return(spectraFigure)
}

# ---------------------------------------------------------------
# plotSpectraDinucleotideSplit takes the spectra as defined above 
#	and plots it individually for each di-nucleotide.
# ---------------------------------------------------------------
plotSpectraDinucleotideSplit <- function(dist, kmerDinucs, distName, spectraDepth){
  
    data <- dist
    colnames(data) <- c("kmer", "occurance")
    data <- subset(data, occurance > 0)
    data <- subset(data, occurance <= spectraDepth)

    dinuc.levels <- data.frame( kmer = kmerDinucs[,1], as.data.frame(sapply(kmerDinucs[,-1], function(x) factor(ifelse(x >= 3, "3+", x), levels = c("0", "1", "2", "3+")))) )
    data.dinuc <- merge(data, dinuc.levels, by = "kmer", all.x = TRUE, all.y = FALSE)

    data.dinuc$kmer <- NULL
    data.melted <- melt(data.dinuc, id.var = "occurance")
    colnames(data.melted) <- c("occurance", "dinuc", "dinucGroup")
    
    dinucKmer <- ddply(data.melted, .(dinuc, dinucGroup, occurance), nrow)
    plot.data <- subset(dinucKmer, occurance <= spectraDepth )
    
    multiDinucSpectra <- ggplot(plot.data, aes(occurance, V1)) + 
                          geom_bar(stat = "identity", aes(fill = factor(dinucGroup))) + 
                          mytheme + labs(x = paste0(params$k, "-mer Occurance"), y = paste0( distName, " Occurance Frequency")) +
                          scale_fill_manual(values = c(as.vector(CRUKcolours[1:4]), as.vector(CRUKcolours[1]))) + facet_wrap(~dinuc, ncol=4)
    
    return(multiDinucSpectra)
}    

# ---------------------------------------------------------------
# plotSpectraDinucleotideSplit takes the spectra as defined above 
# and plots it individually for each di-nucleotide.
# ---------------------------------------------------------------
plotSpectraSplitOneDinuc <- function(dist, dinuc, kmerDinucs, seqType, spectraDepth) {
  
  data <- dist
  colnames(data) <- c("kmer", "occurance")
  data <- subset(data, occurance > 0)
  if(!is.na(spectraDepth)){data <- subset(data, occurance <= spectraDepth)}
  
  kmerDinuc <- kmerDinucs[, c("kmer", dinuc)]
  
  data.dinuc <- merge(data, kmerDinuc, by = "kmer", all.x = TRUE, all.y = FALSE)
  data.dinuc$kmer <- NULL
  colnames(data.dinuc) <- c("occurance", "dinucGroup")
  
  dinucKmer <- ddply(data.dinuc, .(dinucGroup, occurance), nrow)
  
  spectraFigure <- ggplot(dinucKmer, aes(occurance, V1)) + 
    geom_bar(stat = "identity", alpha = 0.4, aes(fill = factor(dinucGroup))) + 
    geom_smooth(method = "loess", se = FALSE, span = 0.2, aes(group = dinucGroup, colour = factor(dinucGroup))) + 
    mytheme + 
    labs(x = paste0(k, "-mer Occurance"), y = paste(seqType, "Occurance Frequency")) +
    scale_x_continuous(expand=c(0.01,0), limits = c(0, spectraDepth)) +
    scale_y_continuous(expand=c(0.01,0), limits = c(0, NA)) + 
    scale_fill_manual(values = as.vector(CRUKcolours[1:4])) + 
    scale_color_manual(values = as.vector(CRUKcolours[1:4]))
  
  return(spectraFigure)
} 

# ---------------------------------------------------------------
# plotSpectraDinucleotideSplit takes the spectra as defined above 
# and plots it individually for each di-nucleotide.
# ---------------------------------------------------------------
plotSpectraAminoAcidSplit <- function(dist, distName, spectraDepth){
  
    data <- dist
    colnames(data) <- c("kmer", "occurance")
    data <- subset(data, occurance > 0)
    if(!is.na(spectraDepth)){data <- subset(data, occurance <= spectraDepth)}

    aa.levels <- data.frame( occurance = data[,2], as.data.frame(sapply(aminoAcids, grepl, x = data[,1])) )
    

    data.melted <- melt(aa.levels, id.var = "occurance")
    colnames(data.melted) <- c("occurance", "aa", "aaGroup")
    
    aaKmer <- ddply(data.melted, .(aa, aaGroup, occurance), nrow)
    plot.data <- subset(aaKmer, occurance <= spectraDepth )
    
    multiAASpectra <- ggplot(plot.data, aes(occurance, V1)) + 
                        geom_bar(stat = "identity", aes(fill = factor(aaGroup))) + 
                        mytheme +  labs(x = paste0("AA ", as.numeric(params$k)/3, "-mer Occurance"), y = paste0( "AA Occurance Frequency")) +
                        scale_fill_manual(values = as.vector(CRUKcolours)) + facet_wrap(~aa, ncol=4)
    
    return(multiAASpectra)
}    

# ---------------------------------------------------------------
# violinPanelFigure takes the distributions and splits the data 
#	similarly to the spectra dinucleotide figures. Here though a 
#	violin plot of the surprisal score for the k-mer is plotted
#	against the number of dinucleotides the k-mer contains.
#
# This takes a long time for large datasets.
# ---------------------------------------------------------------
violinPanelFigure <- function(dist, kmerDinucs){
  
	data <- dist
	colnames(data) <- c("kmer", "occurance")
	data <- subset(data, occurance > 0)

  dinuc.levels <- data.frame( kmer = kmerDinucs[,1], as.data.frame(sapply(kmerDinucs[,-1], function(x) factor(ifelse(x >= 3, "3+", x), levels = c("0", "1", "2", "3+" )))) )
  data.dinuc <- merge(data, dinuc.levels, by = "kmer", all.x = TRUE, all.y = FALSE)

	totalDNA <- sum(data.dinuc$occurance)
	data.dinuc$surprisal <- log2(1/(data.dinuc$occurance/totalDNA))

	data.toMelt <- data.dinuc[,!(names(data.dinuc) %in% c("kmer", "occurance"))]

	data.melted <- melt(data.toMelt, id.vars = "surprisal")
  
	multiViolinPlot <- ggplot(data.melted, aes(x = value, y = surprisal, fill = factor(value)) ) +
						geom_violin() + geom_boxplot(width = .1, outlier.shape = NA, alpha = 0) + 
						mytheme + labs(x = paste0("Di-Nucleotide Occurance in ", params$k, "-mer"), y = "Surprisal (bits)") +
						scale_fill_manual(values = as.vector(CRUKcolours)) + facet_wrap(~variable, ncol=4)
  
	return(multiViolinPlot)
}

#
calcEntropy <- function(occurances){
  
  occurances <- occurances[occurances > 0]
  totalOcc <- sum(occurances)
  surprisal <- log2(1/(occurances/totalOcc))
  entropy <- mean(surprisal)
  
  return(entropy)
} 

#
calcSurprisal <- function(occurances){
  
  totalOcc <- sum(occurances)
  surprisal <- ifelse(occurances == 0, 0, log2(1/(occurances/totalOcc)))
  
  return(surprisal)
} 


# ---------------------------------------------------------------
# plotEIbound plots the exon-intron boundarys in terms of their 
#	surprisal scores. 
#
# This takes a long time for large datasets.
# ---------------------------------------------------------------
plotEIbound.surprisal.test <- function(boundaries){
  
  data <- boundaries

  data.exon5p.df <- data[,grep("surprisal5p_", colnames(data))]
  data.exon5p.df <- data.exon5p.df[complete.cases(data.exon5p.df),]
  data.exon5p.mat <- as.matrix(data.exon5p.df[apply(data.exon5p.df, 1, function(row) {all(row > 0)}), ])
  
  data.exon3p.df <- as.matrix(data[,grep("surprisal3p_", colnames(data))])
  data.exon3p.df <- data.exon3p.df[complete.cases(data.exon3p.df),]
  data.exon3p.mat <- as.matrix(data.exon3p.df[apply(data.exon3p.df, 1, function(row) {all(row > 0)}), ])
  

  # mean scaled by the same amount - not going to scale it for the figure
  data.exon5p.means <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), type = rep('data.5',length(seq(-params$EIdepth, params$EIdepth-1))), 
                                  MeanSuprisal = unname(colMeans(data.exon5p.mat)), 
                                  Entropy = colSums(2^(-data.exon5p.mat)*log2(2^(-data.exon5p.mat))),
                                  SDSuprisal = colSds(data.exon5p.mat))
  
  data.exon3p.means <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), type = rep('data.3',length(seq(-params$EIdepth, params$EIdepth-1))), 
                                  MeanSuprisal = unname(colMeans(data.exon3p.mat)), 
                                  SDSuprisal = colSds(data.exon3p.mat))


  gg.dna.exon5p <- ggplot(data.exon5p.means, aes(x = Position, y = Entropy)) +
  		geom_line(size = 2, colour = as.vector(CRUKcolours[4])) +
  		geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
  		#geom_hline(yintercept = unname(mean(colMeans(data.exon5p.mat))), colour = 'black', alpha = 0.2, size = 1 ) +
      mytheme + theme(plot.margin = margin(r = -0.4, unit = "cm")) + labs( x = "",  y = "DNA entropy (bits)") +
  		scale_x_continuous(expand=c(0.01,0), breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), labels = as.character(c( -params$EIdepth, -params$EIdepth/2, "5' ss", params$EIdepth/2, params$EIdepth) )) #+ 
     # scale_y_continuous(expand=c(0.01,0), limits = c( min(c(lowLimit(data.exon5p.means$MeanSuprisal),lowLimit(data.exon3p.means$MeanSuprisal))), 
      #    max(c(highLimit(data.exon5p.means$MeanSuprisal),highLimit(data.exon3p.means$MeanSuprisal))) ))
  

  gg.dna.exon3p <- ggplot(data.exon3p.means, aes(x = Position, y = MeanSuprisal)) +
  		geom_line(size = 2, colour = as.vector(CRUKcolours[4])) +
  		geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
  		geom_hline(yintercept = unname(mean(colMeans(data.exon5p.mat))), colour = 'black', alpha = 0.2, size = 1 ) +
  	   mytheme + theme(plot.margin = margin(l = -0.4, unit = "cm")) + labs( x = "",  y = "DNA entropy (bits)") +
  		scale_x_continuous(expand=c(0.01,0), breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), labels = c( -params$EIdepth, -params$EIdepth/2, "3' ss", params$EIdepth/2, params$EIdepth) ) + 
      scale_y_continuous(expand=c(0.01,0) , position = "right", limits = c( min(c(lowLimit(data.exon5p.means$MeanSuprisal),lowLimit(data.exon3p.means$MeanSuprisal))), 
        max(c(highLimit(data.exon5p.means$MeanSuprisal),highLimit(data.exon3p.means$MeanSuprisal))) )) 

  x.grob <- textGrob(paste0("Distance from splice site"), gp=gpar(fontfamily = "Helvetica", fontface="bold", fontsize=16))

  ExonIntronMultiplot <- grid.arrange(arrangeGrob( plot_grid(gg.dna.exon5p, gg.dna.exon3p, ncol = 2, align = "v") + 
    theme(plot.margin = unit(c(1,1,1,1), "cm")), bottom = x.grob))
      
  return(ExonIntronMultiplot)
}


# ---------------------------------------------------------------
# plotEIbound plots the exon-intron boundarys in terms of their 
#	surprisal scores. 
#
# This takes a long time for large datasets.
# ---------------------------------------------------------------
plotEIbound.surprisal <- function(boundaries){
  
  data <- boundaries
  
  data.exon5p.df <- data[,grep("surprisal5p_", colnames(data))]
  data.exon5p.df <- data.exon5p.df[complete.cases(data.exon5p.df),]
  data.exon5p.mat <- as.matrix(data.exon5p.df[apply(data.exon5p.df, 1, function(row) {all(row > 0)}), ])
  
  data.exon3p.df <- as.matrix(data[,grep("surprisal3p_", colnames(data))])
  data.exon3p.df <- data.exon3p.df[complete.cases(data.exon3p.df),]
  data.exon3p.mat <- as.matrix(data.exon3p.df[apply(data.exon3p.df, 1, function(row) {all(row > 0)}), ])
  
  # mean scaled by the same amount - not going to scale it for the figure
  data.exon5p.means <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), type = rep('data.5',length(seq(-params$EIdepth, params$EIdepth-1))), 
                                  MeanSuprisal = unname(colMeans(data.exon5p.mat)), 
                                  SDSuprisal = colSds(data.exon5p.mat))
  
  data.exon3p.means <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), type = rep('data.3',length(seq(-params$EIdepth, params$EIdepth-1))), 
                                  MeanSuprisal = unname(colMeans(data.exon3p.mat)), 
                                  SDSuprisal = colSds(data.exon3p.mat))
  
  
  gg.dna.exon5p <- ggplot(data.exon5p.means, aes(x = Position, y = MeanSuprisal)) +
    geom_line(size = 2, colour = as.vector(CRUKcolours[4])) +
    geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
    geom_hline(yintercept = unname(mean(colMeans(data.exon5p.mat))), colour = 'black', alpha = 0.2, size = 1 ) +
    mytheme + theme(plot.margin = margin(r = -0.4, unit = "cm")) + labs( x = "",  y = "DNA entropy (bits)") +
    scale_x_continuous(expand=c(0.01,0), breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), labels = as.character(c( -params$EIdepth, -params$EIdepth/2, "5' ss", params$EIdepth/2, params$EIdepth) )) + 
    scale_y_continuous(expand=c(0.01,0), limits = c( min(c(lowLimit(data.exon5p.means$MeanSuprisal),lowLimit(data.exon3p.means$MeanSuprisal))), 
                                                     max(c(highLimit(data.exon5p.means$MeanSuprisal),highLimit(data.exon3p.means$MeanSuprisal))) ))
  
  
  gg.dna.exon3p <- ggplot(data.exon3p.means, aes(x = Position, y = MeanSuprisal)) +
    geom_line(size = 2, colour = as.vector(CRUKcolours[4])) +
    geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
    geom_hline(yintercept = unname(mean(colMeans(data.exon5p.mat))), colour = 'black', alpha = 0.2, size = 1 ) +
    mytheme + theme(plot.margin = margin(l = -0.4, unit = "cm")) + labs( x = "",  y = "DNA entropy (bits)") +
    scale_x_continuous(expand=c(0.01,0), breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), labels = c( -params$EIdepth, -params$EIdepth/2, "3' ss", params$EIdepth/2, params$EIdepth) ) + 
    scale_y_continuous(expand=c(0.01,0) , position = "right", limits = c( min(c(lowLimit(data.exon5p.means$MeanSuprisal),lowLimit(data.exon3p.means$MeanSuprisal))), 
                                                                          max(c(highLimit(data.exon5p.means$MeanSuprisal),highLimit(data.exon3p.means$MeanSuprisal))) )) 
  
  x.grob <- textGrob(paste0("Distance from splice site"), gp=gpar(fontfamily = "Helvetica", fontface="bold", fontsize=16))
  
  ExonIntronMultiplot <- grid.arrange(arrangeGrob( plot_grid(gg.dna.exon5p, gg.dna.exon3p, ncol = 2, align = "v") + 
                                                     theme(plot.margin = unit(c(1,1,1,1), "cm")), bottom = x.grob))
  
  return(ExonIntronMultiplot)
}

# ---------------------------------------------------------------
#
# ---------------------------------------------------------------
plotEIbound.dinucleotides <- function(data, dinuc){
  
  dinuccontent.exon5p.df <- data[complete.cases(data[, grep("dinuc5p_", colnames(data))]),  grep("dinuc5p_", colnames(data))] 
  totalNumberOfBoundaries.5p <- length(dinuccontent.exon5p.df[,1])
  dinuccontent.exon5p <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), dinuccontent = colSums(dinuccontent.exon5p.df == dinuc))
  
  dinuccontent.exon3p.df <- data[complete.cases(data[, grep("dinuc3p_", colnames(data))]),  grep("dinuc3p_", colnames(data)) ] 
  totalNumberOfBoundaries.3p <- length(dinuccontent.exon3p.df[,1])
  dinuccontent.exon3p <- data.frame(Position = seq(-params$EIdepth, params$EIdepth-1), dinuccontent = colSums(dinuccontent.exon3p.df == dinuc))
  
  topLimit <- max(max(dinuccontent.exon5p$dinuccontent), max(dinuccontent.exon3p$dinuccontent))
  
  gg.dinuc.exon5p <- ggplot(dinuccontent.exon5p, aes(x = Position, y = dinuccontent/totalNumberOfBoundaries.5p)) +
      geom_bar(stat = "identity", colour = as.vector(CRUKcolours[4])) +
      geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
      mytheme + theme(  plot.margin = margin(0, 0, 0, 0, "cm"))  + labs( x = "",  y = "") + 
      scale_x_continuous(expand=c(0.01,0),  breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), 
        labels = c( -params$EIdepth, -params$EIdepth/2, "5'ss", params$EIdepth/2, "") ) + 
      scale_y_continuous(expand=c(0.01,0), limits = c(0,1.1*topLimit/totalNumberOfBoundaries.5p), labels = scales::number_format(accuracy = 0.01))

  gg.dinuc.exon3p <- ggplot(dinuccontent.exon3p, aes(x = Position, y = dinuccontent/totalNumberOfBoundaries.3p)) +
      geom_bar(stat = "identity", colour = as.vector(CRUKcolours[4])) +
      geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
      mytheme + theme(  plot.margin = margin(0, 0, 0, 0, "cm")) + labs( x = "",  y = "") +
      scale_x_continuous(expand=c(0.01,0),  breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), 
        labels = c( "", -params$EIdepth/2, "3'ss", params$EIdepth/2, params$EIdepth) ) + 
      scale_y_continuous(expand=c(0.01,0), position = "right", limits = c(0,1.1*topLimit/totalNumberOfBoundaries.5p), labels = scales::number_format(accuracy = 0.01))

  lab.grob <- textGrob(dinuc, gp=gpar(fontfamily = "Helvetica", fontface="bold", fontsize=12))

  ExonIntronMultiplot.dinuc <-  grid.arrange(gg.dinuc.exon5p, gg.dinuc.exon3p, top = lab.grob, ncol = 2)
      
  return(ExonIntronMultiplot.dinuc)
}

# ---------------------------------------------------------------
#
# ---------------------------------------------------------------
plotExonExonBoundaryEntropy <- function(bound.DNA, bound.CDS, bound.AA){

  DNA.mat <- as.matrix(bound.DNA[,grep("surprisal", colnames(bound.DNA))])
  DNA.mat <- DNA.mat[!(apply(DNA.mat, 1, function(y) any(y == 0))),]
  DNA.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)), 
    type = rep('DNA',length(seq(-(params$EEdepth), (params$EEdepth - 1)))), 
    MeanSurprisal = unname(colMeans(DNA.mat)) - mean(unname(colMeans(DNA.mat))), 
    sds = colSds(DNA.mat))

  CDS.mat <- as.matrix(bound.CDS[,grep("surprisal", colnames(bound.CDS))])
  CDS.mat <- CDS.mat[!(apply(CDS.mat, 1, function(y) any(y == 0))),]
  CDS.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)),
   type = rep('CDS',length(seq(-(params$EEdepth), (params$EEdepth - 1)))),
   MeanSurprisal = unname(colMeans(CDS.mat)) - mean(unname(colMeans(CDS.mat))), 
   sds = colSds(CDS.mat))

  AA.mat <- as.matrix(bound.AA[,grep("surprisal", colnames(bound.AA))])
  AA.mat <- AA.mat[!(apply(AA.mat, 1, function(y) any(y == 0))),]
  AA.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth -1)), 
    type = rep('AA',length(seq(-(params$EEdepth), (params$EEdepth -1)))), 
    MeanSurprisal = unname(colMeans(AA.mat)) - mean(unname(colMeans(AA.mat))), 
    sds = colSds(AA.mat))

  All.means <- rbind(DNA.mean, CDS.mean, AA.mean)

  gg.means <- ggplot(All.means, aes(x = Position, y = MeanSurprisal, group = type, colour = type)) +
      geom_line(size = 1, alpha = 0.7) +
      geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) +
      geom_vline(xintercept = 0, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_vline(xintercept = -1, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_vline(xintercept = 1, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_vline(xintercept = -params$k - 1, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_vline(xintercept = -params$k, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_vline(xintercept = -params$k + 1, colour = 'red', alpha = 0.2, size = 1 ) + 
      geom_hline(yintercept = 0, colour = 'black', alpha = 0.2, size = 1 ) +
      mytheme + labs(x = "",  y = "") + 
      theme(legend.title = element_blank(), legend.position = "right") + scale_color_manual(values = as.vector(CRUKcolours[2:4])) +
      scale_x_continuous(expand=c(0.01,0)) + scale_y_continuous(expand=c(0.01,0), 
        limits = c(-0.30, 0.30))

  return(gg.means)
}

plotIndividualExonExonBoundaryEntropy <- function(bound, label){
  
  mat <- as.matrix(bound[,grep("surprisal", colnames(bound))])
  mat <- mat[!(apply(mat, 1, function(y) any(y == 0))),]
  mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)), 
                     MeanSurprisal = unname(colMeans(mat)) - mean(unname(colMeans(mat))), 
                     sds = colSds(mat))
                     
  gg.means <- ggplot(mean, aes(x = Position, y = MeanSurprisal)) +
    geom_line(size = 1, alpha = 0.7) +
    geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) +
    geom_vline(xintercept = 0, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_vline(xintercept = -1, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_vline(xintercept = 1, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_vline(xintercept = -params$k - 1, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_vline(xintercept = -params$k, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_vline(xintercept = -params$k + 1, colour = 'red', alpha = 0.2, size = 1 ) + 
    geom_hline(yintercept = 0, colour = 'black', alpha = 0.2, size = 1 ) +
    mytheme + labs(x = "",  y = label) +
    scale_x_continuous(expand=c(0.01,0)) + 
    scale_y_continuous(expand=c(0.01,0), limits = c(-0.30, 0.30))
  
  return(gg.means)
}




# ---------------------------------------------------------------
#
# ---------------------------------------------------------------
plotEEbound.dinucleotides <- function(data, dinuc){

  totalNumberOfBoundaries <- length(data[,1])

  dinucCount <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)), 
    dinucCount = colSums(data[,grep("dinuc", colnames(data))] == dinuc))

  gg.EE_dinucFreq <- ggplot(dinucCount, aes(x = Position, y = dinucCount/totalNumberOfBoundaries)) +
      geom_bar(stat = "identity", fill = as.vector(CRUKcolours[4])) +
      geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 2 ) + 
       mytheme + theme(  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + labs(x = "",  y = "") + 
      scale_x_continuous(expand=c(0.01,0), breaks = c(-params$EIdepth, -params$EIdepth/2, -0.5, params$EIdepth/2, params$EIdepth-1), 
        labels = c( -params$EIdepth, -params$EIdepth/2, "5'ss", params$EIdepth/2, params$EIdepth) ) + 
      scale_y_continuous(expand=c(0.01,0), labels = scales::number_format(accuracy = 0.01))

  lab.grob <- textGrob(dinuc, gp=gpar(fontfamily = "Helvetica", fontface="bold", fontsize=12))

  EEplot.dinuc <-  grid.arrange(gg.EE_dinucFreq, top = lab.grob)
      
  return(EEplot.dinuc)

}



# -------------------------------------
# 
# -------------------------------------
plotSurprisaDinucSplit <- function(boundary, dinuc, position, boundType) {
  
    bound <- boundary
  
    bound.dinuc <- bound[bound[, which(colnames(bound) == paste0("dinuc", position))] == dinuc,]
    bound.notdinuc <-  bound[bound[, which(colnames(bound) == paste0("dinuc", position))] != dinuc,]
    
    bound.dinuc.mat <- as.matrix(bound.dinuc[,grep("surprisal", colnames(bound.dinuc))])
    bound.notdinuc.mat <- as.matrix(bound.notdinuc[,grep("surprisal", colnames(bound.notdinuc))])
    
    dinuc.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)), 
      type = rep(dinuc,length(seq(-(params$EEdepth), (params$EEdepth - 1)))), 
      MeanSuprisal = unname(colMeans(bound.dinuc.mat)),
      SDSuprisal = colSds(bound.dinuc.mat))

    notdinuc.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth - 1)), 
      type = rep( paste("not", dinuc, sep = " ") ,length(seq(-(params$EEdepth), (params$EEdepth - 1)))), 
      MeanSuprisal = unname(colMeans(bound.notdinuc.mat)), 
      SDSuprisal = colSds(bound.notdinuc.mat))
  
    plotData <- rbind(dinuc.mean, notdinuc.mean)

    gg.dinuc <- ggplot(plotData, aes(x = Position, y = MeanSuprisal, group = type, colour = type)) +
        geom_line(size = 2) +
        geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
        mytheme + theme(legend.position = "right") +
        labs( x = "",  y = paste(boundType, " entropy (bits)")) + scale_color_manual(values = as.vector(CRUKcolours[2:4])) +
        scale_x_continuous(expand=c(0.01,0)) + scale_y_continuous(expand=c(0.01,0)) + 
        ggtitle(paste0(params$species, ", k = ", params$k, ", D = ", params$EEdepth, ", n(dinuc) = ", 
          length(bound.dinuc$chr), ", n(notdinuc) = ", length(bound.notdinuc$chr)))

    return(gg.dinuc)
}





# ---------------------------------------------------------------
#
# ---------------------------------------------------------------
getInformationGainMatrix <- function(bound.1, bound.2){

tmp.1 <- bound.1[,grep("surprisal", colnames(bound.1))]
tmp.1$ID <- paste(bound.1$transID, bound.1$prevExonID, bound.1$nextExonID, sep = "_")

tmp.2 <- bound.2[,grep("surprisal", colnames(bound.2))]
tmp.2$ID <- paste(bound.2$transID, bound.2$prevExonID, bound.2$nextExonID, sep = "_")

tmp.2 <- tmp.2[match(tmp.1$ID, tmp.2$ID),]

  if(all(tmp.2$ID == tmp.1$ID)){
    infoGain <- tmp.1[,1:(2*params$EEdepth)] - tmp.2[,1:(2*params$EEdepth)]
    infoGain <- cbind(infoGain, ID = tmp.1$ID)
    infoGain <- unique(infoGain, by = "ID")
    rownames(infoGain) <- infoGain$ID
    infoGain$ID <- NULL
    return(infoGain)
  }
  else{
    print("Error, the bounadries do not matchup") 
    return(NULL)
  }
}


# ---------------------------------------------------------------
#
# ---------------------------------------------------------------
plotIndividualBoundary.matrix <- function(boundary, ylab){

boundary.mean <- data.frame(Position = seq(-(params$EEdepth), (params$EEdepth -1)), 
  MeanSurprisal = unname(colMeans(boundary)), 
  sds = colSds(as.matrix(boundary)))

gg.boundary <- ggplot(boundary.mean, aes(x = Position, y = MeanSurprisal)) +
      geom_line(size = 2, alpha = 0.7) +
      geom_vline(xintercept = -0.5, colour = 'black', alpha = 0.2, size = 1 ) + 
      mytheme + labs(x = "",  y = ylab) + 
      theme(legend.title = element_blank()) + scale_color_manual(values = as.vector(CRUKcolours[2])) +
      scale_x_continuous(expand=c(0.01,0)) + scale_y_continuous(expand=c(0.01,0))
}


