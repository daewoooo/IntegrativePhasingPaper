#Required packages to run plotting functions
library("gsubfn")
library("zoo")
library("ggplot2")
library("RColorBrewer")
library("egg")
library("biovizBase")
library("GenomicRanges")

#' Function to preprocess table with comparison data
#' 
#' @param table.data Path to the table to preprocess for plotting
#' @param id.col Column contatinig info about number of SS cells and PB coverage

preprocessTable <- function(table.data, id.col=6) {
  
  data <- read.table(table.data, header = T, stringsAsFactors = F, comment.char = '&') #load the data
  data$ID <- paste(data[,3], data[,4], sep = "_") #merge column 3 and 4 to create ID for various comparisons
  
  ## take info about SS cells and PB coverage and store them in separate columns
  cells <- strapplyc(data[,id.col], "cells(\\d+)", simplify = T)
  cov <- strapplyc(data[,id.col], "cov(\\d+)", simplify = T)
  trials <- strapplyc(data[,id.col], "matrix-(\\d+)", simplify = T)
  technology <- strapplyc(data[,id.col], "^(\\w+)/", simplify = T)
  
  data$SS_sample <- as.numeric(cells)
  data$PB_sample <- as.numeric(cov)
  data$PB_sample[is.na(data$PB_sample)] <- 31
  data$trials <- as.numeric(trials)
  data$technology <- technology
  
  return(data)
}  


#' Function to make lineplot from summary table
#' 
#' @param plt.data Preprocessed data to be plotted.
#' @param title Title used for the plot, if not specified column name is used.
#' @param plt.column Number of the column to be plotted.
#' @param normTotalSNVs Total number of SNVs to calculate % of of covered variants.

plotLines <- function(plt.data, title=NULL, plt.column=NULL, normTotalSNVs=NULL) {
  
  if (is.null(plt.column)) {
    message("Argument plt.column have to be specidfied!!!")
  }  
  
  if (is.null(title)) { 
    title <- names(plt.data)[plt.column]
  }
  
  if (!is.null(normTotalSNVs)) {
    plt.data[,plt.column] <- (plt.data[,plt.column]/normTotalSNVs)*100
  }
  
  column.trials <- split(plt.data[,plt.column], plt.data$trials)
  column.means <- as.data.frame(apply(simplify2array(lapply(column.trials, as.matrix)),1:2,mean, na.rm = TRUE))
  trials <- split(plt.data, plt.data$trials)
  plt.df <- cbind(column.means, trials[[1]][,c('SS_sample','PB_sample','technology')])
  
  plt <- ggplot(plt.df[plt.df$PB_sample != 0,], aes(x=SS_sample, y=V1, group=PB_sample, color=factor(PB_sample))) + facet_grid(technology ~ .) + geom_point() + geom_line() + scale_x_continuous(breaks = c(0,5,10,20,40,60,80,100,120,134)) + scale_colour_brewer(palette = "Set1", name="Depth") + theme_bw() + ylab(title)
  
  return(plt)
}

#' Function to collapse subsequent genomic ranges with the same ID speciefied in id.field
#' 
#' @param gr Genomic ranges object.
#' @param id.field Metadata field specifing unique ID for each genomic range

collapseBins <- function(gr, id.field=0) {
  ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
  ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
  collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
  names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
  return(collapsed.gr)
}

#' Function to collapse subsequent genomic ranges with the same ID speciefied in id.field
#' 
#' @param vcfFile VCF file to load in.
#' @param field Genotype field to load in.
#' @param filterUnphased Specify if unphased variants should be filtered out.

vcf2rangesWH <- function(vcfFile=NULL, field=1, filterUnphased=T) {
  vcf <- read.table(vcfFile, stringsAsFactors = F)
  
  genotypeField <- field + 9
  vcf <- vcf[,c(1:9, genotypeField)]
  
  if (filterUnphased) {
    mask <- grepl(pattern = "\\|", vcf[,10])
    vcf <- vcf[mask,]
  }
  
  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\|")
  alleles <- do.call(rbind, alleles)
  
  #filter only HET pos
  mask <- alleles[,1] != alleles[,2]
  gen.block <- as.matrix(gen.block[mask,])
  vcf <- vcf[mask,]
  
  offset <- round(runif(1, 0, 1000000))
  
  if (dim(gen.block)[2] > 1) {
    hap.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2), blockID=as.numeric(gen.block[,2])+offset)
  } else {
    hap.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2), blockID=rep(1+offset, nrow(gen.block)))  
  }  
  
  return(hap.gr)
} 

#' @param vcfpath
#' @param gapfile
#' @param normTotalSNVs Total number of SNVs to calculate % of of covered variants.
#' @param plotBy ID field use for faceting: SS_cells or Depth

vcfpath = "/home/daewoooo/WORK/Integrative_phasing_paper/Data_analysis/PlatinumHaps_analysis/VCFs/TEST/"
gapfile = "/home/daewoooo/WORK/Integrative_phasing_paper/Data_analysis/LTS_analysis_fixedIllumina/gap_chr1_GRCh37"

plotVCFbyCells <- function(vcfpath=NULL, gapfile=NULL, normTotalSNVs=157131) {

files2plot <- list.files(vcfpath, pattern = "vcf$")
gapfile.chr1 <- read.table(gapfile, header=F)

all.files <- list()
for (infile in files2plot) {
  file <- file.path(vcfpath, infile)
  file.size <- file.info(file)$size
  
  if (file.size) {
    data.tab.gr <- vcf2rangesWH(vcfFile = file, field = 1, filterUnphased = T)
   
    cells <- as.numeric( strapplyc(infile, "cells(\\d+)", simplify = T) )
    cov <- as.numeric( strapplyc(infile, "cov(\\d+)", simplify = T) )
    cov[is.na(cov)] <- 31  
    
    data.tab.gr$SS_cells <- cells
    data.tab.gr$Depth <- cov
    all.files[[1+length(all.files)]] <- as.data.frame(data.tab.gr)  
  }  
}

all.files.df <- do.call(rbind, all.files)
all.files.df$blockID <- as.character(all.files.df$blockID)

perc.cov.count <- sapply(all.files, function(x) max(table(x[,6])))
perc.total <- (perc.cov.count/normTotalSNVs)*100
perc.total.df <- as.data.frame(perc.total)
perc.total.df$ID <- rownames(perc.total.df)
perc.total.df$Tech <- names(all.files)

hg19Ideogram <- getIdeogram("hg19", cytoband = FALSE)
chr.df <- as.data.frame(hg19Ideogram[seqnames(hg19Ideogram) == 'chr1'])

col.levels <- levels(factor(all.files.df$blockID))
col.needed <- length(col.levels)
col.pallete <- rep(brewer.pal(n=9, name="Set1")[-1], ceiling(col.needed/8))
col.pallete <- col.pallete[1:col.needed]
maxBlck <- sapply(all.files, function(x) names(which.max(table(x[,6]))) )
col.pallete[which(col.levels %in% maxBlck)] <- "#E41A1C"

my.theme <- theme(legend.position="none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  strip.text.y = element_text(angle = 180) ) 

wrap <- as.character(unique(all.files.df[,7]))

as.formula(paste(all.files.df[,7], "~ ."))

#add 50000 offset to each phased region to get better visibility of pixels
p1 <- ggplot() + geom_rect(data=all.files.df, aes(xmin=start, xmax=end+50000, ymin=0, ymax=1, fill=factor(blockID))) + scale_fill_manual(values = col.pallete) + my.theme + facet_grid(SS_cells ~ . , switch = "y") + geom_rect(data=gapfile.chr1, aes(xmin=V2, xmax=V3, ymin=0, ymax=1), fill="white", inherit.aes=F) + geom_rect(data=chr.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), color="black", fill=NA, inherit.aes=F)
p2 <- ggplot() + geom_bar(data=perc.total.df, aes(x=reorder(ID, -perc.total) , y=perc.total), stat="identity") + coord_flip() + ylab("") + xlab("") + my.theme + scale_x_discrete(expand=c(0,0)) + scale_y_continuous(limits = c(0,100))

ggarrange(p1, p2, ncol = 2, widths = c(6,1))
}


plotVCFbyDepth <- function(vcfpath=NULL, gapfile=NULL, normTotalSNVs=157131) {
  
  files2plot <- list.files(vcfpath, pattern = "vcf$")
  
  all.files <- list()
  for (infile in files2plot) {
    file <- file.path(vcfpath, infile)
    file.size <- file.info(file)$size
    
    if (file.size) {
      data.tab.gr <- vcf2rangesWH(vcfFile = file, field = 1, filterUnphased = T)
      
      cells <- as.numeric( strapplyc(infile, "cells(\\d+)", simplify = T) )
      cov <- as.numeric( strapplyc(infile, "cov(\\d+)", simplify = T) )
      cov[is.na(cov)] <- 31  
      
      data.tab.gr$SS_cells <- cells
      data.tab.gr$Depth <- cov
      all.files[[1+length(all.files)]] <- as.data.frame(data.tab.gr)  
    }  
  }
  
  all.files.df <- do.call(rbind, all.files)
  all.files.df$blockID <- as.character(all.files.df$blockID)
  
  perc.cov.count <- sapply(all.files, function(x) max(table(x[,6])))
  perc.total <- (perc.cov.count/normTotalSNVs)*100
  perc.total.df <- as.data.frame(perc.total)
  perc.total.df$ID <- rownames(perc.total.df)
  perc.total.df$Tech <- names(all.files)
  
  hg19Ideogram <- getIdeogram("hg19", cytoband = FALSE)
  chr.df <- as.data.frame(hg19Ideogram[seqnames(hg19Ideogram) == 'chr1'])
  
  col.levels <- levels(factor(all.files.df$blockID))
  col.needed <- length(col.levels)
  col.pallete <- rep(brewer.pal(n=9, name="Set1")[-1], ceiling(col.needed/8))
  col.pallete <- col.pallete[1:col.needed]
  maxBlck <- sapply(all.files, function(x) names(which.max(table(x[,6]))) )
  col.pallete[which(col.levels %in% maxBlck)] <- "#E41A1C"
  
  my.theme <- theme(legend.position="none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    strip.text.y = element_text(angle = 180) ) 
  
  wrap <- as.character(unique(all.files.df[,7]))
  
  as.formula(paste(all.files.df[,7], "~ ."))
  
  #add 50000 offset to each phased region to get better visibility of pixels
  p1 <- ggplot() + geom_rect(data=all.files.df, aes(xmin=start, xmax=end+50000, ymin=0, ymax=1, fill=factor(blockID))) + scale_fill_manual(values = col.pallete) + my.theme + facet_grid(Depth ~ . , switch = "y") + geom_rect(data=gapfile.chr1, aes(xmin=V2, xmax=V3, ymin=0, ymax=1), fill="white", inherit.aes=F) + geom_rect(data=chr.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), color="black", fill=NA, inherit.aes=F)
  p2 <- ggplot() + geom_bar(data=perc.total.df, aes(x=reorder(ID, -perc.total) , y=perc.total), stat="identity") + coord_flip() + ylab("") + xlab("") + my.theme + scale_x_discrete(expand=c(0,0)) + scale_y_continuous(limits = c(0,100))
  
  ggarrange(p1, p2, ncol = 2, widths = c(6,1))
}  
