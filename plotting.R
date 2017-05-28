library("gsubfn")
library("zoo")
library("ggplot2")
library("RColorBrewer")

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
#' @param normTotalSNVs Total number of SNVs to calculate % of of covered variants

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

