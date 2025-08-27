#!/usr/bin/env Rscript

#### DOC ####
'Given an .eigen file containing eigenvectors from Plink, plot PCA plots of the specified PCs grouped by species.

Usage:
  plink-eigen2pca.R <EIGEN> [--pcx=INT] [--pcy=INT] [--miss=STR] [--eigenvals=STR] [-v -h]

  Arguments:
   EIGEN                 A .eigenvec file obtained from Plink.

  Options:
    --pcx=INT             The number of the PC that will be plotted in the x-axis [default: 1].
    --pcy=INT             The number of the PC that will be plotted in the y-axis [default: 2].
    --miss=STR            A .imiss file obtained from vcftools, indicating the missigness by sample. If provided, PCA is colored by missingness.
    --eigenvals=STR       A .eigenval file obtained from Plink to add the % variance explained by each PC.
    -v, --verbose         Print the progression of the program execution to the terminal (Standard Error).
    -h, --help            Show this message and exit.

' -> doc

#### LIBS ####
library(docopt)
library(ggplot2)
library(viridis)

#### ARGS ####
args <- docopt(doc, version = 'plink-eigen2pca.R 0.0.1')
print(args)

#### MAIN ####
# Prepare eigenvectors data
data = read.delim(args$EIGEN, sep = " ", header = F)
# data_raw = data
data = data[, -1]
colnames(data)[1] <- "sampleID"

data$species = unlist(lapply(X = as.list(data$sampleID), 
                             FUN = function (x) strsplit(x = x, split = ".", fixed = T)[[1]][1]))

genus = strsplit(x = data$species[1], split = "_")[[1]][1]

mask <- grep(pattern = "V", x = colnames(data), fixed = T)
pcs <- colnames(data)[mask]

# If eigenvalues are provided, make scree plot with % variance explained
if (!is.null(args$eigenvals)){
  eigenvals <- read.delim(args$eigenvals, header = F)
  colnames(eigenvals) <- "Eigenvalue"
  eigenvals$PC <- c(1:nrow(eigenvals))
  eigenvals$var_expl <- 100*eigenvals$Eigenvalue/sum(eigenvals$Eigenvalue)
  plot <- ggplot(data = eigenvals, aes(y = var_expl, x = PC)) +
    geom_point() + 
    geom_line() +
    ggtitle(paste0("Screeplot for ", genus)) +
    ylab("% variance explained")
  ggsave(filename = paste0("screeplot_", genus, ".png"), plot = plot, dpi = 300)
  
  x_label <- paste0("PC", args$pcx, " (", round(eigenvals$var_expl[eigenvals$PC == args$pcx], digits = 2), " %)")
  y_label <- paste0("PC", args$pcy, " (", round(eigenvals$var_expl[eigenvals$PC == args$pcy], digits = 2), " %)")
}

# Make PCA plot
if (!is.null(args$miss)) {
  miss_df <- read.delim(args$miss)
  for (id in data$sampleID){
    data$miss[data$sampleID == id] <- miss_df$F_MISS[miss_df$INDV == id]
  }
  
  plot_pca <- function(pcx, pcy) {
    pca_plot <- ggplot(data = data, aes(x = .data[[pcx]], y = .data[[pcy]],
                                        color = miss)) +
      geom_point(aes(shape = species)) +
      scale_color_viridis() +
      xlab(x_label) +
      ylab(y_label) +
      ggtitle(paste0("PCA for ", genus)) +
      theme(plot.title = element_text(hjust = 0.5))
    return(pca_plot)
  }
  
} else {
  plot_pca <- function(pcx, pcy) {
    pca_plot <- ggplot(data = data, aes(x = .data[[pcx]], y = .data[[pcy]],
                                        color = species)) +
      geom_point() +
      xlab(x_label) +
      ylab(y_label) +
      ggtitle(paste0("PCA for ", genus)) +
      theme(plot.title = element_text(hjust = 0.5))
    return(pca_plot)
  }
}

pca_plot <- plot_pca(pcs[1], pcs[2])

ggsave(filename = paste0("pca_", genus, ".png"), plot = pca_plot, dpi = 300)
ggsave(filename = paste0("pca_", genus, ".svg"), plot = pca_plot, dpi = 300)