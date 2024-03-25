#!/usr/bin/env Rscript

#### LIBS ####
library(tidyverse)
library(ggplot2)

#### RAW ####
data <- list()
data[[1]] <- read.delim2("/home/etd530/Documents/TFM_holobiome/9_nonBfly_blobs/iphiclides_podalirius/flye/iphiclides_podalirius_IP504_flye_nonPolished.blobDB.table.txt",
                    skip = 11)

data[[2]] <- read.delim2("/home/etd530/Documents/TFM_holobiome/9_nonBfly_blobs/iphiclides_podalirius/spades/iphiclides_podalirius_IP504_spades_error_corrected_scaffolds.blobDB.table.txt",
                         skip = 11)

names(data) <- c("Flye", "Spades")

#### MAIN ####
tax_levels <- c("phylum.t.12.s", "order.t.16.s", "family.t.20.s", "genus.t.24.s", "species.t.28.s")

data_no_Arthropoda <- lapply(X = data,
                             FUN = function(df){
                               df %>% 
                                 filter(phylum.t.12.s != "Arthropoda") -> df_no_Arthropoda;
                               return(df_no_Arthropoda)
                               }
                             )


  
Richness <- lapply(X = data_no_Arthropoda, FUN = function(df, tax_levels){
  for (i in 1:length(tax_levels)){
    varname = parse(text=tax_levels[i])
    df %>%
      filter(eval(varname) != "undef" &
               eval(varname) != "unresolved" & 
               eval(varname) != "no-hit") %>% 
      pull(eval(varname)) %>% 
      unique() %>% length() -> tax_num
    tmp <- paste(tax_levels[i], tax_num, collapse = '\t')
    Richness <- if(i == 1){tmp} else {paste(Richness, tmp, sep = "\\n")}
  }
  return(Richness)
},
tax_levels = tax_levels
)

taxons_summary <- lapply(X = data_no_Arthropoda, FUN = function(df){
  df %>%
    group_by(phylum.t.12.s, order.t.16.s, family.t.20.s, genus.t.24.s, species.t.28.s) %>%
    summarise(contig_num = n(),
              span = sum(length),
              median_contig_size = median(length)) -> summary
  return(summary)
})


# Plot span histograms
span_histograms <- mapply(df = taxons_summary, FUN = function(df, names){
    plot <- ggplot(data = df, mapping = aes(x = log(span, base = 10),
                                            fill = order.t.16.s)) +
      geom_histogram() + labs(title = names)
    print(plot)
    return(plot)
},
names = names(taxons_summary)
)

data_no_Arthropoda_sorted <- data_no_Arthropoda[order(data_no_Arthropoda$length, decreasing = T),]
data_no_Arthropoda_sorted$ordering <- as.factor(c(1:nrow(data_no_Arthropoda_sorted)))


ggplot(data = data_no_Arthropoda_sorted, aes(x = ordering, y = length)) +
  geom_bar(stat="identity")

