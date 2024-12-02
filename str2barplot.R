#!/usr/bin/env Rscript

#### LIBS ####
library(ggplot2)
library(tidyverse)

#### VARS ####
setwd("/home/etode/Documents/Plebejus_phylogeny_ddRADSeq/5_Structure/map_figure")
df <- read.csv(file = "p_idas.k3.clumpak.good.csv", header = T)


#### MAIN ####
# Make data long
df %>% pivot_longer(cols = c(P1, P2, P3), names_to = "group", values_to = "proportion") -> df.long

# Turn samples into factor to preserve order
df.long$Sample <- factor(df.long$Sample, levels = unique(df$Sample))

# Plot and save
fig <- ggplot(data=df.long, aes(x=Sample, y=proportion, fill=group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#a2fc3c", "#4686fb", "#30123b")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1),
        axis.title = element_blank(),
        legend.position = "none")

figheight = 3000
figwidth = 3*figheight
for (i in c("png", "pdf")){
  name = paste0("structure.plebejus_idas.k3.", i)
  ggsave(filename = name, plot = fig, width = figwidth, height = figheight, units = "px")
}


  