# Figure 7
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)

# select samples with qPCR data not from day 0
samples.sel <- which(!is.na(meta$qPCR) & meta$days!=0)

# select ASVs that are CPR
asv.sel <- which(taxa$Phylum=="Patescibacteria")

# calculate absolute abundances
asv[,asv.F1A] %>% dim # good ASVs (see Fig 1A)
asv.Fig7.abs <- asv[samples.sel,asv.sel]/rowSums(asv[samples.sel,asv.F1A])*meta[samples.sel,]$qPCR

# add metadata column for oxygen
asv.Fig7.abs$oxygen <- meta[samples.sel,]$oxygen

# change to long format
asv.Fig7.filt <- asv.Fig7.abs %>% pivot_longer(cols=1:(ncol(asv.Fig7.abs)-1), names_to="ASV", values_to="abundance") %>% filter(abundance!=0)

# helper df for taxonomy assignment
asv.taxa <- data.frame(ASV=colnames(asv[,asv.sel]), Class=taxa[asv.sel,]$Class)

# add class information
asv.Fig7.tax <- join(asv.Fig7.filt, asv.taxa, by="ASV")

# filter for CPR of specific classes
asv.Fig7.tax <- asv.Fig7.tax %>% filter(Class %in% c("ABY1","Berkelbacteria","Gracilibacteria","Parcubacteria","Saccharimonadia"))

# Plot Figure 7
ggplot(asv.Fig7.tax, aes(x = Class, y = abundance, fill = oxygen)) +
 geom_boxplot(outlier.shape = NA,position = position_dodge(width = 0.8)) +
 geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "#7b7a7aff") +
 scale_y_log10() +
 scale_fill_manual(values = c("oxic" = "lightblue3", "anoxic" = "#fbd9d3")) +
 theme_minimal() +
 labs(x="Class", y="Abs. Abundance / gene copies liter-1", fill="Oxygen Condition")
