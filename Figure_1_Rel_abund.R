# Figure 1A

library(tidyr)
library(ggplot2)
# required data
dim(asv)

# only use ASVs that are Bacteria and that are classified on Phylum level (good ASVs)
asv.F1A <- which(taxa$Kingdom=="Bacteria" & !is.na(taxa$Phylum))
dim(asv[,asv.F1A])

# calculate relative abundances of CPR in these samples
which(taxa$Phylum=="Patescibacteria") # ASVs of CPR

rowSums(asv[,which(taxa$Phylum=="Patescibacteria")]) # sum of CPR reads per sample
rowSums(asv[,asv.F1A]) # sum of all reads of good ASVs per sample

meta$CPR.rel.abund <- rowSums(asv[,which(taxa$Phylum=="Patescibacteria")])/rowSums(asv[,asv.F1A]) # CPR relative abundances added to metadata

# only use samples from 0.2 µm filter fraction
dim(meta[which(meta$filter.fraction==0.2),])
# select relevant data from metadata
F1A.raw <- meta[which(meta$filter.fraction==0.2),c("Csource","Period","CPR.rel.abund")]
F1A.raw

# add in situ data
in.situ <- read.csv(file="import/patescibacterial_RA-data0.csv", header=TRUE, sep=",")

F1A <- rbind(F1A.raw, in.situ[,c("Csource","Period","CPR.rel.abund")])

# plot Figure
ggplot(F1A, aes(x=Period, y=CPR.rel.abund*100)) + 
 geom_boxplot(color="blue", fill="blue", alpha=0.5, outlier.colour="black", outlier.size=3) + 
 scale_x_discrete(limits=c("in situ", "day0","1week","2weeks","3weeks",">3weeks")) +  coord_cartesian(ylim = c(NA, 55)) #

