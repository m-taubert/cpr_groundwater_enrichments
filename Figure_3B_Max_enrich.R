# Figure 3B

library(tidyr)
library(ggplot2)
library(forcats)
library(plyr)


# calculate absolute abundance for all CPR ASVs in all samples
asv.F3B <- asv[,which(taxa$Phylum=="Patescibacteria")]/rowSums(asv[,which(taxa$Kingdom=="Bacteria" & !is.na(taxa$Phylum))])*meta$qPCR

asv.F3B[1:20,1:5] # will produce NAs in all samples without qPCR data


asv.taxa <- data.frame(ASV=colnames(asv.F3B), Class=taxa[which(taxa$Phylum=="Patescibacteria"),]$Class) # for later Class assignment

# convert data
asv.F3B$sample <- rownames(asv.F3B)
asv.F3B$treatment <- meta$Csource # add required per-sample-information
asv.F3B.long <- pivot_longer(asv.F3B, cols=!c("sample","treatment"), values_to="abundance", names_to="ASV")
asv.F3B.long
asv.F3B.long.tax <- join(asv.F3B.long, asv.taxa, by="ASV")
asv.F3B.long.tax

# filter out all values of NA or 0, filter out all treatments not used
asv.F3B.long.filt <- asv.F3B.long.tax %>% filter(!is.na(abundance) & abundance > 0 & treatment %in% c("autotrophy","methylotrophy"))

# calculate the absolute abundance in situ of ASVs from each of the five classes
# this will give conservative estimates of enrichment
CPR.class <- c("ABY1","Berkelbacteria","Gracilibacteria","Parcubacteria","Saccharimonadia")

CPR.day0 <- data.frame(Class=CPR.class, class.max=sapply(CPR.class, function(c) {
 asv.F3B[which(meta$Period=="day0" & meta$filter.fraction==0.2),which(taxa[which(taxa$Phylum=="Patescibacteria"),]$Class==c)] %>% max(na.rm=TRUE)})) # max abundance of ASVs from CPR classes in day0 samples
CPR.day0

# integrate into filtered data frame
asv.F3B.long.plot <- asv.F3B.long.filt %>% join(CPR.day0, by="Class")
asv.F3B.long.plot
# calculate fold.max enrichment
asv.F3B.long.plot <- asv.F3B.long.plot %>% mutate(fold.max=abundance/class.max)
asv.F3B.long.plot

subset(asv.F3B.long.plot, fold.max >= 1) # to plot only values above 1 in the figure



# plot the figure
ggplot(subset(asv.F3B.long.plot, fold.max >= 1), aes(x=treatment, y = fold.max)) +
 geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
 geom_jitter(aes(color = Class), alpha = 0.7, size = 2.5, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) + 
 scale_color_manual(values = c("Saccharimonadia" = "#0056BE", "Parcubacteria" = "#BF0000", "Gracilibacteria" = "#739900", "Berkelbacteria"="#4D009B", "ABY1" = "#BF7300")) +  
 scale_y_log10(breaks = c(1, 5, 25, 125)) +  
 theme_minimal() + labs(x = "Treatment", y = "X-fold enrichment", color = "Class") +  
 theme(legend.position = "right")
