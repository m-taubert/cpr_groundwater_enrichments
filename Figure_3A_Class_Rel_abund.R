# Figure 3A

library(tidyr)
library(forcats)

# collect existinc combinations of Class and Order for CPR
ClassXOrder <- unique(taxa[which(taxa$Phylum=="Patescibacteria"),c("Class","Order")])
ClassXOrder

# calculate the maximal relative abundance of each of these in each treatment
# example
rowSums(asv[,which(taxa$Class==ClassXOrder[1,1] & ((taxa$Order==ClassXOrder[1, 2]) | ((is.na(taxa$Order) & is.na(ClassXOrder[1, 2])))))])/rowSums(asv[,asv.F1A]) # relative abundances per sample, works also for Order==NA
# asv.F1A is the vector of all good ASVs created in script for Fig 1A, run it before this script

# loop over all of:
1:nrow(ClassXOrder)
# with function
ClassXOrder.abund <- sapply(1:nrow(ClassXOrder), function(i) {
 rowSums(asv[,which(taxa$Class==ClassXOrder[i,1] & ((taxa$Order==ClassXOrder[i, 2]) | ((is.na(taxa$Order) & is.na(ClassXOrder[i, 2]))))), drop = FALSE])/rowSums(asv[,asv.F1A])
})

# transpose
ClassXOrder.abund <- as.data.frame(t(ClassXOrder.abund))
ClassXOrder.abund

# summarize over the samples of each treatment
meta$Csource
treatments <- c("autotrophy","methylotrophy","organic_defined","organic_undefined")


ClassXOrder.treat <- sapply(treatments, function(t) {
 # cycles over each of the four treatments
 apply(ClassXOrder.abund[,which(meta$Csource==t & meta$filter.fraction==0.2)], 1, max, na.rm = TRUE)
 # determines the maximum value for each sample in the respective treatment
})

# add Class and Order information
ClassXOrder.treat <- cbind(ClassXOrder.treat, ClassXOrder)
ClassXOrder.treat

# filter to include abundant orders only
ClassXOrder.treat[,1:4] %>% rowSums
ClassXOrder.filter <- c(1:17,18,22,42)

# transform data for plotting
Fig3A <- pivot_longer(ClassXOrder.treat[ClassXOrder.filter,], cols=1:4, names_to="Csource", values_to="abundance")
Fig3A$Class[is.na(Fig3A$Class)] <- "unknown"
Fig3A$Order[is.na(Fig3A$Order)] <- "unknown"
Fig3A$Taxa <- interaction(Fig3A$Class, Fig3A$Order, drop = TRUE)
Fig3A$Taxa <- factor(Fig3A$Taxa, levels = sort(as.character(unique(Fig3A$Taxa))))
Fig3A$Taxa

# plot Figure 3A
ggplot(Fig3A, aes(x = Csource, y = Taxa, size = abundance*100, color = Class)) +
 geom_point(alpha = 0.7) +
 scale_size_continuous(range = c(0.1, 18), breaks = c(1, 10, 20, 30)) +  # Set size range and specific breaks
 scale_color_manual(values = c("Saccharimonadia" = "#0056BE", "Parcubacteria" = "#BF0000", "Microgenomatia" = "#734D26", "Gracilibacteria" = "#739900", "Berkelbacteria"="#4D009B", "ABY1" = "#BF7300")) +  # Apply custom colors
 labs(x = "Treatment",
      y = "Taxonomy1",
      size = "Abundance",
      color = "Taxonomy") +
 theme_minimal() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1))


