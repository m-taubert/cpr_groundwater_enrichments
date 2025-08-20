# Figure 5A
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)

# define selected samples for analysis (only those with qPCR data)
samples.sel <- which(meta$filter.fraction==0.2 & !is.na(meta$qPCR))
which(meta$filter.fraction==0.2 & !is.na(meta$qPCR)) %>% length 

# define selected ASVs for analysis (all CPR)
asv.sel <- which(taxa$Phylum=="Patescibacteria")
which(taxa$Phylum=="Patescibacteria") %>% length 

# test selection
asv[samples.sel,asv.sel][1:10,1:10]
asv[samples.sel,asv.sel] %>% dim

# determine relative and then directly absolute abundance
asv[,asv.F1A] %>% dim # good ASVs (see Fig 1A)
asv.Fig5A.abs <- asv[samples.sel,asv.sel]/rowSums(asv[samples.sel,asv.F1A])*meta[samples.sel,]$qPCR
asv.Fig5A.abs[1:5,1:10]

# add metadata columns for incubation time and experiment
asv.Fig5A.abs$days <- meta[samples.sel,]$days
asv.Fig5A.abs$folder <- meta[samples.sel,]$folder

# change to long format
asv.Fig5A.filt <- asv.Fig5A.abs %>% pivot_longer(cols=1:(ncol(asv.Fig5A.abs)-2), names_to="ASV", values_to="abundance") %>% filter(abundance!=0)

# calculate growth rates for each ASV within each experiment (folder)
asv.Fig5A.rates <- asv.Fig5A.filt %>%
 group_by(ASV, folder) %>%
 do({
  mod <- lm(log(abundance) ~ days, data = .)
  slope <- coef(mod)[2]
  tibble(growth_rate = slope)
 }) %>%
 ungroup()

# filter to keep only positive growth rates
asv.Fig5A.rates <- asv.Fig5A.rates %>% filter(growth_rate>0)
asv.Fig5A.rates

# calculate doubling times
asv.Fig5A.rates$doubling_time <- log(2)/asv.Fig5A.rates$growth_rate

# helper df for taxonomy assignment
asv.taxa <- data.frame(ASV=colnames(asv[,asv.sel]), Class=taxa[asv.sel,]$Class)
asv.taxa 

# add class information
asv.Fig5A.rates.tax <- join(asv.Fig5A.rates, asv.taxa, by="ASV")
asv.Fig5A.rates.tax 
asv.Fig5A.rates.tax <- asv.Fig5A.rates.tax %>% filter(!is.na(Class)) # filter out CPR with unassigned Class


# plot figure
ggplot(asv.Fig5A.rates.tax, aes(x = Class, y = growth_rate, color = Class)) +
 geom_boxplot(outlier.shape = NA) + 
 geom_jitter(color = "gray48", width = 0.2, height = 0, size = 1.5, alpha = 1) + 
 scale_color_manual(values = c("Saccharimonadia" = "#0056BE", "Parcubacteria" = "#BF0000", "Microgenomatia" = "#734D26", "Gracilibacteria" = "#739900", "Berkelbacteria"="#4D009B", "ABY1" = "#BF7300")) + 
 labs(y = "Growth Rate") +
 theme_minimal() +
 scale_y_continuous(name = "Growth Rate", breaks = seq(0, 0.8, by = 0.2), limits=c(0,1)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()

