# Figure 5B
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)

# define selected samples for analysis (only those with qPCR data)
samples.sel <- which(meta$filter.fraction==0.2 & !is.na(meta$qPCR))
which(meta$filter.fraction==0.2 & !is.na(meta$qPCR)) %>% length 

# define selected ASVs for analysis (abundant non-CPR orders)
non.CPR <- c("Microtrichales","Rhizobiales","Flavobacteriales","Sphingobacteriales","Chlamydiales","Burkholderiales","Pseudomonadales","Nitrospirales","Pirellulales","Opitutales","Pedosphaerales")
# check ASV numbers for these orders
sapply(non.CPR, function(o) {which(taxa$Order==o) %>% length})
# create vector for selection
asv.sel <- which(taxa$Order %in% non.CPR)
asv.sel %>% length 

# test selection
asv[samples.sel,asv.sel][1:10,1:10]
asv[samples.sel,asv.sel] %>% dim

# determine relative and then directly absolute abundance
asv[,asv.F1A] %>% dim # good ASVs (see Fig 1A)
asv.Fig5B.abs <- asv[samples.sel,asv.sel]/rowSums(asv[samples.sel,asv.F1A])*meta[samples.sel,]$qPCR
asv.Fig5B.abs[1:5,1:10]

# add metadata columns for incubation time and experiment
asv.Fig5B.abs$days <- meta[samples.sel,]$days
asv.Fig5B.abs$folder <- meta[samples.sel,]$folder

# change to long format
asv.Fig5B.filt <- asv.Fig5B.abs %>% pivot_longer(cols=1:(ncol(asv.Fig5B.abs)-2), names_to="ASV", values_to="abundance") %>% filter(abundance!=0)

# calculate growth rates for each ASV within each experiment (folder)
asv.Fig5B.rates <- asv.Fig5B.filt %>%
 group_by(ASV, folder) %>%
 do({
  mod <- lm(log(abundance) ~ days, data = .)
  slope <- coef(mod)[2]
  tibble(growth_rate = slope)
 }) %>%
 ungroup()

# filter to keep only positive growth rates
asv.Fig5B.rates <- asv.Fig5B.rates %>% filter(growth_rate>0)
asv.Fig5B.rates

# calculate doubling times
asv.Fig5B.rates$doubling_time <- log(2)/asv.Fig5B.rates$growth_rate

# helper df for taxonomy assignment
asv.taxa <- data.frame(ASV=colnames(asv[,asv.sel]), Order=taxa[asv.sel,]$Order)
asv.taxa 

# add Order information
asv.Fig5B.rates.tax <- join(asv.Fig5B.rates, asv.taxa, by="ASV")
asv.Fig5B.rates.tax 
asv.Fig5B.rates.tax <- asv.Fig5B.rates.tax %>% filter(!is.na(Order)) # filter out ASVs with unassigned Order


# plot figure
ggplot(asv.Fig5B.rates.tax, aes(x = Order, y = growth_rate, color = Order)) +
 geom_boxplot(outlier.shape = NA, color="#000000") + 
 geom_jitter(color = "gray48", width = 0.2, height = 0, size = 1.5, alpha = 1) + 
 labs(y = "Growth Rate") +
 theme_minimal() +
 scale_y_continuous(name = "Growth Rate", breaks = seq(0, 0.8, by = 0.2), limits=c(0,1)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()

