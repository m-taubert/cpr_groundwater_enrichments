library("plyr")
library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")
library("forcats")

###################################################
##########          Plot Figure 8A       ##########
###################################################
# run script Figure_8_prep.R before running this script

# add subnetwork classification to CPR annotations
CPR_annotations.EC <- join(CPR_annotations, subnetwork, by="EC")

# count number of genes for each network category
MAG.network.wide <- CPR_annotations.EC %>% filter(!is.na(network)) %>% group_by(MAG, network) %>% dplyr::summarise(count=n(), .groups="drop") %>% pivot_wider(names_from = network, values_from = count, values_fill = 0)

# add taxonomy information
identical(MAG_taxonomy_CPR$sequence,MAG.network.wide$MAG) # the order is the same in both dfs
MAG.network.wide$Class <- MAG_taxonomy_CPR$Class

# subsetting for taxonomy and converting to long format
MAG.network.long <- pivot_longer(subset(MAG.network.wide, Class %in% c("ABY1","Berkelbacteria","Gracilibacteria","Microgenomatia","Paceibacteria","Saccharimonadia")), cols=-c("MAG","Class"), names_to="network", values_to="EC.count")
MAG.network.long

# plot Figure
ggplot(MAG.network.long, aes(x=fct_rev(Class), y=EC.count, fill=network)) + 
 geom_boxplot(outlier.colour="black", outlier.fill="black", outlier.size=3) +
 scale_fill_manual(breaks=c("oxic","augmented","anoxic"), labels=c("Oxic","Augmented","Anoxic"), values=c("#bbc6e9","#8B8386","#fbd9d3")) +
 labs(x="Class", y="EC count per MAG", fill="Reactions") +
 coord_flip() + theme_minimal()
 




##### statistical analysis for paper
MAG.network.long


# MAGs with at least one gene from "oxic"
median.summary <- MAG.network.long %>% group_by(Class) %>% dplyr::summarise(oxic_pos_count = sum(network=="oxic" & EC.count>0, na.rm = TRUE), .groups = "drop")
median.summary$oxic_pos_count %>% sum
MAG.network.long$MAG %>% unique %>% length
MAG.network.long %>% group_by(Class) %>% dplyr::summarise(oxic_pos_count = sum(network=="oxic" & EC.count==0, na.rm = TRUE), .groups = "drop")

# out of the 577 CPR MAGs, 503 contained genes mapping to the oxic reaction network, 

# median EC counts for each Class and network
MAG.network.long %>% group_by(Class, network) %>% dplyr::summarise(median_EC = median(EC.count, na.rm = TRUE), .groups = "drop")

# with median counts of 2 to 7 genes per MAG
