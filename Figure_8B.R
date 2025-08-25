library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
#library("forcats")


###################################################
##########          Plot Figure 8B       ##########
###################################################
# run script Figure_8_prep.R and Figure_8A.R before running this script


# add class information to the annotations for genes from the oxic network
CPR_annotations.oxic <- left_join(subset(CPR_annotations.EC, network=="oxic"), MAG_taxonomy_CPR[,c("sequence","Class")], by=c("MAG"="sequence"))

# define classes to be kept
CPR.class.select <- (CPR_annotations.oxic$Class %>% unique)[c(2,3,5,6,1,4)]

# count occurrence of kegg_ids in MAGs of these classes
CPR_annotations.counts <- subset(CPR_annotations.oxic, Class %in% CPR.class.select) %>% group_by(Class, kegg_id, across(-c(MAG, Class))) %>% dplyr::summarise(count = n_distinct(MAG), .groups = "drop")

# calculate percentage of MAGs per Class with functions
MAG_taxonomy_CPR.count <- plyr::count(subset(MAG_taxonomy_CPR, Class %in% CPR.class.select)$Class)
colnames(MAG_taxonomy_CPR.count) <- c("Class","total")
CPR_annotations.counts <- CPR_annotations.counts %>% left_join(MAG_taxonomy_CPR.count, by = "Class") %>% mutate(perc = count / total*100)

# pivot to wide format
CPR_annotations.wide <- CPR_annotations.counts %>% select(-"total") %>% pivot_wider(names_from = Class, values_from = c(count,perc), values_fill = list(count=0, perc=0))
CPR_annotations.wide

# define KEGG IDs of interest
KO.Fig8B <- c("K00273","K00278","K00457","K00486","K00544","K00812","K01426","K01692","K07130","K09471","K01091","K01601","K03821","K05349","K05350","K02847","K11936","K12983","K12990","K12992","K14340","K07104","K14292","K22083","K00275","K00588","K01434","K17836","K00395","K00956","K00957","K00958","K01082")
intersect(CPR_annotations.wide$kegg_id,KO.Fig8B) %>% length # all 33 are present

# this is the data shown in table S5
write.table(subset(CPR_annotations.wide, kegg_id %in% KO.Fig8B), file="Table_S5.csv",sep=",",dec=".")

# plot figure 8B
# subset(CPR_annotations.counts, kegg_id %in% KO.Fig8B) # means only plot relevant KOs

ggplot(subset(CPR_annotations.counts, kegg_id %in% KO.Fig8B), aes(x=Class, y=kegg_hit, fill=perc)) +
 geom_tile() + theme_minimal() + theme(panel.grid=element_blank()) + labs(x="", y="") +
 scale_fill_gradientn(name="% of MAGs", colors=c("#f0f4fb","#bbc6e9","#7a8cc2","#5063a1","#2f356b"), values=scales::rescale(c(0,1,10,30,60)), breaks=c(0,1,10,30,60), labels=c("0","1","10","30","60"), oob = scales::squish) +
 theme(axis.text.x=element_text(angle=60, hjust=1), axis.text.y=element_text(angle=0, hjust=1)) +
 coord_fixed(ratio=1)
# y-axis descriptions in the paper have been modified manually. Values can be found in Table S5
