library("plyr")
library("dplyr")
library("stringr")

###################################################
##########          Data import          ##########
###################################################

# import functions table
MAG_annotations <- read.delim2(file="MAGs/dram_all_annotations2019.tsv", sep="\t", header=TRUE)
colnames(MAG_annotations) <- c("gene",colnames(MAG_annotations[2:22]))

# import taxonomy csv
MAG_taxonomy <- read.table(file="MAGs/Taxonomy_GW_MT.csv", sep=",", header=TRUE, fileEncoding='UTF-8-BOM')
MAG_taxonomy$Class[MAG_taxonomy$Class=="UBA1384"] <- "Berkelbacteria"
MAG_taxonomy$Class[MAG_taxonomy$Class=="Paceibacteria_A"] <- "Paceibacteria"
MAG_taxonomy_CPR <- MAG_taxonomy %>% filter(Phylum=="Patescibacteria")

# unify MAG names
MAG_annotations$MAG <- gsub("-contigs", "", MAG_annotations$fasta)

# filter for CPR, select MAG and KEGG ID, and filter for those with KEGG ID
CPR_annotations <- subset(MAG_annotations, MAG %in% subset(MAG_taxonomy, Phylum=="Patescibacteria")$sequence)[,c("MAG","kegg_id","kegg_hit")] %>% subset(kegg_id!="")

# extract EC number from annotations
CPR_annotations$EC <- str_extract(CPR_annotations$kegg_hit,"(?<=\\[EC:)[^]]+(?=\\])")

# import subnetwork database
subnetwork <- readRDS("MAGs/subnetwork.rData")

# remove annotations df, because they are large and we already have them subsetted now
rm("MAG_annotations")
