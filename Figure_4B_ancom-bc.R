# Figure 4B
library(tidyr)
library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(circlize)
library(plyr)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)

#####     selection of data

# remove all 0.1 µm samples
sample.removed <- which(meta$filter.fraction=="0.1" | meta$Period=="day0") 
ASVs.removed <- which(colSums(asv[-sample.removed,])==0)

# Filter for CPR
selected.asvs <- which(taxa$Phylum=="Patescibacteria" & colSums(asv[-sample.removed,])!=0)
selected.asvs # so this is the vector of all "good" ASVs

# the reduced asv data frame
asv.select <- asv[-sample.removed,selected.asvs]
# all Patescibacteria
taxa.select <- taxa[selected.asvs,]
# the reduced metadata table
meta.select <- meta[-sample.removed,]

#########################################
#####     start of the analysis     #####
#########################################
# prepare data for ancombc analysis
asv.phylo <- otu_table(as.matrix(asv.select), taxa_are_rows = FALSE)

meta.phylo <- sample_data(meta.select)

phylo.seq <- phyloseq(asv.phylo, meta.phylo)

# perform ancombc analysis
#####     be adviced that the following code need a substantial time to run
res.ammonium <- ancombc2(phylo.seq, fix_formula="ammonium", prv_cut=0.0)

res.nitrate <- ancombc2(phylo.seq, fix_formula="nitrate", prv_cut=0.0)

res.nitrite <- ancombc2(phylo.seq, fix_formula="nitrite", prv_cut=0.0)

res.thiosulfate <- ancombc2(phylo.seq, fix_formula="thiosulfate", prv_cut=0.0)

res.methylamine <- ancombc2(phylo.seq, fix_formula="methylamine", prv_cut=0.01)

res.oxygen <- ancombc2(phylo.seq, fix_formula="oxygen", prv_cut=0.0)

# the remaining parameters did not yield ASVs with significant correlations
#res.cellulose <- ancombc2(phylo.seq, fix_formula="cellulose", prv_cut=0.0)
#res.starch <- ancombc2(phylo.seq, fix_formula="starch", prv_cut=0.0)
#res.veratric.acid <- ancombc2(phylo.seq, fix_formula="veratric.acid", prv_cut=0.0)
#res.soil.seepage <- ancombc2(phylo.seq, fix_formula="soil.seepage", prv_cut=0.0)
#res.necromass <- ancombc2(phylo.seq, fix_formula="necromass", prv_cut=0.0)
#res.R2A <- ancombc2(phylo.seq, fix_formula="R2A", prv_cut=0.0)
#res.leaf.leachate <- ancombc2(phylo.seq, fix_formula="leaf.leachate", prv_cut=0.0)


########################################
#####     combining of results     #####
########################################

# remove correlations that are not significant
res.ammonium$res[res.ammonium$res[[11]] > 0.05, 3] <- NA
res.nitrate$res[res.nitrate$res[[11]] > 0.05, 3] <- NA
res.nitrite$res[res.nitrite$res[[11]] > 0.05, 3] <- NA
res.thiosulfate$res[res.thiosulfate$res[[11]] > 0.05, 3] <- NA
res.oxygen$res[res.oxygen$res[[11]] > 0.05, 3] <- NA
res.methylamine$res[res.methylamine$res[[11]] > 0.05, 3] <- NA

# merge data (those not included have no significant correlations)
anc.results.temp <- cbind(res.ammonium$res[,c(1,3,11)], res.nitrate$res[,c(3,11)], res.nitrite$res[,c(3,11)], res.thiosulfate$res[,c(3,11)], res.oxygen$res[,c(3,11)])
anc.results <- join(anc.results.temp, res.methylamine$res[,c(1,3,11)], by="taxon")

# add taxonomy to the data frame, the order is the same
anc.results <- cbind(anc.results,Class=taxa.select$Class)

# find rows which only have NAs for removal of ASVs without correlations
anc.results.sub <- anc.results[which(rowSums(!is.na(anc.results[c(2,4,6,8,10,12)]))>0),]

#####     creating heatmap     #####

# sort it based on taxonomy
anc.results.sub <- anc.results.sub[order(anc.results.sub$Class),]


anc.results.sub.long <- anc.results.sub %>% pivot_longer(cols = -c(taxon, Class), names_to = c(".value", "supplement"), names_pattern = "(lfc|q)_(.*)")

write.table(anc.results.sub.long, file="Table_S4_ANCOM-BC.csv")

#define colors
hm.col <- colorRamp2(c(-2, 0, 2), c("red", "white", "#018060")) 
# plot heatmap for figure 4B
pheatmap(anc.results.sub[,c(2,4,6,8,12,10)], color=hm.col, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = TRUE, na_col = "white", annotation_row=anc.results.sub[,c("Class","Class")], annotation_colors=list(Class = c("ABY1" = "#D29D4D", "Berkelbacteria" = "#844FB9", "Gracilibacteria" = "#9DB84D",  "Parcubacteria" = "#D24D4D", "Saccharimonadia" = "#4D89D2")))

# check number of correlating ASVs for paper text
subset(anc.results.sub[,c(2,4,6,8,12,10)], lfc_nitrateyes>0) %>% nrow
subset(anc.results.sub[,c(2,4,6,8,12,10)], lfc_thiosulfateyes>0) %>% nrow
subset(anc.results.sub[,c(2,4,6,8,12,10)], lfc_ammoniumyes>0) %>% nrow
subset(anc.results.sub, Class %in% c("Parcubacteria") & (lfc_nitrateyes>0 | lfc_thiosulfateyes>0 | lfc_ammoniumyes>0)) %>% nrow
subset(anc.results.sub, Class %in% c("ABY1") & (lfc_nitrateyes>0 | lfc_thiosulfateyes>0 | lfc_ammoniumyes>0)) %>% nrow
subset(anc.results.sub, lfc_methylamineyes<0) %>% nrow
subset(anc.results.sub, Class %in% c("ABY1") & lfc_methylamineyes>0) %>% nrow

subset(anc.results.sub, lfc_oxygenoxic>0) %>% nrow

