# script for Figure 9
library("dplyr")
library("tidyr")
library("ggplot2")

# import metabolic results
GW_metabolic <- read.csv(file="MAGs/1-Gw-Metabolic.csv", fileEncoding='UTF-8-BOM') # annotations
FunctionDescribe <- read.csv(file="MAGs/FunctionList.csv", header=TRUE, sep=",") # functional descriptions

# convert Presence/Abscence to 1 and 0 and transpose
GW_metabolic[,-1] <- lapply(GW_metabolic[,-1], function(x) ifelse(x == "Present", 1, 0))
GW_metabolic_t <- as.data.frame(t(GW_metabolic[,-1]))
names(GW_metabolic_t) <- GW_metabolic$FunctionAn

# unify MAG names
rownames(GW_metabolic_t) <- gsub("\\.bin", "-bin", rownames(GW_metabolic_t))

# assign Class and MAG names to columns
GW_metabolic_t$Class <- MAG_taxonomy$Class
GW_metabolic_t$sequence <- rownames(GW_metabolic_t)

# create subset for Saccharimonadia
subset_data_saccharimonadia <- subset(GW_metabolic_t, Class=="Saccharimonadia")
# create subset for Berkelbacteria
subset_data_berkelbacteria <- subset(GW_metabolic_t, Class=="Berkelbacteria")

# create subset for other CPR classes
subset_data_rest <- subset(GW_metabolic_t, Class %in% unique(subset(MAG_taxonomy, Phylum=="Patescibacteria")$Class)[c(-4,-6)])

##### perform Kruskal Wallis test for Saccharimonadia vs. other CPR (except Berkelbacteria) and vice versa ####

# initialize data frames for both groups
results_Sacc <- data.frame(Gene = character(), P_Value = numeric(), Statistic = numeric(), stringsAsFactors = FALSE)
results_Berk <- data.frame(Gene = character(), P_Value = numeric(), Statistic = numeric(), stringsAsFactors = FALSE)

# perform test for Saccharimonadia
for (value in GW_metabolic$FunctionAn) {
 combined_df <- rbind(data.frame(Gene = subset_data_saccharimonadia[[value]], Class = "Saccharimonadia"), data.frame(Gene = subset_data_rest[[value]], Class = "Rest"))
 
 combined_df$Gene <- as.numeric(as.character(combined_df$Gene))
 
 test_result <- kruskal.test(as.formula(paste("Gene ~ Class")), data = combined_df)
 
 temp_df <- data.frame(Gene = value, P_Value = test_result$p.value, Statistic = test_result$statistic)
 
 results_Sacc <- rbind(results_Sacc, temp_df)}

# perform test for Berkelbacteria
for (value in GW_metabolic$FunctionAn) {
 combined_df <- rbind(data.frame(Gene = subset_data_berkelbacteria[[value]], Class = "Berkelbacteria"), data.frame(Gene = subset_data_rest[[value]], Class = "Rest"))
 
 combined_df$Gene <- as.numeric(as.character(combined_df$Gene))
 
 test_result <- kruskal.test(as.formula(paste("Gene ~ Class")), data = combined_df)
 
 temp_df <- data.frame(Gene = value, P_Value = test_result$p.value, Statistic = test_result$statistic)
 
 results_Berk <- rbind(results_Berk, temp_df)}

# filter results
results_Sacc_sig <- subset(results_Sacc, P_Value <= 0.05)
results_Berk_sig <- subset(results_Berk, P_Value <= 0.05)

# combine lists of significant genes
significant_genes <- c(results_Sacc_sig$Gene, results_Berk_sig$Gene) %>% unique

# Filter dataset to keep only relevant columns
CPR_AfterKruskal <- GW_metabolic_t[which(MAG_taxonomy$Phylum=="Patescibacteria"), c("sequence", "Class", significant_genes), drop = FALSE]

# transform into long format
CPR_AfterKruskal_long <- CPR_AfterKruskal %>% pivot_longer(cols = starts_with("Function"), names_to = "Function", values_to = "Presence")

# filter out unwanted Classes
CPR_AfterKruskal_long <- CPR_AfterKruskal_long %>%  filter(!Class %in% c("CPR2", "Andersenbacteria", "WWE3", "Doudnabacteria"))

# combine with descriptions
CPR_AfterKruskal_long1 <- merge(CPR_AfterKruskal_long,FunctionDescribe,by.x = "Function", by.y = "FunctionAn", all.x = TRUE )

# determine number of MAGs per CPR class
class_total <- GW_metabolic_t[which(MAG_taxonomy$Phylum=="Patescibacteria"),c("sequence","Class")] %>% group_by(Class) %>% dplyr::summarize(TotalMags = n(), .groups = 'drop')

# add information to annotation data
filtered_data <- merge(CPR_AfterKruskal_long1,class_total,by = "Class", all.x = TRUE )

# filter data
filtered_data_unique <- filtered_data %>% group_by(Class,sequence) %>% distinct(Gene.abbreviation,Corresponding.KO, .keep_all = TRUE)

# calculate percent of MAGs per class with particular function
Category_Sum <-  filtered_data_unique  %>% group_by(Class,Category,Gene.abbreviation,TotalMags,Corresponding.KO) %>% dplyr::summarize(Sum_Presence = sum(Presence), .groups = 'drop')
Category_Sum$Percent <- Category_Sum$Sum_Presence / Category_Sum$TotalMags *100

# reformat descriptions
Category_Sum <- Category_Sum %>% mutate(Category_Gene = paste(Category, Gene.abbreviation, sep = " - "))
Category_Sum <- Category_Sum %>% mutate(Category_GeneKO = paste(Category_Gene,Corresponding.KO, sep = " - ")) 
Category_Sum_f <- Category_Sum %>% group_by(Class, Category,Category_Gene,Corresponding.KO) %>% dplyr::summarize(Percent1 = max(Percent, na.rm = TRUE))

# Filter the dataset to remove repeated Corresponding.KO rows and select specific functions for plotting
Fig9.filtered <- Category_Sum_f %>% group_by(Category_Gene,Corresponding.KO,Class) %>% filter(row_number() == 1) %>% ungroup()

Fig9.filtered <- Fig9.filtered %>% mutate(Category_GeneKO = paste(Category_Gene,Corresponding.KO, sep = " - "))

# List of functions of interest was manually derived from initial results
Fig9.filtered_forPlot <- Fig9.filtered %>% filter(Category_Gene %in% c('Central carbohydrate metabolism - Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate',
'Aromatic amino acid metabolism - Phenylalanine biosynthesis, chorismate => phenylpyruvate => phenylalanine',
'Nitrogen metabolism - Denitrification, nitrate => nitrogen',
'Cofactor and vitamin metabolism - Cobalamin biosynthesis, cobyrinate a,c-diamide => cobalamin',
'Cofactor and vitamin metabolism - NAD biosynthesis, aspartate => quinolinate => NAD',
'Cofactor and vitamin metabolism - NAD biosynthesis, tryptophan => quinolinate => NAD',
'Polyamine biosynthesis - Polyamine biosynthesis, arginine => agmatine => putrescine => spermidine',
'Cofactor and vitamin metabolism - Molybdenum cofactor biosynthesis, GTP => molybdenum cofactor',
'Carbon fixation - Reductive pentose phosphate cycle, glyceraldehyde-3P => ribulose-5P',
'Branched-chain amino acid metabolism - Isoleucine biosynthesis, threonine => 2-oxobutanoate => isoleucine',
'Central carbohydrate metabolism - Pentose phosphate pathway (Pentose phosphate cycle)',
'Cofactor and vitamin metabolism - Pimeloyl-ACP biosynthesis, BioC-BioH pathway, malonyl-ACP => pimeloyl-ACP',
'Cofactor and vitamin metabolism - Heme biosynthesis, bacteria, glutamyl-tRNA => coproporphyrin III => heme',
'Arginine and proline metabolism - Creatine pathway',
'Sulfur cycling enzymes (detailed) - metE',
'Central carbohydrate metabolism - Pentose phosphate pathway, non-oxidative phase, fructose 6P => ribose 5P',
'Amino acid utilization - ornithine/acetylornithine aminotransferase',
'Amino acid utilization - phosphoserine aminotransferase',
'Pyrimidine metabolism - Pyrimidine deoxyribonucleotide biosynthesis, UDP => dTTP',
'Other carbohydrate metabolism - Nucleotide sugar biosynthesis, galactose => UDP-galactose',
'Lipid metabolism - Phosphatidylethanolamine (PE) biosynthesis, PA => PS => PE',
'Complex carbon degradation - alpha-amylase',
'Complex carbon degradation - glucoamylase',
'Polyamine biosynthesis - Polyamine biosynthesis, arginine => ornithine => putrescine',
'Terpenoid backbone biosynthesis - C5 isoprenoid biosynthesis, mevalonate pathway',
'Arginine and proline metabolism - Arginine biosynthesis, ornithine => arginine',
'Oxidative phosphorylation - atpA (F-type)',
'Oxidative phosphorylation - atpD (F-type)',
'Oxygen metabolism (Oxidative phosphorylation Complex IV) - cyoA',
'Oxygen metabolism (Oxidative phosphorylation Complex IV) - cyoB',
'Oxygen metabolism (Oxidative phosphorylation Complex IV) - cyoC',
'Oxygen metabolism (Oxidative phosphorylation Complex IV) - cyoD',
'Oxygen metabolism (Oxidative phosphorylation Complex IV) - qoxA',
'Lysine metabolism - Lysine degradation, bacteria, L-lysine => succinate',
'Sterol biosynthesis - Cholesterol biosynthesis, squalene 2,3-epoxide => cholesterol',
'Pathogenicity - Vibrio cholerae pathogenicity signature, toxin coregulated pilus',
'Macrolide biosynthesis - Avermectin biosynthesis, 2-methylbutanoyl-CoA/isobutyryl-CoA => 6,8a-Seco-6,8a-deoxy-5-oxoavermectin 1a/1b aglycone => avermectin A1a/B1a/A1b/B1b',
'Carbon fixation - aclA',
'Carbon fixation - aclB',
'Methane metabolism - mmoB',
'Biosynthesis of other antibiotics - Roseoflavin biosynthesis, FMN => roseoflavin',
'Iron cycling - DmkB',
'Iron cycling - Ndh2'))

# Select unique rows based on specific columns
Fig9.plotting <- Fig9.filtered_forPlot %>% dplyr::select(Category_Gene, Class, Corresponding.KO,Percent1) %>% dplyr::distinct()

# prepare and export Supplementary Table S6
Table.S6 <- pivot_wider(Fig9.plotting, names_from=Class,values_from=Percent1)
write.table(Table.S6,file="Table_S6.csv", dec=",",sep=";")
# this is the data shown in Table S6. In the table, duplicates were removed, overlapping functions were combined, and additional annotation data based on KO numbers was collected from KEGG

# select colors for Classes in the plot
custom_colors <- c("WWE3"="grey","Saccharimonadia" = "#0056BE","Paceibacteria" = "#959BA3","NA" = "grey","Microgenomatia" = "#C7C9CB","Kazania"="grey","Gracilibacteria" = "#AEB2B8","CPR2"="grey","UBA1384"="#4D009B","ABY1" = "#D7D7D8") 

# plot figure 9
ggplot(Fig9.plotting, aes(x = Percent1, y = Category_Gene, fill = Class)) +
 geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = custom_colors) +  
 theme_minimal() + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 3)) +
 labs(fill = "Class")

