# Figure 4A
library(vegan)
library("tidyr")
library("ggplot2")
library("vegan")
library("forcats")
library("ggrepel")
library("purrr")
library("dplyr")
library("stringr")
library(pheatmap)
library(ComplexHeatmap)
library(tibble)
library(circlize)

#####     data preparation

# transform qPCR data
qPCRdata <- meta[,c(1,29,31)]
qPCRdata
qPCRdata$qPCR1 <- as.numeric(qPCRdata$qPCR)
qPCRdata$dev <- as.character(qPCRdata$PNK)
qPCRdata$dev1 <-as.numeric(qPCRdata$dev)

# add the genera as column names
taxa <- taxa %>%  unite("TaxaAll", Kingdom:Species, sep = "_", na.rm = FALSE, remove = FALSE)
asv1 <- asv
colnames(asv1) <- taxa$TaxaAll
asv1 <- asv1[,]

# select columns with Bacteria
Asvall <- asv1[, grep("Bacteria_", colnames(asv1))]

# remove ASVs with Bacteria_NA
Asvall_cleaned <- Asvall[, !grepl("Bacteria_NA", colnames(Asvall))]

# calculate relative abundances
Asvall_RA <- Asvall_cleaned / rowSums(Asvall_cleaned)
Asvall_RA <- Asvall_RA[apply(Asvall_RA, 1, function(x) !all(is.na(x))), ]

# add metadata to asv table
Asvall_RA1 <- data.frame(Asvall_RA)
Asvall_RA1$sample <- meta$sample
Asvall_RA1$filter <- meta$filter.fraction
Asvall_RA1$Csource <- meta$CsourceBeta
Asvall_RA1$Period <- meta$Period
Asvall_RA1$days <- meta$days
Asvall_RA1$folder <- meta$folder
Asvall_RA1$oxygen <- meta$oxygen
Asvall_RA1$thiosulfate <- meta$thiosulfate
Asvall_RA1$oxygen <- meta$oxygen
Asvall_RA1$ammonium <- meta$ammonium
Asvall_RA1$nitrite <- meta$nitrite
Asvall_RA1$nitrate <- meta$nitrate
Asvall_RA1$qPCR1 <- qPCRdata$qPCR1
Asvall_RA1$qdev1 <- qPCRdata$dev1
Asvall_RA1$biomass <- meta$biomass.enriched
Asvall_RA1$enrichment <- meta$enrichment

# remove samples from 0.1 fraction
Asvall_RAFilter02 <- Asvall_RA1[Asvall_RA1$filter == 0.2, ] 

# remove columns with 0 abundance after filtering
Asvall_RA02Remove0 <- Asvall_RAFilter02[, colSums(Asvall_RAFilter02 != 0, na.rm = TRUE) > 0, drop = FALSE]
Asvall_RAFilter02 <- Asvall_RA02Remove0



# collect data for ASV abundance of CPR and for supplements
Asvall_RA_Patesci <- Asvall_RAFilter02 %>% dplyr::select(all_of(union(c("sample","enrichment", "Csource", "Period", "days", "folder", "oxygen", "thiosulfate", "ammonium", "nitrite", "nitrate", "qPCR1", "qdev1"), grep("Patescibacteria", names(Asvall_RAFilter02), value = TRUE))))

# MT douple check if this works when I use my own df
# restructure supplement data
Asvall_RA_Patesci <- Asvall_RA_Patesci %>% 
 mutate(
  Anoxic = ifelse(oxygen == "anoxic", 1, 0),
  Oxic = ifelse(oxygen == "oxic", 1, 0),
  VeratricAcid = ifelse(enrichment %in% c("deuterated veratric acid","veratric acid"), 1, 0),
  Methylamine = ifelse(enrichment %in% c("methylamine", "Methylamine"), 1, 0),
  Methanol = ifelse(enrichment == "Methanol", 1, 0),
  SoilSeepage = ifelse(enrichment == "Soil Seepage", 1, 0),
  R2A = ifelse(enrichment == "R2A", 1, 0),
  Leachate = ifelse(enrichment == "Leachate", 1, 0),
  Starch = ifelse(enrichment == "starch", 1, 0),
  Necromass = ifelse(enrichment == "Necromass", 1, 0),
  Cellulose = ifelse(enrichment == "cellulose", 1, 0)
 )

rda_par <- Asvall_RA_Patesci[,c("Anoxic", "Oxic","VeratricAcid", "Methylamine", "Methanol", "SoilSeepage", "R2A", "Leachate", "Starch", "Necromass", "Cellulose", "thiosulfate", "ammonium", "nitrite", "nitrate")]
rda_par$thiosulfate <- as.numeric(rda_par$thiosulfate == "yes")
rda_par$ammonium    <- as.numeric(rda_par$ammonium == "yes")
rda_par$nitrite     <- as.numeric(rda_par$nitrite == "yes")
rda_par$nitrate     <- as.numeric(rda_par$nitrate == "yes")

# create df with only ASV data
asv_Pat_RA_PCoa <- Asvall_RA_Patesci %>% dplyr::select(all_of(grep("Patescibacteria", names(Asvall_RA_Patesci), value = TRUE)))

# combine ASV data with restructured supplement data
asv_Pat_RA_spearman <- cbind(rda_par,asv_Pat_RA_PCoa)

# calculate the correlation matrix
correlation_matrix <- cor(asv_Pat_RA_spearman, method = "spearman")

######################## determine p-values ################################
#Correlations were calculated using Spearman rank correlation coefficients (rS) with acceptance of significance at p < 0.05

# Initialize a matrix to store the p-values
p_value_matrix <- matrix(nrow = ncol(asv_Pat_RA_spearman), ncol = ncol(asv_Pat_RA_spearman))

#####     the following code needs a long time to run
#####     hence it is commented out
#####     rather continue with the stored matrix
# Loop over pairs of columns to test for correlations
#for(i in seq_len(ncol(asv_Pat_RA_spearman))){
#  for(j in seq_len(ncol(asv_Pat_RA_spearman))){
#    # Perform cor.test and store the p-value
#    test_result <- cor.test(asv_Pat_RA_spearman[, i], asv_Pat_RA_spearman[, j], method = "spearman")
#    p_value_matrix[i, j] <- test_result$p.value
#  }
#}

#saveRDS(p_value_matrix, file="Fig4_p_value_matrix.rData")

####     restore saved version of the data

p_value_matrix <- readRDS(file="data/Fig4_p_value_matrix.rData")

# rename columns and rows
colnames(p_value_matrix) <- colnames(correlation_matrix)
rownames(p_value_matrix) <- rownames(correlation_matrix)

# subset so only rows with ASVs are kept and cols with parameters
correlation_matrix <- correlation_matrix[!(rownames(correlation_matrix) %in% c("Anoxic","Oxic","VeratricAcid","Methylamine", "Methanol", "SoilSeepage", "R2A", "Leachate", "Starch", "Necromass", "Cellulose", "thiosulfate", "ammonium", "nitrite", "nitrate")), 1:15]
p_value_matrix <- p_value_matrix[!(rownames(p_value_matrix) %in% c("Anoxic","Oxic","VeratricAcid","Methylamine", "Methanol", "SoilSeepage", "R2A", "Leachate", "Starch", "Necromass", "Cellulose", "thiosulfate", "ammonium", "nitrite", "nitrate")), 1:15]

# mask low correlation values
p_value_matrix[correlation_matrix > -0.3 & correlation_matrix < 0.3] <- NA
correlation_matrix[correlation_matrix > -0.3 & correlation_matrix < 0.3] <- NA

# Keep only rows with any non-NA values
correlation_matrix <- correlation_matrix[rowSums(!is.na(correlation_matrix)) > 0, ]
correlation_matrix <- correlation_matrix[order(rownames(correlation_matrix)), ]
p_value_matrix <- p_value_matrix[rowSums(!is.na(p_value_matrix)) > 0, ]
p_value_matrix <- p_value_matrix[order(rownames(p_value_matrix)), ]

# perform Bonferroni correction
p_value_matrix_bonf <- p_value_matrix * sum(!is.na(p_value_matrix))

# redo removal of values with non-significant p-values
correlation_matrix[p_value_matrix_bonf >= 0.001] <- NA
p_value_matrix_bonf[p_value_matrix_bonf >= 0.001] <- NA

# remove ASVs without significant correlations
correlation_matrix <- correlation_matrix[rowSums(!is.na(correlation_matrix)) > 0, ]
correlation_matrix <- correlation_matrix[order(rownames(correlation_matrix)), ]
p_value_matrix_bonf <- p_value_matrix_bonf[rowSums(!is.na(p_value_matrix_bonf)) > 0, ]
p_value_matrix_bonf <- p_value_matrix_bonf[order(rownames(p_value_matrix_bonf)), ]

# Convert Bonferroni-corrected p-values to long format and remove NAs
p_bonf_long <- p_value_matrix_bonf %>% as.data.frame() %>% rownames_to_column("ASV") %>% pivot_longer(cols = -ASV, names_to = "supplement", values_to = "p_bonf") %>% filter(!is.na(p_bonf))

# Convert correlation matrix to long format and remove NAs
correlation_matrix_long <- correlation_matrix %>% as.data.frame() %>% rownames_to_column("ASV") %>% pivot_longer(cols = -c("ASV"), names_to = "supplement", values_to = "correlation") %>% filter(!is.na(correlation))

# Merge the filtered Bonferroni p-values with corresponding correlations
bonferrony <- merge(p_bonf_long, correlation_matrix_long, by.x = c("ASV", "supplement"), by.y = c("ASV", "supplement"))

bonferrony <- bonferrony %>% mutate(ASV=gsub("\\Candidatus.","Candidatus ", ASV))

#bonferrony$ASV1 <- paste0("X",sapply(bonferrony$ASV, function(x) which(colnames(asv)==x)))
#bonferrony$ASV1[353] <- paste0("X",which(colnames(asv)=="Bacteria_Patescibacteria_Parcubacteria_GWA2-38-13b_NA_NA_NA"))
#write.table(bonferrony, file="Table_S2_Spearmans.csv", sep=",")

#####     plotting
# Create an annotation data frame
annotation_df <- data.frame(rowname = gsub("Bacteria_Patescibacteria_", "", rownames(correlation_matrix)))
annotation_df$first_part <- sapply(strsplit(annotation_df$rowname, "_"), `[`, 1)
annotation_df$cluster <- annotation_df$first_part
rownames(correlation_matrix) <- annotation_df$rowname
annotation_df <- annotation_df[!grepl("NA", annotation_df$first_part), ]
correlation_matrix[is.na(correlation_matrix)] <- 0

# Define color palette
my_palette <- colorRamp2(c(-1, 0, 1), c("red", "white", "#018060"))

# define row order for supplements
custom_row_order <- c("ammonium","nitrate","nitrite","thiosulfate","Methylamine","Cellulose","Starch",  "VeratricAcid","Leachate","Necromass","R2A","SoilSeepage","Anoxic","Oxic")

# produce final data heatmap
FinalHeat <- t(correlation_matrix)    
FinalHeat <- FinalHeat[, !colnames(FinalHeat) %in% "NA_NA_NA_NA_NA"]

# create vector for splitting columns of the heatpam based on taxonomy.
new1_vector <- t(annotation_df$cluster)[1, ]

# Reorder the rows of the correlation matrix
reordered_matrix <- FinalHeat[match(custom_row_order, rownames(FinalHeat)), ]

# plot Figure 4A
Heatmap(reordered_matrix, col = my_palette, na_col = "grey", cluster_rows = FALSE, cluster_columns = FALSE, column_split = new1_vector, column_gap = unit(3, "mm"), show_row_names = TRUE, show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 8), column_title_side = "bottom", column_title_gp = gpar(fontsize = 14, fontface = "bold"), column_names_gp = gpar(fontsize = 14), heatmap_legend_param = list(title_position = "topleft", legend_direction = "vertical", legend_width = unit(6, "cm"), title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8)), top_annotation = HeatmapAnnotation(spacer = anno_empty(), height = unit(2, "mm"))) 
