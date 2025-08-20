# data preparation for Figure 6
# run this script before any other Figure 6 scripts
library(tidyr)
library(dplyr)


# transform qPCR data
qPCRdata <- meta[,c(1,29,31)]
#qPCRdata$qPCR <- gsub("^[^0-9\\.Ee]+|[^0-9\\.Ee]+$", "", qPCRdata$X16S.gene.copies.per.L)
#qPCRdata$qPCR <- as.character(qPCRdata$X16S.gene.copies.per.L)
qPCRdata$qPCR1 <- as.numeric(qPCRdata$qPCR)
qPCRdata$dev <- as.character(qPCRdata$PNK)
qPCRdata$dev1 <-as.numeric(qPCRdata$dev)

# add the genera as column names
taxa <- taxa %>%  unite("TaxaAll", Kingdom:Species, sep = "_", na.rm = FALSE, remove = FALSE)
asv1 <- asv
colnames(asv1) <- taxa$TaxaAll

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

