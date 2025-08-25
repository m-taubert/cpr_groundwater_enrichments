# script for data import

# if anything goes wrong in the following scripts, try rerunning this here first

# import ASV table
asv <- readRDS("import/asv.RData")
nrow(asv)

# import metadata
meta <- read.csv(file="import/metadata.csv", header=TRUE, sep=";")
meta <- meta[,-1]
rownames(meta) <- meta$sample
meta <- meta %>% dplyr::rename(qPCR=X16S.gene.copies.per.L, qPCR.dev=dev)
summary(meta)

# check that metadata aligns with asv table
identical(meta$sample,rownames(asv))

# import taxonomy data
taxa.raw <- readRDS("import/taxa.RData")
taxa <- cbind(sequence=rownames(taxa.raw), data.frame(taxa.raw, row.names=NULL)) # put sequences in own row
taxa[1:10,2:7]
