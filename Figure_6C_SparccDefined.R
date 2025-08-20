# Figure 6C
# make sure to run Figure_6_prep.R first
library(vegan)
library("tidyr")
library("ggplot2")
library("vegan")
library("forcats")
library("purrr")
library("dplyr")
library("stringr")
library(devtools)
library(SpiecEasi)
library(igraph)
library(reshape2)
library("circlize")

############################ Select treatment
AsvSelect <- Asvall_RAFilter02[grep("^organic_defined$", Asvall_RAFilter02$Csource, ignore.case = TRUE), ]
abundance_matrixAsv <- AsvSelect[, grep("Bacteria", names(AsvSelect))]
abundance_matrixAsv[] <- lapply(abundance_matrixAsv, function(x) as.numeric(as.character(x)))
# scale to percentages
selected_rows_multiplied <- abundance_matrixAsv * 100

# Subset the data frame to keep only those columns
filtered_columns <- selected_rows_multiplied[, colSums(selected_rows_multiplied)>0]
abundance_matrixAsv2 <- as.matrix(filtered_columns)
funyG <- apply(abundance_matrixAsv2, 2, as.numeric)

############################# backup and run sparcc analysis
#save(funyG, file = "funyGdefined.RData")

# this code takes a long time to run, please use the provided .RData files instead
#net <-sparcc(funyG)
#save(net, file = "netorganic_defined.RData")

####################################################### Plot
load("sparcc/netorganic_defined.RData")

edg_graph <- net$Cor
colnames(edg_graph) <- colnames(abundance_matrixAsv2)
rownames(edg_graph) <- colnames(abundance_matrixAsv2)

# Find the row and column indices that contain "Patescibacteria"
patescib_indices <- which(grepl("Patescibacteria", colnames(abundance_matrixAsv2)))

# Subset the matrix to include only those rows 
patescib_edges <- edg_graph[patescib_indices, ]

# Remove "Bacteria_" and NAs prefix from row and column names
rownames(patescib_edges) <- gsub("Bacteria_", "", rownames(patescib_edges))
colnames(patescib_edges) <- gsub("Bacteria_", "", colnames(patescib_edges))
rownames(patescib_edges) <- gsub("_NA_NA_NA_NA", "", rownames(patescib_edges))
colnames(patescib_edges) <- gsub("_NA_NA_NA_NA", "", colnames(patescib_edges))

# Filter the edges 
filtered_patescib_edges0 <-  melt(patescib_edges, varnames = c("sourcenode", "targetnode"), value.name = "edgeweight")
filtered_patescib_edges1 <- filtered_patescib_edges0[!(grepl("Patescibacteria", filtered_patescib_edges0$sourcenode) & grepl("Patescibacteria", filtered_patescib_edges0$targetnode)), ]
filtered_patescib_edges <- subset(filtered_patescib_edges1, edgeweight > 0.1)
filtered_patescib_edges$edgeweight <- as.numeric(filtered_patescib_edges$edgeweight)
filtered_patescib_edges <- filtered_patescib_edges %>% mutate(sourcenode = as.character(sourcenode), targetnode = as.character(targetnode)) 

# Identify duplicate edgeweights and keep only the first occurrence
filtered_patescib_edges <- filtered_patescib_edges %>% mutate(edgeweight = as.numeric(edgeweight), sourcenode = as.character(sourcenode), targetnode = as.character(targetnode)) %>% mutate(edgeweight = round(edgeweight, digits = 7))  
filtered_patescib_edges_no_duplicates <-  filtered_patescib_edges %>%  distinct(edgeweight, .keep_all = TRUE)

# Create new columns for taxonomy
cleaned_table <- filtered_patescib_edges_no_duplicates %>% mutate(sourcenode = sub("^[^_]*_([^_]*)_([^_]*)_.*", "\\1_\\2", sourcenode), targetnode = sub("^[^_]*_([^_]*)_([^_]*)_.*", "\\1_\\2", targetnode))

# Remove "Patescibacteria" from the sourcenode
cleaned_table <- cleaned_table %>% mutate(sourcenode = gsub("^Patescibacteria_", "", sourcenode))

# Remove the dot and following number from sourcenode
filtered_plot1 <- cleaned_table %>% mutate(sourcenode = gsub("\\.[0-9]+$", "", sourcenode), targetnode = gsub("\\.[0-9]+$", "", targetnode))
filtered_plot <- filtered_plot1 %>% mutate(sourcenode = gsub("\\.[0-9]+$", "", sourcenode), targetnode = gsub("\\.[0-9]+$", "", targetnode))
filtered_plot$sourcenode <- paste0("CPR_", filtered_plot$sourcenode)

############################ Select taxonomy and assign colors

# Define the list of targetnodes to filter
targetnode_filter <- c("Alphaproteobacteria_Caulobacterales","Alphaproteobacteria_Rhizobiales","Alphaproteobacteria_Rickettsiales","Alphaproteobacteria_Sphingomonadales","Bacteroidia_Chitinophagales","Bacteroidia_Cytophagales","Bacteroidia_Flavobacteriales","Bacteroidia_Sphingobacteriales","Blastocatellia_Blastocatellales","Chlamydiae_Chlamydiales","Chloroflexi_KD4","Gammaproteobacteria_Burkholderiales","Gammaproteobacteria_Diplorickettsiales","Gammaproteobacteria_Pseudomonadales","Gammaproteobacteria_Xanthomonadales","Hydrogenedentia_Hydrogenedentiales","Nitrospiria_Nitrospirales","Planctomycetes_Gemmatales","Planctomycetes_Pirellulales","Planctomycetes_Planctomycetales","Planctomycetota_OM190","Verrucomicrobiae_Verrucomicrobiales")
# Create a named vector of colors for the sectors
color_assignments <- list("Parcubacteria" = "#BF0000","Gracilibacteria" = "#739900","Saccharimonadia" = "#0056BE","ABY1" = "#BF7300","Gammaproteobacteria" = "#538B8B","Alphaproteobacteria" = "#88789b","Bacteroidia" = "#e75480","Verrucomicrobiae" = "#FFB6C1","Chlamydiae" = "#F0E68C","Acidimicrobiia|Blastocatellia" = "darkslateblue","Planctomycetes|Planctomycetota" = "darkgreen","Nitrospirota|Nitrospiria" = "#96a9c6","WPS" = "#7F798F","Chloroflexi" = "bisque4","Berkelbacteria" = "darkslategrey","Cyanobacteria" = "darkseagreen","Actinobacteria|Actinomycetota" = "#8B7D6B","Thermoleophilia"= "#CDC0B0","Bdellovibrionia" ="#8B8970","Omnitrophia" = "bisque4","Brocadiae" = "#8A2BE2","Vampirivibrionia" = "#8B8970","Polyangia" = "#4682B4","Hydrogenedentia" = "#8B8878","Latescibacterota" ="#8B795E","Campylobacteria" ="#8B7D7B")

# Filter rows where targetnode matches the specified values
filtered_plot <- filtered_plot[filtered_plot$targetnode %in% targetnode_filter, ]
filtered_plot <- filtered_plot %>% filter(targetnode != "NA_NA" & sourcenode != "NA_NA")

# Create a named color vector for each unique node
unique_nodes <- unique(c(filtered_plot$sourcenode, filtered_plot$targetnode))
grid_colors <- rep("grey", length(unique_nodes))  # Default color is grey
names(grid_colors) <- unique_nodes

# Assign colors based on the node name
for (group in names(color_assignments)) {
 grid_colors[grepl(group, unique_nodes, ignore.case = TRUE)] <- color_assignments[[group]]
}

# Sort nodes alphabetically
sorted_nodes <- sort(unique(c(filtered_plot$sourcenode, filtered_plot$targetnode)))

# Reorder filtered_plot so that nodes appear alphabetically
filtered_plot <- filtered_plot %>% arrange(factor(sourcenode, levels = sorted_nodes), factor(targetnode, levels = sorted_nodes))

# plot figure 6C
chordDiagram(x = filtered_plot[, c("sourcenode", "targetnode")], grid.col = grid_colors, transparency = 0.5, order = sorted_nodes)

# Add additional labels 
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
 sector_name = get.cell.meta.data("sector.index")
 circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector_name, 
             facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 10/12)  
}, bg.border = NA)  
# comment: in the manuscript figure, only the 200 co-occurrences with the highest correlation value were selected, and only non-CPR orders with more than two co-occurrences were included
