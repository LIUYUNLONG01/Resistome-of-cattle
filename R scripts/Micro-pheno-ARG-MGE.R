# Load required libraries
library(readxl)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(openxlsx)
library(viridis)
library(viridisLite)
library(grid)

# Load the filtered data with significance annotations
data <- read.xlsx("filtered_microbial_phenotype_data_with_significance.xlsx")

# Remove specified columns
data <- data[!(data$Parameter %in% c("Total solid", "istA11", "istA13", "istA15", "istA17")), ]

# Specify the order for the vertical axis
parameter_order <- c("Methane", "insertion_element_IS91", "integrase", "ISCR", "plasmid", 
                     "Tn916", "transposase", "aminoglycoside", "bacitracin", "fusidic_acid", 
                     "multidrug", "mupirocin", "quinolone", "tetracycline", "vancomycin")

# Ensure the specified order is respected
data$Parameter <- factor(data$Parameter, levels = parameter_order)
data <- data[order(data$Parameter), ]

# Create a pivot table for the heatmap with microbial data on the horizontal axis
heatmap_data <- dcast(data, Parameter ~ Microbiome, value.var = "Correlation")

# Replace NA values with a number for masking
heatmap_matrix <- as.matrix(heatmap_data[,-1])
rownames(heatmap_matrix) <- heatmap_data[,1]
heatmap_matrix[is.na(heatmap_matrix)] <- NA

svg("micro-pheno-arg-mge-heatmap.svg", width = 12, height = 8)
tiff("micro-pheno-arg-mge-heatmap.tiff", width = 1200, height = 800)

# Plot the heatmap
pheatmap(heatmap_matrix, 
         color = colorRampPalette(c("#C8DEDE", "#AB8CBB"))(100),
         main = "Heatmap of Correlation (|R| >= 0.5 and PValueAdjusted <= 0.05)",
         border_color = "grey", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         labels_row = rownames(heatmap_matrix),
         na_col = "white",
         fontsize_row = 12,  # Row label font size
         fontsize_col = 12,
         label_col = grid::gpar(col = "black", fontface = "italic"),
         label_row = grid::gpar(col = "black", fontface = "bold"))  

dev.off()

dev.cur()

