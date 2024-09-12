# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(showtext)

showtext.auto()

# Read the combined CSV file
data <- read.csv("combined_arg_mge_data.csv")

# Count the frequency of each genus
genus_frequency <- table(data$Genus)

# Filter data to include only genera with a frequency of 5 or more
frequent_genera <- names(genus_frequency[genus_frequency >= 5])
filtered_data <- subset(data, Genus %in% frequent_genera)

# Sort the data by Family and then by Genus
sorted_data <- filtered_data %>% arrange(Family, Genus)

# Determine the type of gene (ARG or MGE)
sorted_data$GeneType <- ifelse(grepl("ARG", sorted_data$Type), "ARG", "MGE")

# Create a pivot table (similar to the pivot table in Python)
pivot_data <- table(sorted_data$Genus, sorted_data$Gene)

# Convert the pivot table to a data frame for ggplot2
df <- as.data.frame(as.table(pivot_data))

# Merge with the gene type information
gene_types <- unique(sorted_data[, c("Gene", "GeneType")])
df <- merge(df, gene_types, by.x = "Var2", by.y = "Gene")

# Plot the heatmap
ggplot(df, aes(x = Var2, y = Var1, fill = Freq > 0)) +
  geom_tile(color = "gray") +
  scale_fill_manual(values = c("white", "#AB8CBB")) +
  facet_grid(~GeneType, scales = "free_x", space = "free_x") +
  labs(title = "Presence or Absence of Genes in Genera with Frequency ≥ 5 (Grouped by Family)",
       x = "Gene",
       y = "Genus") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, color = "black", family = "Arial"),
        axis.text.y = element_text(size = 12, color = "black", family = "Arial", face = "italic"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 16, face = "bold"))

ggsave("ARG-MGE-GENUS-HOST-HEATMAP.tiff", width = 32, height = 20, dpi = 300 )
ggsave("ARG-MGE-GENUS-HOST-HEATMAP.pdf", width = 32, height = 20, dpi = 300 )
ggsave("ARG-MGE-GENUS-HOST-HEATMAP1.svg", width = 32, height = 20, dpi = 300 )
