# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data
data <- read.csv("unique_genus_per_gene.csv")

# Select the top 20 genes with the most unique genera
top_20_genes <- data %>% 
  arrange(desc(Number.of.Unique.Genus)) %>%
  head(20)

# Create the bar plot
ggplot(top_20_genes, aes(x = reorder(Gene, -Number.of.Unique.Genus), y = Number.of.Unique.Genus)) +
  geom_bar(stat = "identity", fill = "#AB8CBB") +
  labs(x = "Gene", y = "Number of Unique Genus", title = "Top 20 Genes with the Most Unique Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line.x = element_line(color = "darkgray"),  # Add X-axis line
    axis.line.y = element_line(color = "darkgray")   # Add Y-axis line
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjust y-axis to start at zero
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) + # Add space between Y-axis and bars
  theme(
    axis.ticks.length = unit(0.1, "cm"), # Length of ticks
    axis.ticks = element_line(color = "darkgray", size = 0.5) # Color and size of ticks
  )

ggsave("Genes_unique_genus.tiff", width = 12, height = 8)
ggsave("Genes_unique_genus.svg", width = 12, height = 8)
