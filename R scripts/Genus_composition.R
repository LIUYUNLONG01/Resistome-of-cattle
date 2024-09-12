# Load necessary library
library(ggplot2)

# Load the dataset
data <- read.csv("combined_data_G_.csv")

# Calculate the sum of each genus across all samples
genus_abundance <- colSums(data[,-1], na.rm = TRUE)

# Calculate relative abundance
total_abundance <- sum(genus_abundance)
relative_abundance <- genus_abundance / total_abundance

# Create a data frame for plotting
df <- data.frame(
  genus = names(relative_abundance),
  Abundance = relative_abundance
)

# Sort by abundance
df <- df[order(-df$Abundance),]

# Select the top 20 genus
top_20 <- df[1:20,]

# Plotting the top 20 genus by relative abundance
ggplot(top_20, aes(x = reorder(genus, -Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "#AB8CBB", color = "black", width = 0.8) +
  scale_x_discrete(expand = c(0.03, 0.03)) + 
  theme_minimal() +
  labs(title = "Genus Composition (Top 20 by Relative Abundance)",
       x = " ",
       y = "Relative Abundance") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"), # Customize X axis text
    axis.text.y = element_text(size = 10, color = "black"),  # Customize Y axis text
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),  # Customize Y axis title
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add closed box
    panel.grid = element_blank()  # Remove grid lines
  )

# Save the plot as SVG
ggsave("genus_composition.svg", width = 16, height = 8, dpi = 300)

# Save the plot as PDF
ggsave("genus_composition.pdf", width = 16, height = 8, dpi = 300)
