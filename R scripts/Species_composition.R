# Load necessary library
library(ggplot2)

# Load the dataset
data <- read.csv("combined_data_S.csv")

# Calculate the sum of each species across all samples
species_abundance <- colSums(data[,-1], na.rm = TRUE)

# Calculate relative abundance
total_abundance <- sum(species_abundance)
relative_abundance <- species_abundance / total_abundance

# Create a data frame for plotting
df <- data.frame(
  species = names(relative_abundance),
  Abundance = relative_abundance
)

# Sort by abundance
df <- df[order(-df$Abundance),]

# Select the top 20 species
top_20 <- df[1:20,]

# Plotting the top 20 species by relative abundance
ggplot(top_20, aes(x = reorder(species, -Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "#AB8CBB", color = "black", width = 0.8) +
  scale_x_discrete(expand = c(0.03, 0.03)) + 
  theme_minimal() +
  labs(title = "Species Composition (Top 20 by Relative Abundance)",
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
ggsave("species_composition.svg", width = 16, height = 8, dpi = 300)

# Save the plot as PDF
ggsave("species_composition.pdf", width = 16, height = 8, dpi = 300)
