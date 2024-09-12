# Load necessary libraries
library(ggplot2)

# Read the CSV file with the sorted data
data <- read.csv("unique_gene_names_arg_mge_data.csv")

# Select the top 20 genera based on ARG_MGE_Sum
top20_data <- head(data, 20)

# Convert Genus to a factor to maintain the order in the plot
top20_data$Genus <- factor(top20_data$Genus, levels = top20_data$Genus)

# Create a bar plot
ggplot(top20_data, aes(x = Genus)) +
  geom_bar(aes(y = ARG), stat = "identity", position = position_nudge(x = -0.2), width = 0.4, fill = "#AB8CBB") +
  geom_bar(aes(y = MGE), stat = "identity", position = position_nudge(x = 0.2), width = 0.4, fill = "#C8DEDE") +
  labs(title = "Top 20 Genera by ARG and MGE Counts", x = "Genus", y = "Count") +
  theme_minimal() +
  theme(
    # Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Set white background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Rotate x-axis text
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    # Customize axis titles
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    # Customize plot title
    plot.title = element_text(size = 14, face = "bold", color = "black"),
    # Add and customize axis lines
    axis.line = element_line(color = "darkgray"),
    axis.ticks = element_line(color = "darkgray")
  ) +
  scale_y_continuous(expand = c(0, 0))

ggsave("ARG_MGR_GENUS_COUNTS.tiff", width = 12, height = 8)
ggsave("ARG_MGR_GENUS_COUNTS.svg", width = 12, height = 8)

