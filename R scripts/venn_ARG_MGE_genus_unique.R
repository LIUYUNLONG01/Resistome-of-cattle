

# Load the VennDiagram library
library(VennDiagram)

# Load the datasets
arg_data <- read.csv("unique_genus_gene_data_ARG.csv")
mge_data <- read.csv("unique_mge_genus_gene_data.csv")

# Extract unique genera
arg_genera <- unique(arg_data$Genus)
mge_genera <- unique(mge_data$Genus)


tiff("venn_arg_mge_genus.tiff", width = 2400, height = 2400, res = 300)
svg("venn_arg_mge_genus.svg", width = 12, height = 12)
png("venn_diagram_arg.png", width = 800, height = 800)

# Create a Venn diagram
venn.plot <- venn.diagram(
  x = list(
    "ARG Genera" = arg_genera,
    "MGE Genera" = mge_genera
  ),
  category.names = c("ARG Genera", "MGE Genera"),
  fill = c("#AB8CBB", "#C8DEDE"),  # Custom fill colors
  alpha = 0.8,                     # Transparency level
  lty = "solid",                  # Line type (dashed, solid, etc.)
  lwd = 2,                         # Line width
  col = c("darkgray", "darkgray"),   # Border line colors
  cat.col = c("#AB8CBB", "#C8DEDE"), # Category label colors
  cat.cex = 1.5,                   # Category label size
  cat.fontface = "bold",           # Category label font face
  cat.pos = c(0, 180), 
  cat.dist = 0.05,                 # Distance of labels from the diagram
  cex = 2,                         # Size of the numbers
  fontface = "bold",                 # Font face of the numbers
  fontfamily = "sans",               # Font family of the numbers
  filename = NULL,                   # If NULL, the plot is displayed instead of saved
  output = TRUE
)

# Plot the Venn diagram
grid.draw(venn.plot)

dev.off()

dev.cur()

