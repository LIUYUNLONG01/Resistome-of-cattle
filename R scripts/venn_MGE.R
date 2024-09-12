# Install and load necessary packages
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(futile.logger)
library(readxl)

# Load the data
file_path <- "mge.xlsx"
df <- read_excel(file_path, sheet = 1)

# Clean the data by removing NAs and normalizing text
df <- df[!is.na(df$Genus), ]
df$Genus <- tolower(trimws(df$Genus))

# Extract subsets based on direction
positive_Genus <- unique(df[df$direction == "positive", "Genus"])
negative_Genus <- unique(df[df$direction == "negative", "Genus"])

# Convert tibbles to character vectors
positive_Genus <- as.character(positive_Genus$Genus)
negative_Genus <- as.character(negative_Genus$Genus)

# Now print and check the structures again
print("Positive Genera (Vector):")
print(positive_Genus)
print("Negative Genera (Vector):")
print(negative_Genus)

# Calculate the set sizes for the Venn diagram again
length_positive <- length(positive_Genus)
length_negative <- length(negative_Genus)
length_intersection <- length(intersect(positive_Genus, negative_Genus))

# Print the lengths to verify
print(paste("Length of Positive Set: ", length_positive))
print(paste("Length of Negative Set: ", length_negative))
print(paste("Length of Intersection: ", length_intersection))

tiff("venn_diagram_mge.tiff", width = 2500, height = 2500, res = 300)
svg("venn_diagram.svg", width = 8, height = 8)

# Create Venn diagram
venn_plot <- draw.pairwise.venn(
  area1 = length_positive,
  area2 = length_negative,
  cross.area = length_intersection,
  category = c("Positive", "Negative"),
  fill = c("#240046", "#ff6d00"),  
  lty = 1,                        
  lwd = 2,                        
  col = "darkgray",               
  cex = 0,                        
  cat.cex = 1.4,                   
  cat.col = c("darkblue", "red"),
  cat.pos = c(-30, 30),           
  cat.dist = 0.06,                
  rotation.degree = 0,            
  fontfamily = "serif"            
)

grid.draw(venn_plot)

positive_unique_Genus <- setdiff(positive_Genus, negative_Genus)

half_length <- ceiling(length(positive_unique_Genus) / 2)
Genus_first_half <- sort(positive_unique_Genus)[1:half_length]
Genus_second_half <- sort(positive_unique_Genus)[(half_length + 1):length(positive_unique_Genus)]

grid.text(
  label = paste(Genus_first_half, collapse = "\n"), 
  x = 0.1, 
  y = 0.5, 
  just = "left", 
  gp = gpar(fontsize = 7.5, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)

grid.text(
  label = paste(Genus_second_half, collapse = "\n"), 
  x = 0.22, 
  y = 0.5, 
  just = "left", 
  gp = gpar(fontsize = 7.5, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)

negative_unique_Genus <- setdiff(negative_Genus, positive_Genus)

grid.text(
  label = paste(sort(negative_unique_Genus), collapse = "\n"), 
  x = 0.76,  
  y = 0.5, 
  just = "center", 
  gp = gpar(fontsize = 7.5, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)

shared_Genus <- intersect(positive_Genus, negative_Genus)

grid.text(
  label = paste(sort(shared_Genus), collapse = "\n"), 
  x = 0.55,  
  y = 0.5, 
  just = "center", 
  gp = gpar(fontsize = 7.5, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)


grid.draw(venn_plot)

dev.off()

dev.cur()

