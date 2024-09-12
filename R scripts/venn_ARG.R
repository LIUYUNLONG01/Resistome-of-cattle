# Install and load necessary packages
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
library(futile.logger)
library(readxl)

# Load the data
file_path <- "arg.xlsx"
df <- read_excel(file_path, sheet = 1)

# Clean the data by removing NAs and normalizing text
df <- df[!is.na(df$genus), ]
df$genus <- tolower(trimws(df$genus))

# Extract subsets based on direction
positive_genus <- unique(df[df$direction == "positive", "genus"])
negative_genus <- unique(df[df$direction == "negative", "genus"])

# Convert tibbles to character vectors
positive_genus <- as.character(positive_genus$genus)
negative_genus <- as.character(negative_genus$genus)

# Now print and check the structures again
print("Positive Genera (Vector):")
print(positive_genus)
print("Negative Genera (Vector):")
print(negative_genus)

# Calculate the set sizes for the Venn diagram again
length_positive <- length(positive_genus)
length_negative <- length(negative_genus)
length_intersection <- length(intersect(positive_genus, negative_genus))

# Print the lengths to verify
print(paste("Length of Positive Set: ", length_positive))
print(paste("Length of Negative Set: ", length_negative))
print(paste("Length of Intersection: ", length_intersection))

# 打开 TIFF 设备
tiff("venn_diagram_arg.tiff", width = 2500, height = 2500, res = 300)
svg("venn_diagram_arg.svg", width = 8, height = 8)
png("venn_diagram_arg.png", width = 800, height = 800)

# Create Venn diagram
venn_plot <- draw.pairwise.venn(
  area1 = length_positive,
  area2 = length_negative,
  cross.area = length_intersection,
  category = c("Positive", "Negative"),
  fill = c(adjustcolor("#240046", alpha.f = 1), adjustcolor("#ff6d00", alpha.f = 1)), 
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

positive_unique_genus <- setdiff(positive_genus, negative_genus)

grid.text(
  label = paste(sort(positive_unique_genus), collapse = "\n"), 
  x = 0.85, 
  y = 0.5, 
  just = "center", 
  gp = gpar(fontsize = 8, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)

negative_unique_genus <- setdiff(negative_genus, positive_genus)

grid.text(
  label = paste(sort(negative_unique_genus), collapse = "\n"), 
  x = 0.3, 
  y = 0.5, 
  just = "center", 
  gp = gpar(fontsize = 8, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)

shared_genus <- intersect(positive_genus, negative_genus)

grid.text(
  label = paste(sort(shared_genus), collapse = "\n"), 
  x = 0.62,  
  y = 0.5, 
  just = "center", 
  gp = gpar(fontsize = 8, col = "white", fontface = "bold.italic", fontfamily = "Arial")
)


dev.off()

dev.cur()

