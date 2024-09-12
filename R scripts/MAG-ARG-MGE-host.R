
library(readxl)
data <- read_excel("MAG-ARG-MGE.xlsx")

library(dplyr)
library(circlize)

# Create a table of connections
connections <- data %>%
  select(Species = species, gene = gene) %>%
  na.omit() %>%
  count(Species, gene)

# Convert to matrix format for circlize
connection_matrix <- xtabs(n ~ Species + gene, data = connections)

circos.clear()

circos.par(gap.degree = 2, 
           track.margin = c(0.01, 0.01), 
           cell.padding = c(0.02, 0.02, 0.02, 0.02),
           track.height = 0.1)

Species_list <- unique(connections$Species)
gene_list <- unique(connections$gene)

Species_colors <- colorRampPalette(c("#AB8CBB", "#C8DEDE"))(length(Species_list))
gene_colors <- colorRampPalette(c("#FDC563", "#92C6FC"))(length(gene_list))

names(Species_colors) <- Species_list
names(gene_colors) <- gene_list

all_colors <- c(Species_colors, gene_colors)

svg("MAG-ARG-MGE-HOST3.svg", width = 12, height = 12)

tiff("MAG-ARG-MGE-HOST3.tiff", width = 10, height = 10, units = "in", res = 350)

chordDiagram(connection_matrix, 
             grid.col = all_colors, 
             transparency = 0.2, 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.3))

circos.track(track.index = 1, 
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$ylim[1], 
                           CELL_META$sector.index, 
                           facing = "clockwise", 
                           niceFacing = TRUE, 
                           adj = c(0, 0.5),
                           cex = 0.9)
             }, 
             bg.border = NA)

dev.off()

circos.clear()

