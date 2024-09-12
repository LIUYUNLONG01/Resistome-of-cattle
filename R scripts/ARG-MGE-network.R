arg_data <- read.table("ppm_modified.type.txt", header = TRUE, sep="\t")
mge_data <- read.table("corrected_format_abundance.txt", header = TRUE, sep="\t")

merged_data <- merge(arg_data, mge_data, by = "readsID")

library(dplyr)
library(ggplot2)
library(Hmisc)
library(reshape2)

merged_data <- merge(arg_data, mge_data, by = "readsID")

results <- data.frame(ARG = character(), MGE = character(), Correlation = numeric(), P_Value = numeric())

for (arg_col in colnames(arg_data)[-1]) {  
  for (mge_col in colnames(mge_data)[-1]) {  
    test_result <- cor.test(merged_data[[arg_col]], merged_data[[mge_col]], method = "pearson")
    
    results <- rbind(results, data.frame(ARG = arg_col, MGE = mge_col, Correlation = test_result$estimate, P_Value = test_result$p.value))
  }
}

library(igraph)

significant_edges <- results[results$P_Value < 0.05, ]

write.csv(significant_edges,"significant_edges.csv")

g <- graph_from_data_frame(d = significant_edges[, c("ARG", "MGE")], directed = FALSE)

degree_centrality <- degree(g)
betweenness_centrality <- betweenness(g)

fc <- cluster_fast_greedy(g)
communities <- membership(fc)

V(g)$size <- degree_centrality * 2  
V(g)$color <- rainbow(max(communities))[communities] 
E(g)$width <- 1  


par(mfrow = c(1, 1))
svg("Network with Degree Centrality and Visualization.svg", width = 10, height = 7)
plot(g, vertex.size = degree_centrality, layout = layout_with_fr(g), main = "Network with Degree Centrality", vertex.frame.color = NA, 
     vertex.label.cex = 0.6,  
     vertex.label.color = "black",  
     vertex.label.dist = 1.0,  
     vertex.label.font = 1)

plot(fc, g, layout = layout_with_fr(g), main = "Community Structure",vertex.size = V(g)$size, vertex.color = V(g)$color, edge.width = E(g)$width, main = "Network Visualization with Significant Correlations", vertex.frame.color = NA, 
     vertex.label.cex = 0.8,  
     vertex.label.color = "black",
     vertex.label.dist = 1.0, 
     vertex.label.font = 1) 
dev.off()
svg("Network Visualization.svg", width = 10, height = 7)

unique_nodes <- unique(c(significant_edges$ARG, significant_edges$MGE))

nodes <- data.frame(id = unique_nodes)

write.csv(nodes, "nodes.csv", row.names = FALSE, quote = FALSE)
write.csv(significant_edges, "edges_for_cytoscape.csv", row.names = FALSE, quote = FALSE)

