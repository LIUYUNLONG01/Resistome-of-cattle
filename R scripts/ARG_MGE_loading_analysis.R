#ARG
data <- read.csv("/permanova_full_data.csv")

columns_to_scale <- 1:(ncol(data) - 4) 

data_scaled <- data 
data_scaled[, columns_to_scale] <- scale(data[, columns_to_scale])

library(stats)
library(ggplot2) 
library(ggrepel)


pca_result <- prcomp(data_scaled[, columns_to_scale], scale. = TRUE) 

loadings <- pca_result$rotation

print(loadings[, 1:2])

top_loadings_pc1 <- sort(abs(loadings[, "PC1"]), decreasing = TRUE)[1:5]
top_loadings_pc2 <- sort(abs(loadings[, "PC2"]), decreasing = TRUE)[1:5]

cat("Top Loadings for PC1:\n")
print(top_loadings_pc1)

cat("Top Loadings for PC2:\n")
print(top_loadings_pc2)


loadings_df <- as.data.frame(loadings[, 1:2])
names(loadings_df) <- c("PC1", "PC2")


loadings_df$Magnitude <- sqrt(loadings_df$PC1^2 + loadings_df$PC2^2)


ggplot(loadings_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = Magnitude, color = Magnitude), alpha = 1.0) + 
  geom_text_repel(aes(label = rownames(loadings_df)), size = 4, box.padding = 0.35, point.padding = 0.5) + 
  scale_color_gradient(low = "#C8DEDE", high = "#AB8CBB") + 
  scale_size(range = c(1, 10)) +  
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") + 
  theme_minimal() +
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "black"),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines")) +  
  labs(title = "",
       subtitle = " ",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave("ARG_loadings_plot.svg", width = 10, height = 8)



#MGE

data <- read.csv("/merged_data_with_quartiles_and_parity.csv")

columns_to_scale <- 2:(ncol(data) - 4)  


data_scaled <- data 
data_scaled[, columns_to_scale] <- scale(data[, columns_to_scale])


library(stats)
library(ggplot2) 
library(ggrepel)

pca_result <- prcomp(data_scaled[, columns_to_scale], scale. = TRUE) 

loadings <- pca_result$rotation


print(loadings[, 1:2])


top_loadings_pc1 <- sort(abs(loadings[, "PC1"]), decreasing = TRUE)[1:5]
top_loadings_pc2 <- sort(abs(loadings[, "PC2"]), decreasing = TRUE)[1:5]


cat("Top Loadings for PC1:\n")
print(top_loadings_pc1)

cat("Top Loadings for PC2:\n")
print(top_loadings_pc2)


loadings_df <- as.data.frame(loadings[, 1:2]) 
names(loadings_df) <- c("PC1", "PC2")

loadings_df$Magnitude <- sqrt(loadings_df$PC1^2 + loadings_df$PC2^2)

ggplot(loadings_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(size = Magnitude, color = Magnitude), alpha = 1.0) + 
  geom_text_repel(aes(label = rownames(loadings_df)), size = 4, box.padding = 0.35, point.padding = 0.5) +  
  scale_color_gradient(low = "#C8DEDE", high = "#AB8CBB") + 
  scale_size(range = c(1, 10)) + 
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +  
  theme_minimal() +
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "black"),   
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "lines")) +   
  labs(title = " ",
       subtitle = " ",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave("MGE_loadings_plot.svg", width = 10, height = 8)
