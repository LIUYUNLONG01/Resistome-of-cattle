
arg_data <- read.table("ppm_modified.type.txt", header = TRUE, sep="\t")
mge_data <- read.table("corrected_format_abundance.txt", header = TRUE, sep="\t")

merged_data <- merge(arg_data, mge_data, by = "readsID")

library(dplyr)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(viridis)  
library(viridisLite)

merged_data <- merge(arg_data, mge_data, by = "readsID")

results <- data.frame(ARG = character(), MGE = character(), Correlation = numeric(), P_Value = numeric())

for (arg_col in colnames(arg_data)[-1]) { 
  for (mge_col in colnames(mge_data)[-1]) { 
    test_result <- cor.test(merged_data[[arg_col]], merged_data[[mge_col]], method = "pearson")
    
    results <- rbind(results, data.frame(ARG = arg_col, MGE = mge_col, Correlation = test_result$estimate, P_Value = test_result$p.value))
  }
}

results$Significance <- ifelse(results$P_Value <= 0.001, "***",
                               ifelse(results$P_Value < 0.01, "**",
                                      ifelse(results$P_Value < 0.05, "*", "")))

cor_matrix <- dcast(results, ARG ~ MGE, value.var = "Correlation")

sig_matrix <- dcast(results, ARG ~ MGE, value.var = "Significance")

cor_long <- melt(cor_matrix, id.vars = "ARG", variable.name = "MGE", value.name = "Correlation")
sig_long <- melt(sig_matrix, id.vars = "ARG", variable.name = "MGE", value.name = "Significance")

cor_sig_long <- merge(cor_long, sig_long, by = c("ARG", "MGE"))


ggplot(cor_sig_long, aes(x = MGE, y = ARG, fill = Correlation)) +
  geom_tile(color = "gray") +
  scale_fill_gradient2(low = "#36C25F", high = "#731681", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="Correlation") +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = "black", size = 8, face = "bold"),
    axis.text.y = element_text(color = "black", size = 8, face = "bold")
  ) +
  labs(x = " ", y = " ", title = " ")

ggsave("heatmap_ARG_MGE.pdf", width = 10, height = 8, bg = "white")
ggsave("heatmap_ARG_MGE.svg", width = 10, height = 8, bg = "white")
ggsave("heatmap_ARG_MGE.tiff", width = 10, height = 8, bg = "white")

write.csv(results,"ARG-MGE-relationship.csv", row.names = FALSE)
