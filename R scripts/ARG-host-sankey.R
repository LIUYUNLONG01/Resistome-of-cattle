
install.packages("ggplot2")
install.packages("ggtree")
install.packages("tidyverse")

library(ggplot2)
library(ggtree)
library(tidyverse)
library(ggalluvial)
library(showtext)

showtext_auto()


data <- readxl::read_excel("/Filtered_ARG_Data_by_Genus_8_genes.xlsx")

data_subset <- data %>%
  distinct(Gene, Genus, Order, Family, Class, Phylum)

my_colors <- c(
  "Bacillota" = "#AB8CBB",     
  "Bacteroidota" = "#C8DEDE", 
  "Euryarchaeota" = "#92C6FC",
  "Spirochaetota" = "#FDC563" 
)


ggplot(data = data_subset,
       aes(axis1 = Phylum, axis2 = Class, axis3 = Order, axis4 = Family, axis5 = Genus, axis6 = Gene)) +
  geom_alluvium(aes(fill = Phylum), width = 1/12) +
  geom_stratum(width = 1/16, fill = "grey90", color = "grey60", size = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, hjust = 0, vjust = 0.5, nudge_x = 0.05, family = "Arial") +
  scale_x_discrete(limits = c("Phylum", "Class", "Order", "Family", "Genus", "Gene")) +
  scale_fill_manual(values = my_colors) + 
  labs(title = " ") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks = element_blank(), 
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12), 
    legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)  
  )


ggsave("ARG-host-sankey.pdf", width = 12, height = 20)

ggsave("ARG-host-sankey.svg", width = 12, height = 20)
