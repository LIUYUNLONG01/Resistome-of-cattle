
library(ggplot2)
library(readxl)
library(ggrepel)
library(showtext)

showtext_auto()

data <- read_excel("MGE subtypes.xlsx", skip = 1)

data$Detection_count <- as.numeric(data$`Detected count`)
data$Mean_abundance_ppm <- as.numeric(data$`Average abundance (ppm)`)
data$Detection_frequency <- as.numeric(data$`Detection frequency`)

data$Log10_Mean_abundance_ppm <- log10(data$Mean_abundance_ppm + 1e-9)

data$Frequency <- cut(data$Detection_frequency,
                            breaks = c(-Inf, 0.1, 0.5, 0.9, Inf),
                            labels = c("<=10%", ">10%-<=50%", ">50%-<90%", ">=90%"))

top10 <- data[order(-data$Detection_frequency), ][1:10, ]

ggplot(data, aes(x = Detection_count, y = Log10_Mean_abundance_ppm, fill = Frequency)) +
  geom_point(size = 5, shape = 21, color = "gray", stroke = 1.0) +
  labs(x = "Detection count", y = "Log10 Mean abundance ppm", title = " ") +
  scale_fill_manual(values = c("#FDC563", "#92C6FC", "#C8DEDE", "#AB8CBB")) + 
  geom_text_repel(data = top10, aes(label = `MGE subtypes`), size = 3.5, color = "black", 
                                   box.padding = 0.5, point.padding = 0.5, max.overlaps = 20) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    text = element_text(family = "Arial", size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 20, face = "bold", color = "darkred", hjust = 0.5)
  )


ggsave("MGE_plot.pdf", width = 12, height = 10)

ggsave("MGE_plot.svg", width = 12, height = 10)
