
install.packages("DALEX")
install.packages("caret")
install.packages("ggplot2")
install.packages("DALEX")
install.packages("pdp")
library(randomForest) 
library(DALEX) 
library(caret)
library(lattice) 
library(ggplot2) 
library(pdp) 
library(patchwork) 
library(magrittr)
library(extrafont)
font_import()
loadfonts(device = "pdf")


data <- read.csv('mge_methane.csv')

data <- data[, !names(data) %in% c("Sample")] 

data[is.na(data)] <- 0

X <- data[, -which(names(data) == "Methane.intensity")] 
y <- data$Methane.intensity  

set.seed(42) 
trainIndex <- sample(1:nrow(data), 0.8 * nrow(data)) 
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ] 
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

train_control <- trainControl(method = "cv", number = 5)
model_cv <- train(Methane.intensity ~ ., data = data, method = "ranger",
                  importance = 'impurity', trControl = train_control, 
                  tuneLength = 1, num.trees = 500, seed = 42)
print(model_cv)

best_model <- model_cv$finalModel

y_pred <- predict(best_model, X_test)$predictions 

str(y_test)
str(y_pred)

mse <- mean((y_test - y_pred)^2) 
print(paste("Mean Squared Error:", mse))
r2 <- cor(y_test, y_pred)^2 
print(paste("R^2 Score:", r2))

pred_obs <- data.frame(predicted = y_pred, observed = y_test)

ggplot(data = pred_obs, aes(x = predicted, y = observed)) +
  geom_point(size = 8, color = "#440154FF") +
  geom_smooth(method = "lm", se = FALSE, color = "#FDE725FF", linetype = "dashed") + 
  xlab("Predicted Methane Intensity") +
  ylab("Observed Methane Intensity") +
  ggtitle("Predicted vs Observed Methane Intensity") +
  theme_minimal(base_size = 15) + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"), 
    axis.text = element_text(size = 12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotate("text", x = max(pred_obs$predicted), y = min(pred_obs$observed), 
           label = paste("MSE:", round(mse, 2), "\nR^2:", round(r2, 2)), 
           hjust = 1, vjust = 0, size = 5, color = "darkred")


ggsave(filename = "Predicted_vs_Observed_Methane_Intensity_MGE.tiff", width = 8, height = 8, dpi = 300)

ggsave(filename = "Predicted_vs_Observed_Methane_Intensity_MGE.pdf", width = 8, height = 8)

print(best_model)

importance_values <- best_model$variable.importance
importance_df <- data.frame(Feature = names(importance_values), 
                            Importance = as.numeric(importance_values))
top_20_features <- importance_df[order(-importance_df$Importance), ][1:20, ] 
print(top_20_features)

ggplot(top_20_features, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = 'identity', width = 0.7, fill = "#C8DEDE") + 
  coord_flip() +
  xlab("Feature") +
  ylab("Importance") +
  ggtitle("Top 20 Important Features") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", family = "Arial", color = "black"), 
    axis.title.x = element_text(size = 14, face = "bold", family = "Arial", color = "black"), 
    axis.title.y = element_text(size = 14, face = "bold", family = "Arial", color = "black"), 
    axis.text = element_text(size = 12, family = "Arial", color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  geom_text(aes(label = round(Importance, 2)), hjust = -0.2, size = 4, color = "black") 

ggsave(filename = "Top_20_Important_Features_MGE.tiff", width = 20, height = 8, dpi = 300)
ggsave(filename = "Top_20_Important_Features_MGE.svg", width = 20, height = 8)

graphics.off()

top_features <- top_20_features$Feature[1:20]
pd_plots <- list()
for (feature in top_features) {
  pd_plot <- partial(best_model, pred.var = feature, train = X_train, rug = TRUE) %>% autoplot() + 
    geom_hline(yintercept = mean(y_train), linetype = 2, color = "gray") + 
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.background = element_blank())
  pd_plots <- append(pd_plots, list(pd_plot))
  print(paste0("Partial dependence of ", feature))
}

wrap_plots(pd_plots)

ggsave(filename = "Partial_Dependence_Plots_MGE.tiff", width = 20, height = 14, dpi = 300)


ggsave(filename = "Partial_Dependence_Plots_MGE.pdf", width = 20, height = 14, device = cairo_pdf)

pd_plot <- partial(best_model, pred.var = "X1780_tnpA_EF682209.1", train = X_train, rug = TRUE) %>% 
  autoplot() + 
  geom_hline(yintercept = mean(y_train), linetype = 2, color = "gray") + 
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank())
print(pd_plot)


