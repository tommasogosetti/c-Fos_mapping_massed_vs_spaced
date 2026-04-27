library(ggplot2)
library(dplyr)

# I tuoi dati
fos <- c(#fos_z_score_x_each_region_of_group_y)
r2 <- c(#r2_score_x_each_region_of_group_y)
# Combine data in one single dataframe
my_data <- data.frame(r2 = r2, fos = fos)

# Calculate corelation and  p-value
correlation <- cor(my_data$fos, my_data$r2, method = "spearman")
correlation_test <- cor.test(my_data$fos, my_data$r2, method = "spearman")
p_value <- correlation_test$p.value

# Define the desired color RGB(252,124,8)
blue_border <- rgb(210, 48, 0, maxColorValue = 255)

# Create graph
plot <- ggplot(my_data, aes(x = fos, y = r2)) +
  geom_point(size = 5, shape = 21, fill = NA, color = blue_border, stroke = 1.5) + # contorno arancione
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Scatterplot of r² vs Fos Signal",
       x = "Magnitude of Fos Signal",
       y = "Mean Correlation Strength (r²)",
       caption = "Data from Brain Regions") +
  theme_minimal(base_size = 18) +
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(family = "Century Gothic", face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_line(size = 1, color = "black")
  ) +
  annotate("text",
           x = max(my_data$fos) * 0.8,
           y = max(my_data$r2) * 0.9,
           label = paste0("r = ", round(correlation, 3), "\np-value = ", signif(p_value, 3)),
           hjust = 0, vjust = 1,
           color = "black",
           size = 7,
           fontface = "bold",
           family = "Century Gothic")

# Print graph
print(plot)