# Load necessary libraries
library(ggplot2)
library(dplyr)

# Step 1: Load your data from a CSV file
data <- read.csv("C:/Users/tomgo/Desktop/mapping 17 12/PCA/PCA_grezzi.csv", header = TRUE)

# Step 2: Prepare the data for PCA
data_pca <- data %>% select(-animale)  # Remove the 'animal' column for PCA

# Standardize the data
data_scaled <- scale(data_pca)

# Perform PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Create a data frame with PCA results
pca_df <- as.data.frame(pca_result$x)
pca_df$animal <- data$animal

# Step 3: Define colors based on the group
pca_df$group <- ifelse(grepl("^SS", pca_df$animal), "SS",
                       ifelse(grepl("^MS", pca_df$animal), "MS",
                              ifelse(grepl("^SC", pca_df$animal), "SC",
                                     ifelse(grepl("^MC", pca_df$animal), "MC",
                                            ifelse(grepl("^HC", pca_df$animal), "HC", "Other")))))

# Define colors for each group
group_colors <- c("SS" = "darkblue", "MS" = "darkred", "SC" = "darkorange", 
                  "MC" = "darkgreen", "HC" = "gray32")

# Reorder the group factor levels
pca_df$group <- factor(pca_df$group, levels = c("SS", "MS", "SC", "MC", "HC"))

# Calculate variance explained
explained_variance <- summary(pca_result)$importance[2, ] * 100  # Proportion of variance explained

# Create the plot
p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = animal, color = group)) +
  geom_point(size = 7) +  # Increase the size of the dots
  geom_text(vjust = -1, size = 3) +
  stat_ellipse(aes(group = group, color = group), level = 0.95, alpha = 0.2, size = 3.5) +  # Set ellipse border width
  scale_color_manual(values = group_colors) +  # Use the defined colors
  labs(title = "PCA of Brain Regions",
       x = paste("PC1 (", round(explained_variance[1], 2), "% Variance)", sep = ""), 
       y = paste("PC2 (", round(explained_variance[2], 2), "% Variance)", sep = "")) +
  theme_minimal() +
  theme(
    text = element_text(family = "Century Gothic"),   # Set font family for all text elements
    plot.title = element_text(size = 16, face = "bold"),             # Increase title size and make it bold
    axis.title.x = element_text(size = 14, face = "bold"),           # Increase x-axis label size and make it bold
    axis.title.y = element_text(size = 14, face = "bold"),           # Increase y-axis label size and make it bold
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),   # Increase x-axis tick label size and make it bold
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),   # Increase y-axis tick label size and make it bold
    axis.line.x = element_line(size=1, color="black"),                # Make x-axis line bold
    axis.line.y = element_line(size=1, color="black"),                # Make y-axis line bold
    legend.title = element_blank(),                                   # Remove legend title
    legend.text = element_text(size = 14),                           # Increase legend text size
    legend.key.size = unit(1.5, 'cm'),                               # Increase size of legend keys
    legend.background = element_blank(),                               # Remove legend background rectangle
    legend.key = element_blank()                                       # Remove rectangles behind legend points
  )

# Print the plot
print(p)
