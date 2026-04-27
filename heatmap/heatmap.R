# Carica i pacchetti necessari
library(pheatmap)
library(grid)

# Definisci la palette di colori personalizzata
num_colors <- 400
custom_colors <- colorRampPalette(c("gray0", "gray46", "lightgoldenrod2", "yellow2", "gold", "darkgoldenrod1"))(200)

# Calcola il colore di mezzo (per zero e NA)
midpoint_index <- ceiling(num_colors / 2)
midpoint_color <- custom_colors[midpoint_index]
na_color <- midpoint_color

# Genera i break per la legenda (opzionale, puoi usarli se vuoi)
# breaks <- seq(min(heatmap, na.rm = TRUE), max(heatmap, na.rm = TRUE), length.out = num_colors + 1)

# Supponendo che 'heatmap' sia la tua matrice di dati numerici
# Esempio di matrice (sostituisci con i tuoi dati)
# heatmap <- matrix(rnorm(100), nrow = 10)
# rownames(heatmap) <- paste0("Gene", 1:10)
# colnames(heatmap) <- paste0("Sample", 1:10)

# Crea la heatmap con le impostazioni richieste
pheatmap(heatmap,
         cellwidth = NA,
         cellheight = NA,
         fontsize = 10,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = custom_colors,
         na_col = na_color,
         valueslegend = TRUE,
         legend_title = "Legend Title",
         legend_labels = NULL,
         legend_size = 1.5,
         legend_position = "topright",
         legend_x = -2,
         legend_y = 10.5,
         main = "Heatmap Title",
         fontsize_row = 6,                # font più piccolo per asse Y
         fontsize_col = 12,
         angle_col = 45,
         fontface = "bold",               # testo in grassetto
         fontfamily = "Century Gothic")   # font Century Gothic
