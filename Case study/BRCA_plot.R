# Run after GLASD analysis in MATLAB
rm(list=ls())
setwd("U:/GLASD_Robust_corr/Case study")
library(cowplot)
library(ggtext)
library(ggplot2)
library(reshape2)
library(magick)
library(gridExtra)
library(dplyr) 

C_huber <- read.csv("C_huber_optimal.csv", header = FALSE)


### Names of the proteins and pathways #########################################
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA",
                "FOXM1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")
proteins_here <- unique(c(Breast_reactive, Cell_cycle, Hormone_receptor, 
                          Hormone_signaling_Breast))
################################################################################

# Read matrix
C_huber <- read.csv("C_huber_optimal.csv", header = FALSE)
colnames(C_huber) <- proteins_here
rownames(C_huber) <- proteins_here
cor_mat <- as.matrix(C_huber)
cor_melt <- melt(cor_mat)

# Label colors
label_colors <- setNames(rep("black", length(proteins_here)), proteins_here)
label_colors[Breast_reactive] <- "forestgreen"
label_colors[Cell_cycle] <- "firebrick"
label_colors[Hormone_receptor] <- "royalblue"
label_colors[Hormone_signaling_Breast] <- "orange3"

colored_labels <- sapply(names(label_colors), function(prot) {
  sprintf("<span style='color:%s'>%s</span>", label_colors[prot], prot)
})

# Heatmap with built-in correlation color scale
heatmap_plot <- ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1),
                       name = "Correlation") +
  coord_fixed() +
  theme_minimal() +
  #labs(title = "", x = "", y = "") +
  labs(title = "Robust Correlation Heatmap (Huber Loss)", x = "", y = "") +
  scale_x_discrete(labels = colored_labels) +
  scale_y_discrete(labels = colored_labels) +
theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    axis.text.y = element_markdown(),
    legend.position = "right",
    legend.justification = c(0, 1),
    legend.box.just = "left",
    legend.margin = margin(l = -130),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_colorbar(barheight = 8))


# Extract the legend from the heatmap
legend_only <- get_legend(heatmap_plot)

# Create custom pathway legend as vertical boxes + text
legend_df <- data.frame(
  label = c("Breast Reactive", "Cell Cycle", "Hormone Receptor", "Hormone Signaling"),
  color = c("forestgreen", "firebrick", "royalblue", "orange3")
)

pathway_legend <- ggplot(legend_df, aes(x = 1, y = reorder(label, desc(label)), fill = color)) +
  geom_tile(width = 0.2, height = 0.5, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = c(-0.25,-0.38,-0.2,-0.2), size = 4) +
  xlim(0.8, 2.5) +
  scale_fill_identity() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  coord_fixed(ratio = 0.2)

# Combine colorbar + pathway legend vertically
right_stack <- plot_grid(legend_only, pathway_legend +
                           theme(plot.margin = margin(0, 0, 0, -170)),  # <--- shifted closer
                         ncol = 1, rel_heights = c(4, 4))

# Final layout: heatmap left, legends right (tighter)
final_plot <- plot_grid(heatmap_plot + theme(legend.position = "none"),
                        right_stack,
                        ncol = 2,
                        rel_widths = c(1, 0.22),  # <--- narrower gap
                        align = "h")
# Save and show
print(final_plot)
ggsave("C_huber_heatmap.jpg", final_plot,
       width = 12, height = 6, dpi = 300, units = "in")