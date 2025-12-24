################################################################################
# CRC MICROBIOME CORRELATION HEATMAP (ORDER-BLOCKED)
# Standalone plotting script
################################################################################

rm(list = ls())

setwd("U:/GLASD_Robust_corr/Case study high dim")

################################################################################
# 1. READ INPUT FILES
################################################################################

# (A) Read correlation matrix (from GLASD / sample / Tyler / etc.)
R_hat <- as.matrix(
  read.table("C_huber_optimal.csv", sep = ",", header = FALSE)
)

# Safety check: square
stopifnot(nrow(R_hat) == ncol(R_hat))

# Enforce exact diagonal = 1
diag(R_hat) <- 1

# (B) Read order-block vector (length should be = number of columns)
Order_vec <- scan("Real data/Order_block_vector.csv",
                  what = "character",
                  sep  = ",")

stopifnot(length(Order_vec) == ncol(R_hat))

################################################################################
# 2. Block sizes & label locations
################################################################################

library(dplyr)

# Frequencies by order-block
block_sizes <- table(Order_vec)

# Ensure "Other" is last (if present)
if ("Other" %in% names(block_sizes)) {
  block_sizes <- c(
    block_sizes[names(block_sizes) != "Other"],
    block_sizes["Other"]
  )
}

block_edges   <- cumsum(block_sizes)
block_centers <- block_edges - block_sizes / 2
block_labels  <- names(block_sizes)

################################################################################
# 3. Convert correlation matrix to long format
################################################################################

R_long <- expand.grid(
  row = seq_len(nrow(R_hat)),
  col = seq_len(ncol(R_hat))
)

R_long$corr <- as.vector(R_hat)

################################################################################
# 4. Nature-style heatmap
################################################################################

library(ggplot2)

p <- ggplot(R_long, aes(x = col, y = row, fill = corr)) +
  geom_raster() +
  
  scale_fill_gradient2(
    low   = "#3B4CC0",
    mid   = "white",
    high  = "#B40426",
    limits = c(-1, 1),
    oob = scales::squish,   # <- CRUCIAL LINE
    name   = "Correlation"
  ) +
  
  scale_x_continuous(
    breaks = block_centers,
    labels = block_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = block_centers,
    labels = block_labels,
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(face = "bold",
                                angle = 90,        # vertical labels
                                hjust = 0.5,
                                vjust = 0.5,
                                size  = 8
    ),
    axis.text.y  = element_text(face = "bold", size = 8),
    
    #axis.title.x = element_text(face = "bold", size = 10),
    #axis.title.y = element_text(face = "bold", size = 10),
    
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )


################################################################################
# 5. Add block boundaries
################################################################################

p <- p +
  geom_vline(
    xintercept = block_edges + 0.5,
    color = "black",
    linewidth = 0.25
  ) +
  geom_hline(
    yintercept = block_edges + 0.5,
    color = "black",
    linewidth = 0.25
  )

print(p)

################################################################################
# 6. Save publication-quality figures
################################################################################

ggsave(
  filename = "CRC_Correlation_OrderLevel.pdf",
  plot     = p,
  width    = 9,
  height   = 9,
  dpi      = 600,
  useDingbats = FALSE
)

ggsave(
  filename = "CRC_Correlation_OrderLevel.png",
  plot     = p,
  width    = 9,
  height   = 9,
  dpi      = 600
)
