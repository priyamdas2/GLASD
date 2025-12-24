# =========================================================
#  ORDERED EIGENVALUE SPECTRAL PLOT (FINAL – REVIEWER SAFE)
#  Methods: Naive, Tyler, GLASD-Huber
#  Vertical lines ONLY at FIRST NEGATIVE eigenvalue index
#  NO line for GLASD if all eigenvalues are positive
# =========================================================

rm(list = ls())

setwd("U:/GLASD_Robust_corr/Case study high dim")

library(tidyverse)

# ------------------- INPUT FILES -------------------
files <- c(
  "C_naive.csv",
  "C_tyler.csv",
  "C_huber_optimal.csv"
)

method_names <- c(
  "Naive",
  "Tyler",
  "GLASD-Huber"
)

# ------------------- SAFE MATRIX READER -------------------
safe_read_matrix <- function(f) {
  
  M <- as.matrix(read.csv(f, header = FALSE, check.names = FALSE))
  
  # Drop row-name column if present
  if (!is.numeric(M[1, 1])) {
    M <- M[, -1, drop = FALSE]
  }
  
  # Force numeric
  M <- apply(M, 2, as.numeric)
  
  # Squareness check
  if (nrow(M) != ncol(M)) {
    stop(paste("Matrix in", f, "is NOT square:", nrow(M), "x", ncol(M)))
  }
  
  return(M)
}

# ------------------- READ + EIGEN EXTRACTION -------------------
eig_df <- map2_dfr(files, method_names, function(f, method) {
  
  R <- safe_read_matrix(f)
  
  # Enforce symmetry & unit diagonal
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  # True eigenvalues (NO clipping)
  eigvals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  
  tibble(
    method = method,
    eigval = eigvals
  )
})

# ------------------- ORDER EIGENVALUES -------------------
eig_df <- eig_df %>%
  group_by(method) %>%
  arrange(desc(eigval), .by_group = TRUE) %>%
  mutate(index = row_number()) %>%
  ungroup()

# ------------------- LOG SCALE (POSITIVE ONLY, NO WARNINGS) -------------------
eig_df <- eig_df %>%
  mutate(
    eig_plot = if_else(eigval > 0,
                       log10(eigval),
                       as.numeric(NA))
  )

# ------------------- FIND FIRST NEGATIVE EIGEN INDEX ONLY -------------------
neg_start_df <- eig_df %>%
  group_by(method) %>%
  summarise(
    first_neg_index = ifelse(any(eigval < 0),
                             min(index[eigval < 0]),
                             NA_integer_),
    num_negative = sum(eigval < 0),
    rank_est = sum(eigval > 1e-10),
    .groups = "drop"
  ) %>%
  filter(!is.na(first_neg_index))    # GLASD automatically removed if no negatives

print(neg_start_df)

# ------------------- PROFESSIONAL FRAMED PLOT -------------------
p <- ggplot(eig_df, aes(x = index, y = eig_plot, color = method)) +
  
  geom_line(linewidth = 1.15, na.rm = TRUE) +
  
  # Vertical dashed lines ONLY at FIRST NEGATIVE eigenvalue
  geom_vline(
    data = neg_start_df,
    aes(xintercept = first_neg_index, color = method),
    linetype = "dashed",
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Naive"       = "#4E79A7",
      "Tyler"       = "#F28E2B",
      "GLASD-Huber" = "#E15759"
    )
  ) +
  
  labs(
    x = "Ordered Index",
    y = expression(log[10](Eigenvalue)),
    color = "Method"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    # ---- LEGEND POSITION ----
    legend.position = c(0.02, 0.05),     
    legend.justification = c("left", "bottom"),
    
    # ---- LARGER AXIS LABELS ----
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    
    # ---- LARGER AXIS TICKS (optional but recommended) ----
    axis.text.x  = element_text(size = 13),
    axis.text.y  = element_text(size = 13),
    
    # ---- LARGER LEGEND TEXT ----
    legend.title = element_text(size = 17),#, face = "bold"),
    legend.text  = element_text(size = 13),
    
    plot.title = element_blank(),
    panel.grid.minor = element_blank(),
    
    # ---- WHITE LEGEND BOX ----
    legend.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.4
    ),
    legend.key = element_rect(
      fill = "white",
      colour = NA
    ),
    
    # ---- CLEAN PROFESSIONAL FRAME ----
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.85
    )
  )


# ------------------- SAVE FIGURE -------------------
ggsave(
  filename = "Ordered_Eigenvalue_Spectral_Comparison_Final.jpg",
  plot = p,
  width = 8.5, height = 6, dpi = 600
)

p
cat("\nSaved: Ordered_Eigenvalue_Spectral_Comparison_Final.jpg\n")
