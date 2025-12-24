## ------------------- Setup -------------------
# You may need:
# install.packages(c("tidyverse", "patchwork"))
setwd("U:/GLASD_Robust_corr/Benchmark/ALL_Outputs")
## ------------------- Setup -------------------
library(tidyverse)
library(patchwork)

## ===== UPDATED: NOW 4 FUNCTIONS =====
funcNames <- c("ackley", "griewank", "rastrigin", "rosenbrock")
funcNames.display <- c("Ackley", "Griewank", "Rastrigin", "Rosenbrock")

M_values  <- c(10)
num_reps  <- 100

methods <- c(
  "GLASD",
  "Active-set",
  "Interior-point",
  "Sqp",
  "Barzilai-borwein",
  "Conjugate-gradient",
  "Steepest-descent",
  "Trust-region"
)

fun_plots  <- list()
time_plots <- list()
plot_idx <- 1

for (ii in seq_along(funcNames)) {
  f <- funcNames[ii]
  f_disp <- funcNames.display[ii]
  for (M in M_values) {
    
    ## ---------- Read matrices ----------
    funvals_others  <- as.matrix(read.csv(
      sprintf("Funvals_100_reps_others_%s_M_%d_reps_100.csv", f, M),
      header = FALSE))
    
    comp_times_others <- as.matrix(read.csv(
      sprintf("Comp_times_100_reps_others_%s_M_%d_reps_100.csv", f, M),
      header = FALSE))
    
    funvals_GLASD <- as.matrix(read.csv(
      sprintf("Funvals_100_reps_GLASD_%s_M_%d_reps_100.csv", f, M),
      header = FALSE))
    
    comp_times_GLASD <- as.matrix(read.csv(
      sprintf("Comp_times_100_reps_GLASD_%s_M_%d_reps_100.csv", f, M),
      header = FALSE))
    
    funvals    <- cbind(funvals_GLASD, funvals_others)
    comp_times <- cbind(comp_times_GLASD, comp_times_others)
    
    colnames(funvals)    <- methods
    colnames(comp_times) <- methods
    
    ## ---------- Tidy data ----------
    fun_df <- as.data.frame(funvals)
    fun_df$rep <- seq_len(num_reps)
    
    fun_long <- fun_df |>
      pivot_longer(cols = all_of(methods),
                   names_to = "Method",
                   values_to = "FunVal") |>
      mutate(
        Method   = factor(Method, levels = methods),  
        Function = f,
        M        = M
      )
    
    time_df <- as.data.frame(comp_times)
    time_df$rep <- seq_len(num_reps)
    
    time_long <- time_df |>
      pivot_longer(cols = all_of(methods),
                   names_to = "Method",
                   values_to = "CompTime") |>
      mutate(
        Method   = factor(Method, levels = methods),   
        Function = f,
        M        = M
      )
    
    
    base_theme <- theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1, vjust = 1),
        plot.title   = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    ## ---------- Function value boxplot ----------
    p_fun <- ggplot(fun_long, aes(x = Method, y = FunVal)) +
      geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
      scale_y_log10() +
      labs(
        x = "",
        y = expression(log[10]*"(obj. func. value)"),   
        title = sprintf("%s: Objective value", f_disp)
      ) +
      base_theme
    
    ## ---------- Computation time boxplot ----------
    p_time <- ggplot(time_long, aes(x = Method, y = CompTime)) +
      geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
      scale_y_log10() +
      labs(
        x = "",
        y = expression(log[10]*"(comp. time)"),  
        title = sprintf("%s: Computation time", f_disp)
      ) +
      base_theme
    
    fun_plots[[plot_idx]]  <- p_fun
    time_plots[[plot_idx]] <- p_time
    plot_idx <- plot_idx + 1
  }
}

## ==========================================================
## Combine into one 4×2 panel (DYNAMIC & SAFE)
## ==========================================================
combined_plot <- wrap_plots(
  map(seq_along(fun_plots), 
      ~ fun_plots[[.x]] | time_plots[[.x]]),
  ncol = 1
)

## ---------- Save final figure ----------
ggsave(
  filename = "Boxplots_Fun_and_Time_all_functions.pdf",
  plot     = combined_plot,
  width    = 10,
  height   = 12
)

ggsave(
  filename = "Boxplots_Fun_and_Time_all_functions.png",
  plot     = combined_plot,
  width    = 10,
  height   = 12,
  dpi      = 300
)
