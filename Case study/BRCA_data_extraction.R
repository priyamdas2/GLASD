# List of cancers to work with: BRCA
rm(list=ls())
setwd("U:/GLASD_Robust_corr/Case study")
library(cowplot)
library(ggtext)
library(ggplot2)
library(reshape2)
library(magick)
library(gridExtra)
library(dplyr) 
################################################################################
### Reading dataset ############################################################
################################################################################

names_vector <- read.csv("NExUS data/cancer_serial_no.csv", header = FALSE)[[1]]
target_terms <- "BRCA"
match_indices <- which(names_vector %in% target_terms)
matched_names <- names_vector[match_indices]
matrix_list <- lapply(match_indices, function(i) {
  file_path <- file.path("NExUS data", paste0(i, ".csv"))
  as.matrix(read.csv(file_path, header = FALSE))
})
names(matrix_list) <- matched_names


BRCA <- matrix_list$BRCA


################################################################################
### Selecting variables ########################################################
################################################################################

### Full pathway information ###################################################
# load("U:/MOOP_Christine/Github_vignette/Case study/NExUS data/RPPA_12_pathway.rda")
name_list <- vector("list", 12)
array_names_FULL <- c("APOPTOSIS", "CELL CYCLE", "DNA DMG RSPNS", "EMT", 
                      "HORMONE RECPTR", "HORMONE SIG BRST", "PI3K/AKT",
                      "RAS/MAPK", "RTK", "TSC/mTOR", "BREAST REACTIVE", 
                      "CORE REACTIVE") 
Apoptosis <- c("BAK", "BAX", "BID", "BIM", "CASPASE7CLEAVEDD198", "BAD_pS112", 
               "BCL2", "BCLXL", "CIAP")
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA",
                "FOXM1")
DNA_damage_response <- c("X53BP1", "ATM", "CHK1_pS345", "CHK2_pT68", "KU80", 
                         "MRE11", "P53", "RAD50", "RAD51", "XRCC1")
EMT <- c("FIBRONECTIN", "NCADHERIN", "COLLAGENVI", "CLAUDIN7", "ECADHERIN", 
         "BETACATENIN", "PAI1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")
PI3K_AKT <- c("P27_pT157", "P27_pT198", "INPP4B", "AKT_pS473", "AKT_pT308", 
              "GSK3ALPHABETA_pS21S9",
              "GSK3_pS9", "PRAS40_pT246", "TUBERIN_pT1462", "PTEN")
RAS_MAPK <- c("ARAF_pS299", "CJUN_pS73", "CRAF_pS338", "JNK_pT183Y185", 
              "MAPK_pT202Y204", "MEK1_pS217S221", "P38_pT180Y182", 
              "P90RSK_pT359S363", "YB1_pS102")
RTK <- c("EGFR_pY1068", "EGFR_pY1173", "HER2_pY1248", "HER3_pY1289", "SHC_pY317",
         "SRC_pY416", "SRC_pY527")
TSC_mTOR <- c("X4EBP1_pS65", "X4EBP1_pT37T46", "X4EBP1_pT70", "P70S6K_pT389", 
              "MTOR_pS2448", "S6_pS235S236", "S6_pS240S244", "RB_pS807S811")
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")
Core_reactive <- c("CLAUDIN7", "ECADHERIN", "BETACATENIN", "CAVEOLIN1", "RBM15")
################################################################################

proteins_here <- unique(c(Breast_reactive, Cell_cycle, Hormone_receptor, Hormone_signaling_Breast))
selected_variables <- read.csv("NExUS data/selected_variables.csv", header = FALSE)[[1]]
match_indices <- match(proteins_here,selected_variables)


BRCA_subset <- BRCA[, match_indices]
p <- dim(BRCA_subset)[2]
sample_sizes <- dim(BRCA_subset)[1]

write.csv(BRCA_subset, file = "BRCA_subset.csv", row.names = FALSE)



################################################################################
# Define protein groups
Breast_reactive <- c("BETACATENIN", "CAVEOLIN1", "MYH11", "RAB11", "GAPDH", "RBM15")
Cell_cycle <- c("CDK1", "CYCLINB1", "CYCLINE2", "P27_pT157", "P27_pT198", "PCNA", "FOXM1")
Hormone_receptor <- c("ERALPHA", "ERALPHA_pS118", "PR", "AR")
Hormone_signaling_Breast <- c("BCL2", "INPP4B", "GATA3")

# Assign group labels to proteins
protein_group <- sapply(proteins_here, function(p) {
  if (p %in% Breast_reactive) return("Breast Reactive")
  if (p %in% Cell_cycle) return("Cell Cycle")
  if (p %in% Hormone_receptor) return("Hormone Receptor")
  if (p %in% Hormone_signaling_Breast) return("Hormone Signaling")
  return("Unknown")
})

# Compute outlier counts per protein (IQR method)
iqr_counts <- apply(BRCA_subset, 2, function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  sum(x < (q1 - 1.5 * iqr) | x > (q3 + 1.5 * iqr))
})

# Create data frame
bar_df <- data.frame(
  Protein = proteins_here,
  OutlierCount = iqr_counts,
  Pathway = factor(protein_group, levels = c("Breast Reactive", "Cell Cycle", "Hormone Receptor", "Hormone Signaling"))
)

# Use custom colors
custom_colors <- c("Breast Reactive" = "forestgreen", 
                   "Cell Cycle" = "firebrick", 
                   "Hormone Receptor" = "royalblue", 
                   "Hormone Signaling" = "orange3")

# Plot
library(ggplot2)

outlier_plot <- ggplot(bar_df, aes(x = reorder(Protein, OutlierCount), y = OutlierCount, fill = Pathway)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "IQR-Based Outlier Counts per Protein",
    x = "Protein", y = "Outlier Count", fill = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(outlier_plot)

ggsave("Outlier_plot.jpg", outlier_plot,
       width = 10, height = 6, dpi = 300, units = "in")