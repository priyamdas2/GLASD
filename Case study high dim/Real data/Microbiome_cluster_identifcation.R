rm(list = ls())

setwd("U:/GLASD_Robust_corr/Case study high dim/Real data")

library(stringr)
library(tibble)
library(dplyr)

X <- read.table("X.csv", sep = ",", header = FALSE)
X_names <- scan("X_names.csv", what = "character", sep = ",", quiet = TRUE)


tax_df <- tibble(raw = X_names) %>%
  mutate(
    Domain  = str_extract(raw, "d__[^\\.]+") %>% str_remove("d__"),
    Phylum  = str_extract(raw, "p__[^\\.]+") %>% str_remove("p__"),
    Class   = str_extract(raw, "c__[^\\.]+") %>% str_remove("c__"),
    Order   = str_extract(raw, "o__[^\\.]+") %>% str_remove("o__"),
    Family  = str_extract(raw, "f__[^\\.]+") %>% str_remove("f__"),
    Genus   = str_extract(raw, "g__[^\\.]+") %>% str_remove("g__")
  )

length(unique(tax_df$Phylum))  # 13
length(unique(tax_df$Class))   # 16
length(unique(tax_df$Order))   # 36
length(unique(tax_df$Family))  # 67
length(unique(tax_df$Genus))   # 257

# Order frequencies (descending)
order_order <- names(sort(table(tax_df$Order), decreasing = TRUE))

# Create ordered factor
tax_df$Order_ordered <- factor(tax_df$Order, levels = order_order)

# Get column permutation
ord <- order(tax_df$Order_ordered)

# Reorder matrix and names
X_blocked        <- X[, ord]
X_names_blocked  <- X_names[ord]

# Verify the new order sequence
tax_df$Order[ord]


# ---- STEP 1: Frequency table ----
ord_tab <- table(tax_df$Order)

# ---- STEP 2: Keep only orders with >= 5 taxa ----
keep_orders <- names(ord_tab[ord_tab >= 5])

# ---- STEP 3: Sort kept orders by frequency (descending) ----
keep_orders_sorted <- names(sort(ord_tab[keep_orders], decreasing = TRUE))

# ---- STEP 4: Collapse others ----
tax_df$Order_block <- ifelse(
  tax_df$Order %in% keep_orders_sorted,
  tax_df$Order,
  "Other"
)

# ---- STEP 5: FORCE FACTOR LEVEL ORDER (Other LAST) ----
tax_df$Order_block <- factor(
  tax_df$Order_block,
  levels = c(keep_orders_sorted, "Other")
)

# ---- STEP 6: Reorder columns ----
ord <- order(tax_df$Order_block)

X_blocked        <- X[, ord]
X_names_blocked  <- X_names[ord]
tax_df_blocked   <- tax_df[ord, ]

# ---- STEP 7: Block sizes IN DISPLAY ORDER ----
block_sizes <- table(tax_df$Order_block)
block_sizes


# ---- STEP 8: Save outputs ----
# ---- (1) Blocked data matrix ----
write.table(
  X_blocked,
  file = "X_blocked_order.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE
)

# ---- (2) Blocked variable (genus) names ----
write.table(
  X_names_blocked,
  file = "X_names_blocked_order.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# ---- (3) Block labels (ORDER LEVEL) ----
Order_block_vec <- tax_df_blocked$Order_block 

write.table(
  Order_block_vec,
  file = "Order_block_vector.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

