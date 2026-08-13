cat("===== BUILDING MONOCLE 3 BANJO INPUTS =====\n")

# ------------------------------------------------------------
# 1. Load original RNA data and Monocle ordering
# ------------------------------------------------------------
rna <- read.csv(
    "rna.csv",
    check.names = FALSE
)

rna$Cell <- seq_len(nrow(rna))

ord <- read.csv(
    "data/pseudotime/monocle3_order.csv",
    check.names = FALSE
)

# ------------------------------------------------------------
# 2. Validate Monocle ordering
# ------------------------------------------------------------
stopifnot(nrow(ord) == 959)
stopifnot(!anyDuplicated(ord$Cell))
stopifnot(!591 %in% ord$Cell)
stopifnot(all(is.finite(ord$monocle3_pseudotime)))

idx <- match(ord$Cell, rna$Cell)

stopifnot(!anyNA(idx))

ordered_rna <- rna[idx, , drop = FALSE]

stopifnot(all(ordered_rna$Cell == ord$Cell))

cat("Cells matched successfully: ", nrow(ordered_rna), "\n", sep = "")

# ------------------------------------------------------------
# 3. Read EXACT SAME gene sets used for Slingshot
# ------------------------------------------------------------
top10_file <- "gene_sets/top10_genes.txt"

top15_file <- "gene_sets/top15_genes.txt"

top10 <- scan(
    top10_file,
    what = character(),
    quiet = TRUE
)

top15 <- scan(
    top15_file,
    what = character(),
    quiet = TRUE
)

cat("\nTop 10 genes:\n")
print(top10)

cat("\nTop 15 genes:\n")
print(top15)

# Validate gene sets
stopifnot(length(top10) == 10)
stopifnot(length(top15) == 15)

stopifnot(all(top10 %in% colnames(rna)))
stopifnot(all(top15 %in% colnames(rna)))

stopifnot(all(top10 == top15[1:10]))

# ------------------------------------------------------------
# 4. Create output directory
# ------------------------------------------------------------
dir.create(
    "data/banjo_inputs",
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 5. Create Top-10 and Top-15 matrices
# ------------------------------------------------------------
top10_matrix <- ordered_rna[, top10, drop = FALSE]
top15_matrix <- ordered_rna[, top15, drop = FALSE]

write.table(
    top10_matrix,
    file = "data/banjo_inputs/banjo_monocle3_top10.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

write.table(
    top15_matrix,
    file = "data/banjo_inputs/banjo_monocle3_top15.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# ------------------------------------------------------------
# 6. Save ordering + metadata for reproducibility
# ------------------------------------------------------------
metadata <- data.frame(
    Cell = ord$Cell,
    h = ordered_rna$h,
    monocle3_pseudotime = ord$monocle3_pseudotime
)

write.csv(
    metadata,
    "data/pseudotime/monocle3_cell_order.csv",
    row.names = FALSE
)

writeLines(
    top10,
    "gene_sets/top10_genes.txt"
)

writeLines(
    top15,
    "gene_sets/top15_genes.txt"
)

# ------------------------------------------------------------
# 7. Report
# ------------------------------------------------------------
cat("\n===== OUTPUT CHECK =====\n")

cat(
    "Top-10 observations: ",
    nrow(top10_matrix),
    "\n",
    sep = ""
)

cat(
    "Top-10 genes: ",
    ncol(top10_matrix),
    "\n",
    sep = ""
)

cat(
    "Top-15 observations: ",
    nrow(top15_matrix),
    "\n",
    sep = ""
)

cat(
    "Top-15 genes: ",
    ncol(top15_matrix),
    "\n",
    sep = ""
)

cat("\nFirst 10 ordered cells:\n")
print(
    metadata[1:10, ]
)

cat("\nLast 10 ordered cells:\n")
print(
    tail(metadata, 10)
)

cat("\nFiles created:\n")
cat("monocle3_banjo/banjo_monocle3_top10.txt\n")
cat("monocle3_banjo/banjo_monocle3_top15.txt\n")
cat("monocle3_banjo/monocle3_cell_order.csv\n")
cat("monocle3_banjo/top10_genes.txt\n")
cat("monocle3_banjo/top15_genes.txt\n")

cat("\nMONOCLE 3 BANJO INPUTS COMPLETE\n")
