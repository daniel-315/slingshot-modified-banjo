# Generic gene subset generator for Banjo
#
# Usage:
# Rscript scripts/gene_selection/make_banjo_subset.R \
#   <ordered_expression.txt> \
#   <gene_list.txt> \
#   <output.txt>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript make_banjo_subset.R ",
      "<ordered_expression.txt> ",
      "<gene_list.txt> ",
      "<output.txt>\n"
    )
  )
}

expression_file <- args[1]
gene_file <- args[2]
output_file <- args[3]

# Read pseudotime-ordered expression matrix
expr <- read.delim(
  expression_file,
  header = TRUE,
  check.names = FALSE
)

# Read one gene name per line
genes <- readLines(gene_file)

# Remove whitespace and blank lines
genes <- trimws(genes)
genes <- genes[genes != ""]

# Remove duplicate gene names while preserving order
genes <- unique(genes)

cat("Input expression file:", expression_file, "\n")
cat("Observations:", nrow(expr), "\n")
cat("Available genes:", ncol(expr), "\n")
cat("Requested genes:", length(genes), "\n\n")

# Check for genes that do not exist
missing_genes <- setdiff(genes, colnames(expr))

if (length(missing_genes) > 0) {
  cat("ERROR: These genes were not found:\n")
  cat(paste(missing_genes, collapse = "\n"), "\n")
  quit(status = 1)
}

# Select requested genes
subset_expr <- expr[, genes, drop = FALSE]

# Write TAB-delimited Banjo input
write.table(
  subset_expr,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Selected genes:\n")
for (i in seq_along(genes)) {
  cat(sprintf("%2d. %s\n", i, genes[i]))
}

cat("\nBANJO INPUT CREATED\n")
cat("Output:", output_file, "\n")
cat("Observations:", nrow(subset_expr), "\n")
cat("Genes:", ncol(subset_expr), "\n")
