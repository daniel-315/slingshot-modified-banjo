suppressPackageStartupMessages({
    library(monocle3)
    library(SingleCellExperiment)
    library(ggplot2)
})


# ------------------------------------------------------------
# Monocle 3 v1.3.1 partition compatibility fix
#
# v1.3.1 can generate NaN values in compute_partitions()
# when there are zero edges between detected clusters.
# Current upstream Monocle 3 explicitly converts those
# NaN values to zero before constructing the partition graph.
#
# This runtime patch reproduces that upstream fix while
# keeping Monocle 3 v1.3.1 and the rest of the analysis intact.
# ------------------------------------------------------------

compute_partitions_fixed <- function(g,
                                     optim_res,
                                     qval_thresh = 0.05,
                                     verbose = FALSE) {

    cell_membership <- as.factor(igraph::membership(optim_res))

    membership_matrix <-
        Matrix::sparse.model.matrix(~ cell_membership + 0)

    num_links <-
        Matrix::t(membership_matrix) %*%
        igraph::as_adjacency_matrix(g) %*%
        membership_matrix

    diag(num_links) <- 0

    louvain_modules <- levels(cell_membership)

    edges_per_module <- Matrix::rowSums(num_links)
    total_edges <- sum(num_links)

    if (verbose) {
        message("Between-cluster total edges: ", total_edges)
    }

    theta <-
        (as.matrix(edges_per_module) / total_edges) %*%
        Matrix::t(edges_per_module / total_edges)

    var_null_num_links <-
        theta * (1 - theta) / total_edges

    num_links_ij <-
        num_links / total_edges - theta

    cluster_mat <-
        monocle3:::pnorm_over_mat(
            as.matrix(num_links_ij),
            var_null_num_links
        )

    num_links <- num_links_ij / total_edges

    # Upstream Monocle 3 fix for zero total edges
    num_links[is.nan(num_links)] <- 0
    cluster_mat[is.nan(cluster_mat)] <- 0

    cluster_mat <-
        matrix(
            stats::p.adjust(cluster_mat),
            nrow = length(louvain_modules),
            ncol = length(louvain_modules)
        )

    sig_links <- as.matrix(num_links)

    row.names(sig_links) <-
        colnames(sig_links) <-
        louvain_modules

    sig_links[cluster_mat > qval_thresh] <- 0
    diag(sig_links) <- 0

    cluster_g <-
        igraph::graph_from_adjacency_matrix(
            sig_links,
            weighted = TRUE,
            mode = "undirected"
        )

    list(
        cluster_g = cluster_g,
        num_links = num_links,
        cluster_mat = cluster_mat
    )
}

assignInNamespace(
    "compute_partitions",
    compute_partitions_fixed,
    ns = "monocle3"
)

cat("Applied Monocle 3 zero-between-cluster-edge fix.\n")


set.seed(2016)

cat("===== MONOCLE 3 TRAJECTORY =====\n")

dir.create(
    "data/pseudotime",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "results/monocle3/trajectory",
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Load RNA data
# ------------------------------------------------------------
rna <- read.csv("rna.csv", check.names = FALSE)

rna$Cell <- seq_len(nrow(rna))

# Match primary Slingshot analysis
rna <- rna[rna$Cell != 591, , drop = FALSE]

gene_names <- setdiff(colnames(rna), c("h", "Cell"))

cat("Cells used: ", nrow(rna), "\n", sep = "")
cat("Genes used: ", length(gene_names), "\n", sep = "")

# Monocle expects genes x cells
expr <- t(as.matrix(rna[, gene_names, drop = FALSE]))

rownames(expr) <- gene_names
colnames(expr) <- paste0("Cell_", rna$Cell)

# ------------------------------------------------------------
# 2. Metadata
# ------------------------------------------------------------
cell_metadata <- data.frame(
    h = rna$h,
    original_cell = rna$Cell,
    row.names = colnames(expr)
)

gene_metadata <- data.frame(
    gene_short_name = gene_names,
    row.names = gene_names
)

# ------------------------------------------------------------
# 3. Create Monocle 3 CDS
# ------------------------------------------------------------
cds <- new_cell_data_set(
    expression_data = expr,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
)

cat("CDS created successfully.\n")

# ------------------------------------------------------------
# 4. PCA preprocessing
#
# rna.csv contains continuous normalized values rather than
# raw integer counts, so do NOT normalize a second time.
#
# scaling=TRUE parallels the scaled PCA used in Slingshot.
# ------------------------------------------------------------
cds <- preprocess_cds(
    cds,
    method = "PCA",
    num_dim = 10,
    norm_method = "none",
    scaling = TRUE
)

cat("PCA preprocessing complete.\n")

# ------------------------------------------------------------
# 5. UMAP
#
# Keep deterministic settings.
# ------------------------------------------------------------
cds <- reduce_dimension(
    cds,
    reduction_method = "UMAP",
    preprocess_method = "PCA",
    umap.fast_sgd = FALSE,
    cores = 1
)

cat("UMAP complete.\n")

# ------------------------------------------------------------
# 6. Clustering
# ------------------------------------------------------------
cds <- cluster_cells(
    cds,
    reduction_method = "UMAP"
)

cat("\nClusters:\n")
print(table(clusters(cds)))

cat("\nPartitions before graph learning:\n")
print(table(partitions(cds)))

# ------------------------------------------------------------
# 7. Learn ONE principal graph across all cells
#
# This gives us one pseudotime coordinate for the same 959 cells,
# which is appropriate for comparison with the single Slingshot
# ordering used for Banjo.
# ------------------------------------------------------------
cds <- learn_graph(
    cds,
    use_partition = FALSE
)

cat("\nPrincipal graph learned.\n")

# ------------------------------------------------------------
# 8. Select trajectory root programmatically from h = 0 cells
#
# Find the principal graph node with the largest number of
# earliest experimental cells.
# ------------------------------------------------------------
early_cells <- which(colData(cds)$h == 0)

closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

closest_vertex <-
    as.matrix(closest_vertex[colnames(cds), , drop = FALSE])

root_index <- as.numeric(
    names(
        which.max(
            table(closest_vertex[early_cells, ])
        )
    )
)

root_pr_node <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[root_index]

cat("Selected root principal node: ", root_pr_node, "\n", sep = "")
cat("Number of h=0 cells: ", length(early_cells), "\n", sep = "")

# ------------------------------------------------------------
# 9. Order cells in pseudotime
# ------------------------------------------------------------
cds <- order_cells(
    cds,
    root_pr_nodes = root_pr_node
)

pt <- pseudotime(cds)

cat("\n===== PSEUDOTIME CHECK =====\n")
cat("Total cells: ", length(pt), "\n", sep = "")
cat("Finite pseudotime: ", sum(is.finite(pt)), "\n", sep = "")
cat("Infinite pseudotime: ", sum(is.infinite(pt)), "\n", sep = "")
cat("NA pseudotime: ", sum(is.na(pt)), "\n", sep = "")

finite <- is.finite(pt)

rho <- cor(
    pt[finite],
    colData(cds)$h[finite],
    method = "spearman"
)

cat("Spearman pseudotime vs experimental h: ",
    round(rho, 4),
    "\n",
    sep = "")

cat("Pseudotime minimum: ", min(pt[finite]), "\n", sep = "")
cat("Pseudotime maximum: ", max(pt[finite]), "\n", sep = "")

# ------------------------------------------------------------
# 10. Export cell ordering
# ------------------------------------------------------------
umap <- reducedDims(cds)$UMAP

results <- data.frame(
    Cell = colData(cds)$original_cell,
    h = colData(cds)$h,
    monocle3_pseudotime = pt,
    UMAP1 = umap[, 1],
    UMAP2 = umap[, 2],
    cluster = as.character(clusters(cds)),
    partition = as.character(partitions(cds))
)

results <- results[
    order(results$monocle3_pseudotime, results$Cell),
]

write.csv(
    results,
    "data/pseudotime/monocle3_order.csv",
    row.names = FALSE
)

saveRDS(
    cds,
    "data/pseudotime/monocle3_model.rds"
)

# ------------------------------------------------------------
# 11. Save trajectory figures
# ------------------------------------------------------------
p_pt <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    label_roots = TRUE
)

ggsave(
    "results/monocle3/trajectory/monocle3_trajectory_pseudotime.png",
    p_pt,
    width = 8,
    height = 6,
    dpi = 200
)

p_time <- plot_cells(
    cds,
    color_cells_by = "h",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    label_roots = FALSE
)

ggsave(
    "results/monocle3/trajectory/monocle3_trajectory_experimental_time.png",
    p_time,
    width = 8,
    height = 6,
    dpi = 200
)

cat("\n===== FILES CREATED =====\n")
cat("monocle3_order.csv\n")
cat("monocle3_model.rds\n")
cat("monocle3_trajectory_pseudotime.png\n")
cat("monocle3_trajectory_experimental_time.png\n")

cat("\nMONOCLE 3 TRAJECTORY COMPLETE\n")
