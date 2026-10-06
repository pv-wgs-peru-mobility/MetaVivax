#!/usr/bin/env Rscript

#===============================================================================
# Minimum spanning tree from pairwise SNP distance
# Author: Mahdi Safarpour
# Purpose: Build an MST from the 1 - PS distance matrix
#===============================================================================

# --------------------------- Load libraries -----------------------------------
suppressPackageStartupMessages({
  library(igraph)
})

# --------------------------- Define paths -------------------------------------
distance_file <- "results/pairwise_allele_sharing/pairwise_distance_matrix.tsv"
metadata_file <- "data/analysis_ready/pv_56samples_metadata.tsv"
outdir <- "results/minimum_spanning_tree"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- Read distance matrix -----------------------------
cat("Loading pairwise distance matrix...\n")

dist_tab <- read.delim(distance_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
sample_ids <- as.character(dist_tab$sample)

dist_mat <- as.matrix(dist_tab[, -1])
rownames(dist_mat) <- sample_ids
storage.mode(dist_mat) <- "numeric"

dist_mat <- dist_mat[sample_ids, sample_ids]
rownames(dist_mat) <- sub("\\.GT$", "", rownames(dist_mat))
colnames(dist_mat) <- sub("\\.GT$", "", colnames(dist_mat))

diag(dist_mat) <- 0

if (any(is.na(dist_mat))) stop("Distance matrix contains NA values.")
cat("Samples:", nrow(dist_mat), "\n")

# --------------------------- Read metadata ------------------------------------
cat("Loading sample metadata...\n")

if (grepl("\\.csv$", metadata_file, ignore.case = TRUE)) {
  meta <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  meta <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

sample_col <- intersect(c("Sample", "sample", "id"), names(meta))[1]
if (is.na(sample_col)) stop("Metadata must contain Sample, sample, or id column.")

meta$Sample <- sub("\\.GT$", "", as.character(meta[[sample_col]]))
meta <- meta[match(rownames(dist_mat), meta$Sample), ]
if (any(is.na(meta$Sample))) stop("Some samples in the distance matrix are missing from metadata.")
rownames(meta) <- meta$Sample

meta$community <- as.character(meta$community)
meta$community[meta$community == "UrcoMirano"] <- "Urco Miraño"

travel_text <- tolower(as.character(meta$travel))
meta$traveler <- travel_text %in% c("1", "yes", "y", "true", "traveler", "traveller")

# --------------------------- Build MST ----------------------------------------
cat("Building minimum spanning tree...\n")

pairs <- which(upper.tri(dist_mat), arr.ind = TRUE)
edges <- data.frame(
  from = rownames(dist_mat)[pairs[, 1]],
  to = colnames(dist_mat)[pairs[, 2]],
  weight = dist_mat[pairs],
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(edges, directed = FALSE)
mst_net <- mst(g, weights = E(g)$weight)

V(mst_net)$community <- meta[V(mst_net)$name, "community"]
V(mst_net)$traveler <- meta[V(mst_net)$name, "traveler"]

# --------------------------- Save tables --------------------------------------
mst_edges <- as_data_frame(mst_net, what = "edges")
names(mst_edges)[names(mst_edges) == "weight"] <- "distance_1_minus_PS"
write.table(mst_edges, file.path(outdir, "mst_edges.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write_graph(mst_net, file.path(outdir, "mst_graph.graphml"), format = "graphml")

# --------------------------- Plot MST -----------------------------------------
cat("Creating MST plot...\n")

community_colors <- c(
  "Libertad" = "#3D7BFF",
  "Gamitanacocha" = "#00C853",
  "Urco Miraño" = "#FF7A00"
)

node_colors <- community_colors[V(mst_net)$community]
node_colors[is.na(node_colors)] <- "grey80"
