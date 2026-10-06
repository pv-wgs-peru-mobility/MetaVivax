#!/usr/bin/env Rscript

#===============================================================================
# Hierarchical clustering from pairwise SNP distance
# Author: Mahdi Safarpour
# Purpose: Run average-linkage clustering from the 1 - PS distance matrix
#===============================================================================

# --------------------------- Define paths -------------------------------------
distance_file <- "results/pairwise_allele_sharing/pairwise_distance_matrix.tsv"
metadata_file <- "data/analysis_ready/pv_56samples_metadata.tsv"
outdir <- "results/hierarchical_clustering"
k_clusters <- 4

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

# --------------------------- Hierarchical clustering --------------------------
cat("Running average-linkage hierarchical clustering...\n")

hc <- hclust(as.dist(dist_mat), method = "average")
clusters <- cutree(hc, k = k_clusters)

cluster_table <- data.frame(
  Sample = names(clusters),
  cluster = as.integer(clusters),
  community = meta[names(clusters), "community"],
  traveler = meta[names(clusters), "traveler"],
  stringsAsFactors = FALSE
)

write.table(cluster_table, file.path(outdir, "hierarchical_cluster_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --------------------------- Plot dendrogram ----------------------------------
cat("Creating dendrogram...\n")

community_colors <- c(
  "Libertad" = "#3D7BFF",
  "Gamitanacocha" = "#00C853",
  "Urco Miraño" = "#FF7A00"
)
cluster_colors <- c("#00C853", "#3D7BFF", "#FF5252", "#FFB300")

hc_plot <- hc
hc_plot$height <- hc$height^(1 / 4)

ordered_samples <- hc$labels[hc$order]
ordered_meta <- meta[ordered_samples, ]
tip_colors <- community_colors[ordered_meta$community]
tip_colors[is.na(tip_colors)] <- "grey80"

plot_dendrogram <- function() {
  par(mar = c(4, 4, 3, 8), xpd = NA)
  plot(
    hc_plot,
    labels = FALSE,
    hang = -1,
    main = paste0("Hierarchical clustering (K = ", k_clusters, ")"),
    xlab = "Samples",
    ylab = "1 - PS distance"
  )
  rect.hclust(hc_plot, k = k_clusters, border = cluster_colors[seq_len(k_clusters)])
  points(seq_along(ordered_samples), rep(0, length(ordered_samples)),
         pch = 21, bg = tip_colors, col = "black", cex = 0.9)
  traveler_pos <- which(ordered_meta$traveler)
  points(traveler_pos, rep(0, length(traveler_pos)), pch = 4, col = "black", cex = 0.6, lwd = 1.1)
  legend("topright", inset = c(-0.33, 0), legend = names(community_colors),
         pch = 21, pt.bg = community_colors, pt.cex = 1.3, bty = "n", title = "Community")
}

pdf(file.path(outdir, "hierarchical_clustering_k4.pdf"), width = 11, height = 7)
plot_dendrogram()
dev.off()

png(file.path(outdir, "hierarchical_clustering_k4.png"), width = 3300, height = 2100, res = 300)
plot_dendrogram()
dev.off()

cat("Done. Output directory:", outdir, "\n")
