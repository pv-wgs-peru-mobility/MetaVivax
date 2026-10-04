#!/usr/bin/env Rscript

# Minimum spanning tree from the shared pairwise SNP allele-sharing distance.
# The distance matrix is produced by pairwise_snp_allele_sharing_distance.R.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

options <- list(
  make_option("--distance_matrix", type = "character",
              default = file.path("results", "pairwise_allele_sharing",
                                  "pairwise_distance_matrix.tsv"),
              help = "Pairwise 1-PS distance matrix [default: %default]."),
  make_option("--metadata", type = "character",
              default = file.path("data", "analysis_ready",
                                  "pv_56samples_metadata.tsv"),
              help = "Sample metadata with Sample, community, travel [default: %default]."),
  make_option("--out_dir", type = "character",
              default = file.path("results", "minimum_spanning_tree"),
              help = "Output directory [default: %default]."),
  make_option("--expected_samples", type = "integer", default = 56L,
              help = "Expected number of samples; 0 disables the check [default: %default]."),
  make_option("--expected_travelers", type = "integer", default = 16L,
              help = "Expected number of travelers; -1 disables the check [default: %default].")
)
opt <- parse_args(OptionParser(option_list = options))

if (opt$expected_samples < 0L || opt$expected_travelers < -1L) {
  stop("Invalid expected sample or traveler count.")
}
if (!file.exists(opt$distance_matrix)) stop("Missing distance matrix: ", opt$distance_matrix)
if (!file.exists(opt$metadata)) stop("Missing metadata: ", opt$metadata)

read_distance_matrix <- function(path, expected_samples) {
  tab <- read.delim(path, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, comment.char = "")
  if (ncol(tab) < 4L || names(tab)[1L] != "sample") {
    stop("Distance matrix must have a first column named 'sample' and at least 3 samples.")
  }

  ids <- trimws(as.character(tab[[1L]]))
  col_ids <- trimws(names(tab)[-1L])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      any(!nzchar(col_ids)) || anyDuplicated(col_ids) ||
      !setequal(ids, col_ids)) {
    stop("Distance matrix row and column sample IDs must be unique and identical.")
  }
  if (nrow(tab) != length(col_ids)) stop("Distance matrix is not square.")
  if (expected_samples > 0L && length(ids) != expected_samples) {
    stop("Expected ", expected_samples, " samples; found ", length(ids), ".")
  }

  d <- as.matrix(tab[, -1L, drop = FALSE])
  suppressWarnings(storage.mode(d) <- "double")
  dimnames(d) <- list(ids, col_ids)
  d <- d[sort(ids), sort(ids), drop = FALSE]

  if (any(!is.finite(d))) stop("Distance matrix contains missing or nonnumeric values.")
  if (any(d < -1e-10 | d > 1 + 1e-10)) stop("1-PS distances must be between 0 and 1.")
  if (any(abs(diag(d)) > 1e-10)) stop("Distance matrix diagonal must be zero.")
  if (any(abs(d - t(d)) > 1e-10)) stop("Distance matrix must be symmetric.")
  d
}

read_metadata <- function(path, d, expected_travelers) {
  sep <- if (grepl("\\.csv$", tolower(path))) "," else "\t"
  meta <- read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                     stringsAsFactors = FALSE, comment.char = "", quote = "\"")
  required <- c("Sample", "community", "travel")
  if (!all(required %in% names(meta))) {
    stop("Metadata must contain: ", paste(required, collapse = ", "), ".")
  }
  meta <- meta[, required, drop = FALSE]
  meta$Sample <- trimws(as.character(meta$Sample))
  if (anyNA(meta$Sample) || any(!nzchar(meta$Sample)) || anyDuplicated(meta$Sample)) {
    stop("Metadata Sample IDs must be present and unique.")
  }

  ids <- rownames(d)
  if (!setequal(ids, meta$Sample)) {
    without_gt <- sub("\\.GT$", "", ids)
    if (all(grepl("\\.GT$", ids)) && !anyDuplicated(without_gt) &&
        setequal(without_gt, meta$Sample)) {
      dimnames(d) <- list(without_gt, without_gt)
      ids <- without_gt
    } else {
      stop("Metadata Sample IDs do not match distance matrix IDs (with or without .GT).")
    }
  }
  if (nrow(meta) != length(ids)) stop("Metadata must have exactly one row per sample.")
  meta <- meta[match(ids, meta$Sample), , drop = FALSE]
  rownames(meta) <- meta$Sample

  key <- iconv(trimws(as.character(meta$community)), to = "ASCII//TRANSLIT")
  key <- gsub("[^a-z]", "", tolower(key))
  community_names <- c(libertad = "Libertad",
                       gamitanacocha = "Gamitanacocha",
                       urcomirano = "Urco Miraño")
  if (anyNA(key) || any(!key %in% names(community_names))) {
    stop("Community must be Libertad, Gamitanacocha, or Urco Miraño.")
  }
  meta$community <- unname(community_names[key])

  travel <- tolower(trimws(as.character(meta$travel)))
  yes <- travel %in% c("1", "traveler", "traveller", "yes", "true")
  no <- travel %in% c("0", "non-traveler", "non-traveller", "nontraveler",
                     "nontraveller", "no", "false")
  if (anyNA(travel) || any(!(yes | no))) {
    stop("Metadata travel must encode traveler as 1 and non-traveler as 0.")
  }
  meta$travel <- as.integer(yes)
  if (expected_travelers >= 0L && sum(meta$travel) != expected_travelers) {
    stop("Expected ", expected_travelers, " travelers; found ", sum(meta$travel), ".")
  }

  list(distance = d, metadata = meta)
}

d <- read_distance_matrix(opt$distance_matrix, opt$expected_samples)
input <- read_metadata(opt$metadata, d, opt$expected_travelers)
d <- input$distance
meta <- input$metadata

# Explicitly include every pair, including pairs with distance zero. A weighted
# adjacency matrix would drop zero-distance edges before constructing the MST.
ij <- which(upper.tri(d), arr.ind = TRUE)
edges <- data.frame(from = rownames(d)[ij[, 1L]],
                    to = colnames(d)[ij[, 2L]],
                    weight = d[ij], stringsAsFactors = FALSE)
g <- graph_from_data_frame(edges, directed = FALSE,
                           vertices = data.frame(name = rownames(d)))
tree <- mst(g, weights = E(g)$weight)
if (vcount(tree) != nrow(d) || ecount(tree) != nrow(d) - 1L ||
    !is_connected(tree)) {
  stop("The MST did not connect every sample.")
}

V(tree)$community <- meta[V(tree)$name, "community"]
V(tree)$travel <- meta[V(tree)$name, "travel"]

mst_edges <- as_data_frame(tree, what = "edges")
names(mst_edges)[names(mst_edges) == "weight"] <- "distance_1_minus_PS"
mst_edges$community_1 <- meta[mst_edges$from, "community"]
mst_edges$community_2 <- meta[mst_edges$to, "community"]
mst_edges$travel_1 <- meta[mst_edges$from, "travel"]
mst_edges$travel_2 <- meta[mst_edges$to, "travel"]
mst_edges <- mst_edges[order(mst_edges$from, mst_edges$to), , drop = FALSE]

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(mst_edges, file.path(opt$out_dir, "mst_edges.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write_graph(tree, file.path(opt$out_dir, "mst_graph.graphml"), format = "graphml")

community_colors <- c("Libertad" = "#3D7BFF",
                      "Gamitanacocha" = "#00C853",
                      "Urco Miraño" = "#FF7A00")
set.seed(323)
lay <- layout_with_fr(tree, niter = 3000L, weights = rep(1, ecount(tree)))

plot_mst <- function() {
  par(mar = c(1, 1, 2, 1))
  plot(tree, layout = lay, vertex.size = 11,
       vertex.color = unname(community_colors[V(tree)$community]),
       vertex.frame.color = "black",
       vertex.label = ifelse(V(tree)$travel == 1L, "x", ""),
       vertex.label.dist = 0, vertex.label.cex = 0.8,
       vertex.label.color = "black", vertex.label.font = 2,
       edge.color = "#444444", edge.width = 1.2, margin = 0.13)
  title(main = "Minimum spanning tree (1 - PS)")
  legend("topright", legend = names(community_colors), pch = 21,
         pt.bg = unname(community_colors), col = "black", pt.cex = 1.3,
         bty = "n", cex = 0.85, title = "Community")
  legend("topleft", legend = "x inside node: traveler-associated infection",
         bty = "n", cex = 0.75)
}

pdf(file.path(opt$out_dir, "mst_allele_sharing.pdf"),
    width = 9, height = 9, useDingbats = FALSE)
plot_mst()
dev.off()
png(file.path(opt$out_dir, "mst_allele_sharing.png"),
    width = 2700, height = 2700, res = 300)
plot_mst()
dev.off()

message("MST completed: ", vcount(tree), " samples, ", ecount(tree), " edges.")
message("Output: ", opt$out_dir)
