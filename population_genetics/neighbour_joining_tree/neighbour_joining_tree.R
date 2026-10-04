#!/usr/bin/env Rscript

# Exploratory neighbour-joining tree from the shared pairwise SNP allele-sharing
# distance.

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
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
              default = file.path("results", "neighbour_joining_tree"),
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

tree <- nj(as.dist(d))
if (length(tree$tip.label) != nrow(d) ||
    !setequal(tree$tip.label, rownames(d)) ||
    any(!is.finite(tree$edge.length))) {
  stop("Neighbour-joining returned an invalid tree.")
}
if (any(tree$edge.length < 0)) {
  warning("Neighbour-joining produced negative branch lengths; inspect the Newick file. ",
          "The topology-only figures do not display branch lengths.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
write.tree(tree, file = file.path(opt$out_dir, "neighbour_joining_tree.newick"))

community_colors <- c("Libertad" = "#3D7BFF",
                      "Gamitanacocha" = "#00C853",
                      "Urco Miraño" = "#FF7A00")
tip_meta <- meta[match(tree$tip.label, meta$Sample), , drop = FALSE]
tip_colors <- unname(community_colors[tip_meta$community])
traveler_tips <- which(tip_meta$travel == 1L)

plot_tree <- function(type) {
  par(mar = c(1, 1, 2, 1))
  plot.phylo(tree, type = type, use.edge.length = FALSE,
             show.tip.label = FALSE, edge.width = 1, no.margin = FALSE)
  tiplabels(tip = seq_along(tree$tip.label), pch = 21,
            bg = tip_colors, col = "black", cex = 1.1)
  if (length(traveler_tips) > 0L) {
    tiplabels(tip = traveler_tips, pch = 4,
              col = "black", cex = 0.65, lwd = 1.5)
  }
  title(main = paste("Neighbour-joining tree (", type, "; topology only)", sep = ""))
  legend("topright", legend = names(community_colors), pch = 21,
         pt.bg = unname(community_colors), col = "black", pt.cex = 1.3,
         bty = "n", cex = 0.85, title = "Community")
  legend("topleft", legend = "x inside tip: traveler-associated infection",
         bty = "n", cex = 0.75)
}

for (type in c("unrooted", "fan")) {
  pdf(file.path(opt$out_dir, paste0("neighbour_joining_", type, ".pdf")),
      width = 10, height = 10, useDingbats = FALSE)
  plot_tree(type)
  dev.off()
  png(file.path(opt$out_dir, paste0("neighbour_joining_", type, ".png")),
      width = 3000, height = 3000, res = 300)
  plot_tree(type)
  dev.off()
}

message("Neighbour-joining completed: ", length(tree$tip.label), " samples.")
message("Output: ", opt$out_dir)
