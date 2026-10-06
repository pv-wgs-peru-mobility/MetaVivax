#!/usr/bin/env Rscript

# Identity-by-descent (IBD) segment inference and threshold network analysis.
# Run with run_IBD.slurm. Paths and selected settings can be overridden with
# the IBD_* environment variables documented in README.md.

suppressPackageStartupMessages({
  library(isoRelate)
  library(igraph)
})

options(stringsAsFactors = FALSE)

ibd_root <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD"
input_dir <- Sys.getenv("IBD_INPUT_DIR", file.path(ibd_root, "New", "77"))
metadata_file <- Sys.getenv("IBD_METADATA", file.path(ibd_root, "SNP_matrix_with_metadata.csv"))
out_dir <- Sys.getenv("IBD_OUT_DIR", file.path(ibd_root, "results"))
prefix <- Sys.getenv("IBD_PREFIX", "77.19298")
cores <- as.integer(Sys.getenv("IBD_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "3")))
# The manuscript's 40% cutoff was upstream QC; IBD itself uses 30%.
# An upstream count of 56 does not establish the count after IBD filtering.
expected_samples <- as.integer(Sys.getenv("IBD_EXPECTED_SAMPLES", "0"))
isolate_max_missing <- as.numeric(Sys.getenv("IBD_ISOLATE_MAX_MISSING", "0.3"))
snp_max_missing <- 0.6
maf <- 0.01
minimum_snps <- 450L
minimum_length_bp <- 700000L
genotyping_error <- 0.001
bp_per_cm <- as.numeric(Sys.getenv("IBD_BP_PER_CM", "13700"))
map_mode <- Sys.getenv("IBD_MAP_MODE", "from_bp")
# PED column 5 must contain validated per-isolate isoRelate model labels:
# 1 = single infection; 2 = multiple infection (not necessarily two clones).
moi_mode <- Sys.getenv("IBD_MOI_MODE", "ped")
thresholds_pct <- c(5L, 20L, 40L, 60L, 80L)

if (is.na(cores) || cores < 1L || is.na(expected_samples) || expected_samples < 0L ||
    is.na(isolate_max_missing) || isolate_max_missing < 0 || isolate_max_missing > 1 ||
    is.na(bp_per_cm) || bp_per_cm <= 0) {
  stop("Invalid IBD_CORES, IBD_EXPECTED_SAMPLES, IBD_ISOLATE_MAX_MISSING, or IBD_BP_PER_CM.")
}
if (!map_mode %in% c("from_bp", "input")) stop("IBD_MAP_MODE must be 'from_bp' or 'input'.")
if (!moi_mode %in% c("all_two", "ped")) stop("IBD_MOI_MODE must be 'all_two' or 'ped'.")

ped_path <- file.path(input_dir, paste0(prefix, ".ped"))
map_path <- file.path(input_dir, paste0(prefix, ".map"))
for (path in c(ped_path, map_path, metadata_file)) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
save_pdf <- function(path, expr, width = 10, height = 8) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
clean_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.(TGT|GT)$", "", x, ignore.case = TRUE)
}
clean_travel <- function(x) {
  value <- tolower(trimws(as.character(x)))
  yes <- value %in% c("1", "yes", "y", "true", "traveler", "traveller")
  no <- value %in% c("0", "no", "n", "false", "non-traveler", "non-traveller",
                   "nontraveler", "nontraveller")
  if (anyNA(value) || any(!(yes | no))) {
    stop("Every retained sample needs a valid travel value (0/1, no/yes, or non-traveler/traveler).")
  }
  as.integer(yes)
}
finite_or_na <- function(value) if (is.finite(value)) value else NA_real_

message("Reading PED, MAP, and metadata...")
ped <- read.table(ped_path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
map <- read.table(map_path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(ped) < 8L || ncol(map) != 4L || ncol(ped) != 6L + 2L * nrow(map)) {
  stop("PED must have six pedigree fields and two allele columns per SNP in the four-column MAP.")
}
names(ped)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")
names(map) <- c("chr", "snp_id", "pos_cM", "pos_bp")
ped$fid <- as.character(ped$fid)
ped$iid <- as.character(ped$iid)
if (anyNA(ped$fid) || anyNA(ped$iid) || any(!nzchar(ped$fid)) ||
    any(!nzchar(ped$iid)) || anyDuplicated(paste(ped$fid, ped$iid, sep = "/"))) {
  stop("PED family/isolate combinations must be present and unique.")
}
if (moi_mode == "all_two") {
  ped$moi <- rep(2L, nrow(ped))
  warning("MOI mode all_two assumes every isolate is mixed; this reproduces the supplied scripts' assumption, not an estimated MOI.")
} else {
  moi_labels <- trimws(as.character(ped$moi))
  if (anyNA(moi_labels) || any(!moi_labels %in% c("1", "2"))) {
    stop("With IBD_MOI_MODE=ped, PED column 5 must contain curated 1 (single) or 2 (mixed) labels; ordinary PLINK sex codes are not MOI.")
  }
  ped$moi <- as.integer(moi_labels)
  message("MOI mode ped: using PED column 5 as curated single/mixed labels; verify that these are MOI classifications, not PLINK sex codes.")
  if (all(ped$moi == 2L)) {
    warning("PED column 5 labels every isolate mixed (2). Verify this was established for each isolate.")
  }
}

map$pos_bp <- suppressWarnings(as.numeric(as.character(map$pos_bp)))
if (anyNA(map$pos_bp) || any(!is.finite(map$pos_bp)) || any(map$pos_bp <= 0)) {
  stop("MAP base-pair positions must be positive numbers.")
}
if (map_mode == "from_bp") {
  map$pos_cM <- map$pos_bp / bp_per_cm
} else {
  map$pos_cM <- suppressWarnings(as.numeric(as.character(map$pos_cM)))
}
if (anyNA(map$pos_cM) || any(!is.finite(map$pos_cM)) || any(map$pos_cM <= 0)) {
  stop("MAP genetic distances must be positive numbers in cM.")
}
for (idx in split(seq_len(nrow(map)), as.character(map$chr))) {
  if (any(diff(map$pos_bp[idx]) < 0) || any(diff(map$pos_cM[idx]) < 0)) {
    stop("MAP markers must be ordered by increasing position within each chromosome.")
  }
}

message("Filtering genotypes with isoRelate...")
geno <- isoRelate::getGenotypes(
  ped.map = list(ped, map), reference.ped.map = NULL, maf = maf,
  isolate.max.missing = isolate_max_missing,
  snp.max.missing = snp_max_missing,
  input.map.distance = "cM", reference.map.distance = "cM"
)
if (!all(c("pedigree", "genotypes") %in% names(geno))) {
  stop("Unexpected isoRelate genotype result; expected pedigree and genotypes.")
}
filtered_ped <- geno$pedigree
if (nrow(filtered_ped) < 2L) stop("Fewer than two isolates remain after filtering.")
if (expected_samples > 0L && nrow(filtered_ped) != expected_samples) {
  stop("Expected ", expected_samples, " isolates after filtering; found ", nrow(filtered_ped),
       ". Check input/QC or set IBD_EXPECTED_SAMPLES=0 to disable this check.")
}

# Resolve the sample names in metadata against the filtered PED. Do not silently
# drop missing records or collapse family/isolate identifiers.
id_col <- if ("Sample" %in% names(meta)) "Sample" else if ("id" %in% names(meta)) "id" else NA_character_
if (is.na(id_col) || !"travel" %in% names(meta)) {
  stop("Metadata must have Sample (or id) and travel columns.")
}
meta_ids <- clean_id(meta[[id_col]])
if (anyNA(meta_ids) || any(!nzchar(meta_ids)) || anyDuplicated(meta_ids)) {
  stop("Metadata sample IDs must be nonempty and unique after removing .TGT/.GT suffixes.")
}
iso_id <- paste(filtered_ped$fid, filtered_ped$iid, sep = "/")
id_candidates <- list(
  full = clean_id(iso_id),
  iid = clean_id(filtered_ped$iid),
  fid = clean_id(filtered_ped$fid)
)
valid <- vapply(id_candidates, function(ids) !anyDuplicated(ids) && all(ids %in% meta_ids), logical(1))
if (!any(valid)) {
  stop("Metadata IDs do not match all filtered PED isolates by FID/IID, IID, or FID. Examples: ",
       paste(head(iso_id, 4L), collapse = ", "))
}
matches <- lapply(id_candidates[valid], match, table = meta_ids)
if (length(matches) > 1L && any(vapply(matches[-1L], function(x) !identical(x, matches[[1L]]), logical(1)))) {
  stop("Metadata IDs match more than one PED identifier scheme with different sample assignments.")
}
meta_aligned <- meta[matches[[1L]], , drop = FALSE]
sample_map <- data.frame(
  iso_id = iso_id,
  sample_id = meta_ids[matches[[1L]]],
  community = if ("community" %in% names(meta_aligned)) trimws(as.character(meta_aligned$community)) else NA_character_,
  travel = clean_travel(meta_aligned$travel),
  stringsAsFactors = FALSE
)
message("Matched ", nrow(sample_map), " isolates to metadata using ", names(matches)[1L], ".")
save_csv(sample_map, file.path(out_dir, paste0(prefix, "_sample_map.csv")))
saveRDS(geno, file.path(out_dir, paste0(prefix, "_genotypes.rds")))

message("Estimating IBD parameters on ", cores, " cores...")
params <- isoRelate::getIBDparameters(ped.genotypes = geno, number.cores = cores)
saveRDS(params, file.path(out_dir, paste0(prefix, "_parameters.rds")))

message("Inferring IBD segments...")
segments <- isoRelate::getIBDsegments(
  ped.genotypes = geno, parameters = params, number.cores = cores,
  minimum.snps = minimum_snps, minimum.length.bp = minimum_length_bp,
  error = genotyping_error
)
if (!is.data.frame(segments)) stop("isoRelate did not return a segment data frame.")
saveRDS(segments, file.path(out_dir, paste0(prefix, "_segments.rds")))
save_csv(segments, file.path(out_dir, paste0(prefix, "_segments.csv")))

if (nrow(segments) > 0L) {
  summary_lines <- capture.output(isoRelate::getIBDsummary(ped.genotypes = geno, ibd.segments = segments))
  writeLines(summary_lines, file.path(out_dir, paste0(prefix, "_segment_summary.txt")))
  # These two outputs are SNP-level summaries. They are not pairwise genome fractions.
  ibd_matrix <- isoRelate::getIBDmatrix(geno, segments)
  snp_proportion <- isoRelate::getIBDproportion(geno, ibd_matrix)
  saveRDS(ibd_matrix, file.path(out_dir, paste0(prefix, "_snp_ibd_matrix.rds")))
  save_csv(snp_proportion, file.path(out_dir, paste0(prefix, "_snp_ibd_proportion.csv")))
  save_pdf(file.path(out_dir, paste0(prefix, "_snp_ibd_proportion.pdf")),
           isoRelate::plotIBDproportions(snp_proportion, plot.title = "Pairs IBD at each SNP"))
  save_pdf(file.path(out_dir, paste0(prefix, "_ibd_segments.pdf")),
           isoRelate::plotIBDsegments(geno, segments, plot.title = "Detected IBD segments"))
} else {
  writeLines("No IBD segments met the segment filters.",
             file.path(out_dir, paste0(prefix, "_segment_summary.txt")))
}

# Match isoRelate's genome-fraction definition used by getIBDpclusters:
# sum of each pair's segment lengths / sum of filtered SNP spans by chromosome.
marker_map <- geno$genotypes[, c("chr", "pos_bp")]
chr_spans <- vapply(split(marker_map$pos_bp, as.character(marker_map$chr)),
                    function(x) max(x) - min(x), numeric(1))
genome_span_bp <- sum(chr_spans)
if (!is.finite(genome_span_bp) || genome_span_bp <= 0) {
  stop("Filtered markers do not define a positive genome span.")
}

indices <- utils::combn(seq_len(nrow(sample_map)), 2L)
pairs <- data.frame(
  sample_1 = sample_map$sample_id[indices[1L, ]],
  sample_2 = sample_map$sample_id[indices[2L, ]],
  iso_1 = sample_map$iso_id[indices[1L, ]],
  iso_2 = sample_map$iso_id[indices[2L, ]],
  n_segments = integer(ncol(indices)),
  ibd_length_bp = numeric(ncol(indices)),
  stringsAsFactors = FALSE
)
if (nrow(segments) > 0L) {
  required <- c("fid1", "iid1", "fid2", "iid2", "length_bp")
  if (!all(required %in% names(segments))) stop("IBD segments have unexpected columns.")
  first <- match(paste(segments$fid1, segments$iid1, sep = "/"), sample_map$iso_id)
  second <- match(paste(segments$fid2, segments$iid2, sep = "/"), sample_map$iso_id)
  if (anyNA(first) || anyNA(second) || any(first == second)) {
    stop("IBD segment sample identifiers do not match the filtered pedigree.")
  }
  lengths <- suppressWarnings(as.numeric(as.character(segments$length_bp)))
  if (anyNA(lengths) || any(!is.finite(lengths)) || any(lengths < 0)) {
    stop("IBD segment lengths must be finite and nonnegative.")
  }
  pair_key <- paste(indices[1L, ], indices[2L, ], sep = "_")
  segment_key <- paste(pmin(first, second), pmax(first, second), sep = "_")
  pair_index <- match(segment_key, pair_key)
  if (anyNA(pair_index)) stop("An IBD segment could not be matched to a sample pair.")
  sums <- tapply(lengths, pair_index, sum)
  pairs$ibd_length_bp[as.integer(names(sums))] <- as.numeric(sums)
  pairs$n_segments <- tabulate(pair_index, nbins = nrow(pairs))
}
pairs$ibd_fraction <- pairs$ibd_length_bp / genome_span_bp
if (any(pairs$ibd_fraction > 1 + 1e-8)) {
  stop("At least one pair has IBD fraction >1: overlapping segments need review.")
}
pairs$ibd_fraction[pairs$ibd_fraction > 1] <- 1
save_csv(pairs, file.path(out_dir, paste0(prefix, "_pairwise_ibd_fraction.csv")))

ibd_fraction_matrix <- matrix(0, nrow(sample_map), nrow(sample_map),
                              dimnames = list(sample_map$sample_id, sample_map$sample_id))
ibd_fraction_matrix[indices[1L, ] + (indices[2L, ] - 1L) * nrow(sample_map)] <- pairs$ibd_fraction
ibd_fraction_matrix[indices[2L, ] + (indices[1L, ] - 1L) * nrow(sample_map)] <- pairs$ibd_fraction
diag(ibd_fraction_matrix) <- 1
save_csv(data.frame(sample = rownames(ibd_fraction_matrix), ibd_fraction_matrix,
                    check.names = FALSE),
         file.path(out_dir, paste0(prefix, "_pairwise_ibd_fraction_matrix.csv")))

community_colors <- c("Libertad" = "#3D7BFF", "Gamitanacocha" = "#00C853",
                      "UrcoMirano" = "#FF7A00", "Urco Miraño" = "#FF7A00")

save_pdf(file.path(out_dir, paste0(prefix, "_pairwise_ibd_heatmap.pdf")), {
  ord <- order(sample_map$community, sample_map$sample_id, na.last = TRUE)
  heatmap(ibd_fraction_matrix[ord, ord, drop = FALSE], Rowv = NA, Colv = NA,
          scale = "none", col = colorRampPalette(c("white", "#7B3294"))(100),
          breaks = seq(0, 1, length.out = 101), margins = c(9, 9),
          main = "Pairwise genome fraction IBD")
}, width = 12, height = 11)

network_rows <- list()
mixing_rows <- list()
node_rows <- list()
component_rows <- list()
for (threshold_pct in thresholds_pct) {
  threshold <- threshold_pct / 100
  tag <- sprintf("thr%02d", threshold_pct)
  threshold_dir <- file.path(out_dir, tag)
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)
  message("Building IBD network at ", threshold_pct, "%...")

  selected <- pairs[pairs$ibd_fraction >= threshold,
                    c("sample_1", "sample_2", "ibd_fraction", "ibd_length_bp", "n_segments")]
  names(selected)[1:2] <- c("from", "to")
  selected$community_1 <- sample_map$community[match(selected$from, sample_map$sample_id)]
  selected$community_2 <- sample_map$community[match(selected$to, sample_map$sample_id)]
  selected$travel_1 <- sample_map$travel[match(selected$from, sample_map$sample_id)]
  selected$travel_2 <- sample_map$travel[match(selected$to, sample_map$sample_id)]
  vertices <- data.frame(name = sample_map$sample_id, community = sample_map$community,
                         travel = sample_map$travel)
  graph <- igraph::graph_from_data_frame(selected, directed = FALSE, vertices = vertices)
  if (igraph::vcount(graph) != nrow(sample_map) || igraph::ecount(graph) != nrow(selected)) {
    stop("Network vertex or edge count is inconsistent at ", tag, ".")
  }
  save_csv(selected, file.path(threshold_dir, paste0(prefix, "_", tag, "_edges.csv")))
  saveRDS(graph, file.path(threshold_dir, paste0(prefix, "_", tag, "_network.rds")))
  igraph::write_graph(graph, file.path(threshold_dir, paste0(prefix, "_", tag, "_network.graphml")),
                      format = "graphml")

  vertex_names <- igraph::V(graph)$name
  vertex_travel <- sample_map$travel[match(vertex_names, sample_map$sample_id)]
  vertex_community <- sample_map$community[match(vertex_names, sample_map$sample_id)]
  degree <- igraph::degree(graph)
  components <- igraph::components(graph)
  local_clustering <- igraph::transitivity(graph, type = "localundirected", isolates = "zero")
  local_clustering[!is.finite(local_clustering)] <- 0
  node_data <- data.frame(
    threshold_pct = threshold_pct, sample_id = vertex_names,
    community = vertex_community, travel = vertex_travel,
    degree = as.integer(degree), betweenness = as.numeric(igraph::betweenness(graph, directed = FALSE)),
    local_clustering = as.numeric(local_clustering),
    component_id = as.integer(components$membership),
    component_size = as.integer(components$csize[components$membership]),
    is_isolate = degree == 0,
    stringsAsFactors = FALSE
  )
  node_rows[[tag]] <- node_data
  save_csv(node_data, file.path(threshold_dir, paste0(prefix, "_", tag, "_nodes.csv")))

  component_data <- do.call(rbind, lapply(seq_len(components$no), function(k) {
    members <- which(components$membership == k)
    data.frame(threshold_pct = threshold_pct, component_id = k,
               n_samples = length(members),
               n_travelers = sum(vertex_travel[members] == 1L),
               n_nontravelers = sum(vertex_travel[members] == 0L),
               members = paste(vertex_names[members], collapse = ";"))
  }))
  component_rows[[tag]] <- component_data
  save_csv(component_data, file.path(threshold_dir, paste0(prefix, "_", tag, "_components.csv")))

  no_no <- sum(selected$travel_1 == 0L & selected$travel_2 == 0L)
  yes_yes <- sum(selected$travel_1 == 1L & selected$travel_2 == 1L)
  no_yes <- sum(selected$travel_1 != selected$travel_2)
  assortativity <- if (igraph::ecount(graph) > 0L && length(unique(vertex_travel)) > 1L) {
    finite_or_na(igraph::assortativity_nominal(graph, vertex_travel, directed = FALSE))
  } else NA_real_
  mixing_rows[[tag]] <- data.frame(
    threshold_pct = threshold_pct, no_no = no_no, yes_yes = yes_yes, no_yes = no_yes,
    pct_no_no = if (nrow(selected)) 100 * no_no / nrow(selected) else NA_real_,
    pct_yes_yes = if (nrow(selected)) 100 * yes_yes / nrow(selected) else NA_real_,
    pct_no_yes = if (nrow(selected)) 100 * no_yes / nrow(selected) else NA_real_,
    assortativity_travel = assortativity
  )
  save_csv(mixing_rows[[tag]], file.path(threshold_dir, paste0(prefix, "_", tag, "_travel_mixing.csv")))

  transitivity <- igraph::transitivity(graph, type = "globalundirected")
  network_rows[[tag]] <- data.frame(
    threshold_pct = threshold_pct, n_samples = igraph::vcount(graph),
    n_total_pairs = nrow(pairs), n_edges = igraph::ecount(graph),
    density = igraph::edge_density(graph, loops = FALSE),
    n_connected_samples = sum(degree > 0), n_isolates = sum(degree == 0),
    n_components = components$no,
    n_clusters_2plus = sum(components$csize >= 2L),
    n_singletons = sum(components$csize == 1L),
    largest_component = max(components$csize),
    mean_degree = mean(degree),
    transitivity = finite_or_na(transitivity),
    triangle_count = length(igraph::triangles(graph)) / 3,
    stringsAsFactors = FALSE
  )
  save_csv(network_rows[[tag]], file.path(threshold_dir, paste0(prefix, "_", tag, "_summary.csv")))

  set.seed(123L + threshold_pct)
  layout <- if (igraph::ecount(graph) > 0L) igraph::layout_with_fr(graph) else igraph::layout_in_circle(graph)
  colors <- unname(community_colors[vertex_community])
  colors[is.na(colors)] <- "grey75"
  save_csv(data.frame(sample_id = vertex_names, x = layout[, 1L], y = layout[, 2L]),
           file.path(threshold_dir, paste0(prefix, "_", tag, "_layout.csv")))
  save_pdf(file.path(threshold_dir, paste0(prefix, "_", tag, "_network.pdf")), {
    par(mar = c(1, 1, 2, 1))
    plot(graph, layout = layout, vertex.size = 9, vertex.color = colors,
         vertex.frame.color = "black", vertex.label = ifelse(vertex_travel == 1L, "x", ""),
         vertex.label.color = "black", vertex.label.dist = 0, vertex.label.cex = 0.75,
         edge.color = "grey65", edge.width = 0.8,
         main = paste0("IBD network: ", threshold_pct, "% threshold"))
    legend("topleft", legend = c("Libertad", "Gamitanacocha", "Urco Miraño", "Traveler: x"),
           pch = c(21, 21, 21, NA),
           pt.bg = c("#3D7BFF", "#00C853", "#FF7A00", NA),
           bty = "n", cex = 0.8)
  })
}

network_summary <- do.call(rbind, network_rows)
travel_mixing <- do.call(rbind, mixing_rows)
node_metrics <- do.call(rbind, node_rows)
component_summary <- do.call(rbind, component_rows)
if (any(diff(network_summary$n_edges) > 0L)) {
  stop("Edge counts increased with a stricter IBD threshold; check the pairwise data.")
}
save_csv(network_summary, file.path(out_dir, paste0(prefix, "_threshold_network_summary.csv")))
save_csv(travel_mixing, file.path(out_dir, paste0(prefix, "_threshold_travel_mixing.csv")))
save_csv(node_metrics, file.path(out_dir, paste0(prefix, "_threshold_node_metrics.csv")))
save_csv(component_summary, file.path(out_dir, paste0(prefix, "_threshold_component_summary.csv")))

save_pdf(file.path(out_dir, paste0(prefix, "_threshold_overview.pdf")), {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(network_summary$threshold_pct, network_summary$n_edges, type = "b", pch = 19,
       xlab = "IBD threshold (%)", ylab = "Edges", main = "IBD links")
  plot(network_summary$threshold_pct, network_summary$n_connected_samples, type = "b", pch = 19,
       xlab = "IBD threshold (%)", ylab = "Samples", main = "Connected samples")
  plot(network_summary$threshold_pct, network_summary$n_components, type = "b", pch = 19,
       xlab = "IBD threshold (%)", ylab = "Components", main = "Including singletons")
  plot(network_summary$threshold_pct, network_summary$density, type = "b", pch = 19,
       xlab = "IBD threshold (%)", ylab = "Density", main = "Network density")
  par(op)
})

writeLines(c(
  paste("Input PED:", ped_path), paste("Input MAP:", map_path),
  paste("Metadata:", metadata_file), paste("Filtered samples:", nrow(sample_map)),
  paste("Filtered SNPs:", nrow(geno$genotypes)), paste("Detected segments:", nrow(segments)),
  paste("Travelers:", sum(sample_map$travel == 1L)),
  paste("Genome span from filtered markers (bp):", genome_span_bp),
  paste("MOI mode:", moi_mode), paste("Map mode:", map_mode),
  paste("Input PED labeled single infection:", sum(ped$moi == 1L)),
  paste("Input PED labeled multiple infection:", sum(ped$moi == 2L)),
  paste("BP per cM:", if (map_mode == "from_bp") bp_per_cm else "input MAP distances"),
  paste("Isolate max missing:", isolate_max_missing),
  paste("SNP max missing:", snp_max_missing), paste("MAF:", maf),
  paste("Minimum segment SNPs:", minimum_snps),
  paste("Minimum segment length (bp):", minimum_length_bp),
  paste("Genotyping error:", genotyping_error),
  paste("Thresholds (%):", paste(thresholds_pct, collapse = ", "))
), file.path(out_dir, paste0(prefix, "_run_settings.txt")))
writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0(prefix, "_session_info.txt")))
message("IBD analysis complete: ", nrow(sample_map), " samples; outputs in ", out_dir)
