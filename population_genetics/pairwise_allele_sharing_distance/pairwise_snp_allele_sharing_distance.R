#!/usr/bin/env Rscript

# ============================================================
# Pairwise SNP allele-sharing distance
#
# Purpose:
#   Calculate pairwise genetic similarity (PS) and genetic distance
#   (1 - PS) for minimum spanning tree and hierarchical clustering.
#
# Input:
#   Analysis-ready SNP genotype table containing the final 56 samples
#   and 16,485 polymorphic biallelic SNPs after upstream quality control.
#   The common sample- and SNP-level QC filters are not repeated here.
#
# Allele-sharing calculation:
#   - Missing genotype calls are excluded pairwise at each SNP.
#   - Identical allele composition contributes 1.0 to PS.
#   - Genotypes sharing one allele contribute 0.5 to PS.
#   - Genotypes sharing no alleles contribute 0 to PS.
#   - PS is the mean allele-sharing score across comparable SNPs.
#   - Genetic distance is defined as 1 - PS.
#
# Required input format:
#   Tab-delimited table with columns:
#     CHROM, POS, TYPE, sample1.GT, sample2.GT, ...
#   Optional REF and ALT columns are allowed.
#
# Default project paths:
#   Input : data/analysis_ready/pv_56samples_16485snps_gt.tsv
#   Output: results/pairwise_allele_sharing/
#
# Outputs:
#   1) pairwise_allele_sharing_long.tsv
#   2) pairwise_distance_matrix.tsv
#   3) pairwise_PS_matrix.tsv
#   4) pairwise_shared_loci_matrix.tsv
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(
    "--gt_table",
    type = "character",
    default = file.path("data", "analysis_ready", "pv_56samples_16485snps_gt.tsv"),
    help = "Analysis-ready tab-delimited GT table [default: %default]."
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = file.path("results", "pairwise_allele_sharing"),
    help = "Output directory [default: %default]."
  ),
  make_option(
    "--min_shared_loci",
    type = "integer",
    default = 1L,
    help = "Minimum number of comparable SNPs required to estimate PS [default: %default]."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$gt_table)) {
  stop("Input genotype table not found: ", opt$gt_table)
}
if (opt$min_shared_loci < 1L) {
  stop("--min_shared_loci must be >= 1.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Genotype parsing
# ------------------------------------------------------------
parse_gt <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", ".", "NA", "NaN", "nan", "./.", ".|.")] <- NA_character_

  allele1 <- rep(NA_character_, length(x))
  allele2 <- rep(NA_character_, length(x))

  idx <- which(!is.na(x))
  if (length(idx) > 0L) {
    parts <- strsplit(x[idx], "[/|]")
    valid <- lengths(parts) == 2L &
      vapply(parts, function(z) !any(z %in% c("", ".", "NA", "NaN", "nan")), logical(1))

    valid_idx <- idx[valid]
    valid_parts <- parts[valid]

    if (length(valid_idx) > 0L) {
      allele1[valid_idx] <- vapply(valid_parts, `[[`, character(1), 1L)
      allele2[valid_idx] <- vapply(valid_parts, `[[`, character(1), 2L)
    }
  }

  list(a1 = allele1, a2 = allele2)
}

# ------------------------------------------------------------
# Read analysis-ready genotype table
# ------------------------------------------------------------
message("Reading genotype table: ", opt$gt_table)
gt <- fread(opt$gt_table, data.table = TRUE, na.strings = NULL)

required_cols <- c("CHROM", "POS", "TYPE")
if (!all(required_cols %in% names(gt))) {
  stop("Input table must contain CHROM, POS, and TYPE columns.")
}

metadata_cols <- intersect(c("CHROM", "POS", "TYPE", "REF", "ALT"), names(gt))
sample_cols <- setdiff(names(gt), metadata_cols)

if (length(sample_cols) < 2L) {
  stop("At least two sample genotype columns are required.")
}

message("Samples: ", length(sample_cols))
message("SNP loci: ", nrow(gt))

# ------------------------------------------------------------
# Parse genotype alleles once for efficient pairwise comparison
# ------------------------------------------------------------
message("Parsing genotype calls...")
parsed <- lapply(gt[, ..sample_cols], parse_gt)

allele1 <- as.data.table(lapply(parsed, `[[`, "a1"))
allele2 <- as.data.table(lapply(parsed, `[[`, "a2"))
setnames(allele1, sample_cols)
setnames(allele2, sample_cols)

samples <- sample_cols
n_samples <- length(samples)

ps_mat <- matrix(
  NA_real_, nrow = n_samples, ncol = n_samples,
  dimnames = list(samples, samples)
)
dist_raw_mat <- ps_mat
shared_mat <- matrix(
  0L, nrow = n_samples, ncol = n_samples,
  dimnames = list(samples, samples)
)

long_out <- vector("list", n_samples * (n_samples - 1L) / 2L)
out_idx <- 1L

# ------------------------------------------------------------
# Pairwise allele-sharing similarity and distance
# ------------------------------------------------------------
message("Calculating pairwise allele-sharing distances...")

for (i in seq_len(n_samples)) {
  a1_i <- allele1[[samples[i]]]
  a2_i <- allele2[[samples[i]]]

  for (j in i:n_samples) {
    a1_j <- allele1[[samples[j]]]
    a2_j <- allele2[[samples[j]]]

    ok <- !is.na(a1_i) & !is.na(a2_i) & !is.na(a1_j) & !is.na(a2_j)
    n_shared <- sum(ok)

    if (n_shared >= opt$min_shared_loci) {
      same_order <- as.integer(a1_i[ok] == a1_j[ok]) +
                    as.integer(a2_i[ok] == a2_j[ok])
      swapped_order <- as.integer(a1_i[ok] == a2_j[ok]) +
                       as.integer(a2_i[ok] == a1_j[ok])

      shared_alleles <- pmax(same_order, swapped_order)
      locus_score <- shared_alleles / 2

      ps <- mean(locus_score)
      distance <- 1 - ps

      n_full <- sum(locus_score == 1)
      n_partial <- sum(locus_score == 0.5)
      n_zero <- sum(locus_score == 0)
    } else {
      ps <- NA_real_
      distance <- NA_real_
      n_full <- NA_integer_
      n_partial <- NA_integer_
      n_zero <- NA_integer_
    }

    ps_mat[i, j] <- ps_mat[j, i] <- ps
    dist_raw_mat[i, j] <- dist_raw_mat[j, i] <- distance
    shared_mat[i, j] <- shared_mat[j, i] <- n_shared

    if (i < j) {
      long_out[[out_idx]] <- data.table(
        sample_1 = samples[i],
        sample_2 = samples[j],
        n_shared_nonmissing_loci = n_shared,
        n_full_sharing_loci = n_full,
        n_partial_sharing_loci = n_partial,
        n_zero_sharing_loci = n_zero,
        PS = ps,
        distance_1_minus_PS = distance
      )
      out_idx <- out_idx + 1L
    }
  }
}

long_out <- rbindlist(long_out, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------
# Complete distance matrix for MST / hierarchical clustering
# ------------------------------------------------------------
dist_complete_mat <- dist_raw_mat

off_diag <- row(dist_raw_mat) != col(dist_raw_mat)
observed_distances <- dist_raw_mat[off_diag & is.finite(dist_raw_mat)]

if (length(observed_distances) == 0L) {
  stop("No pairwise distances could be estimated.")
}

max_observed_distance <- max(observed_distances)
missing_pairs <- off_diag & is.na(dist_complete_mat)
dist_complete_mat[missing_pairs] <- max_observed_distance

diag(dist_complete_mat) <- 0
diag(ps_mat) <- 1

if (nrow(long_out) > 0L) {
  long_out[, distance_for_mst_clustering := distance_1_minus_PS]
  long_out[is.na(distance_for_mst_clustering),
           distance_for_mst_clustering := max_observed_distance]
  long_out[, distance_replaced := is.na(distance_1_minus_PS)]
}

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------
fwrite(
  long_out,
  file.path(opt$out_dir, "pairwise_allele_sharing_long.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(dist_complete_mat, keep.rownames = "sample"),
  file.path(opt$out_dir, "pairwise_distance_matrix.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(ps_mat, keep.rownames = "sample"),
  file.path(opt$out_dir, "pairwise_PS_matrix.tsv"),
  sep = "\t"
)

fwrite(
  as.data.table(shared_mat, keep.rownames = "sample"),
  file.path(opt$out_dir, "pairwise_shared_loci_matrix.tsv"),
  sep = "\t"
)

message("Done.")
message("Output directory: ", opt$out_dir)

