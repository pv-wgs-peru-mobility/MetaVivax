#!/usr/bin/env Rscript

#===============================================================================
# Genetic differentiation: Weir-Cockerham Fst
# Author: Mahdi Safarpour
#
# Steps:
#   1. Read the analysis-ready, two-allele GT table and sample metadata.
#      Common quality control is upstream; no extra MAF filter is applied.
#   2. Preserve both alleles: C/C -> two C; C/T -> one C and one T.
#      Preserve mixed-call counts for the Weir-Cockerham estimator.
#   3. Calculate pairwise community Fst for all infections, non-travelers
#      and travelers, plus traveler/non-traveler Fst within each community.
#   4. Use diveRsity to compute the Weir-Cockerham variance components.
#      Pool components across SNPs, then calculate sum(a)/sum(a+b+c).
#   5. Obtain 95% percentile confidence intervals from 1,000 locus bootstraps.
#   6. For within-community tests, permute whole-infection travel labels
#      1,000 times, keeping group sizes fixed; use Bonferroni adjustment.
#   7. Save summary and per-locus tables for plot_fst.R.
#===============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(diveRsity)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gt_table", default = "data/analysis_ready/pv_56samples_16485snps_gt.tsv"),
  make_option("--metadata", default = "data/analysis_ready/pv_56samples_metadata.tsv"),
  make_option("--out_dir", default = "results/fst"),
  make_option("--n_boot", type = "integer", default = 1000L),
  make_option("--n_perm", type = "integer", default = 1000L),
  make_option("--seed", type = "integer", default = 1L)
)))
if (opt$n_boot < 2L || opt$n_perm < 0L) stop("Use n_boot >= 2 and n_perm >= 0.")
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

# --------------------------- Read genotypes and metadata ----------------------
gt <- fread(opt$gt_table, colClasses = "character", na.strings = NULL)
meta <- fread(opt$metadata, colClasses = "character")

if (!all(c("CHROM", "POS", "TYPE") %in% names(gt))) {
  stop("Genotype table must contain CHROM, POS, TYPE and sample GT columns.")
}
sample_cols <- setdiff(names(gt), c("CHROM", "POS", "TYPE", "REF", "ALT"))
sample_ids <- sub("\\.GT$", "", sample_cols)
if (length(sample_ids) < 2L || anyDuplicated(sample_ids)) {
  stop("At least two uniquely named samples are required.")
}
if (nrow(gt) == 0L || any(tolower(gt$TYPE) != "snp") ||
    anyDuplicated(paste(gt$CHROM, gt$POS, sep = ":"))) {
  stop("Use the analysis-ready table of uniquely identified biallelic SNPs.")
}

id_col <- intersect(c("id", "Sample", "sample", "sample_id"), names(meta))[1]
if (is.na(id_col) || !all(c("community", "travel") %in% names(meta))) {
  stop("Metadata must contain a sample ID column, community and travel (0/1).")
}
setnames(meta, id_col, "id")
meta[, id := sub("\\.GT$", "", trimws(id))]
if (anyDuplicated(meta$id) || any(!sample_ids %in% meta$id)) {
  stop("Every genotype sample must have one matching metadata row.")
}
meta <- meta[match(sample_ids, id)]
key <- tolower(gsub("[[:space:]_-]", "",
                    iconv(trimws(meta$community), to = "ASCII//TRANSLIT")))
meta[, community := c(libertad = "Libertad", gamitanacocha = "Gamitanacocha",
                       urcomirano = "Urco Mirano")[key]]
if (anyNA(meta$community) || anyNA(meta$travel) ||
    any(!meta$travel %in% c("0", "1"))) {
  stop("Check community names and travel coding (0 = non-traveler, 1 = traveler).")
}
meta[, travel := as.integer(travel)]

# Preserve both alleles; phase and allele order do not change these summaries.
calls <- as.matrix(gt[, ..sample_cols])
calls[] <- toupper(trimws(gsub("|", "/", calls, fixed = TRUE)))
missing <- is.na(calls) | calls %in% c("", ".", "NA", "NAN", "./.") |
  grepl("^[01ACGT]/\\.$|^\\./[01ACGT]$", calls)
calls[missing] <- NA_character_
if (any(!is.na(calls) & !grepl("^([01]/[01]|[ACGT]/[ACGT])$", calls))) {
  stop("Use complete GT calls (0/0, 0/1, 1/1 or A/A, A/G, etc.), not single-allele numeric codes.")
}

allele1 <- allele2 <- calls
allele1[] <- substr(calls, 1, 1)
allele2[] <- substr(calls, 3, 3)
alleles <- lapply(seq_len(nrow(gt)), function(i) {
  sort(unique(c(allele1[i, ], allele2[i, ])), na.last = NA)
})
if (any(lengths(alleles) > 2L) ||
    any(vapply(alleles, function(z) any(z %in% c("0", "1")) &&
                 any(z %in% c("A", "C", "G", "T")), logical(1)))) {
  stop("Each SNP must use one consistent biallelic genotype encoding.")
}
allele_A <- vapply(alleles, function(z) if (length(z)) z[1] else NA_character_,
                   character(1))
allele_B <- vapply(alleles, function(z) if (length(z) == 2L) z[2] else NA_character_,
                   character(1))
dosage <- (allele1 != allele_A) + (allele2 != allele_A)
colnames(dosage) <- sample_ids
loci <- gt[, .(CHROM, POS)]
loci[, `:=`(allele_A = allele_A, allele_B = allele_B)]
community_order <- c("Libertad", "Gamitanacocha", "Urco Mirano")
category_order <- c("Overall", "Non-travelers", "Travelers")
message("Input: ", ncol(dosage), " infections; ", nrow(dosage), " SNPs.")

# --------------------------- Weir-Cockerham components -------------------------
# Both allele frequencies AND mixed-call counts are needed by this estimator.
wc_components <- function(mat1, mat2) {
  n1 <- rowSums(!is.na(mat1))
  n2 <- rowSums(!is.na(mat2))
  p1 <- rowSums(mat1, na.rm = TRUE) / (2 * n1)
  p2 <- rowSums(mat2, na.rm = TRUE) / (2 * n2)
  n_het <- rowSums(mat1 == 1, na.rm = TRUE) + rowSums(mat2 == 1, na.rm = TRUE)
  indices <- which(n1 >= 2L & n2 >= 2L)
  rows <- lapply(indices, function(i) {
    af <- rbind(A = c(1 - p1[i], 1 - p2[i]), B = c(p1[i], p2[i]))
    hsum <- c(A = n_het[i], B = n_het[i])
    vc <- glbWCcpp(hsum = hsum, af = af, indtyp = c(n1[i], n2[i]))
    data.table(locus_index = i, n1 = n1[i], n2 = n2[i],
                a = vc$a, b = vc$b, c = vc$c)
  })
  if (!length(rows)) return(data.table(locus_index = integer(), n1 = integer(),
                                       n2 = integer(), a = numeric(), b = numeric(),
                                       c = numeric(), denominator = numeric()))
  result <- rbindlist(rows)
  result[, denominator := a + b + c]
  result[is.finite(a) & is.finite(denominator) & denominator > 0]
}
pooled_fst <- function(components) {
  if (!nrow(components)) return(NA_real_)
  sum(components$a) / sum(components$denominator)
}

# --------------------------- Define comparisons -------------------------------
comparisons <- list()
for (category_name in category_order) {
  for (pair in combn(community_order, 2, simplify = FALSE)) {
    comparisons[[length(comparisons) + 1L]] <- list(
      category = category_name, group1 = pair[1], group2 = pair[2],
      ids1 = meta[community == pair[1] &
                    (category_name == "Overall" |
                       travel == as.integer(category_name == "Travelers")), id],
      ids2 = meta[community == pair[2] &
                    (category_name == "Overall" |
                       travel == as.integer(category_name == "Travelers")), id]
    )
  }
}
for (comm in community_order) {
  comparisons[[length(comparisons) + 1L]] <- list(
    category = "Within-community", group1 = paste(comm, "Non-travelers"),
    group2 = paste(comm, "Travelers"),
    ids1 = meta[community == comm & travel == 0L, id],
    ids2 = meta[community == comm & travel == 1L, id]
  )
}

# --------------------------- Estimates, confidence intervals and tests --------
summaries <- per_locus <- list()
for (comparison in comparisons) {
  message(comparison$category, ": ", comparison$group1, " vs ", comparison$group2)
  ids1 <- comparison$ids1
  ids2 <- comparison$ids2
  components <- wc_components(dosage[, ids1, drop = FALSE], dosage[, ids2, drop = FALSE])
  observed <- pooled_fst(components)
  ci <- c(NA_real_, NA_real_)
  p_perm <- NA_real_
  if (nrow(components)) {
    boot <- replicate(opt$n_boot, {
      rows <- sample.int(nrow(components), replace = TRUE)
      pooled_fst(components[rows])
    })
    ci <- quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

    if (comparison$category == "Within-community" && opt$n_perm > 0L) {
      all_ids <- c(ids1, ids2)
      permuted <- replicate(opt$n_perm, {
        shuffled <- sample(all_ids)
        perm_ids1 <- shuffled[seq_along(ids1)]
        perm_ids2 <- shuffled[seq.int(length(ids1) + 1L, length(shuffled))]
        pooled_fst(wc_components(dosage[, perm_ids1, drop = FALSE],
                                 dosage[, perm_ids2, drop = FALSE]))
      })
      valid <- is.finite(permuted)
      p_perm <- (1 + sum(permuted[valid] >= observed)) / (1 + sum(valid))
      if (!any(valid)) p_perm <- NA_real_
    }
    locus_rows <- copy(loci[components$locus_index])
    locus_rows[, `:=`(category = comparison$category, group1 = comparison$group1,
                       group2 = comparison$group2, n1 = components$n1, n2 = components$n2,
                       a = components$a, b = components$b, c = components$c,
                       Fst_locus = components$a / components$denominator)]
    per_locus[[length(per_locus) + 1L]] <- locus_rows
  }
  summaries[[length(summaries) + 1L]] <- data.table(
    category = comparison$category, group1 = comparison$group1, group2 = comparison$group2,
    n_group1 = length(ids1), n_group2 = length(ids2), n_loci = nrow(components),
    Fst = observed, CI_lower = ci[1], CI_upper = ci[2], p_permutation = p_perm
  )
}
summary <- rbindlist(summaries)
summary[, p_adjusted_bonferroni := NA_real_]
summary[category == "Within-community",
        p_adjusted_bonferroni := p.adjust(p_permutation, method = "bonferroni")]
fwrite(summary, file.path(opt$out_dir, "fst_summary.tsv"), sep = "\t")
fwrite(rbindlist(per_locus), file.path(opt$out_dir, "fst_per_locus.tsv"), sep = "\t")
fwrite(meta[, .(n_total = .N, n_nontravelers = sum(travel == 0L),
                n_travelers = sum(travel == 1L)), by = community],
       file.path(opt$out_dir, "fst_sample_counts.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(opt$out_dir, "sessionInfo.txt"))
message("Done. Outputs: ", opt$out_dir)

