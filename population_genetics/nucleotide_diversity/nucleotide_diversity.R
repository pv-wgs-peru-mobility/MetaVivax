#!/usr/bin/env Rscript

#===============================================================================
# SNP-based nucleotide diversity
# Author: Mahdi Safarpour
#
# Steps:
#   1. Read the analysis-ready, two-allele GT table and sample metadata.
#      Use all 16,485 polymorphic biallelic SNPs; common QC is upstream.
#   2. Count both observed alleles: C/C -> two C; C/T -> one C and one T.
#      Missing or incomplete calls contribute neither allele.
#   3. At each SNP and within each community/travel group, calculate
#      pi = 1 - sum(p_i^2). Require at least two non-missing infections.
#      Retain zero values when a SNP is monomorphic within a group.
#   4. Save locus values and group summaries.
#   5. Compare communities using Friedman tests matched by SNP, followed
#      by paired Wilcoxon tests. Use rank-sum tests for within-community
#      traveler/non-traveler comparisons and Bonferroni adjustment.
#   6. Draw the nine-group violin plot.
#===============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gt_table", default = "data/analysis_ready/pv_56samples_16485snps_gt.tsv"),
  make_option("--metadata", default = "data/analysis_ready/pv_56samples_metadata.tsv"),
  make_option("--out_dir", default = "results/nucleotide_diversity")
)))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
min_calls <- 2L

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

# --------------------------- SNP-level diversity ------------------------------
pi_rows <- list()
for (comm in community_order) {
  for (category_name in category_order) {
    ids <- meta[community == comm &
                  (category_name == "Overall" |
                     travel == as.integer(category_name == "Travelers")), id]
    group_dosage <- dosage[, ids, drop = FALSE]
    n <- rowSums(!is.na(group_dosage))
    p_B <- rowSums(group_dosage, na.rm = TRUE) / (2 * n)
    pi_locus <- 2 * p_B * (1 - p_B)
    pi_locus[n < min_calls] <- NA_real_
    group_rows <- copy(loci)
    group_rows[, `:=`(community = comm, category = category_name,
                       n_group_samples = length(ids), n_nonmissing = n,
                       p_A = 1 - p_B, p_B = p_B, pi_locus = pi_locus)]
    group_rows[!is.finite(p_B), `:=`(p_A = NA_real_, p_B = NA_real_)]
    pi_rows[[length(pi_rows) + 1L]] <- group_rows
  }
}
pi_values <- rbindlist(pi_rows)
summary <- pi_values[is.finite(pi_locus), .(
  n_group_samples = n_group_samples[1], n_loci = .N,
  n_zero_loci = sum(pi_locus == 0), mean_pi = mean(pi_locus),
  median_pi = median(pi_locus), Q1_pi = quantile(pi_locus, 0.25),
  Q3_pi = quantile(pi_locus, 0.75)
), by = .(community, category)]
fwrite(pi_values, file.path(opt$out_dir, "pi_per_locus.tsv"), sep = "\t")
fwrite(summary, file.path(opt$out_dir, "pi_summary.tsv"), sep = "\t")

# --------------------------- Matched community comparisons --------------------
global_tests <- pairwise_tests <- list()
community_pairs <- combn(community_order, 2, simplify = FALSE)
for (category_name in category_order) {
  wide <- dcast(pi_values[category == category_name], CHROM + POS ~ community,
                value.var = "pi_locus")
  paired <- wide[complete.cases(wide[, ..community_order])]
  if (nrow(paired) < 2L) next
  mat <- as.matrix(paired[, ..community_order])
  fried <- friedman.test(mat)
  global_tests[[length(global_tests) + 1L]] <- data.table(
    category = category_name, test = "Friedman",
    n_matched_loci = nrow(paired), statistic = unname(fried$statistic),
    p_value = fried$p.value
  )
  tests <- rbindlist(lapply(community_pairs, function(pair) {
    wt <- wilcox.test(paired[[pair[1]]], paired[[pair[2]]],
                      paired = TRUE, exact = FALSE)
    data.table(category = category_name, community1 = pair[1], community2 = pair[2],
               test = "Paired Wilcoxon signed-rank", n_matched_loci = nrow(paired),
               p_value = wt$p.value)
  }))
  tests[!is.finite(p_value), p_value := NA_real_]
  tests[, p_adjusted_bonferroni := p.adjust(p_value, method = "bonferroni")]
  pairwise_tests[[length(pairwise_tests) + 1L]] <- tests
}
fwrite(rbindlist(global_tests), file.path(opt$out_dir, "pi_community_tests.tsv"), sep = "\t")
fwrite(rbindlist(pairwise_tests), file.path(opt$out_dir, "pi_community_pairwise_tests.tsv"), sep = "\t")

# --------------------------- Within-community travel comparisons --------------
travel_tests <- list()
for (comm in community_order) {
  x <- pi_values[community == comm & category == "Non-travelers" & is.finite(pi_locus), pi_locus]
  y <- pi_values[community == comm & category == "Travelers" & is.finite(pi_locus), pi_locus]
  p <- if (length(x) && length(y)) {
    wilcox.test(x, y, paired = FALSE, exact = FALSE)$p.value
  } else NA_real_
  travel_tests[[length(travel_tests) + 1L]] <- data.table(
    community = comm, test = "Wilcoxon rank-sum", p_value = p
  )
}
travel_tests <- rbindlist(travel_tests)
travel_tests[!is.finite(p_value), p_value := NA_real_]
travel_tests[, p_adjusted_bonferroni := p.adjust(p_value, method = "bonferroni")]
fwrite(travel_tests, file.path(opt$out_dir, "pi_travel_tests.tsv"), sep = "\t")

# --------------------------- Violin plot --------------------------------------
plot_data <- copy(pi_values[is.finite(pi_locus)])
plot_data[, community := factor(community, levels = community_order)]
plot_data[, category := factor(category, levels = category_order)]
p <- ggplot(plot_data, aes(category, pi_locus, fill = category)) +
  geom_violin(trim = TRUE, scale = "width", adjust = 2, color = "black") +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 1.5) +
  facet_grid(. ~ community, drop = FALSE) +
  scale_fill_manual(values = c("Overall" = "#F4A582", "Non-travelers" = "#D6604D",
                               "Travelers" = "#9E3D3D"), guide = "none") +
  scale_x_discrete(labels = c("Overall", "Non-\ntravelers", "Travelers"), drop = FALSE) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = NULL, y = expression("SNP-based diversity (" * pi * ")")) +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA))
ggsave(file.path(opt$out_dir, "pi_violin.pdf"), p, width = 11.5, height = 5.8)
ggsave(file.path(opt$out_dir, "pi_violin.png"), p, width = 11.5, height = 5.8, dpi = 300)
writeLines(capture.output(sessionInfo()), file.path(opt$out_dir, "sessionInfo.txt"))
message("Done. Outputs: ", opt$out_dir)

