#!/usr/bin/env Rscript

#===============================================================================
# Expected heterozygosity
# Author: Mahdi Safarpour
#
# Steps:
#   1. Read the analysis-ready, two-allele SNP GT table and sample metadata.
#      Common genotype, sample and SNP quality control is applied upstream.
#   2. Retain SNPs with global MAF >= 0.05, as specified in the manuscript.
#   3. Count both observed alleles: C/C -> two C; C/T -> one C and one T.
#      Each infection has equal weight; incomplete calls are excluded.
#   4. Within each community and travel group, calculate SNP-level
#      He = n / (n - 1) * (1 - sum(p_i^2)).
#      n is the number of non-missing infections at that SNP, not 2n alleles.
#   5. Save locus values and summaries, including group-specific zero values.
#   6. Run Kruskal-Wallis, Dunn and within-community Wilcoxon rank-sum tests,
#      using Bonferroni correction for multiple pairwise comparisons.
#   7. Draw the nine-group boxplot.
#===============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dunn.test)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gt_table", default = "data/analysis_ready/pv_56samples_16485snps_gt.tsv"),
  make_option("--metadata", default = "data/analysis_ready/pv_56samples_metadata.tsv"),
  make_option("--out_dir", default = "results/expected_heterozygosity")
)))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
min_maf <- 0.05
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

# --------------------------- Global MAF filter --------------------------------
n_global <- rowSums(!is.na(dosage))
p_B_global <- rowSums(dosage, na.rm = TRUE) / (2 * n_global)
maf <- pmin(p_B_global, 1 - p_B_global)
keep <- n_global >= min_calls & is.finite(maf) & maf >= min_maf
qc <- copy(loci)
qc[, `:=`(n_nonmissing = n_global, global_MAF = maf, retained = keep)]
fwrite(qc, file.path(opt$out_dir, "He_global_MAF.tsv"), sep = "\t")
if (!any(keep)) stop("No SNPs pass the global MAF >= 0.05 filter.")
dosage <- dosage[keep, , drop = FALSE]
loci <- loci[keep]
message("SNPs retained at global MAF >= 0.05: ", nrow(dosage))

# --------------------------- SNP-level He -------------------------------------
he_rows <- list()
for (comm in community_order) {
  for (category_name in category_order) {
    ids <- meta[community == comm &
                  (category_name == "Overall" |
                     travel == as.integer(category_name == "Travelers")), id]
    group_dosage <- dosage[, ids, drop = FALSE]
    n <- rowSums(!is.na(group_dosage))
    p_B <- rowSums(group_dosage, na.rm = TRUE) / (2 * n)
    diversity <- 2 * p_B * (1 - p_B)
    He <- n / (n - 1) * diversity
    He[n < min_calls] <- NA_real_
    group_rows <- copy(loci)
    group_rows[, `:=`(community = comm, category = category_name,
                       n_group_samples = length(ids), n_nonmissing = n,
                       p_A = 1 - p_B, p_B = p_B, He = He)]
    group_rows[!is.finite(p_B), `:=`(p_A = NA_real_, p_B = NA_real_)]
    he_rows[[length(he_rows) + 1L]] <- group_rows
  }
}
he <- rbindlist(he_rows)
summary <- he[is.finite(He), .(
  n_group_samples = n_group_samples[1], n_loci = .N,
  n_zero_loci = sum(He == 0), mean_He = mean(He), median_He = median(He),
  Q1_He = quantile(He, 0.25), Q3_He = quantile(He, 0.75)
), by = .(community, category)]
fwrite(he, file.path(opt$out_dir, "He_per_locus.tsv"), sep = "\t")
fwrite(summary, file.path(opt$out_dir, "He_summary.tsv"), sep = "\t")

# --------------------------- Community comparisons ----------------------------
global_tests <- pairwise_tests <- list()
for (category_name in category_order) {
  sub <- he[category == category_name & is.finite(He)]
  if (uniqueN(sub$community) < 2L || uniqueN(sub$He) < 2L) next
  kw <- kruskal.test(He ~ community, data = sub)
  global_tests[[length(global_tests) + 1L]] <- data.table(
    category = category_name, test = "Kruskal-Wallis",
    statistic = unname(kw$statistic), p_value = kw$p.value
  )
  dt <- dunn.test(sub$He, sub$community, method = "bonferroni",
                  kw = FALSE, altp = TRUE, list = TRUE)
  pairwise_tests[[length(pairwise_tests) + 1L]] <- data.table(
    category = category_name, comparison = dt$comparisons,
    test = "Dunn", p_value = dt$altP, p_adjusted_bonferroni = dt$altP.adjusted
  )
}
fwrite(rbindlist(global_tests), file.path(opt$out_dir, "He_community_tests.tsv"), sep = "\t")
fwrite(rbindlist(pairwise_tests), file.path(opt$out_dir, "He_community_pairwise_tests.tsv"), sep = "\t")

# --------------------------- Travel comparisons -------------------------------
travel_tests <- list()
for (comm in community_order) {
  x <- he[community == comm & category == "Non-travelers" & is.finite(He), He]
  y <- he[community == comm & category == "Travelers" & is.finite(He), He]
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
fwrite(travel_tests, file.path(opt$out_dir, "He_travel_tests.tsv"), sep = "\t")

# --------------------------- Boxplot ------------------------------------------
plot_data <- copy(he[is.finite(He)])
plot_data[, community := factor(community, levels = community_order)]
plot_data[, category := factor(category, levels = category_order)]
p <- ggplot(plot_data, aes(category, He)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA, width = 0.65) +
  stat_summary(fun = median, geom = "point", color = "red", size = 1.8) +
  facet_grid(. ~ community, drop = FALSE) +
  scale_x_discrete(labels = c("Overall", "Non-\ntravelers", "Travelers"), drop = FALSE) +
  labs(x = NULL, y = "Expected heterozygosity (He)") +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(color = "black", fill = NA))
ggsave(file.path(opt$out_dir, "He_boxplot.pdf"), p, width = 11.5, height = 5.8)
ggsave(file.path(opt$out_dir, "He_boxplot.png"), p, width = 11.5, height = 5.8, dpi = 300)
writeLines(capture.output(sessionInfo()), file.path(opt$out_dir, "sessionInfo.txt"))
message("Done. Outputs: ", opt$out_dir)

