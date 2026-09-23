############################################################
# Expected heterozygosity analysis for SNP data
#
# Alleles: 1, 2, 3, 4, NA
# He formula: He = 1 - sum(p_i^2)
#
# QC:
# 1. SNP missingness <= 30%
# 2. Global He > 0
# 3. Global MAF >= 0.10
#
# For the boxplot:
# Plot only loci with He > 0 within each group
############################################################

library(tidyverse)

############################################################
# STEP 0 — Paths
############################################################

input_dir  <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/"
output_dir <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/Expected_He/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

geno_csv <- "77.19298.TGT_nucleotide_numeric.csv"
meta_csv <- "metadata.csv"

geno_path <- file.path(input_dir, geno_csv)
meta_path <- file.path(input_dir, meta_csv)

out_sample_missingness <- file.path(output_dir, "sample_missingness_QC_table.csv")
out_snp_qc             <- file.path(output_dir, "SNP_QC_missingness_He_MAF_table.csv")
out_he_values_all      <- file.path(output_dir, "He_values_all_loci_after_QC.csv")
out_he_values_plot     <- file.path(output_dir, "He_values_polymorphic_loci_for_boxplot.csv")
out_he_summary_all     <- file.path(output_dir, "He_summary_all_loci_after_QC.csv")
out_he_summary_plot    <- file.path(output_dir, "He_summary_polymorphic_loci_for_boxplot.csv")
out_boxplot            <- file.path(output_dir, "He_boxplot_polymorphic_loci.jpeg")

############################################################
# STEP 1 — Read files
############################################################

geno <- read.csv(geno_path, check.names = FALSE, na.strings = c("NA", "", "NaN"))
meta <- read.csv(meta_path, check.names = FALSE)

meta <- meta %>%
  mutate(
    id = as.character(id),
    community = str_trim(community),
    community = case_when(
      str_to_lower(community) %in% c("libertad", "libertda", "librtad") ~ "Libertad",
      str_to_lower(community) %in% c("gamitanacocha", "gamitanachocha") ~ "Gamitanacocha",
      str_to_lower(community) %in% c("urcomirano", "urcomarino") ~ "Urcomirano",
      TRUE ~ community
    ),
    travel = as.integer(travel),
    travel_group = ifelse(travel == 1, "Travelers", "Non-travelers")
  )

############################################################
# STEP 2 — Convert genotype data to long format
############################################################

geno <- geno %>%
  mutate(locus = paste(CHROM, POS, TYPE, sep = "_"))

sample_cols <- setdiff(colnames(geno), c("CHROM", "POS", "TYPE", "locus"))

geno_long <- geno %>%
  select(locus, all_of(sample_cols)) %>%
  pivot_longer(
    cols = -locus,
    names_to = "id",
    values_to = "allele"
  )

merged_long <- geno_long %>%
  left_join(meta, by = "id")

############################################################
# STEP 3 — Sample missingness QC table
############################################################

sample_missingness <- merged_long %>%
  group_by(id, community, travel, travel_group) %>%
  summarise(
    n_loci_total = n(),
    n_loci_missing = sum(is.na(allele)),
    n_loci_non_missing = sum(!is.na(allele)),
    sample_missingness = n_loci_missing / n_loci_total,
    .groups = "drop"
  ) %>%
  arrange(desc(sample_missingness))

write.csv(sample_missingness, out_sample_missingness, row.names = FALSE)

############################################################
# STEP 4 — Functions
############################################################

calc_he <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) return(NA_real_)

  freqs <- table(x) / length(x)
  He <- 1 - sum(freqs^2)

  return(as.numeric(He))
}

calc_maf <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) return(NA_real_)

  freqs <- table(x) / length(x)

  if (length(freqs) < 2) return(0)

  return(as.numeric(min(freqs)))
}

############################################################
# STEP 5 — SNP QC
############################################################

snp_qc <- merged_long %>%
  group_by(locus) %>%
  summarise(
    n_samples_total = n_distinct(id),
    n_missing = sum(is.na(allele)),
    n_non_missing = sum(!is.na(allele)),
    missingness = n_missing / n_samples_total,
    n_alleles_global = n_distinct(allele[!is.na(allele)]),
    global_He = calc_he(allele),
    global_MAF = calc_maf(allele),
    keep_missingness30 = missingness <= 0.30,
    keep_polymorphic_global = !is.na(global_He) & global_He > 0,
    keep_MAF10 = !is.na(global_MAF) & global_MAF >= 0.10,
    keep_final = keep_missingness30 & keep_polymorphic_global & keep_MAF10,
    .groups = "drop"
  )

write.csv(snp_qc, out_snp_qc, row.names = FALSE)

loci_keep <- snp_qc %>%
  filter(keep_final) %>%
  pull(locus)

filtered_long <- merged_long %>%
  filter(locus %in% loci_keep)

cat("\nQC summary:\n")
cat("Total SNPs:", n_distinct(merged_long$locus), "\n")
cat("SNPs retained after missingness <=30%, global He>0, global MAF>=0.10:",
    length(loci_keep), "\n")

############################################################
# STEP 6 — Calculate He for 9 groups
############################################################

community_order <- c("Libertad", "Gamitanacocha", "Urcomirano")
group_order <- c("Overall", "Non-travelers", "Travelers")

he_list <- list()

for (comm in community_order) {

  dat_comm <- filtered_long %>%
    filter(community == comm)

  he_list[[paste0(comm, "_Overall")]] <-
    dat_comm %>%
    group_by(locus) %>%
    summarise(
      He = calc_he(allele),
      n_samples = n_distinct(id[!is.na(allele)]),
      .groups = "drop"
    ) %>%
    mutate(community = comm, group = "Overall")

  he_list[[paste0(comm, "_NonTravelers")]] <-
    dat_comm %>%
    filter(travel == 0) %>%
    group_by(locus) %>%
    summarise(
      He = calc_he(allele),
      n_samples = n_distinct(id[!is.na(allele)]),
      .groups = "drop"
    ) %>%
    mutate(community = comm, group = "Non-travelers")

  he_list[[paste0(comm, "_Travelers")]] <-
    dat_comm %>%
    filter(travel == 1) %>%
    group_by(locus) %>%
    summarise(
      He = calc_he(allele),
      n_samples = n_distinct(id[!is.na(allele)]),
      .groups = "drop"
    ) %>%
    mutate(community = comm, group = "Travelers")
}

he_df <- bind_rows(he_list) %>%
  filter(!is.na(He), He >= 0, He <= 0.5) %>%
  mutate(
    community = factor(community, levels = community_order),
    group = factor(group, levels = group_order),
    x_group = case_when(
      community == "Libertad" & group == "Overall" ~ 1,
      community == "Libertad" & group == "Non-travelers" ~ 2,
      community == "Libertad" & group == "Travelers" ~ 3,
      community == "Gamitanacocha" & group == "Overall" ~ 4,
      community == "Gamitanacocha" & group == "Non-travelers" ~ 5,
      community == "Gamitanacocha" & group == "Travelers" ~ 6,
      community == "Urcomirano" & group == "Overall" ~ 7,
      community == "Urcomirano" & group == "Non-travelers" ~ 8,
      community == "Urcomirano" & group == "Travelers" ~ 9
    )
  ) %>%
  filter(!is.na(x_group))

write.csv(he_df, out_he_values_all, row.names = FALSE)

############################################################
# STEP 7 — Create plotting dataset
# Only loci with He > 0 within each group are plotted.
############################################################

he_df_plot <- he_df %>%
  filter(He > 0)

write.csv(he_df_plot, out_he_values_plot, row.names = FALSE)

############################################################
# STEP 8 — Summary tables
############################################################

he_summary_all <- he_df %>%
  group_by(community, group, x_group) %>%
  summarise(
    n_samples = max(n_samples, na.rm = TRUE),
    n_loci_all = n(),
    n_polymorphic_loci = sum(He > 0),
    percent_polymorphic = 100 * mean(He > 0),
    mean_He_all_loci = mean(He, na.rm = TRUE),
    median_He_all_loci = median(He, na.rm = TRUE),
    Q1_He_all_loci = quantile(He, 0.25, na.rm = TRUE),
    Q3_He_all_loci = quantile(He, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(x_group)

write.csv(he_summary_all, out_he_summary_all, row.names = FALSE)

he_summary_plot <- he_df_plot %>%
  group_by(community, group, x_group) %>%
  summarise(
    n_samples = max(n_samples, na.rm = TRUE),
    n_loci_plotted = n(),
    mean_He = mean(He, na.rm = TRUE),
    median_He = median(He, na.rm = TRUE),
    Q1_He = quantile(He, 0.25, na.rm = TRUE),
    Q3_He = quantile(He, 0.75, na.rm = TRUE),
    min_He = min(He, na.rm = TRUE),
    max_He = max(He, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(x_group)

write.csv(he_summary_plot, out_he_summary_plot, row.names = FALSE)

cat("\nSummary using all retained loci:\n")
print(he_summary_all)

cat("\nSummary using polymorphic loci only for boxplot:\n")
print(he_summary_plot)

############################################################
# STEP 9 — Boxplot with 9 groups
############################################################

x_labels <- c(
  "Overall", "Non-travelers", "Travelers",
  "Overall", "Non-travelers", "Travelers",
  "Overall", "Non-travelers", "Travelers"
)

p_box <- ggplot(he_df_plot, aes(x = x_group, y = He, group = x_group)) +
  geom_boxplot(
    fill = "white",
    color = "black",
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.55,
    color = "red",
    linewidth = 0.5
  ) +
  geom_text(
    data = he_summary_plot,
    aes(x = x_group, y = -0.025, label = paste0("n=", n_samples)),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  annotate("text", x = 2, y = -0.33, label = "Libertad", size = 4.5) +
  annotate("text", x = 5, y = -0.33, label = "Gamitanacocha", size = 4.5) +
  annotate("text", x = 8, y = -0.33, label = "Urcomirano", size = 4.5) +
  scale_x_continuous(
    breaks = 1:9,
    labels = x_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.5, 0.1),
    expand = expansion(mult = c(0.18, 0.03))
  ) +
  coord_cartesian(ylim = c(0, 0.5), clip = "off") +
  labs(
    x = "",
    y = "Expected heterozygosity",
    title = "Expected heterozygosity among polymorphic SNPs"
  ) +
  theme_classic() +
  theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
  axis.title.y = element_text(size = 13),
  plot.title = element_text(hjust = 0.5, size = 16),
  plot.margin = margin(10, 10, 140, 10)
)

ggsave(out_boxplot, p_box, width = 11, height = 6.5, dpi = 300)

############################################################
# Final message
############################################################

cat("\nAnalysis completed successfully.\n")
cat("Sample missingness table:", out_sample_missingness, "\n")
cat("SNP QC table:", out_snp_qc, "\n")
cat("All He values:", out_he_values_all, "\n")
cat("Plotted He values:", out_he_values_plot, "\n")
cat("All-loci He summary:", out_he_summary_all, "\n")
cat("Boxplot He summary:", out_he_summary_plot, "\n")
cat("Boxplot:", out_boxplot, "\n")