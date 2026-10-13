#!/usr/bin/env Rscript

#===============================================================================
# Fst figures
# Author: Mahdi Safarpour
#
# Steps:
#   1. Read fst_summary.tsv produced by FST.R; no genotypes are recalculated.
#   2. Create symmetric community heatmaps for all infections, non-travelers
#      and travelers, corresponding to the main differentiation figure.
#   3. Plot within-community traveler/non-traveler Fst with 95% confidence
#      intervals and Bonferroni-adjusted permutation p-values.
#===============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--out_dir", default = "results/fst")
)))
fst <- fread(file.path(opt$out_dir, "fst_summary.tsv"))
communities <- c("Libertad", "Gamitanacocha", "Urco Mirano")
categories <- c("Overall", "Non-travelers", "Travelers")

# --------------------------- Community heatmaps -------------------------------
heat <- fst[category %in% categories, .(category, row = group1, col = group2, Fst)]
heat <- rbind(heat, heat[, .(category, row = col, col = row, Fst)])
diagonal <- CJ(category = categories, row = communities)
diagonal[, `:=`(col = row, Fst = 0)]
heat <- rbind(heat, diagonal)
heat[, category := factor(category, levels = categories)]
heat[, row := factor(row, levels = rev(communities))]
heat[, col := factor(col, levels = communities)]
heat[, label := ifelse(is.finite(Fst), sprintf("%.3f", Fst), "NA")]

p_heat <- ggplot(heat, aes(col, row, fill = Fst)) +
  geom_tile(color = "black", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3.5) +
  facet_wrap(~category, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey90", name = "Fst") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(opt$out_dir, "fst_community_heatmaps.pdf"), p_heat, width = 11.5, height = 4.7, bg = "white")
ggsave(file.path(opt$out_dir, "fst_community_heatmaps.png"), p_heat, width = 11.5, height = 4.7, dpi = 300, bg = "white")

# --------------------------- Within-community confidence intervals ------------
within <- fst[category == "Within-community"]
within[, community := factor(sub(" Non-travelers$", "", group1), levels = communities)]
within[, label := ifelse(is.finite(Fst), sprintf("Fst = %.3f", Fst), "NA")]
within[is.finite(p_adjusted_bonferroni), label := paste0(
  label, "\nAdjusted p = ", format.pval(p_adjusted_bonferroni, digits = 3, eps = 0.001)
)]
range_y <- range(c(within$CI_lower, within$CI_upper, 0), na.rm = TRUE)
padding <- max(diff(range_y), 0.1) * 0.1
within[, label_y := CI_upper + padding]

p_within <- ggplot(within, aes(community, Fst)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.12, linewidth = 0.8) +
  geom_point(color = "red", size = 2.5) +
  geom_text(aes(y = label_y, label = label), size = 3.5, vjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  labs(x = NULL, y = "Weir-Cockerham Fst",
       title = "Non-travelers versus travelers within each community") +
  theme_classic(base_size = 12)
ggsave(file.path(opt$out_dir, "fst_within_community.pdf"), p_within, width = 8.5, height = 5.2)
ggsave(file.path(opt$out_dir, "fst_within_community.png"), p_within, width = 8.5, height = 5.2, dpi = 300)
message("Fst figures saved to: ", opt$out_dir)
