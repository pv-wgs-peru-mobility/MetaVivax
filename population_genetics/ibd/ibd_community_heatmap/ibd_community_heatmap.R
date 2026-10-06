#!/usr/bin/env Rscript

#===============================================================================
# Final professional community-by-community IBD heatmap
# Recalculates directly from pairwise_annotated file
# Threshold: IBD >= 0.50
#===============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

# --------------------------- Paths --------------------------------------------

root_outdir <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis"

pairwise_file <- file.path(
  root_outdir,
  "03_pairwise_ibd_dataset",
  "77.19298_pairwise_ibd_annotated.csv"
)

outdir <- file.path(root_outdir, "12_final_selected_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

main_threshold <- 0.50

# --------------------------- Read pairwise data --------------------------------

pairwise <- read_csv(pairwise_file, show_col_types = FALSE)

required_cols <- c("sample1", "sample2", "ibd_prop", "community1", "community2")
missing_cols <- setdiff(required_cols, names(pairwise))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

community_order <- c("Libertad", "Gamitanacocha", "UrcoMirano")

# --------------------------- Recalculate unordered pair summaries --------------

community_summary <- pairwise %>%
  mutate(
    community1 = as.character(community1),
    community2 = as.character(community2),
    community_low = if_else(
      match(community1, community_order) <= match(community2, community_order),
      community1,
      community2
    ),
    community_high = if_else(
      match(community1, community_order) <= match(community2, community_order),
      community2,
      community1
    ),
    related = ibd_prop >= main_threshold
  ) %>%
  group_by(community_low, community_high) %>%
  summarise(
    n_pairs = n(),
    n_related = sum(related, na.rm = TRUE),
    R = n_related / n_pairs,
    mean_ibd = mean(ibd_prop, na.rm = TRUE),
    median_ibd = median(ibd_prop, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  community_summary,
  file.path(outdir, "final_community_heatmap_threshold_0_50_unordered_summary.csv")
)

# --------------------------- Make symmetric plotting table ---------------------

heat_sym <- bind_rows(
  community_summary %>%
    transmute(
      x_community = community_low,
      y_community = community_high,
      n_pairs,
      n_related,
      R,
      mean_ibd,
      median_ibd
    ),
  community_summary %>%
    filter(community_low != community_high) %>%
    transmute(
      x_community = community_high,
      y_community = community_low,
      n_pairs,
      n_related,
      R,
      mean_ibd,
      median_ibd
    )
) %>%
  distinct(x_community, y_community, .keep_all = TRUE)

heat_complete <- expand_grid(
  x_community = community_order,
  y_community = community_order
) %>%
  left_join(heat_sym, by = c("x_community", "y_community")) %>%
  mutate(
    n_pairs = if_else(is.na(n_pairs), 0L, as.integer(n_pairs)),
    n_related = if_else(is.na(n_related), 0L, as.integer(n_related)),
    R = if_else(is.na(R), 0, R),
    x_community = factor(x_community, levels = community_order),
    y_community = factor(y_community, levels = rev(community_order)),
    R_label = if_else(R == 0, "0", sprintf("%.3f", R)),
    label = paste0(n_related, "/", n_pairs, "\nR = ", R_label),
    label_colour = if_else(R >= 0.55, "white", "#5A0000")
  )

# Check symmetry
symmetry_check <- heat_complete %>%
  select(x_community, y_community, R) %>%
  mutate(
    x_community = as.character(x_community),
    y_community = as.character(y_community)
  ) %>%
  left_join(
    heat_complete %>%
      select(x_community, y_community, R) %>%
      mutate(
        x_community = as.character(x_community),
        y_community = as.character(y_community)
      ) %>%
      rename(
        x_reverse = y_community,
        y_reverse = x_community,
        R_reverse = R
      ),
    by = c("x_community" = "x_reverse", "y_community" = "y_reverse")
  ) %>%
  mutate(R_difference = R - R_reverse)

write_csv(
  symmetry_check,
  file.path(outdir, "final_community_heatmap_threshold_0_50_symmetry_check.csv")
)

write_csv(
  heat_complete,
  file.path(outdir, "final_community_heatmap_threshold_0_50_plot_data.csv")
)

# --------------------------- Plot ---------------------------------------------

p_heat <- ggplot(
  heat_complete,
  aes(x = x_community, y = y_community, fill = R)
) +
  geom_tile(
    colour = "white",
    linewidth = 1.3,
    width = 0.97,
    height = 0.97
  ) +
  geom_text(
    aes(label = label, colour = label_colour),
    size = 4.4,
    lineheight = 0.9,
    fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#fff5f0", "#fee0d2", "#fcbba1", "#fb6a4a", "#cb181d"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    name = "Related-pair\nproportion (R)"
  ) +
  coord_fixed() +
  labs(
    title = "Community-level parasite relatedness",
    subtitle = "Proportion of sample pairs classified as related using IBD = 0.50",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle = element_text(size = 13, hjust = 0, margin = margin(b = 12)),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 13,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 20, 15, 15)
  )

# --------------------------- Save outputs -------------------------------------

ggsave(
  filename = file.path(outdir, "final_community_heatmap_threshold_0_50_corrected.pdf"),
  plot = p_heat,
  width = 8.5,
  height = 7.2,
  units = "in"
)

ggsave(
  filename = file.path(outdir, "final_community_heatmap_threshold_0_50_corrected.png"),
  plot = p_heat,
  width = 8.5,
  height = 7.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "final_community_heatmap_threshold_0_50_corrected.tiff"),
  plot = p_heat,
  width = 8.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

cat("Finished.\n")
cat("Corrected final figure saved in:\n")
cat(outdir, "\n")