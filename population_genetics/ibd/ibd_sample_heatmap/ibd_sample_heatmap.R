#!/usr/bin/env Rscript

#===============================================================================
# Final sample-level pairwise IBD heatmap
# Samples ordered by community and travel status
# With community labels and travel-status symbols
# Updated: separator lines are restricted to the heatmap body
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

metadata_file <- file.path(
  root_outdir,
  "02_metadata_clean",
  "metadata_clean.csv"
)

outdir <- file.path(root_outdir, "12_final_selected_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- Read data ----------------------------------------

pairwise <- read_csv(pairwise_file, show_col_types = FALSE)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

required_pairwise_cols <- c("sample1", "sample2", "ibd_prop")
required_metadata_cols <- c("sample_id", "community", "travel_status")

missing_pairwise_cols <- setdiff(required_pairwise_cols, names(pairwise))
missing_metadata_cols <- setdiff(required_metadata_cols, names(metadata))

if (length(missing_pairwise_cols) > 0) {
  stop("Missing columns in pairwise file: ", paste(missing_pairwise_cols, collapse = ", "))
}

if (length(missing_metadata_cols) > 0) {
  stop("Missing columns in metadata file: ", paste(missing_metadata_cols, collapse = ", "))
}

# --------------------------- Define sample order -------------------------------

community_order <- c("Libertad", "Gamitanacocha", "UrcoMirano")
travel_order <- c("Non-traveler", "Traveler")

metadata_ordered <- metadata %>%
  mutate(
    community = factor(community, levels = community_order),
    travel_status = factor(travel_status, levels = travel_order)
  ) %>%
  arrange(community, travel_status, sample_id) %>%
  mutate(order_index = row_number())

sample_order <- metadata_ordered$sample_id
n_samples <- length(sample_order)

# --------------------------- Build symmetric heatmap matrix --------------------

heat_df <- pairwise %>%
  select(sample1, sample2, ibd_prop) %>%
  bind_rows(
    pairwise %>%
      transmute(
        sample1 = sample2,
        sample2 = sample1,
        ibd_prop = ibd_prop
      )
  ) %>%
  bind_rows(
    tibble(
      sample1 = sample_order,
      sample2 = sample_order,
      ibd_prop = 1
    )
  ) %>%
  left_join(
    metadata_ordered %>%
      select(sample1 = sample_id, x_index = order_index),
    by = "sample1"
  ) %>%
  left_join(
    metadata_ordered %>%
      select(sample2 = sample_id, y_index_raw = order_index),
    by = "sample2"
  ) %>%
  filter(!is.na(x_index), !is.na(y_index_raw)) %>%
  mutate(
    y_index = n_samples - y_index_raw + 1
  )

write_csv(
  heat_df,
  file.path(outdir, "final_sample_level_ibd_heatmap_metadata_ordered_plot_data.csv")
)

# --------------------------- Community label positions -------------------------

community_ranges <- metadata_ordered %>%
  group_by(community) %>%
  summarise(
    x_min = min(order_index),
    x_max = max(order_index),
    x_mid = mean(c(x_min, x_max)),
    y_min_raw = min(order_index),
    y_max_raw = max(order_index),
    y_mid_raw = mean(c(y_min_raw, y_max_raw)),
    .groups = "drop"
  ) %>%
  mutate(
    y_mid = n_samples - y_mid_raw + 1,

    # Top labels
    x_label = case_when(
      community == "UrcoMirano" ~ x_mid + 1.4,
      TRUE ~ x_mid
    ),
    y_label_top = n_samples + 4.2,

    # Left labels
    x_label_left = -4.8,
    y_label_left = case_when(
      community == "UrcoMirano" ~ y_mid - 1.7,
      TRUE ~ y_mid
    ),
    left_label_size = case_when(
      community == "UrcoMirano" ~ 4.0,
      TRUE ~ 4.5
    )
  )

boundary_positions <- community_ranges %>%
  arrange(x_min) %>%
  filter(row_number() < n()) %>%
  transmute(
    boundary_x = x_max + 0.5,
    boundary_y = n_samples - x_max + 0.5
  )

# --------------------------- Travel-status symbols -----------------------------

travel_symbols_top <- metadata_ordered %>%
  mutate(
    x = order_index,
    y = n_samples + 1.5
  ) %>%
  select(sample_id, community, travel_status, x, y)

travel_symbols_left <- metadata_ordered %>%
  mutate(
    x = -1.2,
    y = n_samples - order_index + 1
  ) %>%
  select(sample_id, community, travel_status, x, y)

travel_symbols <- bind_rows(travel_symbols_top, travel_symbols_left)

# --------------------------- Plot ---------------------------------------------

p_heat <- ggplot() +

  geom_tile(
    data = heat_df,
    aes(x = x_index, y = y_index, fill = ibd_prop),
    colour = "white",
    linewidth = 0.15,
    width = 0.98,
    height = 0.98
  ) +

  # Separator lines restricted to heatmap body only
  geom_segment(
    data = boundary_positions,
    aes(
      x = boundary_x,
      xend = boundary_x,
      y = 0.5,
      yend = n_samples + 0.5
    ),
    colour = "grey25",
    linewidth = 0.75
  ) +
  geom_segment(
    data = boundary_positions,
    aes(
      x = 0.5,
      xend = n_samples + 0.5,
      y = boundary_y,
      yend = boundary_y
    ),
    colour = "grey25",
    linewidth = 0.75
  ) +

  geom_point(
    data = travel_symbols,
    aes(x = x, y = y, shape = travel_status),
    size = 2.4,
    colour = "#4A0000",
    stroke = 0.8
  ) +

  geom_text(
    data = community_ranges,
    aes(x = x_label, y = y_label_top, label = community),
    fontface = "bold",
    size = 4.5,
    colour = "black"
  ) +

  geom_text(
    data = community_ranges,
    aes(
      x = x_label_left,
      y = y_label_left,
      label = community,
      size = left_label_size
    ),
    fontface = "bold",
    colour = "black",
    angle = 90
  ) +

  scale_size_identity() +

  scale_fill_gradientn(
    colours = c("#fff5f0", "#fee0d2", "#fcbba1", "#fb6a4a", "#cb181d"),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    name = "IBD"
  ) +

  scale_shape_manual(
    values = c(
      "Non-traveler" = 16,
      "Traveler" = 17
    ),
    name = "Travel status"
  ) +

  coord_fixed(
    xlim = c(-5.8, n_samples + 0.8),
    ylim = c(0.5, n_samples + 5.2),
    clip = "off"
  ) +

  labs(
    title = "Sample-level pairwise IBD heatmap",
    subtitle = "Samples ordered by community and travel status; symbols indicate travel status",
    x = NULL,
    y = NULL
  ) +

  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0),
    plot.subtitle = element_text(size = 13, hjust = 0, margin = margin(b = 12)),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(20, 25, 20, 50)
  )

# --------------------------- Save outputs -------------------------------------

ggsave(
  filename = file.path(outdir, "final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.pdf"),
  plot = p_heat,
  width = 10.5,
  height = 9.2,
  units = "in"
)

ggsave(
  filename = file.path(outdir, "final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.png"),
  plot = p_heat,
  width = 10.5,
  height = 9.2,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(outdir, "final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.tiff"),
  plot = p_heat,
  width = 10.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

cat("Finished.\n")
cat("Corrected annotated sample-level heatmap saved in:\n")
cat(outdir, "\n")