#!/usr/bin/env Rscript

#===============================================================================
# IBD relatedness analysis for community and travel comparisons
# Author: Mahdi Safarpour
#
# Purpose:
#   1. Start IBD analysis from PED/MAP files using IsoRelate
#   2. Build a pairwise sample-level IBD table from the IsoRelate IBD matrix
#   3. Merge pairwise IBD with final metadata
#   4. Compare within- and between-community relatedness
#   5. Compare traveler and non-traveler relatedness
#   6. Generate heatmaps, distributions, map-style plots, sensitivity analysis,
#      and networks
#
# Main threshold:
#   IBD >= 0.50
#
# Sensitivity thresholds:
#   IBD >= 0.05 and IBD >= 0.80
#===============================================================================

# --------------------------- Load libraries -----------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(igraph)
  library(isoRelate)
})

# --------------------------- Define paths -------------------------------------

input_dir <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77"
prefix <- "77.19298"

root_outdir <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis"
metadata_path <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/SNP_matrix_with_metadata.csv"

main_threshold <- 0.50
all_thresholds <- c(0.05, 0.50, 0.80)

# --------------------------- Create output folders ----------------------------

dirs <- list(
  logs               = file.path(root_outdir, "00_logs"),
  ibd_from_scratch   = file.path(root_outdir, "01_ibd_from_scratch"),
  metadata_clean     = file.path(root_outdir, "02_metadata_clean"),
  pairwise_data      = file.path(root_outdir, "03_pairwise_ibd_dataset"),
  community_analysis = file.path(root_outdir, "04_community_relatedness"),
  travel_analysis    = file.path(root_outdir, "05_travel_relatedness"),
  within_between     = file.path(root_outdir, "06_within_between_community"),
  heatmaps           = file.path(root_outdir, "07_heatmaps"),
  distributions      = file.path(root_outdir, "08_ibd_distributions"),
  map_connectivity   = file.path(root_outdir, "09_map_connectivity"),
  sensitivity        = file.path(root_outdir, "10_sensitivity_analysis"),
  network            = file.path(root_outdir, "11_networks")
)

walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dirs$logs, paste0(prefix, "_ibd_pairwise_relatedness_analysis_log.txt"))
sink(log_file, split = TRUE)

on.exit({
  while (sink.number() > 0) {
    sink()
  }
}, add = TRUE)

cat("Analysis started:", as.character(Sys.time()), "\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", root_outdir, "\n")
cat("Main threshold:", main_threshold, "\n")
cat("All thresholds:", paste(all_thresholds, collapse = ", "), "\n\n")

#===============================================================================
# Helper functions
#===============================================================================

clean_id <- function(x) {
  x0 <- as.character(x)
  x0 <- stringr::str_trim(x0)
  x0 <- stringr::str_remove(x0, "\\.TGT$")
  x0 <- stringr::str_remove(x0, "\\.tgt$")

  # Keep two-part and three-part sample IDs:
  # 904016-20
  # 904021-16-19
  # Also handles dot or underscore versions:
  # 904016.20
  # 904021.16.19
  m <- stringr::str_match(
    x0,
    "(90[456][0-9]{3})[-._]([0-9]{2})(?:[-._]([0-9]{2}))?"
  )

  out <- ifelse(
    !is.na(m[, 1]) & !is.na(m[, 4]),
    paste0(m[, 2], "-", m[, 3], "-", m[, 4]),
    ifelse(
      !is.na(m[, 1]),
      paste0(m[, 2], "-", m[, 3]),
      x0
    )
  )

  return(out)
}

threshold_label <- function(x) {
  stringr::str_replace(sprintf("%.2f", x), "\\.", "_")
}

save_plot <- function(plot_obj, filename_base, outdir, width = 10, height = 7) {
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height
  )

  ggplot2::ggsave(
    filename = file.path(outdir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

make_pairwise_from_ibd_matrix <- function(ibd_matrix, metadata_ids, output_dir) {
  df <- as.data.frame(ibd_matrix, check.names = FALSE)

  cat("Columns in IBD matrix:\n")
  print(names(df)[1:min(10, ncol(df))])

  if (ncol(df) <= 4) {
    stop("IBD matrix has no pairwise columns. Expected columns 5 onward to be pairwise sample columns.")
  }

  # In IsoRelate, columns 1-4 are SNP information.
  # Columns 5 onward are binary IBD indicators for sample pairs.
  pair_cols <- names(df)[5:ncol(df)]

  cat("Number of pairwise columns in IBD matrix:", length(pair_cols), "\n")

  metadata_ids <- unique(clean_id(metadata_ids))

  # Sort by decreasing length so IDs such as 904021-16-19 are matched
  # before shorter IDs such as 904021-16.
  metadata_ids_sorted <- metadata_ids[order(nchar(metadata_ids), decreasing = TRUE)]

  make_id_variants <- function(id) {
    c(
      id,
      stringr::str_replace_all(id, "-", "."),
      stringr::str_replace_all(id, "-", "_")
    )
  }

  extract_pair_ids <- function(pair_col) {
    work_col <- pair_col
    hits <- character(0)

    for (id in metadata_ids_sorted) {
      variants <- make_id_variants(id)

      found <- any(vapply(
        variants,
        function(v) stringr::str_detect(work_col, stringr::fixed(v)),
        logical(1)
      ))

      if (found) {
        hits <- c(hits, id)

        # Remove all variants after matching.
        # This prevents shorter IDs being falsely matched inside longer IDs.
        for (v in variants) {
          work_col <- stringr::str_replace_all(work_col, stringr::fixed(v), "")
        }
      }
    }

    hits <- unique(hits)

    if (length(hits) == 2) {
      return(tibble::tibble(
        pair_col = pair_col,
        sample1 = hits[1],
        sample2 = hits[2]
      ))
    }

    return(tibble::tibble(
      pair_col = pair_col,
      sample1 = NA_character_,
      sample2 = NA_character_,
      n_hits = length(hits),
      hits_found = paste(hits, collapse = ";")
    ))
  }

  pair_key <- purrr::map_dfr(pair_cols, extract_pair_ids)

  failed_cols <- pair_key %>%
    filter(is.na(sample1) | is.na(sample2))

  readr::write_csv(
    failed_cols,
    file.path(output_dir, "failed_to_parse_ibd_matrix_pair_columns.csv")
  )

  if (nrow(failed_cols) > 0) {
    cat("ERROR: Some IBD matrix pair columns could not be parsed.\n")
    cat("Number of failed columns:", nrow(failed_cols), "\n")
    cat("See file: failed_to_parse_ibd_matrix_pair_columns.csv\n")
    print(head(failed_cols, 10))
    stop("Stopping because some pair columns could not be converted to sample IDs.")
  }

  pair_props <- purrr::map_dbl(pair_cols, function(col) {
    v <- df[[col]]

    if (is.factor(v)) {
      v <- as.character(v)
    }

    v <- as.numeric(v)

    mean(v, na.rm = TRUE)
  })

  out <- pair_key %>%
    select(pair_col, sample1, sample2) %>%
    mutate(
      ibd_prop = pair_props,
      sample1 = clean_id(sample1),
      sample2 = clean_id(sample2)
    ) %>%
    filter(!is.na(sample1), !is.na(sample2)) %>%
    filter(sample1 != sample2) %>%
    mutate(
      pair_id = if_else(
        sample1 < sample2,
        paste(sample1, sample2, sep = "__"),
        paste(sample2, sample1, sep = "__")
      )
    ) %>%
    distinct(pair_id, .keep_all = TRUE)

  cat("Number of unique pairwise comparisons created:", nrow(out), "\n")

  return(out)
}


summarise_related_pairs <- function(data, group_vars, thresholds) {
  purrr::map_dfr(thresholds, function(thr) {
    data %>%
      mutate(
        threshold = thr,
        related = ibd_prop >= thr
      ) %>%
      group_by(across(all_of(group_vars)), threshold) %>%
      summarise(
        n_pairs = n(),
        n_related = sum(related, na.rm = TRUE),
        R = n_related / n_pairs,
        mean_ibd = mean(ibd_prop, na.rm = TRUE),
        median_ibd = median(ibd_prop, na.rm = TRUE),
        min_ibd = min(ibd_prop, na.rm = TRUE),
        max_ibd = max(ibd_prop, na.rm = TRUE),
        .groups = "drop"
      )
  })
}


pairwise_prop_tests <- function(data, group_col, threshold) {
  df <- data %>%
    mutate(
      related = ibd_prop >= threshold,
      group = as.character(.data[[group_col]])
    ) %>%
    filter(!is.na(group))

  groups <- sort(unique(df$group))

  if (length(groups) < 2) {
    return(tibble::tibble())
  }

  combs <- combn(groups, 2, simplify = FALSE)

  out <- purrr::map_dfr(combs, function(g) {
    sub <- df %>%
      filter(group %in% g)

    tab <- table(sub$group, sub$related)

    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble::tibble(
        group_1 = g[1],
        group_2 = g[2],
        threshold = threshold,
        test = NA_character_,
        p_value = NA_real_
      ))
    }

    expected <- suppressWarnings(chisq.test(tab)$expected)

    if (any(expected < 5)) {
      test_res <- fisher.test(tab)
      test_name <- "Fisher exact test"
    } else {
      test_res <- chisq.test(tab)
      test_name <- "Chi-square test"
    }

    tibble::tibble(
      group_1 = g[1],
      group_2 = g[2],
      threshold = threshold,
      test = test_name,
      p_value = test_res$p.value
    )
  }) %>%
    mutate(
      p_adjusted_bonferroni = p.adjust(p_value, method = "bonferroni")
    )

  return(out)
}

#===============================================================================
# Step 1: Run IsoRelate from scratch
#===============================================================================

cat("\n--- Step 1: Running IsoRelate from scratch ---\n")

ped_path <- file.path(input_dir, paste0(prefix, ".ped"))
map_path <- file.path(input_dir, paste0(prefix, ".map"))

if (!file.exists(ped_path)) {
  stop("PED file not found: ", ped_path)
}

if (!file.exists(map_path)) {
  stop("MAP file not found: ", map_path)
}

if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}

cat("Loading PED file:", ped_path, "\n")

ped <- read.table(
  ped_path,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

names(ped)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")

# Keep the convention from the previous IsoRelate script
ped$moi <- rep(2, nrow(ped))

cat("PED dimensions:", paste(dim(ped), collapse = " x "), "\n")

cat("Loading MAP file:", map_path, "\n")

map <- read.table(
  map_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

names(map)[1:4] <- c("chr", "snp_id", "pos_cM", "pos_bp")

map <- map %>%
  mutate(
    snp_id = as.character(snp_id),
    pos_bp = as.numeric(pos_bp),
    pos_cM = pos_bp / 13700
  )

cat("MAP dimensions:", paste(dim(map), collapse = " x "), "\n")

pedmap <- list(ped, map)

cat("Filtering genotypes with getGenotypes...\n")

geno <- getGenotypes(
  ped.map = pedmap,
  reference.ped.map = NULL,
  maf = 0.01,
  isolate.max.missing = 0.30,
  snp.max.missing = 0.6,
  input.map.distance = "cM",
  reference.map.distance = "cM"
)

saveRDS(
  geno,
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_geno.rds"))
)

cat("Estimating IBD parameters...\n")

param <- getIBDparameters(
  ped.genotypes = geno,
  number.cores = 3
)

saveRDS(
  param,
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_params.rds"))
)

cat("Inferring IBD segments...\n")

ibd_segments <- getIBDsegments(
  ped.genotypes = geno,
  parameters = param,
  number.cores = 3,
  minimum.snps = 450,
  minimum.length.bp = 700000,
  error = 0.01
)

saveRDS(
  ibd_segments,
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_ibd_segments.rds"))
)

readr::write_csv(
  as.data.frame(ibd_segments),
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_ibd_segments.csv"))
)

cat("Number of IBD segments detected:", nrow(ibd_segments), "\n")

cat("Creating IsoRelate summaries...\n")

ibd_summary <- getIBDsummary(
  ped.genotypes = geno,
  ibd.segments = ibd_segments
)

readr::write_csv(
  as.data.frame(ibd_summary),
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_ibd_summary.csv"))
)

ibd_matrix <- getIBDmatrix(
  ped.genotypes = geno,
  ibd.segments = ibd_segments
)

saveRDS(
  ibd_matrix,
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_ibd_matrix.rds"))
)

# This is SNP-level proportion of pairs IBD, not the pairwise sample table.
# It is saved for diagnostic plotting only.
ibd_prop_snp_level <- getIBDproportion(
  ped.genotypes = geno,
  ibd.matrix = ibd_matrix
)

readr::write_csv(
  as.data.frame(ibd_prop_snp_level),
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_ibd_proportion_snp_level.csv"))
)

pdf(
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_IsoRelate_IBD_proportions_snp_level.pdf")),
  width = 10,
  height = 6
)

plotIBDproportions(
  ibd_prop_snp_level,
  plot.title = "Proportion of pairs IBD across SNPs"
)

dev.off()

pdf(
  file.path(dirs$ibd_from_scratch, paste0(prefix, "_IsoRelate_IBD_segments.pdf")),
  width = 12,
  height = 8
)

plotIBDsegments(
  geno,
  ibd_segments,
  plot.title = "Distribution of IBD segments"
)

dev.off()

#===============================================================================
# Step 2: Clean metadata
#===============================================================================

cat("\n--- Step 2: Cleaning metadata ---\n")

if (grepl("\\.csv$", metadata_path, ignore.case = TRUE)) {
  metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)
} else {
  metadata <- readxl::read_excel(metadata_path)
}

cat("Metadata columns:\n")
print(names(metadata))

id_col <- if ("id" %in% names(metadata)) {
  "id"
} else if ("Sample" %in% names(metadata)) {
  "Sample"
} else if ("sample" %in% names(metadata)) {
  "sample"
} else {
  NA_character_
}

required_cols <- c("community", "travel")
missing_cols <- setdiff(required_cols, names(metadata))

if (is.na(id_col)) {
  stop("Metadata is missing a sample identifier column: id, Sample, or sample")
}

if (length(missing_cols) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

metadata_clean <- metadata %>%
  transmute(
    sample_id_raw = .data[[id_col]],
    sample_id = clean_id(.data[[id_col]]),
    community = as.character(community),
    travel = as.integer(travel)
  ) %>%
  mutate(
    community = case_when(
      str_to_lower(community) == "libertad" ~ "Libertad",
      str_to_lower(community) == "gamitanacocha" ~ "Gamitanacocha",
      str_to_lower(community) %in% c("urcomirano", "urco mirano") ~ "UrcoMirano",
      TRUE ~ community
    ),
    travel_status = case_when(
      travel == 1 ~ "Traveler",
      travel == 0 ~ "Non-traveler",
      TRUE ~ NA_character_
    )
  )

duplicated_ids <- metadata_clean %>%
  count(sample_id) %>%
  filter(n > 1)

readr::write_csv(
  duplicated_ids,
  file.path(dirs$metadata_clean, "duplicated_metadata_ids.csv")
)

if (nrow(duplicated_ids) > 0) {
  cat("Duplicated metadata IDs after cleaning:\n")
  print(duplicated_ids)
  stop("Duplicated sample IDs found in metadata after cleaning. See duplicated_metadata_ids.csv")
}

readr::write_csv(
  metadata_clean,
  file.path(dirs$metadata_clean, "metadata_clean.csv")
)

cat("Metadata sample count:", nrow(metadata_clean), "\n")
cat("Community counts:\n")
print(table(metadata_clean$community))
cat("Travel counts:\n")
print(table(metadata_clean$travel_status))

# Community centroids based on average house coordinates.
# In the house-level GPS file, gps_latitude contained longitude-like values
# around -73, and gps_longitude contained latitude-like values around -3.

community_centroids <- tibble::tibble(
  community = c("Libertad", "Gamitanacocha", "UrcoMirano"),
  longitude = c(-73.231521, -73.317568, -73.063149),
  latitude  = c(-3.491239, -3.426290, -3.359659)
)

readr::write_csv(
  community_centroids,
  file.path(dirs$metadata_clean, "community_centroids.csv")
)

#===============================================================================
# Step 3: Build pairwise IBD dataset from IBD matrix
#===============================================================================

cat("\n--- Step 3: Creating pairwise IBD dataset from IBD matrix ---\n")

pairwise_ibd <- make_pairwise_from_ibd_matrix(
  ibd_matrix = ibd_matrix,
  metadata_ids = metadata_clean$sample_id,
  output_dir = dirs$pairwise_data
)

readr::write_csv(
  pairwise_ibd,
  file.path(dirs$pairwise_data, paste0(prefix, "_pairwise_ibd_from_matrix.csv"))
)

matrix_samples <- sort(unique(c(pairwise_ibd$sample1, pairwise_ibd$sample2)))
metadata_samples <- sort(metadata_clean$sample_id)

samples_in_matrix_not_metadata <- setdiff(matrix_samples, metadata_samples)
samples_in_metadata_not_matrix <- setdiff(metadata_samples, matrix_samples)

readr::write_csv(
  tibble::tibble(sample_id = samples_in_matrix_not_metadata),
  file.path(dirs$pairwise_data, "samples_in_matrix_not_metadata.csv")
)

readr::write_csv(
  tibble::tibble(sample_id = samples_in_metadata_not_matrix),
  file.path(dirs$pairwise_data, "samples_in_metadata_not_matrix.csv")
)

cat("Samples in pairwise IBD matrix:", length(matrix_samples), "\n")
cat("Samples in metadata:", length(metadata_samples), "\n")
cat("Samples in matrix but not metadata:", length(samples_in_matrix_not_metadata), "\n")
cat("Samples in metadata but not matrix:", length(samples_in_metadata_not_matrix), "\n")

if (length(samples_in_matrix_not_metadata) > 0) {
  stop("Some samples in the IBD matrix are not present in metadata. See samples_in_matrix_not_metadata.csv")
}

metadata_1 <- metadata_clean %>%
  rename(
    sample1 = sample_id,
    sample1_raw_meta = sample_id_raw,
    community1 = community,
    travel1 = travel,
    travel_status1 = travel_status
  )

metadata_2 <- metadata_clean %>%
  rename(
    sample2 = sample_id,
    sample2_raw_meta = sample_id_raw,
    community2 = community,
    travel2 = travel,
    travel_status2 = travel_status
  )

pairwise_annotated <- pairwise_ibd %>%
  left_join(metadata_1, by = "sample1") %>%
  left_join(metadata_2, by = "sample2")

unmatched <- pairwise_annotated %>%
  filter(
    is.na(community1) |
      is.na(community2) |
      is.na(travel_status1) |
      is.na(travel_status2)
  )

readr::write_csv(
  unmatched,
  file.path(dirs$pairwise_data, paste0(prefix, "_unmatched_pairs_check.csv"))
)

if (nrow(unmatched) > 0) {
  cat("WARNING: Some pairwise IDs did not match metadata.\n")
  cat("See unmatched file in:", dirs$pairwise_data, "\n")
  stop("Stopping because some IBD pairs could not be matched to metadata.")
}

pairwise_annotated <- pairwise_annotated %>%
  rowwise() %>%
  mutate(
    community_a = sort(c(community1, community2))[1],
    community_b = sort(c(community1, community2))[2],
    community_pair = paste(community_a, community_b, sep = " - "),
    within_between = if_else(
      community1 == community2,
      "Within community",
      "Between communities"
    ),
    travel_pair = case_when(
      travel_status1 == "Non-traveler" & travel_status2 == "Non-traveler" ~
        "Non-traveler - Non-traveler",
      travel_status1 == "Traveler" & travel_status2 == "Traveler" ~
        "Traveler - Traveler",
      TRUE ~ "Traveler - Non-traveler"
    ),
    related_main = ibd_prop >= main_threshold
  ) %>%
  ungroup() %>%
  mutate(
    community_pair = factor(
      community_pair,
      levels = c(
        "Gamitanacocha - Gamitanacocha",
        "Gamitanacocha - Libertad",
        "Gamitanacocha - UrcoMirano",
        "Libertad - Libertad",
        "Libertad - UrcoMirano",
        "UrcoMirano - UrcoMirano"
      )
    ),
    within_between = factor(
      within_between,
      levels = c("Within community", "Between communities")
    ),
    travel_pair = factor(
      travel_pair,
      levels = c(
        "Non-traveler - Non-traveler",
        "Traveler - Non-traveler",
        "Traveler - Traveler"
      )
    )
  )

readr::write_csv(
  pairwise_annotated,
  file.path(dirs$pairwise_data, paste0(prefix, "_pairwise_ibd_annotated.csv"))
)

cat("Number of annotated pairwise comparisons:", nrow(pairwise_annotated), "\n")

expected_pairs <- length(matrix_samples) * (length(matrix_samples) - 1) / 2
cat("Expected number of pairs from matrix samples:", expected_pairs, "\n")

#===============================================================================
# Step 4: Community relatedness summaries
#===============================================================================

cat("\n--- Step 4: Community relatedness summaries ---\n")

community_summary <- summarise_related_pairs(
  data = pairwise_annotated,
  group_vars = c("community_pair"),
  thresholds = all_thresholds
)

readr::write_csv(
  community_summary,
  file.path(dirs$community_analysis, "community_pair_relatedness_summary_all_thresholds.csv")
)

community_tests <- purrr::map_dfr(all_thresholds, function(thr) {
  pairwise_prop_tests(pairwise_annotated, "community_pair", thr)
})

readr::write_csv(
  community_tests,
  file.path(dirs$community_analysis, "community_pair_relatedness_tests_all_thresholds.csv")
)

#===============================================================================
# Step 5: Within vs between community analysis
#===============================================================================

cat("\n--- Step 5: Within vs between community analysis ---\n")

within_between_summary <- summarise_related_pairs(
  data = pairwise_annotated,
  group_vars = c("within_between"),
  thresholds = all_thresholds
)

readr::write_csv(
  within_between_summary,
  file.path(dirs$within_between, "within_between_relatedness_summary_all_thresholds.csv")
)

within_between_tests <- purrr::map_dfr(all_thresholds, function(thr) {
  pairwise_prop_tests(pairwise_annotated, "within_between", thr)
})

readr::write_csv(
  within_between_tests,
  file.path(dirs$within_between, "within_between_relatedness_tests_all_thresholds.csv")
)

p_within_between <- ggplot(pairwise_annotated, aes(x = within_between, y = ibd_prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  geom_hline(yintercept = main_threshold, linetype = "dashed") +
  labs(
    title = "Pairwise IBD: within vs between communities",
    subtitle = "Dashed line = main relatedness threshold, IBD >= 0.50",
    x = NULL,
    y = "Pairwise IBD proportion"
  ) +
  theme_bw()

save_plot(
  p_within_between,
  "within_between_community_ibd_boxplot",
  dirs$within_between,
  width = 8,
  height = 6
)

#===============================================================================
# Step 6: Travel-status relatedness analysis
#===============================================================================

cat("\n--- Step 6: Travel-status relatedness analysis ---\n")

travel_summary <- summarise_related_pairs(
  data = pairwise_annotated,
  group_vars = c("travel_pair"),
  thresholds = all_thresholds
)

readr::write_csv(
  travel_summary,
  file.path(dirs$travel_analysis, "travel_pair_relatedness_summary_all_thresholds.csv")
)

travel_tests <- purrr::map_dfr(all_thresholds, function(thr) {
  pairwise_prop_tests(pairwise_annotated, "travel_pair", thr)
})

readr::write_csv(
  travel_tests,
  file.path(dirs$travel_analysis, "travel_pair_relatedness_tests_all_thresholds.csv")
)

p_travel <- ggplot(pairwise_annotated, aes(x = travel_pair, y = ibd_prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  geom_hline(yintercept = main_threshold, linetype = "dashed") +
  labs(
    title = "Pairwise IBD by travel-status comparison",
    subtitle = "Dashed line = main relatedness threshold, IBD >= 0.50",
    x = NULL,
    y = "Pairwise IBD proportion"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(
  p_travel,
  "travel_status_ibd_boxplot",
  dirs$travel_analysis,
  width = 9,
  height = 6
)

#===============================================================================
# Step 7: Community heatmaps
#===============================================================================

cat("\n--- Step 7: Creating community heatmaps ---\n")

community_matrix_summary <- pairwise_annotated %>%
  group_by(community_a, community_b) %>%
  summarise(
    n_pairs = n(),
    mean_ibd = mean(ibd_prop, na.rm = TRUE),
    median_ibd = median(ibd_prop, na.rm = TRUE),
    min_ibd = min(ibd_prop, na.rm = TRUE),
    max_ibd = max(ibd_prop, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  community_matrix_summary,
  file.path(dirs$heatmaps, "community_pair_mean_ibd_matrix_summary.csv")
)

heatmap_data_all <- purrr::map_dfr(all_thresholds, function(thr) {
  pairwise_annotated %>%
    mutate(
      threshold = thr,
      related = ibd_prop >= thr
    ) %>%
    group_by(community_a, community_b, threshold) %>%
    summarise(
      n_pairs = n(),
      n_related = sum(related, na.rm = TRUE),
      R = n_related / n_pairs,
      .groups = "drop"
    )
})

readr::write_csv(
  heatmap_data_all,
  file.path(dirs$heatmaps, "community_heatmap_R_all_thresholds.csv")
)

heatmap_symmetric <- bind_rows(
  heatmap_data_all,
  heatmap_data_all %>%
    filter(community_a != community_b) %>%
    transmute(
      community_a = .data$community_b,
      community_b = .data$community_a,
      threshold = threshold,
      n_pairs = n_pairs,
      n_related = n_related,
      R = R
    )
)

for (thr in all_thresholds) {
  p_heat <- heatmap_symmetric %>%
    filter(threshold == thr) %>%
    ggplot(aes(x = community_a, y = community_b, fill = R)) +
    geom_tile() +
    geom_text(
      aes(label = paste0(n_related, "/", n_pairs, "\nR=", round(R, 3))),
      size = 3
    ) +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "black") +
    labs(
      title = paste0("Community-by-community related-pair proportion, IBD >= ", thr),
      x = NULL,
      y = NULL,
      fill = "R"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  save_plot(
    p_heat,
    paste0("community_heatmap_R_threshold_", threshold_label(thr)),
    dirs$heatmaps,
    width = 8,
    height = 7
  )
}

#===============================================================================
# Step 8: IBD distribution plots
#===============================================================================

cat("\n--- Step 8: Creating IBD distribution plots ---\n")

p_dist_community <- ggplot(pairwise_annotated, aes(x = ibd_prop)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = main_threshold, linetype = "dashed") +
  facet_wrap(~ community_pair, scales = "free_y") +
  labs(
    title = "Pairwise IBD distributions by community comparison",
    subtitle = "Dashed line = main relatedness threshold, IBD >= 0.50",
    x = "Pairwise IBD proportion",
    y = "Number of pairs"
  ) +
  theme_bw()

save_plot(
  p_dist_community,
  "ibd_distribution_by_community_pair",
  dirs$distributions,
  width = 12,
  height = 8
)

p_dist_travel <- ggplot(pairwise_annotated, aes(x = ibd_prop)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = main_threshold, linetype = "dashed") +
  facet_wrap(~ travel_pair, scales = "free_y") +
  labs(
    title = "Pairwise IBD distributions by travel-status comparison",
    subtitle = "Dashed line = main relatedness threshold, IBD >= 0.50",
    x = "Pairwise IBD proportion",
    y = "Number of pairs"
  ) +
  theme_bw()

save_plot(
  p_dist_travel,
  "ibd_distribution_by_travel_pair",
  dirs$distributions,
  width = 11,
  height = 6
)

p_dist_within_between <- ggplot(pairwise_annotated, aes(x = ibd_prop)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = main_threshold, linetype = "dashed") +
  facet_wrap(~ within_between, scales = "free_y") +
  labs(
    title = "Pairwise IBD distributions: within vs between communities",
    subtitle = "Dashed line = main relatedness threshold, IBD >= 0.50",
    x = "Pairwise IBD proportion",
    y = "Number of pairs"
  ) +
  theme_bw()

save_plot(
  p_dist_within_between,
  "ibd_distribution_within_between_community",
  dirs$distributions,
  width = 10,
  height = 6
)

#===============================================================================
# Step 9: Sample-level pairwise IBD heatmaps
#===============================================================================

cat("\n--- Step 9: Creating sample-level IBD heatmaps ---\n")

samples_used <- sort(unique(c(pairwise_annotated$sample1, pairwise_annotated$sample2)))

ibd_mat_sample <- matrix(
  0,
  nrow = length(samples_used),
  ncol = length(samples_used),
  dimnames = list(samples_used, samples_used)
)

diag(ibd_mat_sample) <- 1

for (i in seq_len(nrow(pairwise_annotated))) {
  s1 <- pairwise_annotated$sample1[i]
  s2 <- pairwise_annotated$sample2[i]
  val <- pairwise_annotated$ibd_prop[i]

  ibd_mat_sample[s1, s2] <- val
  ibd_mat_sample[s2, s1] <- val
}

saveRDS(
  ibd_mat_sample,
  file.path(dirs$heatmaps, "sample_level_pairwise_ibd_matrix.rds")
)

readr::write_csv(
  as.data.frame(ibd_mat_sample) %>%
    tibble::rownames_to_column("sample_id"),
  file.path(dirs$heatmaps, "sample_level_pairwise_ibd_matrix.csv")
)

sample_order_metadata <- metadata_clean %>%
  filter(sample_id %in% samples_used) %>%
  arrange(community, travel_status, sample_id) %>%
  pull(sample_id)

sample_heatmap_metadata <- pairwise_annotated %>%
  select(sample1, sample2, ibd_prop) %>%
  bind_rows(
    pairwise_annotated %>%
      transmute(sample1 = sample2, sample2 = sample1, ibd_prop = ibd_prop)
  ) %>%
  bind_rows(
    tibble::tibble(
      sample1 = sample_order_metadata,
      sample2 = sample_order_metadata,
      ibd_prop = 1
    )
  ) %>%
  mutate(
    sample1 = factor(sample1, levels = sample_order_metadata),
    sample2 = factor(sample2, levels = sample_order_metadata)
  )

p_sample_heatmap_metadata <- ggplot(sample_heatmap_metadata, aes(x = sample1, y = sample2, fill = ibd_prop)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "black") +
  labs(
    title = "Sample-level pairwise IBD heatmap",
    subtitle = "Samples ordered by community and travel status",
    x = NULL,
    y = NULL,
    fill = "IBD"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

save_plot(
  p_sample_heatmap_metadata,
  "sample_level_pairwise_ibd_heatmap_metadata_ordered",
  dirs$heatmaps,
  width = 9,
  height = 8
)

hc <- hclust(as.dist(1 - ibd_mat_sample), method = "average")
sample_order_clustered <- hc$labels[hc$order]

sample_heatmap_clustered <- pairwise_annotated %>%
  select(sample1, sample2, ibd_prop) %>%
  bind_rows(
    pairwise_annotated %>%
      transmute(sample1 = sample2, sample2 = sample1, ibd_prop = ibd_prop)
  ) %>%
  bind_rows(
    tibble::tibble(
      sample1 = sample_order_clustered,
      sample2 = sample_order_clustered,
      ibd_prop = 1
    )
  ) %>%
  mutate(
    sample1 = factor(sample1, levels = sample_order_clustered),
    sample2 = factor(sample2, levels = sample_order_clustered)
  )

p_sample_heatmap_clustered <- ggplot(sample_heatmap_clustered, aes(x = sample1, y = sample2, fill = ibd_prop)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "black") +
  labs(
    title = "Sample-level pairwise IBD heatmap",
    subtitle = "Samples clustered using hierarchical clustering on 1 - IBD",
    x = NULL,
    y = NULL,
    fill = "IBD"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

save_plot(
  p_sample_heatmap_clustered,
  "sample_level_pairwise_ibd_heatmap_clustered",
  dirs$heatmaps,
  width = 9,
  height = 8
)

#===============================================================================
# Step 10: Map-style community connectivity
#===============================================================================

cat("\n--- Step 10: Creating map-style community connectivity plot ---\n")

map_edges <- heatmap_data_all %>%
  filter(threshold == main_threshold) %>%
  filter(community_a != community_b) %>%
  left_join(community_centroids, by = c("community_a" = "community")) %>%
  rename(lon_a = longitude, lat_a = latitude) %>%
  left_join(community_centroids, by = c("community_b" = "community")) %>%
  rename(lon_b = longitude, lat_b = latitude)

readr::write_csv(
  map_edges,
  file.path(dirs$map_connectivity, "community_connectivity_edges_threshold_0_50.csv")
)

p_map <- ggplot() +
  geom_segment(
    data = map_edges,
    aes(
      x = lon_a,
      y = lat_a,
      xend = lon_b,
      yend = lat_b,
      linewidth = R,
      alpha = R
    )
  ) +
  geom_point(
    data = community_centroids,
    aes(x = longitude, y = latitude),
    size = 4
  ) +
  geom_text(
    data = community_centroids,
    aes(x = longitude, y = latitude, label = community),
    nudge_y = 0.01,
    size = 4
  ) +
  scale_linewidth_continuous(range = c(0.2, 3)) +
  coord_equal() +
  labs(
    title = "Community genetic connectivity based on IBD",
    subtitle = "Edges weighted by related-pair proportion at IBD >= 0.50",
    x = "Longitude",
    y = "Latitude",
    linewidth = "R",
    alpha = "R"
  ) +
  theme_bw()

save_plot(
  p_map,
  "community_map_connectivity_threshold_0_50",
  dirs$map_connectivity,
  width = 8,
  height = 7
)

#===============================================================================
# Step 11: Sensitivity analysis
#===============================================================================

cat("\n--- Step 11: Sensitivity analysis ---\n")

sensitivity_summary_all <- bind_rows(
  community_summary %>%
    mutate(
      analysis = "Community pair",
      group = as.character(community_pair)
    ) %>%
    select(analysis, group, threshold, n_pairs, n_related, R, mean_ibd, median_ibd),

  within_between_summary %>%
    mutate(
      analysis = "Within/between community",
      group = as.character(within_between)
    ) %>%
    select(analysis, group, threshold, n_pairs, n_related, R, mean_ibd, median_ibd),

  travel_summary %>%
    mutate(
      analysis = "Travel-status pair",
      group = as.character(travel_pair)
    ) %>%
    select(analysis, group, threshold, n_pairs, n_related, R, mean_ibd, median_ibd)
)

readr::write_csv(
  sensitivity_summary_all,
  file.path(dirs$sensitivity, "sensitivity_summary_all_analyses.csv")
)

p_sensitivity <- sensitivity_summary_all %>%
  ggplot(aes(x = threshold, y = R, group = group)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ analysis, scales = "free_y") +
  labs(
    title = "Sensitivity analysis across IBD thresholds",
    subtitle = "Thresholds: 0.05, 0.50, and 0.80",
    x = "IBD threshold",
    y = "Related-pair proportion, R"
  ) +
  theme_bw()

save_plot(
  p_sensitivity,
  "sensitivity_R_across_thresholds",
  dirs$sensitivity,
  width = 11,
  height = 7
)

#===============================================================================
# Step 12: IBD networks at selected thresholds
#===============================================================================

cat("\n--- Step 12: Creating IBD networks at selected thresholds ---\n")

for (thr in all_thresholds) {
  cat("Creating network at threshold:", thr, "\n")

  edge_df <- pairwise_annotated %>%
    filter(ibd_prop >= thr) %>%
    transmute(
      from = sample1,
      to = sample2,
      ibd_prop = ibd_prop,
      community_pair = as.character(community_pair),
      travel_pair = as.character(travel_pair)
    )

  node_df <- metadata_clean %>%
    filter(sample_id %in% samples_used) %>%
    transmute(
      name = sample_id,
      community = community,
      travel_status = travel_status
    )

  g <- igraph::graph_from_data_frame(
    d = edge_df,
    vertices = node_df,
    directed = FALSE
  )

  saveRDS(
    g,
    file.path(
      dirs$network,
      paste0(prefix, "_network_threshold_", threshold_label(thr), ".rds")
    )
  )

  network_summary <- tibble::tibble(
    threshold = thr,
    n_nodes = igraph::vcount(g),
    n_edges = igraph::ecount(g),
    n_components = igraph::components(g)$no,
    n_isolates = sum(igraph::degree(g) == 0)
  )

  readr::write_csv(
    network_summary,
    file.path(
      dirs$network,
      paste0(prefix, "_network_summary_threshold_", threshold_label(thr), ".csv")
    )
  )

  pdf(
    file.path(
      dirs$network,
      paste0(prefix, "_network_threshold_", threshold_label(thr), ".pdf")
    ),
    width = 10,
    height = 10
  )

  plot(
    g,
    vertex.size = 5,
    vertex.label = NA,
    vertex.color = as.numeric(as.factor(igraph::V(g)$community)),
    edge.width = ifelse(igraph::ecount(g) > 0, igraph::E(g)$ibd_prop * 3, 1),
    main = paste0("IBD network, threshold >= ", thr)
  )

  dev.off()
}

#===============================================================================
# Finish
#===============================================================================

cat("\nAnalysis finished:", as.character(Sys.time()), "\n")
cat("All outputs saved in:", root_outdir, "\n")