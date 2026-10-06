#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=admixture
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mahdi.safarpour@uantwerpen.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

module load ADMIXTURE/1.3.0-x86_64
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

BASE=/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv
OUT_DIR=$BASE/ADMIXTURE
PREFIX=77.19298
META=$BASE/IBD/SNP_matrix_with_metadata.csv

K_MIN=2
K_MAX=10
N_REPS=5
BOOTSTRAPS=100
THREADS=${SLURM_CPUS_PER_TASK:-8}

cd "$OUT_DIR"

need_bfile() {
  local stem="$1"
  if [[ ! -f "${stem}.bed" || ! -f "${stem}.bim" || ! -f "${stem}.fam" ]]; then
    echo "ERROR: Missing ${stem}.{bed,bim,fam} in $OUT_DIR"
    exit 1
  fi
}

make_chr_fixed() {
  local stem="$1"
  awk 'BEGIN{OFS="\t"}{$1=1; print}' "${stem}.bim" > "${stem}.chr1.bim"
  cp -f "${stem}.bed" "${stem}.chr1.bed"
  cp -f "${stem}.fam" "${stem}.chr1.fam"
}

run_admixture_series() {
  local stem="$1"
  local cv_table="${stem}.CV_errors.tsv"

  echo -e "K\trep\tseed\tCV" > "$cv_table"

  for K in $(seq "$K_MIN" "$K_MAX"); do
    for rep in $(seq 1 "$N_REPS"); do
      seed=$((1000 + K * 100 + rep))
      log="${stem}.K${K}.rep${rep}.log"

      echo "Running $stem K=$K replicate=$rep seed=$seed"
      admixture --cv -B${BOOTSTRAPS} -j"$THREADS" -s "$seed" "${stem}.bed" "$K" | tee "$log"

      if [[ -f "${stem}.${K}.Q" ]]; then
        mv -f "${stem}.${K}.Q" "${stem}.K${K}.rep${rep}.Q"
      fi
      if [[ -f "${stem}.${K}.P" ]]; then
        mv -f "${stem}.${K}.P" "${stem}.K${K}.rep${rep}.P"
      fi

      cv=$(awk '/CV error/ {print $NF}' "$log" | tail -n 1)
      if [[ -z "$cv" ]]; then
        echo "ERROR: Could not parse CV error from $log"
        exit 1
      fi
      echo -e "${K}\t${rep}\t${seed}\t${cv}" >> "$cv_table"
    done

    best_rep=$(awk -v k="$K" 'BEGIN{best=1; bestcv=1e99} NR>1 && $1==k {if ($4 < bestcv) {bestcv=$4; best=$2}} END{print best}' "$cv_table")
    cp -f "${stem}.K${K}.rep${best_rep}.Q" "${stem}.${K}.Q"
    cp -f "${stem}.K${K}.rep${best_rep}.P" "${stem}.${K}.P"
  done
}

for stem in "${PREFIX}.unpruned" "${PREFIX}.pruned_r2_0.1" "${PREFIX}.pruned_r2_0.2"; do
  need_bfile "$stem"
  make_chr_fixed "$stem"
done

run_admixture_series "${PREFIX}.unpruned.chr1"
run_admixture_series "${PREFIX}.pruned_r2_0.1.chr1"
run_admixture_series "${PREFIX}.pruned_r2_0.2.chr1"

cat > plot_admixture_all.R <<'RSCRIPT'
suppressPackageStartupMessages({
  library(tidyverse)
  library(Cairo)
})

prefix <- Sys.getenv("PREFIX")
dataset <- Sys.getenv("DATASET")
meta_file <- Sys.getenv("META")

stem <- paste0(prefix, ".", dataset)
cvfile <- paste0(stem, ".CV_errors.tsv")
famfile <- paste0(stem, ".fam")

if (!file.exists(cvfile)) stop("Missing CV file: ", cvfile)
if (!file.exists(famfile)) stop("Missing FAM file: ", famfile)
if (!file.exists(meta_file)) stop("Missing metadata file: ", meta_file)

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
id_col <- if ("id" %in% names(meta)) "id" else if ("Sample" %in% names(meta)) "Sample" else if ("sample" %in% names(meta)) "sample" else NA_character_
if (is.na(id_col)) stop("Metadata must contain id, Sample, or sample column")
if (!all(c("community", "travel") %in% names(meta))) stop("Metadata must contain community and travel columns")

meta <- meta %>%
  mutate(
    id_clean = gsub("\\.TGT$", "", .data[[id_col]]),
    travel_text = tolower(as.character(travel)),
    travel = case_when(
      travel_text %in% c("1", "yes", "y", "true", "traveler", "traveller") ~ "Traveler",
      travel_text %in% c("0", "no", "n", "false", "non-traveler", "non-traveller", "nontraveler", "nontraveller") ~ "Non-traveler",
      TRUE ~ "Unknown"
    ),
    community = as.character(community),
    community = if_else(community == "UrcoMirano", "Urco Mirano", community)
  )

cv_raw <- readr::read_tsv(cvfile, show_col_types = FALSE)
cv_df <- cv_raw %>%
  group_by(K) %>%
  summarise(CV = mean(CV), CV_sd = sd(CV), .groups = "drop") %>%
  arrange(K)

Kopt <- cv_df$K[which.min(cv_df$CV)]
write.csv(cv_df, paste0(stem, ".CV_mean_by_K.csv"), row.names = FALSE)

p_cv <- ggplot(cv_df, aes(K, CV)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = CV - CV_sd, ymax = CV + CV_sd), width = 0.15, na.rm = TRUE) +
  theme_bw() +
  ggtitle(paste0("CV curve: ", dataset))

CairoJPEG(paste0(stem, ".CV_curve.jpeg"), width = 2000, height = 1600, res = 300)
print(p_cv)
dev.off()

manuscriptK <- if (dataset == "unpruned.chr1") 4 else 5
Kplot <- sort(unique(c(3, 4, 5, Kopt, manuscriptK)))
Kplot <- Kplot[Kplot >= 2 & Kplot <= 10]

fam <- read.table(famfile, header = FALSE, stringsAsFactors = FALSE)
samples <- tibble(Sample = fam$V2) %>%
  mutate(Sample_clean = gsub("\\.TGT$", "", Sample)) %>%
  left_join(meta %>% transmute(Sample_clean = id_clean, community, travel), by = "Sample_clean") %>%
  mutate(
    community = replace_na(community, "Unknown"),
    travel = replace_na(travel, "Unknown"),
    Sample_label = if_else(travel == "Traveler", paste0(Sample_clean, "*"), Sample_clean)
  )

samples$community <- factor(samples$community, levels = sort(unique(samples$community)))

get_cluster_cols <- function(k) {
  base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
  cols <- if (k > length(base)) c(base, scales::hue_pal()(k - length(base))) else base[1:k]
  setNames(cols, paste0("Cluster", seq_len(k)))
}

for (k in Kplot) {
  qfile <- paste0(stem, ".", k, ".Q")
  if (!file.exists(qfile)) next

  Q <- read.table(qfile, header = FALSE)
  colnames(Q) <- paste0("Cluster", seq_len(k))
  Q$Sample <- fam$V2

  Q_long <- as_tibble(Q) %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry") %>%
    left_join(samples, by = "Sample")

  Q_long$Sample_label <- factor(Q_long$Sample_label, levels = unique(samples$Sample_label))

  p <- ggplot(Q_long, aes(x = Sample_label, y = Ancestry, fill = Cluster)) +
    geom_col(width = 1) +
    facet_grid(~community, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = get_cluster_cols(k)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    ) +
    ylab("Ancestry proportion") +
    xlab("Samples (* traveler-associated infection)") +
    ggtitle(paste0("ADMIXTURE (", dataset, "), K=", k))

  outname <- paste0(stem, ".K", k, ".barplot.community_split.jpeg")
  CairoJPEG(outname, width = 3600, height = 1400, res = 300)
  print(p)
  dev.off()
}
RSCRIPT

PREFIX="$PREFIX" DATASET="unpruned.chr1" META="$META" Rscript plot_admixture_all.R
PREFIX="$PREFIX" DATASET="pruned_r2_0.1.chr1" META="$META" Rscript plot_admixture_all.R
PREFIX="$PREFIX" DATASET="pruned_r2_0.2.chr1" META="$META" Rscript plot_admixture_all.R

echo "ADMIXTURE analysis complete."
