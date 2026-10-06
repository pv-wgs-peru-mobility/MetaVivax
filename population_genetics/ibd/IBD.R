#!/usr/bin/env Rscript

# Identity-by-descent (IBD) segment inference and threshold network analysis.
# Run with run_IBD.slurm. Paths and selected settings can be overridden with
# the IBD_* environment variables documented in README.md.

suppressPackageStartupMessages({
  library(isoRelate)
  library(igraph)
})

options(stringsAsFactors = FALSE)

ibd_root <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD"
input_dir <- Sys.getenv("IBD_INPUT_DIR", file.path(ibd_root, "New", "77"))
metadata_file <- Sys.getenv("IBD_METADATA", file.path(ibd_root, "SNP_matrix_with_metadata.csv"))
out_dir <- Sys.getenv("IBD_OUT_DIR", file.path(ibd_root, "results"))
prefix <- Sys.getenv("IBD_PREFIX", "77.19298")
cores <- as.integer(Sys.getenv("IBD_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "3")))
expected_samples <- as.integer(Sys.getenv("IBD_EXPECTED_SAMPLES", "56"))
isolate_max_missing <- as.numeric(Sys.getenv("IBD_ISOLATE_MAX_MISSING", "0.4"))
snp_max_missing <- 0.6
maf <- 0.01
minimum_snps <- 450L
minimum_length_bp <- 700000L
genotyping_error <- 0.001
bp_per_cm <- as.numeric(Sys.getenv("IBD_BP_PER_CM", "13700"))
map_mode <- Sys.getenv("IBD_MAP_MODE", "from_bp")
moi_mode <- Sys.getenv("IBD_MOI_MODE", "all_two")
thresholds_pct <- c(5L, 20L, 40L, 60L, 80L)

if (is.na(cores) || cores < 1L || is.na(expected_samples) || expected_samples < 0L ||
    is.na(isolate_max_missing) || isolate_max_missing < 0 || isolate_max_missing > 1 ||
    is.na(bp_per_cm) || bp_per_cm <= 0) {
  stop("Invalid IBD_CORES, IBD_EXPECTED_SAMPLES, IBD_ISOLATE_MAX_MISSING, or IBD_BP_PER_CM.")
}
if (!map_mode %in% c("from_bp", "input")) stop("IBD_MAP_MODE must be 'from_bp' or 'input'.")
if (!moi_mode %in% c("all_two", "ped")) stop("IBD_MOI_MODE must be 'all_two' or 'ped'.")

ped_path <- file.path(input_dir, paste0(prefix, ".ped"))
map_path <- file.path(input_dir, paste0(prefix, ".map"))
for (path in c(ped_path, map_path, metadata_file)) {
  if (!file.exists(path)) stop("Missing input file: ", path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
save_pdf <- function(path, expr, width = 10, height = 8) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}
clean_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.(TGT|GT)$", "", x, ignore.case = TRUE)
}
clean_travel <- function(x) {
  value <- tolower(trimws(as.character(x)))
  yes <- value %in% c("1", "yes", "y", "true", "traveler", "traveller")
  no <- value %in% c("0", "no", "n", "false", "non-traveler", "non-traveller",
                   "nontraveler", "nontraveller")
  if (anyNA(value) || any(!(yes | no))) {
    stop("Every retained sample needs a valid travel value (0/1, no/yes, or non-traveler/traveler).")
  }
  as.integer(yes)
}
finite_or_na <- function(value) if (is.finite(value)) value else NA_real_

message("Reading PED, MAP, and metadata...")
ped <- read.table(ped_path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
map <- read.table(map_path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(ped) < 8L || ncol(map) != 4L || ncol(ped) != 6L + 2L * nrow(map)) {
  stop("PED must have six pedigree fields and two allele columns per SNP in the four-column MAP.")
}
names(ped)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")
names(map) <- c("chr", "snp_id", "pos_cM", "pos_bp")
ped$fid <- as.character(ped$fid)
ped$iid <- as.character(ped$iid)
if (anyNA(ped$fid) || anyNA(ped$iid) || any(!nzchar(ped$fid)) ||
    any(!nzchar(ped$iid)) || anyDuplicated(paste(ped$fid, ped$iid, sep = "/"))) {
  stop("PED family/isolate combinations must be present and unique.")
}
if (moi_mode == "all_two") {
  ped$moi <- rep(2L, nrow(ped))
  message("MOI mode all_two: assigning 2 to all PED isolates, as in the supplied scripts.")
} else {
  ped$moi <- suppressWarnings(as.integer(as.character(ped$moi)))
  if (anyNA(ped$moi) || any(!ped$moi %in% c(1L, 2L))) {
    stop("With IBD_MOI_MODE=ped, PED column 5 must contain only 1 or 2.")
  }
}

map$pos_bp <- suppressWarnings(as.numeric(as.character(map$pos_bp)))
if (anyNA(map$pos_bp) || any(!is.finite(map$pos_bp)) || any(map$pos_bp <= 0)) {
