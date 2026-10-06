#!/usr/bin/env Rscript

#===============================================================================
# IBD analysis with isoRelate
# Author: Mahdi Safarpour
# Purpose: Estimate identity-by-descent sharing from PLINK PED/MAP data
#===============================================================================

# --------------------------- Load libraries -----------------------------------
suppressPackageStartupMessages({
  library(isoRelate)
  library(igraph)
})

# --------------------------- Define paths and settings ------------------------
outdir <- "/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77"
prefix <- "77.19298"
number_cores <- 3
cluster_prop <- 0.40

ped_path <- file.path(outdir, paste0(prefix, ".ped"))
map_path <- file.path(outdir, paste0(prefix, ".map"))

# --------------------------- Load PED and MAP ---------------------------------
cat("Loading PED and MAP files...\n")

if (!file.exists(ped_path)) stop("PED file not found: ", ped_path)
if (!file.exists(map_path)) stop("MAP file not found: ", map_path)

ped <- read.table(ped_path, header = FALSE, sep = "", stringsAsFactors = FALSE)
names(ped)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")
ped$moi <- 2L
cat("PED dimensions:", paste(dim(ped), collapse = " x "), "\n")

map <- read.table(map_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(map)[1:4] <- c("chr", "snp_id", "pos_cM", "pos_bp")
map$snp_id <- as.character(map$snp_id)
map$pos_bp <- as.numeric(map$pos_bp)
map$pos_cM <- map$pos_bp / 13700
cat("MAP dimensions:", paste(dim(map), collapse = " x "), "\n")

pedmap <- list(ped, map)

# --------------------------- Genotype filtering -------------------------------
cat("Filtering genotypes...\n")

geno <- getGenotypes(
  ped.map = pedmap,
  reference.ped.map = NULL,
  maf = 0.01,
  isolate.max.missing = 0.30,
  snp.max.missing = 0.60,
  input.map.distance = "cM",
  reference.map.distance = "cM"
)

saveRDS(geno, file.path(outdir, paste0(prefix, "_geno.rds")))
cat("Filtered genotypes saved.\n")

# --------------------------- Estimate IBD parameters --------------------------
cat("Estimating IBD parameters...\n")

param <- getIBDparameters(
  ped.genotypes = geno,
  number.cores = number_cores
)

saveRDS(param, file.path(outdir, paste0(prefix, "_params.rds")))
cat("IBD parameters estimated and saved.\n")

# --------------------------- Infer IBD segments -------------------------------
cat("Inferring IBD segments...\n")

ibd <- getIBDsegments(
  ped.genotypes = geno,
  parameters = param,
  number.cores = number_cores,
  minimum.snps = 450,
  minimum.length.bp = 700000,
  error = 0.01
)

saveRDS(ibd, file.path(outdir, paste0(prefix, "_ibd_segments.rds")))
write.csv(ibd, file.path(outdir, paste0(prefix, "_ibd_segments.csv")), row.names = FALSE)
cat("IBD segments detected:", nrow(ibd), "\n")

# --------------------------- Summaries ----------------------------------------
cat("Creating IBD summaries...\n")

ibd_summary <- capture.output(getIBDsummary(ped.genotypes = geno, ibd.segments = ibd))
writeLines(ibd_summary, file.path(outdir, paste0(prefix, "_ibd_summary.txt")))

ibd_matrix <- getIBDmatrix(geno, ibd)
saveRDS(ibd_matrix, file.path(outdir, paste0(prefix, "_ibd_matrix.rds")))

ibd_prop <- getIBDproportion(geno, ibd_matrix)
write.csv(ibd_prop, file.path(outdir, paste0(prefix, "_ibd_proportion.csv")), row.names = FALSE)

# --------------------------- Plots --------------------------------------------
cat("Creating IBD plots...\n")

pdf(file.path(outdir, paste0(prefix, "_IBD_proportions.pdf")), width = 10, height = 6)
plotIBDproportions(ibd_prop, plot.title = "Proportion of pairs IBD")
dev.off()

pdf(file.path(outdir, paste0(prefix, "_IBD_segments.pdf")), width = 12, height = 8)
plotIBDsegments(geno, ibd, plot.title = "Distribution of IBD segments")
dev.off()

# --------------------------- IBD clustering -----------------------------------
cat("Performing IBD clustering at prop =", cluster_prop, "...\n")

clusters <- getIBDpclusters(
  ped.genotypes = geno,
  ibd.segments = ibd,
  prop = cluster_prop,
  hi.clust = FALSE
)

saveRDS(clusters, file.path(outdir, paste0(prefix, "_clusters_prop0.40.rds")))

pdf(file.path(outdir, paste0(prefix, "_IBD_clusters_network.pdf")), width = 14, height = 10)
plot(
  clusters$i.network,
  vertex.size = 4,
  vertex.label = NA,
  main = "IBD clusters (prop = 0.40)"
)
dev.off()

# --------------------------- Full network -------------------------------------
cat("Building full IBD network including isolated samples...\n")

net <- clusters$i.network
V(net)$name <- as.character(V(net)$name)

ids_iid <- unique(as.character(geno$pedigree$iid))
ids_fid_iid <- unique(paste(geno$pedigree$fid, geno$pedigree$iid, sep = "_"))

if (sum(ids_iid %in% V(net)$name) >= sum(ids_fid_iid %in% V(net)$name)) {
  all_samples <- ids_iid
} else {
  all_samples <- ids_fid_iid
}

missing_samples <- setdiff(all_samples, V(net)$name)
cat("Samples without IBD links added as isolated nodes:", length(missing_samples), "\n")

if (length(missing_samples) > 0) {
  net <- add_vertices(net, nv = length(missing_samples), name = missing_samples)
}

deg <- degree(net)
V(net)$color <- ifelse(deg == 0, "grey80", "tomato")
layout_full <- layout_with_fr(net)

saveRDS(net, file.path(outdir, paste0(prefix, "_IBD_full_network.rds")))

pdf(file.path(outdir, paste0(prefix, "_IBD_network_full_colored.pdf")), width = 10, height = 10)
plot(
  net,
  layout = layout_full,
  vertex.size = 3,
  vertex.label = NA,
  vertex.frame.color = NA,
  edge.color = "grey80",
  main = "Full IBD network"
)
dev.off()

cat("Done.\n")
