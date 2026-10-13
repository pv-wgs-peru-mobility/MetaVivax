#!/usr/bin/env bash

#SBATCH --job-name=pv-fst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mahdi.safarpour@uantwerpen.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

# Submit from the repository root so that the relative input paths resolve.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

Rscript population_genetics/fst/FST.R \
  --gt_table data/analysis_ready/pv_56samples_16485snps_gt.tsv \
  --metadata data/analysis_ready/pv_56samples_metadata.tsv \
  --out_dir results/fst \
  --n_boot 1000 \
  --n_perm 1000 \
  --seed 1

Rscript population_genetics/fst/plot_fst.R --out_dir results/fst
