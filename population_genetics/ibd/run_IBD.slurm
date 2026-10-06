#!/usr/bin/env bash
#SBATCH --job-name=pv-ibd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mahdi.safarpour@uantwerpen.be
#SBATCH --output=%j-%x-stdout.out
#SBATCH --error=%j-%x-stderr.out

set -euo pipefail

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

ibd_script="/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/IBD.R"
if [[ ! -f "$ibd_script" ]]; then
  echo "IBD script not found: $ibd_script" >&2
  exit 1
fi

export IBD_CORES="${SLURM_CPUS_PER_TASK:-3}"
Rscript "$ibd_script"
