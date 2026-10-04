#!/usr/bin/env bash

#SBATCH --job-name=wgs-variant-calling-malformedvcf  # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks=1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=28          # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=72:00:00             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL  # send email when job begins, ends or fails
#SBATCH --mail-user=pmoris@itg.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

printf "Script execution started at $(date).\n"

# define required software

# if using modules
# module load ...

# if using conda:
# source /data/antwerpen/203/vsc20380/miniforge3/etc/profile.d/conda.sh
# conda activate wgs

# if using hpc-container-wrapper
export PATH="/scratch/antwerpen/203/vsc20380/containers/wgs/wgs-container-lock/bin:$PATH"

module load parallel

# Get file path of script and set project root.
# Otherwise, you would need to make sure to call the script from within the
# directory where it is stored.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")
else
    SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
fi

echo "Script dir = ${SCRIPT_DIR}"
PROJECT_ROOT=$(realpath "${SCRIPT_DIR}/../")
echo "Project root = ${PROJECT_ROOT}"

# NOTE: a samplesheet must be provided in the form of a csv file containing the sample name (part of the fastq filename) and the species (pf,pv,pm,po,pk)
# e.g.
# sample,species
# FCHTVJYCCXY_L1_WHRDMALtbmRAAAAA-23_1.fq.gz,pv
# FCHTVJYCCXY_L1_WHRDMALtbmRAABAA-33_1.fq.gz,pv

# NOTE: the -r1 -r2 read pair suffix, the -e fastq file extension and the -fastq_identifier file name format must be set correctly. Check the alignment script for details on the different options.

# qc
bash "${PROJECT_ROOT}/scripts/qc.sh" \
    --input_dir "${PROJECT_ROOT}/data/fastq/" \
    --output_dir "${PROJECT_ROOT}/results/" \
    -r1 "_1" -r2 "_2" -e ".fastq.gz"

# alignment
# NOTE: remove --extension flag, because qc.sh automatically sets it to .fastq.gz
# NOTE: --competitive enables competitive mapping, omiting it first maps against human and then maps all unmapped reads to plasmodium
bash "${PROJECT_ROOT}/scripts/alignment.sh" \
    --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
    --output_dir "${PROJECT_ROOT}/results/" \
    -r1 "_1" -r2 "_R2" \
    --fastq_identifier name_first

# variant calling
bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
    --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
    --output_dir "${PROJECT_ROOT}/results/"

# annotation
bash "${PROJECT_ROOT}/scripts/annotation-combine.sh" \
    --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
    --output_dir "${PROJECT_ROOT}/results/"

printf "Script finished at $(date).\n"
