#!/usr/bin/env bash

# set bash strict mode
set -euo pipefail

# allow debug mode by running `TRACE=1 ./script.sh` - equivalent to `set -x`
if [[ "${TRACE-0}" == "1" ]]; then set -o xtrace; fi

# get file path of project root to allow it to be run from any working directory
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
PROJECT_ROOT=$(realpath "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)/")
echo "Project root = ${PROJECT_ROOT}"

declare -a genome_array=( 
    "${PROJECT_ROOT}/Pfalciparum/3D7/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
    "${PROJECT_ROOT}/Pvivax/PvPAM/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxPAM_Genome.fasta"
    "${PROJECT_ROOT}/Pvivax/P01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PvivaxP01_Genome.fasta"
    "${PROJECT_ROOT}/Povale/curtisi/PocGH01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PovalecurtisiGH01_Genome.fasta"
    "${PROJECT_ROOT}/Povale/wallikeri/PowCR01/PlasmoDB-release-68/PlasmoDB-68_PovalewallikeriPowCR01_Genome.fasta"
    "${PROJECT_ROOT}/Pmalariae/UG01/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PmalariaeUG01_Genome.fasta"
    "${PROJECT_ROOT}/Pknowlesi/H/PlasmoDB/PlasmoDB-release-68/PlasmoDB-68_PknowlesiH_Genome.fasta"
    # "${PROJECT_ROOT}/human/GRCh38.p14/GENCODE/47/GRCh38.primary_assembly.genome.fa.gz"
    # "${PROJECT_ROOT}/human/GRCh38.p14/RefSeq/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
)

for genome in "${genome_array[@]}"; do
    echo "Creating snpEff database for ${genome}..."

    # fetch genome
    genome_dir="${PROJECT_ROOT}/snpEff/snpEff_database/$(basename ${genome})"
    genome_dir="${genome_dir%_Genome.fasta}"
    mkdir -p "${genome_dir}"

    # if [[ "${genome}" = *.gz ]]; then
    if file --mime-type "${genome}" | grep -q gzip$; then
        ln -sf "${genome}" "${genome_dir}/sequences.fa.gz"
    else
        ln -sf "${genome}" "${genome_dir}/sequences.fa"
    fi

    # fetch annotation
    annotation="$(realpath $(dirname ${genome})/*.gff*)"
    if [ -f "${annotation}" ]; then
        type="gff"
        if file --mime-type "${annotation}" | grep -q gzip$; then
            ln -sf "${annotation}" "${genome_dir}/genes.gff.gz"
        else
            ln -sf "${annotation}" "${genome_dir}/genes.gff"
        fi
    else
        type="gtf"
        annotation="$(realpath $(dirname ${genome})/*.gtf*)"
        if file --mime-type "${annotation}" | grep -q gzip$; then
            ln -sf "${annotation}" "${genome_dir}/genes.gtf.gz"
        else
            ln -sf "${annotation}" "${genome_dir}/genes.gtf"
        fi
    fi
    # annotation_gtf="$(basename ${genome})/*.gtf*"

    # fetch CDS
    cds="${genome%_Genome.fasta}_AnnotatedCDSs.fasta"
    ln -sf ${cds} "${genome_dir}/cds.fa"

    # build database
    if [ ! "${type}" == "gtf" ]; then
        snpEff build -c "${PROJECT_ROOT}/snpEff/snpEff.config" -gff3 -v -noCheckProtein $(basename ${genome_dir})
    else
        snpEff build -d -c "${PROJECT_ROOT}/snpEff/snpEff.config" -gtf22 -v -noCheckProtein $(basename ${genome_dir})
    fi
done

declare -a genome_array=( 
    "${PROJECT_ROOT}/Povale/curtisi/Poc221/1a8e7da061e635fcf1fb214d88279e2d41ad2fad/Poc221.fasta"
    "${PROJECT_ROOT}/Povale/wallikeri/Pow222/1a8e7da061e635fcf1fb214d88279e2d41ad2fad/Pow222.fasta"
)

for genome in "${genome_array[@]}"; do
    echo "Creating snpEff database for ${genome}..."

    # fetch genome
    genome_dir="${PROJECT_ROOT}/snpEff/snpEff_database/$(basename ${genome})"
    genome_dir="${genome_dir%.fasta}"
    mkdir -p "${genome_dir}"

    # if [[ "${genome}" = *.gz ]]; then
    if file --mime-type "${genome}" | grep -q gzip$; then
        ln -sf "${genome}" "${genome_dir}/sequences.fa.gz"
    else
        ln -sf "${genome}" "${genome_dir}/sequences.fa"
    fi

    # fetch annotation
    annotation="$(realpath $(dirname ${genome})/*.gff*)"
    if [ -f "${annotation}" ]; then
        type="gff"
        if file --mime-type "${annotation}" | grep -q gzip$; then
            ln -sf "${annotation}" "${genome_dir}/genes.gff.gz"
        else
            ln -sf "${annotation}" "${genome_dir}/genes.gff"
        fi
    else
        type="gtf"
        annotation="$(realpath $(dirname ${genome})/*.gtf*)"
        if file --mime-type "${annotation}" | grep -q gzip$; then
            ln -sf "${annotation}" "${genome_dir}/genes.gtf.gz"
        else
            ln -sf "${annotation}" "${genome_dir}/genes.gtf"
        fi
    fi

    # fetch protein
    protein="${genome%.fasta}_Proteins.fasta"
    ln -sf ${protein} "${genome_dir}/protein.fa"

    # build database
    if [ ! "${type}" == "gtf" ]; then
        snpEff build -c "${PROJECT_ROOT}/snpEff/snpEff.config" -gff3 -v -noCheckCds $(basename ${genome_dir})
    else
        snpEff build -d -c "${PROJECT_ROOT}/snpEff/snpEff.config" -gtf22 -v -noCheckCds $(basename ${genome_dir})
    fi
done
