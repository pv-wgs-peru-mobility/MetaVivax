#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=admix-ld-prune
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mahdi.safarpour@uantwerpen.be
#SBATCH -o %j-%x-stdout.out
#SBATCH -e %j-%x-stderr.out

module load PLINK/1.9b_6.21-x86_64

BASE=/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv
VCF_PATH=$BASE/IBD/New/77/77.19298.vcf.gz
OUT_DIR=$BASE/ADMIXTURE
PREFIX=77.19298

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

PLINK_CONTIG_OPTS="--allow-extra-chr --chr-set 95 no-xy no-mt"
SET_IDS='--set-missing-var-ids @:#:\$1:\$2'

GENO=0.60
MIND=0.30
MAF=0.001
WIN=50
STEP=10
R2_VALUES=(0.1 0.2)

EXPECTED_SAMPLES=56
EXPECTED_UNPRUNED_SNPS=16485
EXPECTED_PRUNED_R2_0_1=1745
EXPECTED_PRUNED_R2_0_2=1941

rm_prefix() {
  local p="$1"
  rm -f "${p}.bed" "${p}.bim" "${p}.fam" "${p}.log" "${p}.nosex"
}

check_count() {
  local observed="$1"
  local expected="$2"
  local label="$3"
  if [[ "$observed" -ne "$expected" ]]; then
    echo "ERROR: $label count is $observed; expected $expected"
    exit 1
  fi
}

echo "Input VCF: $VCF_PATH"
echo "Output directory: $OUT_DIR"
echo "Filters: geno=$GENO mind=$MIND maf=$MAF"
echo "LD pruning: window=$WIN step=$STEP r2=${R2_VALUES[*]}"

if [[ ! -f "$VCF_PATH" ]]; then
  echo "ERROR: VCF not found: $VCF_PATH"
  exit 1
fi

rm_prefix "${PREFIX}.bedset"
rm_prefix "${PREFIX}.qc"
rm_prefix "${PREFIX}.unpruned"

plink $PLINK_CONTIG_OPTS \
  --vcf "$VCF_PATH" \
  --double-id \
  $SET_IDS \
  --make-bed \
  --out "${PREFIX}.bedset"

plink $PLINK_CONTIG_OPTS \
  --bfile "${PREFIX}.bedset" \
  --geno "$GENO" \
  --mind "$MIND" \
  --maf "$MAF" \
  --make-bed \
  --out "${PREFIX}.qc"

cp -f "${PREFIX}.qc.bed" "${PREFIX}.unpruned.bed"
cp -f "${PREFIX}.qc.bim" "${PREFIX}.unpruned.bim"
cp -f "${PREFIX}.qc.fam" "${PREFIX}.unpruned.fam"

n_unpruned=$(wc -l < "${PREFIX}.unpruned.bim")
n_samples=$(wc -l < "${PREFIX}.unpruned.fam")

echo "Unpruned SNPs: $n_unpruned"
echo "Samples: $n_samples"

check_count "$n_samples" "$EXPECTED_SAMPLES" "sample"
check_count "$n_unpruned" "$EXPECTED_UNPRUNED_SNPS" "unpruned SNP"

prune_make_bed() {
  local r2="$1"
  local ld="${PREFIX}.ld_r2_${r2}"
  local out="${PREFIX}.pruned_r2_${r2}"

  rm -f "${ld}.prune.in" "${ld}.prune.out" "${ld}.log" "${ld}.nosex"
  rm_prefix "$out"

  plink $PLINK_CONTIG_OPTS \
    --bfile "${PREFIX}.qc" \
    --indep-pairwise "$WIN" "$STEP" "$r2" \
    --out "$ld"

  plink $PLINK_CONTIG_OPTS \
    --bfile "${PREFIX}.qc" \
    --extract "${ld}.prune.in" \
    --make-bed \
    --out "$out"

  local n_prune n_bim
  n_prune=$(wc -l < "${ld}.prune.in")
  n_bim=$(wc -l < "${out}.bim")

  echo "Pruned r2=$r2 SNPs: $n_bim"

  if [[ "$n_prune" -ne "$n_bim" ]]; then
    echo "ERROR: prune.in count does not match pruned BIM count for r2=$r2"
    exit 1
  fi

  if [[ "$r2" == "0.1" ]]; then
    check_count "$n_bim" "$EXPECTED_PRUNED_R2_0_1" "r2=0.1 pruned SNP"
  fi
  if [[ "$r2" == "0.2" ]]; then
    check_count "$n_bim" "$EXPECTED_PRUNED_R2_0_2" "r2=0.2 pruned SNP"
  fi
}

for r2 in "${R2_VALUES[@]}"; do
  prune_make_bed "$r2"
done

cat > "${PREFIX}.admixture_input_counts.txt" <<EOF
samples	$n_samples
unpruned_snps	$n_unpruned
pruned_r2_0.1_snps	$(wc -l < "${PREFIX}.pruned_r2_0.1.bim")
pruned_r2_0.2_snps	$(wc -l < "${PREFIX}.pruned_r2_0.2.bim")
EOF

echo "ADMIXTURE input files are ready in: $OUT_DIR"
