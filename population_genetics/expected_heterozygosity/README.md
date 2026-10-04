PAIRWISE SNP ALLELE-SHARING DISTANCE
====================================

Script
------
pairwise_snp_allele_sharing_distance.R

Purpose
-------
Calculate pairwise SNP similarity (PS) and genetic distance (1 - PS) for use in
minimum spanning tree and hierarchical clustering analyses.

Required R packages
-------------------
- optparse
- data.table

Input
-----
Default input file:

  data/analysis_ready/pv_56samples_16485snps_gt.tsv

The script expects the common analysis-ready genotype dataset used for the
population-genetic analyses. Upstream quality control has already been applied,
including genotype depth >5, sample missingness <=40%, SNP missingness <=50%,
and retention of polymorphic biallelic SNPs. The resulting dataset contains
56 samples and 16,485 SNPs.

The input is a tab-delimited genotype table with columns:

  CHROM, POS, TYPE, sample1.GT, sample2.GT, ...

REF and ALT columns may also be present.

Allele-sharing calculation
--------------------------
For each pair of samples, only SNPs with non-missing genotype calls in both
samples are compared.

Per-locus similarity is scored as:

  1.0   identical allele composition
  0.5   one shared allele
  0.0   no shared alleles

Examples:

  A/A vs A/A   -> 1.0
  A/G vs A/G   -> 1.0
  A/G vs A/A   -> 0.5
  A/G vs G/G   -> 0.5
  A/A vs G/G   -> 0.0

Pairwise similarity is the mean score across comparable SNPs:

  PS = mean(per-locus allele-sharing score)

Genetic distance is:

  distance = 1 - PS

No additional SNP- or sample-level filtering is performed in this script.

Missing pairwise comparisons
----------------------------
If a sample pair has no SNPs with non-missing calls in both samples, the direct
pairwise distance is unavailable. For the complete distance matrix used for
MST and hierarchical clustering, this value is replaced with the maximum
observed pairwise distance.

Output
------
Default output directory:

  results/pairwise_allele_sharing/

Files produced:

1. pairwise_allele_sharing_long.tsv
   Pairwise similarity and distance results in long format.

2. pairwise_distance_matrix.tsv
   Complete 1 - PS distance matrix for MST and hierarchical clustering.

3. pairwise_PS_matrix.tsv
   Pairwise PS similarity matrix.

4. pairwise_shared_loci_matrix.tsv
   Number of comparable SNPs for each sample pair.

Example usage
-------------
Using the default project structure:

  Rscript pairwise_snp_allele_sharing_distance.R

Using explicit paths:

  Rscript pairwise_snp_allele_sharing_distance.R \
    --gt_table data/analysis_ready/pv_56samples_16485snps_gt.tsv \
    --out_dir results/pairwise_allele_sharing

