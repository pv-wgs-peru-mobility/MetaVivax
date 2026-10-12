# Pairwise SNP allele-sharing distance

This script calculates pairwise SNP allele-sharing similarity (`PS`) and genetic distance (`1 - PS`) between the 56 *Plasmodium vivax* infections. The resulting distance matrix is used as input for the minimum spanning tree (MST) and hierarchical clustering analyses.

## Script

```text
pairwise_snp_allele_sharing_distance.R
```

## Required R packages

```text
optparse
data.table
```

## Input file

The default input file is:

```text
data/analysis_ready/pv_56samples_16485snps_gt.tsv
```

This is the analysis-ready genotype table containing the final 56 samples and 16,485 polymorphic biallelic SNPs after upstream quality control. The input table is tab-delimited and must contain:

```text
CHROM, POS, TYPE, sample1.GT, sample2.GT, ...
```

Optional `REF` and `ALT` columns may also be present.

## Allele-sharing calculation

For each pair of samples, the script compares only SNPs where both samples have non-missing genotype calls. Missing calls are excluded pairwise, so each sample pair is compared using the SNPs that are informative for that pair.

At each SNP, the genotype is treated as two observed alleles. Homozygous calls contain two copies of the same allele, while heterozygous calls contain two different alleles. The script then counts how many alleles are shared between the two samples at that SNP.

The per-SNP allele-sharing score is:

```text
1.0  = both alleles are shared
0.5  = one allele is shared
0.0  = no alleles are shared
```

For example:

```text
A/A vs A/A  -> both alleles match          -> score 1.0
A/G vs A/G  -> both alleles match          -> score 1.0
A/G vs G/A  -> both alleles match          -> score 1.0
A/G vs A/A  -> one allele is shared        -> score 0.5
A/G vs G/G  -> one allele is shared        -> score 0.5
A/A vs G/G  -> no alleles are shared       -> score 0.0
```

This means that heterozygous calls are not discarded and are not treated as a single separate category. Instead, they contribute according to their allele composition. A heterozygous call can show full sharing with another heterozygous call containing the same two alleles, partial sharing with a homozygous call containing one of those alleles, or no sharing with a genotype containing neither allele.

Pairwise similarity (`PS`) is the mean allele-sharing score across all comparable SNPs for that sample pair:

```text
PS = mean(per-SNP allele-sharing scores)
```

Genetic distance is then calculated as:

```text
distance = 1 - PS
```

Therefore, pairs with higher allele sharing have smaller genetic distances.

## Simple example

The example below shows how the calculation works for two samples across five SNPs.

```text
SNP   Sample 1   Sample 2   Score
1     A/A        A/A        1.0
2     A/G        A/A        0.5
3     C/T        T/C        1.0
4     G/G        A/A        0.0
5     missing    A/G        excluded
```

The missing SNP is excluded for this pair. The pairwise similarity is:

```text
PS = (1.0 + 0.5 + 1.0 + 0.0) / 4 = 0.625
```

The genetic distance is:

```text
distance = 1 - 0.625 = 0.375
```

## Missing pairwise comparisons

For each sample pair, the script counts the number of SNPs where both samples have non-missing genotype calls. This value is reported as `n_shared_nonmissing_loci`.

A pairwise `PS` value is calculated only when the number of comparable SNPs is at least the value set by `--min_shared_loci`. In this script, the default value is:

```text
--min_shared_loci = 1
```

With this default, a pairwise comparison is calculated whenever the two samples share at least one SNP with non-missing genotype calls. This is suitable here because the input is already an analysis-ready dataset with 56 samples and 16,485 SNPs after upstream quality control, so sample pairs are expected to have many comparable SNPs.
If a sample pair has fewer comparable SNPs than --min_shared_loci, the script keeps the pair in the long-format output but reports the direct PS and 1 - PS distance as NA.
MST and hierarchical clustering require a complete numeric distance matrix. Therefore, any unavailable off-diagonal distance is replaced with the maximum observed pairwise distance. These replacements are flagged in pairwise_allele_sharing_long.tsv using:
distance_replaced = TRUE

The matrix file pairwise_shared_loci_matrix.tsv records the number of comparable SNPs for every sample pair, so the amount of data used for each pairwise estimate can be checked.

## Outputs

The default output directory is:

```text
results/pairwise_allele_sharing/
```

The script writes four output files:

```text
pairwise_allele_sharing_long.tsv
pairwise_distance_matrix.tsv
pairwise_PS_matrix.tsv
pairwise_shared_loci_matrix.tsv
```

`pairwise_allele_sharing_long.tsv` contains one row per sample pair, including the number of comparable SNPs, the number of fully shared, partially shared, and non-shared loci, the pairwise similarity (`PS`), and the genetic distance (`1 - PS`).

`pairwise_distance_matrix.tsv` contains the complete `1 - PS` distance matrix used for MST and hierarchical clustering.

`pairwise_PS_matrix.tsv` contains the pairwise allele-sharing similarity matrix.

`pairwise_shared_loci_matrix.tsv` contains the number of non-missing comparable SNPs for each sample pair.

## Example usage

Using the default project structure:

```bash
Rscript pairwise_snp_allele_sharing_distance.R
```

Using explicit paths:

```bash
Rscript pairwise_snp_allele_sharing_distance.R \
  --gt_table data/analysis_ready/pv_56samples_16485snps_gt.tsv \
  --out_dir results/pairwise_allele_sharing
```

## Recommended repository location

Keep this script and README in:

```text
population_genetics/pairwise_allele_sharing_distance/
```

Recommended folder contents:

```text
population_genetics/pairwise_allele_sharing_distance/pairwise_snp_allele_sharing_distance.R
population_genetics/pairwise_allele_sharing_distance/README.md
```
