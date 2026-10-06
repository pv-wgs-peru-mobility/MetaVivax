# Hierarchical clustering

This analysis uses the pairwise genetic distances calculated by [`pairwise_snp_allele_sharing_distance.R`](../pairwise_allele_sharing_distance/pairwise_snp_allele_sharing_distance.R). See that script for the full genotype-parsing and distance-calculation method.

## Distance input

The upstream script reads an analysis-ready table of 16,485 polymorphic biallelic SNPs from 56 samples. For each pair of samples, it compares SNPs with genotype calls in both samples. Identical allele composition scores 1, one shared allele scores 0.5, and no shared alleles scores 0. Pairwise similarity (PS) is the mean score across comparable SNPs; genetic distance is `1 - PS`.

It writes `results/pairwise_allele_sharing/pairwise_distance_matrix.tsv`, which is the input to `hierarchical_clustering.R`. It also writes `pairwise_allele_sharing_long.tsv` (per-pair scores, distances, and sharing counts), `pairwise_PS_matrix.tsv`, and `pairwise_shared_loci_matrix.tsv`. If a sample pair has too few comparable SNPs, the completed distance matrix substitutes the maximum observed distance and flags that substitution in the long table.

## Heterozygous calls and polyclonality

Heterozygous genotype calls are retained as two alleles and can contribute partial sharing (0.5) when a pair shares one allele. This incorporates the observed alleles in mixed calls into the distance calculation. The script does **not** estimate the number of parasite clones.

## Clustering

`hierarchical_clustering.R` reads the distance matrix and sample metadata, then applies average-linkage hierarchical clustering (`hclust(method = "average")`). The script writes cluster assignment and summary tables plus a dendrogram under `results/hierarchical_clustering/`.
