# Hierarchical clustering

This analysis uses the pairwise SNP distance matrix produced by `pairwise_snp_allele_sharing_distance.R`.

## Input

```text
results/pairwise_allele_sharing/pairwise_distance_matrix.tsv
data/analysis_ready/pv_56samples_metadata.tsv
```

The distance matrix contains pairwise genetic distances calculated as `1 - PS`.

## Method

`hierarchical_clustering.R` reads the distance matrix and performs average-linkage hierarchical clustering using:

```r
hclust(as.dist(dist_mat), method = "average")
```

The dendrogram cut is `K = 4`. Samples are colored by community, and traveler-associated infections are marked with a cross.

## Outputs

```text
results/hierarchical_clustering/hierarchical_cluster_assignments.tsv
results/hierarchical_clustering/hierarchical_clustering_k4.pdf
results/hierarchical_clustering/hierarchical_clustering_k4.png
```
