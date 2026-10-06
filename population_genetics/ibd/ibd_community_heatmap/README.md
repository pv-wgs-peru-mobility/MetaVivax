# IBD community heatmap

This script creates the community-by-community IBD heatmap using the annotated pairwise IBD table from the IBD pairwise relatedness analysis.

## Script

```text
ibd_community_heatmap.R
```

## Input

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv
```

Required columns:

```text
sample1, sample2, ibd_prop, community1, community2
```

## Method

Pairs are classified as related when:

```text
IBD ≥ 0.50
```

For each community pair, the script calculates the number of total pairs, the number of related pairs, and the related-pair proportion `R`.

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/12_final_selected_figures/
```

Main outputs:

```text
final_community_heatmap_threshold_0_50_unordered_summary.csv
final_community_heatmap_threshold_0_50_plot_data.csv
final_community_heatmap_threshold_0_50_corrected.pdf
final_community_heatmap_threshold_0_50_corrected.png
final_community_heatmap_threshold_0_50_corrected.tiff
```
