# IBD community heatmap

This script creates the community-by-community identity-by-descent (IBD) heatmap used to summarize parasite relatedness between the three study communities. It is a downstream plotting and summary script: it does not re-run `isoRelate`, but uses the annotated pairwise IBD table generated from the IBD analysis.


## Script

```text
ibd_community_heatmap.R
```

## Input file

The script uses the annotated pairwise IBD table:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv
```

This file contains one row per isolate pair and includes the pairwise IBD proportion plus the community assignment of both isolates.

The required columns are:

```text
sample1, sample2, ibd_prop, community1, community2
```

`ibd_prop` is the pairwise proportion of the genome inferred to be IBD between two isolates. Larger values indicate greater IBD sharing.

## Method

Isolate pairs were classified as related when their pairwise IBD proportion was `>= 0.50` and as unrelated when their pairwise IBD proportion was `< 0.50`. This threshold was used to summarize strong pairwise relatedness in the community-level heatmap. The continuous pairwise IBD values were retained in the input table, while the threshold was used only to calculate the related-pair proportion `R` for each community pair.

For each community pair, the script counts the total number of isolate pairs and the number of pairs with `ibd_prop >= 0.50`. It then calculates the related-pair proportion `R` as:

```text
R = number of pairs with IBD >= 0.50 / total number of pairs

The script treats each community comparison as an unordered pair. For example, Libertad - Gamitanacocha and Gamitanacocha - Libertad are summarized as the same comparison. The unordered summary is then converted into a symmetric table for plotting, so that the heatmap can be read in both directions.
The community order used in the plot is:
Libertad, Gamitanacocha, UrcoMirano

Within-community cells show IBD relatedness among pairs of isolates sampled from the same community. Between-community cells show IBD relatedness among isolate pairs from two different communities. Each heatmap cell reports the number of related pairs, the total number of pairs, and the value of R.
The color scale ranges from 0 to 1, where darker red indicates a larger proportion of related pairs at the IBD >= 0.50 threshold.
```

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/12_final_selected_figures/
```

The main outputs are:

```text
final_community_heatmap_threshold_0_50_unordered_summary.csv
final_community_heatmap_threshold_0_50_plot_data.csv
final_community_heatmap_threshold_0_50_corrected.pdf
final_community_heatmap_threshold_0_50_corrected.png
final_community_heatmap_threshold_0_50_corrected.tiff
```

The script also writes:

```text
final_community_heatmap_threshold_0_50_symmetry_check.csv
```

The PDF is used for vector-quality inspection, while the PNG and TIFF files are publication-ready raster versions of the same heatmap.

## Recommended repository location

Keep this script and README in a subfolder of the IBD analysis:

```text
population_genetics/IBD/ibd_community_heatmap/
```

Recommended folder contents:

```text
population_genetics/IBD/ibd_community_heatmap/ibd_community_heatmap.R
population_genetics/IBD/ibd_community_heatmap/README.md
```
