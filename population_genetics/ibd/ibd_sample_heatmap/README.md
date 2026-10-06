# IBD sample heatmap

This script creates the sample-level pairwise IBD heatmap. Samples are ordered by community and travel status.

## Script

```text
ibd_sample_heatmap.R
```

## Inputs

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/02_metadata_clean/metadata_clean.csv
```

Required pairwise columns:

```text
sample1, sample2, ibd_prop
```

Required metadata columns:

```text
sample_id, community, travel_status
```

## Method

The script builds a symmetric sample-by-sample IBD matrix, adds the diagonal, orders samples by community and travel status, and marks travel status with symbols.

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/12_final_selected_figures/
```

Main outputs:

```text
final_sample_level_ibd_heatmap_metadata_ordered_plot_data.csv
final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.pdf
final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.png
final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.tiff
```
