# IBD pairwise relatedness analysis

This script estimates IBD sharing with `isoRelate` and prepares the pairwise relatedness tables used for community, travel, heatmap, and network summaries.

## Script

```text
ibd_pairwise_relatedness_analysis.R
```

## Inputs

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.ped
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.map
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/SNP_matrix_with_metadata.csv
```

The metadata file must contain sample ID, community, and travel status columns.

## Main settings

| Step | Setting |
| --- | --- |
| MAF filter | MAF ≥ 0.01 |
| Per-isolate missingness | ≤ 30% |
| Per-SNP missingness | ≤ 60% |
| Genetic map conversion | bp / 13,700 cM |
| Minimum IBD segment size | 450 SNPs |
| Minimum IBD segment length | 700,000 bp |
| Main relatedness threshold | IBD ≥ 0.50 |
| Sensitivity thresholds | IBD ≥ 0.05 and IBD ≥ 0.80 |

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/
```

Main outputs include:

```text
01_ibd_from_scratch/
02_metadata_clean/metadata_clean.csv
03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv
04_community_relatedness/
05_travel_relatedness/
06_within_between_community/
07_heatmaps/
08_ibd_distributions/
09_map_connectivity/
10_sensitivity_analysis/
11_networks/
```

The file `03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv` is the input for the final community-level and sample-level heatmap scripts.
