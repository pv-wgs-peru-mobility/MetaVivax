# ADMIXTURE analysis

This folder contains the scripts used to prepare ADMIXTURE input files, perform LD pruning, run ADMIXTURE, and generate ancestry bar plots.

## Scripts

| Script | Purpose |
| --- | --- |
| `ld_pruning_admixture.sh` | Converts the analysis VCF to PLINK BED format, applies the manuscript QC thresholds, and creates unpruned and LD-pruned datasets. |
| `run_admixture.sh` | Runs ADMIXTURE for K = 2 to K = 10 and creates CV and ancestry plots. |

## Input

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.vcf.gz
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/SNP_matrix_with_metadata.csv
```

## Settings checked against the manuscript

| Item | Setting |
| --- | --- |
| Samples expected after filtering | 56 |
| Unpruned SNPs expected | 16,485 |
| LD pruning threshold 1 | r2 = 0.1 |
| LD pruning threshold 2 | r2 = 0.2 |
| SNPs expected after r2 = 0.1 pruning | 1,745 |
| SNPs expected after r2 = 0.2 pruning | 1,941 |
| K values tested | 2 to 10 |
| Independent ADMIXTURE runs per K | 5 |
| Bootstrap replicates per run | 100 |
| Manuscript K for unpruned data | K = 4 |
| Manuscript K for LD-pruned data | K = 5 |

The LD-pruning script stops if the observed sample or SNP counts do not match the manuscript values.

## Run order

```bash
sbatch ld_pruning_admixture.sh
sbatch run_admixture.sh
```

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/ADMIXTURE/
```

Main outputs include:

```text
77.19298.unpruned.*
77.19298.pruned_r2_0.1.*
77.19298.pruned_r2_0.2.*
77.19298.*.CV_errors.tsv
77.19298.*.CV_mean_by_K.csv
77.19298.*.CV_curve.jpeg
77.19298.*.K4.barplot.community_split.jpeg
77.19298.*.K5.barplot.community_split.jpeg
```
