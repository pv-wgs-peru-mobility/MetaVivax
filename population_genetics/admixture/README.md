# ADMIXTURE analysis

This folder contains the scripts used to prepare genotype datasets for ADMIXTURE, run ancestry-inference analyses, and generate ancestry bar plots for the *Plasmodium vivax* WGS dataset.

## Scripts

`ld_pruning_admixture.sh` prepares the ADMIXTURE input files. It converts the VCF to PLINK BED format, applies the filtering settings used for the ADMIXTURE input, creates the unpruned dataset, and generates two LD-pruned datasets.

`run_admixture.sh` runs ADMIXTURE on each prepared dataset for `K = 2` to `K = 10`, summarizes cross-validation error, and creates ancestry bar plots using the sample metadata.

## Input files

The main genotype input is:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.vcf.gz
```

The metadata file used for community labels and traveler annotations is:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/SNP_matrix_with_metadata.csv
```

The metadata file must contain a sample identifier column (`id`, `Sample`, or `sample`) and the columns `community` and `travel`.

## Dataset preparation

The VCF is converted to PLINK BED format using PLINK v1.9. The preparation script applies the missingness and allele-frequency filters before ADMIXTURE input files are generated:

```text
per-SNP missingness <= 60%
per-sample missingness <= 30%
MAF >= 0.01
```

The script checks that the prepared data match the manuscript dataset. After filtering, the expected number of samples is `56`, and the unpruned dataset is expected to contain `16,485` SNPs.

Three ADMIXTURE input datasets are generated.

The unpruned dataset contains all retained SNPs after the filtering step. This dataset is used to capture the full genome-wide ancestry signal.

The first LD-pruned dataset is generated with `r2 = 0.1`. This more stringent LD-pruning threshold removes more correlated SNPs and retains a marker set with stronger independence between nearby variants. The expected number of SNPs after this pruning step is `1,745`.

The second LD-pruned dataset is generated with `r2 = 0.2`. This is a less stringent LD-pruning threshold than `r2 = 0.1`, so it retains more SNPs while still reducing the influence of linked markers. The expected number of SNPs after this pruning step is `1,941`.

Using both LD-pruned datasets allows the ADMIXTURE results to be checked under two levels of marker independence, while the unpruned dataset preserves the full SNP signal.

The preparation script stops if the observed sample count or SNP counts do not match the manuscript values.

## ADMIXTURE run

ADMIXTURE v1.3.0 is run separately for each prepared dataset:

```text
77.19298.unpruned
77.19298.pruned_r2_0.1
77.19298.pruned_r2_0.2
```

For each dataset, ADMIXTURE is run for:

```text
K = 2, 3, 4, 5, 6, 7, 8, 9, 10
```

Here, `K` represents the number of inferred ancestry components tested by ADMIXTURE. Testing `K = 2` to `K = 10` allows the analysis to compare models with different numbers of inferred components rather than selecting a single value in advance. 

For each value of `K`, the script performs `5` independent ADMIXTURE runs using different random seeds. Each run uses `100` bootstrap replicates. The independent runs are used to check that the ancestry estimates and cross-validation results are consistent across repeated runs.

The script extracts the cross-validation error from each ADMIXTURE log file. For each dataset, the mean cross-validation error is calculated for every tested `K`. The best-supported `K` is selected as the value with the lowest mean cross-validation error for that dataset.

ADMIXTURE components are interpreted as inferred genetic structure within the analyzed dataset. They are not interpreted as definitive geographic ancestry or as direct assignment to external source populations.

## Run order

Run the scripts in this order:

```bash
sbatch ld_pruning_admixture.sh
sbatch run_admixture.sh
```

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/ADMIXTURE/
```

The main prepared datasets are:

```text
77.19298.unpruned.*
77.19298.pruned_r2_0.1.*
77.19298.pruned_r2_0.2.*
```

For each dataset, the ADMIXTURE run produces ancestry proportion files (`.Q`), allele-frequency files (`.P`), log files, and cross-validation summaries.

The main summary outputs include:

```text
77.19298.*.CV_errors.tsv
77.19298.*.CV_mean_by_K.csv
77.19298.*.CV_curve.jpeg
77.19298.*.K*.barplot.community_split.jpeg
77.19298.admixture_input_counts.txt
```

The CV files are used to select the best-supported `K` for each dataset. The bar plots show ancestry proportions for each sample, split by community, with traveler-associated infections marked in the sample labels.

## Recommended repository location

Keep these scripts and this README in:

```text
population_genetics/admixture/
```

Recommended folder contents:

```text
population_genetics/admixture/ld_pruning_admixture.sh
population_genetics/admixture/run_admixture.sh
population_genetics/admixture/README.md
```
