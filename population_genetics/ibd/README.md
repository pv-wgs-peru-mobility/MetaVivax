# Identity-by-descent (IBD) analysis

This analysis estimates identity-by-descent (IBD) sharing between *Plasmodium vivax* isolates using `isoRelate`. The script starts from PLINK PED/MAP files prepared for the WGS dataset and produces IBD segments, pairwise IBD summaries, IBD proportion plots, and an IBD network.

The main script is `IBD.R`. The file `run_IBD.slurm` is the SLURM launcher used to run the analysis on the VSC cluster.

## Input files

The default input files are:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.ped
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.map
```

The PED file contains the isolate genotypes and the MAP file contains the SNP positions. Before running `isoRelate`, the script recalculates the genetic map positions from physical positions using:

```text
position in cM = position in bp / 13,700
```

PED column 5 is used by isoRelate as the infection/MOI coding field. It is set to 2 for all isolates before the isoRelate steps so that mixed genotype states can be retained in the IBD analysis. This coding allows the analysis to use the observed allele information from infections with more than one allele at a SNP, rather than forcing all infections to be treated as strictly monoclonal.

## Analysis workflow

`IBD.R` follows the full IBD workflow used for the analysis.

First, the script loads the PED and MAP files and assigns standard names to the first six PED columns. The MAP file is then reformatted so that the SNP ID, physical position, and recalculated genetic position are available in the format required by `isoRelate`.

Next, `getGenotypes()` is used to create the filtered genotype object. The filtering settings are:

- minor allele frequency: MAF >= 0.01
- per-isolate missingness: <= 30%
- per-SNP missingness: <= 60%
- input map distance: cM

The filtered genotype object is saved as `77.19298_geno.rds` and is used for all downstream IBD steps.

The script then estimates the IBD model parameters with `getIBDparameters()`. These parameters are saved as `77.19298_params.rds`.

IBD segments are inferred with `getIBDsegments()`. The segment-calling settings are:

- minimum number of SNPs per IBD segment: 450
- minimum segment length: 700,000 bp
- genotype error rate: 0.01

The detected IBD segments are saved as both an R object and a CSV file.

After segment inference, the script creates the main IBD summaries. It writes the segment summary, the IBD matrix, and the pairwise IBD proportion file. These outputs are used to describe the amount of IBD sharing between isolates and to support the downstream IBD network analysis.

The script also produces two summary plots: one showing the distribution of IBD proportions and one showing the distribution of detected IBD segments.

Finally, the script performs IBD clustering using getIBDpclusters(). In the reported analysis, this clustering step was repeated across a range of IBD proportion thresholds: 0.05, 0.20, 0.40, 0.50, 0.60, 0.80, and 0.90. These thresholds were used to evaluate how the inferred IBD network changed from permissive to more stringent definitions of relatedness.

Networks were retained and summarized up to the 0.80 threshold. At 0.90, no reportable IBD network result was obtained, so the 0.90 threshold was not included in the final reported outputs.
A second full-network object is also created by adding isolates with no IBD links at the selected threshold. This keeps connected and isolated samples in the same network object and allows the final network plot to show both isolates with IBD links and isolates without detectable links at that threshold.

## Run

The analysis can be run on the VSC cluster with:

```bash
sbatch /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/run_IBD.slurm
```

## Outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/
```

The main output files are:

- `77.19298_geno.rds`: filtered genotype object from `getGenotypes()`
- `77.19298_params.rds`: IBD model parameters from `getIBDparameters()`
- `77.19298_ibd_segments.rds`: detected IBD segments
- `77.19298_ibd_segments.csv`: detected IBD segments in table format
- `77.19298_ibd_summary.txt`: IBD segment summary
- `77.19298_ibd_matrix.rds`: pairwise IBD matrix
- `77.19298_ibd_proportion.csv`: pairwise IBD proportions
- `77.19298_IBD_proportions.pdf`: plot of IBD proportions
- `77.19298_IBD_segments.pdf`: plot of detected IBD segments
- `77.19298_clusters_prop0.40.rds`: IBD clustering result at `prop = 0.40`
- `77.19298_IBD_clusters_network.pdf`: IBD cluster network
- `77.19298_IBD_full_network.rds`: full IBD network including isolates without IBD links
- `77.19298_IBD_network_full_colored.pdf`: full IBD network plot

## Recommended repository location

Keep these files together in:

```text
population_genetics/IBD/
```

The recommended folder contents are:

```text
population_genetics/IBD/IBD.R
population_genetics/IBD/run_IBD.slurm
population_genetics/IBD/README.md
```
