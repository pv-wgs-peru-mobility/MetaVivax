# Identity-by-descent (IBD) analysis

This folder contains the R script used to estimate identity-by-descent (IBD) sharing between *Plasmodium vivax* isolates using `isoRelate`.

## Files

| File | Description |
| --- | --- |
| `IBD.R` | Runs IBD filtering, parameter estimation, segment inference, summaries, plots, and IBD clustering. |
| `run_IBD.slurm` | SLURM launcher for running `IBD.R` on the VSC cluster. |

## Input files

Default inputs:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.ped
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.map
```

The MAP genetic positions are recalculated from physical positions using 13.7 kb/cM.

## Analysis settings

| Step | Setting |
| --- | --- |
| Minor allele frequency filter | MAF ≥ 0.01 |
| Per-isolate missingness filter | ≤ 30% |
| Per-SNP missingness filter | ≤ 60% |
| Genetic map conversion | bp / 13,700 cM |
| Minimum IBD segment size | 450 SNPs |
| Minimum IBD segment length | 700,000 bp |
| IBD clustering threshold | prop = 0.40 |

PED column 5 is set to `2` before running `isoRelate`, following the infection coding used for this IBD analysis.

## Run

```bash
sbatch /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/run_IBD.slurm
```

## Main outputs

Outputs are written to:

```text
/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/
```

| Output | Description |
| --- | --- |
| `77.19298_geno.rds` | Filtered genotype object from `getGenotypes()`. |
| `77.19298_params.rds` | IBD model parameters from `getIBDparameters()`. |
| `77.19298_ibd_segments.rds` / `.csv` | IBD segments detected by `getIBDsegments()`. |
| `77.19298_ibd_summary.txt` | IBD summary from `getIBDsummary()`. |
| `77.19298_ibd_matrix.rds` | IBD matrix from `getIBDmatrix()`. |
| `77.19298_ibd_proportion.csv` | SNP-level proportion of sample pairs inferred to be IBD. |
| `77.19298_IBD_proportions.pdf` | Plot of IBD proportions. |
| `77.19298_IBD_segments.pdf` | Plot of detected IBD segments. |
| `77.19298_IBD_clusters_network.pdf` | IBD cluster network. |
| `77.19298_IBD_full_network.rds` | IBD network including isolated samples. |
| `77.19298_IBD_network_full_colored.pdf` | Full IBD network plot. |
