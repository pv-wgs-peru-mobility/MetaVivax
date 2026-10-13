# Genetic differentiation

`FST.R` estimates pairwise genetic differentiation using the Weir-Cockerham estimator implemented in `diveRsity`, with 95% confidence intervals. It covers both the community comparisons in the main manuscript and the within-community traveler/non-traveler comparisons.

## Scripts

`FST.R` calculates differentiation, confidence intervals and within-community permutation tests. `plot_fst.R` reads its summary table to create the community heatmaps and within-community confidence-interval figure. `run_FST.slurm` runs both scripts on the VSC cluster.

## Input

```text
data/analysis_ready/pv_56samples_16485snps_gt.tsv
data/analysis_ready/pv_56samples_metadata.tsv
```

The analysis-ready SNP table must retain both observed alleles in every call, for example `0/0`, `0/1`, `1/1` or `C/C`, `C/T`, `T/T`. It contains `CHROM`, `POS`, `TYPE` and sample columns, with optional `REF` and `ALT`. Sample columns may end in `.GT`. Metadata contains a sample identifier (`id`, `Sample`, `sample` or `sample_id`), `community` and `travel` (`0` = non-traveler, `1` = traveler).

The script uses the final 56-infection SNP dataset after upstream quality control, not an allele-sharing distance matrix or a single-allele numeric table. No additional MAF filter is applied for Fst.

## Heterozygous calls

Homozygous calls contribute two copies of the observed allele; heterozygous calls contribute one copy of each allele. For example, `C/T` contributes half C and half T to allele frequencies and is retained as a mixed call for the observed heterozygosity required by Weir-Cockerham. It is not converted to a third allele, reduced to its first allele or discarded.

For these infection-level estimates, the two observed alleles in a mixed call receive equal weight. This retains the information in polyclonal infections without estimating relative clone abundance or treating the two allele entries as two independently sampled infections.

## Method

Community Fst is calculated separately for all infections, non-travelers and travelers. Each dataset includes the three unordered community comparisons: Libertad-Gamitanacocha, Libertad-Urco Mirano and Gamitanacocha-Urco Mirano. Within each community, a further comparison contrasts non-travelers with travelers.

For each comparison, complete calls supply allele frequencies and mixed-call counts. SNPs require at least two non-missing infections in each group and an estimable variance-component denominator. Missing and partially missing calls do not contribute to allele counts.

The script uses `diveRsity::glbWCcpp()` to calculate the Weir-Cockerham components: `a` represents variation between groups, `b` variation among infections within groups, and `c` variation within their two-allele calls. Multilocus Fst is calculated by pooling the components:

```text
Fst = sum(a) / sum(a + b + c)
```

This is a ratio of summed components, not the mean of per-SNP Fst values. Negative estimates are retained in tables and figures; they are interpreted as no detectable differentiation rather than being silently replaced with zero.

The script obtains 95% percentile confidence intervals by resampling SNP loci with replacement 1,000 times and recalculating the pooled estimate. Within-community significance is assessed with 1,000 permutations of whole-infection travel labels, preserving group sizes and each infection's complete multilocus genotype. Bonferroni-adjusted permutation p-values are also saved. The random seed is `1`.

## Worked example

Consider two groups of four infections at one SNP:

```text
Group A
S1: C/C -> two C alleles         -> C = 2, T = 0
S2: C/T -> one C and one T       -> C = 1, T = 1
S3: C/C -> two C alleles         -> C = 2, T = 0
S4: T/T -> two T alleles         -> C = 0, T = 2
Total                            -> C = 5, T = 3

Group B
S1: T/T -> two T alleles         -> C = 0, T = 2
S2: C/T -> one C and one T       -> C = 1, T = 1
S3: T/T -> two T alleles         -> C = 0, T = 2
S4: T/T -> two T alleles         -> C = 0, T = 2
Total                            -> C = 1, T = 7
```text

```text
Group A: C count = 5, T count = 3; p_C = 5/8 = 0.625
Group B: C count = 1, T count = 7; p_C = 1/8 = 0.125
```

Each group has one mixed `C/T` call among four infections, giving an observed mixed-call frequency of `1/4 = 0.25`. Both the allele frequencies and this mixed-call information enter the Weir-Cockerham calculation. Summing its components over the two alleles at this SNP gives:

```text
a = 0.17708333
b = 0.16666667
c = 0.25000000
Fst = 0.17708333 / (0.17708333 + 0.16666667 + 0.25000000)
    = 0.29824561
```

For the full dataset, the same components are summed across all estimable SNPs before taking the ratio. The simple diversity expression `1 - sum(p_i^2)` explains the allele counting, but is not a replacement for the Weir-Cockerham Fst estimator.

## Run

Run from the repository root. Required packages are `optparse`, `data.table`, `diveRsity` and `ggplot2`.

```bash
Rscript population_genetics/fst/FST.R
Rscript population_genetics/fst/plot_fst.R
```

On the cluster:

```bash
sbatch population_genetics/fst/run_FST.slurm
```

The genotype, metadata and output paths can be set with `--gt_table`, `--metadata` and `--out_dir`. The default resampling settings are `--n_boot 1000`, `--n_perm 1000` and `--seed 1`.

## Main outputs

Results are written to `results/fst/`. `fst_summary.tsv` contains multilocus Fst, confidence intervals, group sizes and within-community permutation p-values. `fst_community_heatmaps.pdf` and `.png` show the main community comparisons. `fst_within_community.pdf` and `.png` show traveler/non-traveler differentiation with confidence intervals. The per-locus component table and sample-count table are also saved.

## References

The estimator follows [Weir and Cockerham (1984)](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x). The calculation uses the documented [Weir-Cockerham component function in diveRsity](https://rdrr.io/cran/diveRsity/man/glbWCcpp.html).

## Repository location

Keep this README, `FST.R`, `plot_fst.R` and `run_FST.slurm` in `population_genetics/fst/`.
