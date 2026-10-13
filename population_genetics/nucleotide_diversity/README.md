# Nucleotide diversity

`nucleotide_diversity.R` calculates SNP-based nucleotide diversity (pi) at each of the 16,485 polymorphic biallelic SNPs, separately by community and travel status.
## Input

```text
data/analysis_ready/pv_56samples_16485snps_gt.tsv
data/analysis_ready/pv_56samples_metadata.tsv
```

The genotype table retains both observed alleles in each call, for example `0/0`, `0/1`, `1/1` or `C/C`, `C/T`, `T/T`. It contains `CHROM`, `POS`, `TYPE` and sample columns, with optional `REF` and `ALT`. Sample names may end in `.GT`. Metadata contains a sample identifier (`id`, `Sample`, `sample` or `sample_id`), `community` and `travel` (`0` = non-traveler, `1` = traveler).

## Method 

At each SNP, allele frequencies are calculated from the complete calls within the selected group. A homozygous call contributes two copies of one allele, and a heterozygous call contributes one copy of each allele. Consequently, each infection has equal total weight and a `C/T` call contributes half C and half T to the allele-frequency calculation. Phased calls and reversed allele order are treated equivalently.

The manuscript's per-SNP diversity measure is:

```text
pi = 1 - sum(p_i^2)
```

For a biallelic SNP, this is equivalent to `2 * p * (1 - p)`. At least two non-missing infections are required in the group at that SNP. Missing or partially missing calls are excluded; fewer than two complete calls give a missing estimate. If the available calls all contain the same allele, `pi = 0`. Zero values are retained in summaries and violin plots.

## Heterozygote example

```text
S1: C/C -> two C alleles         -> C = 2, T = 0
S2: C/T -> one C and one T       -> C = 1, T = 1
S3: C/C -> two C alleles         -> C = 2, T = 0
S4: T/T -> two T alleles         -> C = 0, T = 2
Total                            -> C = 5, T = 3
```

```text
p_C = 5 / 8 = 0.625
p_T = 3 / 8 = 0.375
pi = 1 - (0.625^2 + 0.375^2) = 0.46875
```

The `C/T` call is retained and contributes to both allele frequencies. This is the uncorrected per-SNP measure specified for pi; unlike the He analysis, no `n/(n-1)` factor is applied. For mixed infections, equal weighting represents the observed allele composition, not relative clone abundance.

## Comparisons and figure

Calculations are performed for all infections, non-travelers, and travelers within Libertad, Gamitanacocha and Urco Mirano. The nine distributions are displayed as violin plots with embedded boxplots and median points.

For each travel-status category, community differences are assessed using a Friedman test with the same SNP as the matched unit. Only loci with an estimate in all three communities enter that test. Pairwise community comparisons use paired Wilcoxon signed-rank tests on those matched loci. Within-community traveler/non-traveler comparisons use Wilcoxon rank-sum tests, as specified in the manuscript. Pairwise p-values use Bonferroni correction and `p < 0.05` as the significance threshold.

These results summarize diversity among the retained SNP loci, as described in the Methods and figure caption; they are not a diversity estimate averaged over all genomic bases.

## Run

Run from the repository root. Required packages are `optparse`, `data.table` and `ggplot2`.

```bash
Rscript population_genetics/nucleotide_diversity/nucleotide_diversity.R
```

The options `--gt_table`, `--metadata` and `--out_dir` allow explicit paths.

## Main outputs

Results are written to `results/nucleotide_diversity/`. `pi_per_locus.tsv` contains allele frequencies, non-missing infection counts and pi by SNP/group. `pi_summary.tsv` contains group summaries. `pi_violin.pdf` and `pi_violin.png` contain the diversity figure. Community and travel-comparison statistics are saved separately.

## Repository location

Keep this README and `nucleotide_diversity.R` in `population_genetics/nucleotide_diversity/`.
