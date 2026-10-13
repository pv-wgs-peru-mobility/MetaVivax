# Expected heterozygosity

`expected_heterozygosity.R` calculates SNP-level expected heterozygosity (He) within each community, separately for all infections, non-travelers and traveler-associated infections. 

## Input

The script uses the genotype table:

```text
data/analysis_ready/pv_56samples_16485snps_gt.tsv
data/analysis_ready/pv_56samples_metadata.tsv
```

The genotype table contains `CHROM`, `POS`, `TYPE` and sample genotype columns, with optional `REF` and `ALT`. Calls must retain both alleles, for example `0/0`, `0/1`, `1/1` or `C/C`, `C/T`, `T/T`. Sample columns may end in `.GT`. Metadata must contain a sample identifier (`id`, `Sample`, `sample` or `sample_id`), `community` and `travel`, coded `0` for non-travelers and `1` for travelers.

Common quality control has already been applied upstream. For He, the script additionally retains globally polymorphic biallelic SNPs with global minor allele frequency (MAF) `>= 0.05`, as specified in the Methods. SNP selection is performed across all infections once; the same selected loci are then used for every community and travel-status group.

## Method

At each SNP, allele frequencies are calculated from both alleles in every non-missing genotype. A homozygous call contributes two copies of its allele; a heterozygous call contributes one copy of each allele. Thus, `C/T` has the same frequency contribution as half C and half T. Each infection has equal total weight, and allele order or phasing does not affect the estimate.

He is calculated using the manuscript's formula:

```text
He = n / (n - 1) * (1 - sum(p_i^2))
```

Here, `n` is the number of infections with complete genotype calls at that SNP within the group, and `p_i` is the frequency of allele `i`. The `n/(n-1)` factor is the sample-size correction. `n`is the number of isolates (infections) with complete genotype calls at that SNP within the group being analysed.

For each SNP, He is calculated only when at least two infections in the group have complete genotype calls, such as C/C or C/T. Fully missing calls (./.) and partially missing calls (C/.) are excluded entirely. If fewer than two complete calls remain, He is recorded as NA. SNPs that are monomorphic within a group retain He = 0 in both summaries and plots. Identical heterozygous calls such as C/T still contain two alleles and therefore do not give He = 0.
 
## Heterozygote example

For four infections at one SNP:

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

1 - (p_C^2 + p_T^2) = 1 - (0.390625 + 0.140625) = 0.46875
He = 4 / (4 - 1) * 0.46875 = 0.625
```

The heterozygous infection is retained and contributes to both allele frequencies. For mixed infections, equal weighting describes the observed allele composition; it does not estimate relative clone abundance. The results are infection-level diversity summaries.

## Comparisons and figure

Community distributions are compared using Kruskal-Wallis tests followed by Dunn's pairwise tests, separately for the overall, non-traveler and traveler datasets. Traveler/non-traveler differences within each community are assessed using Wilcoxon rank-sum tests. Pairwise p-values are adjusted using Bonferroni correction, with significance assessed at `p < 0.05`.

The boxplot shows nine distributions: the three travel-status categories within Libertad, Gamitanacocha and Urco Mirano. All estimable values, including zeros and sample-size-corrected values above 0.5, are retained.

## Run

Run from the repository root. Required packages are `optparse`, `data.table`, `ggplot2` and `dunn.test`.

```bash
Rscript population_genetics/expected_heterozygosity/expected_heterozygosity.R
```

The options `--gt_table`, `--metadata` and `--out_dir` allow explicit paths.

## Main outputs

Results are written to `results/expected_heterozygosity/`. `He_per_locus.tsv` contains allele frequencies, the non-missing infection count and He for each SNP/group. `He_summary.tsv` contains group-level summaries. `He_boxplot.pdf` and `He_boxplot.png` contain the diversity figure. The script also saves the global MAF selection table and statistical-test tables.

## Repository location

Keep this README and `expected_heterozygosity.R` in `population_genetics/expected_heterozygosity/`.
