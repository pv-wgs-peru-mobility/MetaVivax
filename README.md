# MetaVivax

Whole-genome sequencing and population genetics of *Plasmodium vivax* infections in the Peruvian Amazon.

This repository contains the bioinformatics workflows used to characterize *P. vivax* infections in individuals with and without reported local travel from Libertad, Gamitanacocha and Urco Mirano, in the Mazan district of Loreto, Peru. The analyses combine whole-genome sequencing data with community and travel-history metadata to compare genetic diversity, differentiation, population structure and relatedness among infections.

The final population-genetic dataset comprises 56 infections and 16,485 polymorphic biallelic SNPs after upstream quality control. Infections are analysed by community and travel status, with analysis-specific SNP selection where required.

## Repository structure

The repository is organised into two main sections:

- [wgs_pipeline](wgs_pipeline/): processing sequencing reads, calling and filtering variants, and preparing genotype data.
- [population_genetics](population_genetics/): analysing genetic diversity, differentiation, population structure and relatedness using the prepared SNP data and sample metadata.

## Whole-genome sequencing pipeline

The [wgs_pipeline](wgs_pipeline/) section contains Bash scripts and SLURM job launchers for processing sequencing data on the VSC computing cluster. The workflow covers read-quality assessment and trimming, removal of human-derived reads, alignment to the *P. vivax* reference genome, alignment processing and duplicate marking, variant calling and joint genotyping, variant filtering and annotation, and preparation of genotype tables for downstream analyses.

Quality control is applied to variants, genotype calls and samples before preparing the final SNP dataset. The principal outputs are filtered VCF files and genotype data used by the population-genetic workflows.

## Population genetics

The [population_genetics](population_genetics/) section contains R scripts and accompanying analysis documentation, together with Bash or SLURM launchers where needed. The workflows examine the following aspects of parasite genetic variation:

- **Genetic diversity:** expected heterozygosity and SNP-based nucleotide diversity, calculated within communities and travel-status groups.
- **Genetic differentiation:** Weir-Cockerham Fst comparisons between communities and between traveler-associated and non-traveler infections within each community.
- **Population structure:** PCA, DAPC and ADMIXTURE analyses to examine genetic clustering and ancestry-component profiles.
- **Genetic relationships:** pairwise SNP allele-sharing distances used for minimum spanning trees, average-linkage hierarchical clustering and neighbour-joining trees.
- **Identity by descent:** IBD estimation with `isoRelate`, followed by pairwise relatedness summaries, sample- and community-level heatmaps, and network analyses.

These analyses compare genomic patterns among traveler-associated and non-traveler infections and assess their relationships to parasite populations sampled in the three communities. Each analysis has its own README describing the input files, calculation method, handling of heterozygous genotype calls, parameters and main outputs.

## Using the workflows

Start with [wgs_pipeline](wgs_pipeline/) to process sequencing reads and prepare the filtered genotype data. The scripts in [population_genetics](population_genetics/) then use the relevant genotype or variant files together with sample metadata. Where an analysis depends on an upstream result, its README identifies that input; for example, MST and hierarchical clustering use the distance matrix produced by the pairwise SNP allele-sharing script.

Refer to the individual scripts and analysis READMEs for software requirements, execution order and output locations. Cluster-specific paths and module settings should be adjusted to the computing environment in which the scripts are run.
