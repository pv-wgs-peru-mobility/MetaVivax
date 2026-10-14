# WGS pipeline

## Overview

This section contains the scripts used to process paired-end whole-genome sequencing reads, align parasite reads, call and filter variants, and annotate their predicted functional consequences. The workflow is organised as a modular Bash pipeline: a SLURM batch launcher runs separate scripts and supplies their settings through command-line options.

The processing order is read quality control and trimming, alignment and BAM preparation, variant calling and filtering, and variant annotation.

## Before running the pipeline

### Project paths

The launcher uses `PROJECT_ROOT` to locate the pipeline files. It determines its own directory, `SCRIPT_DIR`, and defines the project root as that directory's parent. The runtime project therefore needs these paths:

```text
scripts/qc.sh
scripts/alignment.sh
scripts/variant-calling.sh
scripts/annotation-combine.sh
data/fastq/
data/samplesheet.csv
config/snpEff.config
data/snpEff_database/
results/
```

Place the launcher in the project's `scripts/` directory when using its existing project-root calculation. These paths describe the required runtime layout; preserve that layout or update the paths when installing the pipeline elsewhere.

### FASTQ files and samplesheet

Place the raw paired-end FASTQ files in `data/fastq/`. Each sample has a read-1 file and its matching read-2 file, for example:

```text
Sample01_1.fastq.gz
Sample01_2.fastq.gz
Sample02_1.fastq.gz
Sample02_2.fastq.gz
```

The samplesheet is a CSV file with `sample` and `species` columns. Use sample identifiers that match those extracted by `alignment.sh`. For the simple naming convention above:

```csv
sample,species
Sample01,pv
Sample02,pv
```

`pv` identifies *Plasmodium vivax*. If processing other species, use the species codes supported by the pipeline scripts.

### Software environment and SLURM launcher

The `.slurm` file is a **SLURM batch script**. Its `#!/usr/bin/env bash` line selects Bash, while the `#SBATCH` lines tell the scheduler which resources and logging settings to use.

The supplied launcher requests one node, one task and 28 CPUs for that task, with a maximum runtime of 72 hours. The processing scripts read `SLURM_CPUS_PER_TASK` to determine their thread allocation. The launcher sets an email recipient and requests notifications at job start, completion and failure.

Standard output and standard error are written separately using `%j-%x-stdout.out` and `%j-%x-stderr.out`; `%j` is the job ID and `%x` is the job name. The job name is a scheduler label and can be changed to describe the pipeline. These settings are documented in the official [SLURM sbatch manual](https://slurm.schedmd.com/sbatch.html).

The active software setup adds the container-wrapper executable directory to `PATH` and loads the `parallel` module. The commented module and Conda examples describe alternative environment setups. Update the executable path and email address for the account running the job.

The launcher runs the four stages described below. Prepare the SnpEff database described in the annotation step before submitting the full job.

## How the Bash commands work

The stage commands below are **Bash script invocations with command-line options**.

`bash` runs the named shell script. `${PROJECT_ROOT}` expands to the project-directory path, and the double quotes keep each complete path together as one argument. A backslash at the end of a line continues the same command on the following line; it must be the final character on that line.

`--input_dir` specifies the directory containing the raw reads, and `--output_dir` specifies where results are written. These are options defined by `qc.sh`, followed by their values. The shorter options `-r1`, `-r2` and `-e` are also defined by the script. These options tell the script which filename suffix identifies each mate of a paired-end read pair:

```text
Sample01_1.fastq.gz -> suffix _1 -> read 1
Sample01_2.fastq.gz -> suffix _2 -> read 2
```

The shared prefix, `Sample01`, connects the two files as a pair. The `-e ".fastq.gz"` option specifies the raw-file extension separately from the read suffix.

## Pipeline stages

### 1. Read quality control and trimming

`qc.sh` assesses the raw reads with FastQC, screens reads against the configured references with FastQ Screen, trims paired reads using fastp, and assesses the trimmed reads again. Use:

```bash
bash "${PROJECT_ROOT}/scripts/qc.sh" \
  --input_dir "${PROJECT_ROOT}/data/fastq/" \
  --output_dir "${PROJECT_ROOT}/results/" \
  -r1 "_1" -r2 "_2" -e ".fastq.gz"
```

With `-r1 "_1" -r2 "_2"`, trimmed reads are written as:

```text
results/fastp/Sample01_1.trim.fastq.gz
results/fastp/Sample01_2.trim.fastq.gz
```

The raw extension can be `.fq.gz` or `.fastq.gz`, but this QC script writes compressed trimmed files ending in `.trim.fastq.gz`.

### 2. Alignment and BAM preparation

`alignment.sh` reads the trimmed pairs, obtains the species assignment from the samplesheet, aligns reads, and prepares sorted, duplicate-marked BAM files. Without `--competitive`, it first aligns reads to the human reference and then aligns the retained unaligned read pairs to the parasite reference. Adding `--competitive` selects alignment against the configured combined reference instead.

For the `_1`/`_2` convention and sample-first identifiers, use:

```bash
bash "${PROJECT_ROOT}/scripts/alignment.sh" \
  --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
  --output_dir "${PROJECT_ROOT}/results/" \
  -r1 "_1" -r2 "_2" \
  --fastq_identifier name_first
```

The read suffixes must match the trimmed filenames produced by QC. Thus, the `_1`/`_2` example uses the same suffix values in both commands. The supplied alignment implementation explicitly reads `.trim.fastq.gz` files.

### 3. Variant calling and filtering

`variant-calling.sh` uses GATK HaplotypeCaller to generate sample-level GVCFs, combines them with GenomicsDBImport, performs joint genotyping with GenotypeGVCFs, and applies the variant filters defined in the script.

```bash
bash "${PROJECT_ROOT}/scripts/variant-calling.sh" \
  --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
  --output_dir "${PROJECT_ROOT}/results/"
```

### 4. Variant annotation and combination

This stage annotates the filtered variants using a matching SnpEff database, then combines outputs through `annotation-combine.sh`.

#### SnpEff database preparation

SnpEff uses a reference-specific database to annotate variants against gene and transcript models. In this project, the database inputs are stored under:

```text
data/snpEff_database/<reference-name>/
```

The `<reference-name>` must match the genome identifier configured in `config/snpEff.config` and the identifier used for annotation. The files in this directory are:

- **`sequences.fa`:** the reference genome sequence in FASTA format. A symbolic link to the reference FASTA can be used instead of a duplicate copy.
- **`genes.gff`:** the gene annotation in GFF format, describing genes, transcripts and their genomic features.
- **`cds.fa`:** coding DNA sequences used to check the coding sequences reconstructed from the genome and annotation.
- **`protein.fa`:** optional reference protein sequences for checking translated proteins against the reference annotation.

For example, a database directory containing all four files has these paths:

```text
data/snpEff_database/<reference-name>/sequences.fa
data/snpEff_database/<reference-name>/genes.gff
data/snpEff_database/<reference-name>/cds.fa
data/snpEff_database/<reference-name>/protein.fa
```

Use reference and annotation files from the same genome assembly and annotation release. Chromosome names must agree between the genome FASTA, the GFF and the VCF being annotated. Transcript identifiers in the validation FASTA files must correspond to the annotation. See the official [SnpEff database-building documentation](https://pcingola.github.io/SnpEff/snpeff/build_db/) and [GFF input documentation](https://pcingola.github.io/SnpEff/snpeff/build_db_gff_gtf/).

`generate_snpEff.sh` prepares the files in the layout described above. Prepare and build this database before submitting the full SLURM job, so it is available when annotation starts. The supplied launcher does not call the database-preparation script.

For a GFF3 annotation, an example build command is:

```bash
snpEff build -v -gff3 \
  -c "${PROJECT_ROOT}/config/snpEff.config" \
  -dataDir "${PROJECT_ROOT}/data/snpEff_database" \
  "${REFERENCE_NAME}"
```

For an interactive build, set `PROJECT_ROOT` to the runtime project directory and `REFERENCE_NAME` to the exact genome identifier in the configuration first. Here, `-c` selects the configuration file, `-dataDir` explicitly selects the database directory, `-gff3` selects the annotation format, and `-v` requests detailed output. Check the build report before using the database for annotation.

#### Run annotation

After the preceding processing stages have completed and the database is available, the launcher runs:

```bash
bash "${PROJECT_ROOT}/scripts/annotation-combine.sh" \
  --samplesheet "${PROJECT_ROOT}/data/samplesheet.csv" \
  --output_dir "${PROJECT_ROOT}/results/"
```

## Run the complete workflow

After preparing the project inputs, software environment and SnpEff database, submit the launcher. For a launcher named `run_wgs.slurm`, use:

```bash
sbatch scripts/run_wgs.slurm
```

Replace that example filename with the actual launcher filename. The start and finish timestamps and the scheduler log files help track the run.

## Main outputs

The pipeline produces read-quality reports, trimmed paired-end FASTQ files, aligned and duplicate-marked BAM files, sample-level GVCFs, jointly genotyped and filtered VCFs, and annotated variant outputs. Results are written beneath the `results/` directory supplied to the scripts.

## Practical tips

### Match read suffixes and extensions

For a different naming convention, change the values to match the actual filenames:

```text
Sample01_1.fastq.gz / Sample01_2.fastq.gz
  -> -r1 "_1" -r2 "_2" -e ".fastq.gz"

Sample01_R1.fastq.gz / Sample01_R2.fastq.gz
  -> -r1 "_R1" -r2 "_R2" -e ".fastq.gz"

Sample01_R1_001.fastq.gz / Sample01_R2_001.fastq.gz
  -> -r1 "_R1_001" -r2 "_R2_001" -e ".fastq.gz"

Sample01_1.fq.gz / Sample01_2.fq.gz
  -> -r1 "_1" -r2 "_2" -e ".fq.gz"
```

The options identify existing filenames; they do not rename the raw files. Use `"_1"` and `"_2"` in the Bash script, without backslashes before the underscores.

Keep the read-1 and read-2 suffixes consistent between QC and alignment. For example, QC outputs named `_1.trim.fastq.gz` and `_2.trim.fastq.gz` require `-r1 "_1" -r2 "_2"` in the alignment command.

### Choose the correct sample-name parser

The `--fastq_identifier` option controls how the sample identifier is extracted from a filename. In the supplied alignment script, `name_first` uses the first underscore-separated field. For a filename such as `FCHTVJYCCXY_L1_WHRDMALtbmRAAAAA-23_1.fq.gz`, the sample identifier is the third field, `WHRDMALtbmRAAAAA-23`; that convention uses `--fastq_identifier name_middle` instead. The samplesheet should use the corresponding sample identifier consistently across alignment and variant calling.

### Check paths, software and logs

Confirm that the launcher can locate `scripts/`, `data/`, `config/` and `results/` under `PROJECT_ROOT`. Update the container-wrapper path and notification email for the computing account being used. Check the output and error logs for each run, including the SnpEff database-build report.

Adding `set -euo pipefail` immediately after the `#SBATCH` directives makes the launcher stop when a stage fails.

## Downstream analyses

The filtered and annotated variant data support the workflows in [population_genetics](../population_genetics/). Each analysis README describes its required input format, sample metadata and any analysis-specific processing.
