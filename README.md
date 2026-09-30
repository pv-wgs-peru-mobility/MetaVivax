# snpEff database structure

SnpEff expects that all the files it needs to be stored in a specific location and with specific filenames, matching the names used in `./config/snpEff.config`. Briefly,

1. A copy of the reference genome should be placed in `./data/snpEff_database/<reference-name>/sequences.fa` (a symlink also works)
2. A gene annotation file should be placed in `./data/snpEff_database/<reference-name>/genes.gff`.
3. Add a CDS fasta file: `./data/snpEff_database/<reference-name>/cds.fa`.
4. Optionally a protein fasta file can also be added: `./data/snpEff_database/<reference-name>/protein.fa`.

The generate_snpEff.sh script takes care of all of these steps.
