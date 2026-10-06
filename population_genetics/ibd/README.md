# Identity-by-descent (IBD) analysis

`IBD.R` estimates IBD segments from PLINK PED/MAP genotypes with isoRelate, summarizes pairwise genome sharing, and compares IBD networks at **5%, 20%, 40%, 60%, and 80%** sharing thresholds. This analysis starts from the PED/MAP files; it does not use the allele-sharing distance matrix made for the minimum spanning tree and hierarchical clustering.

## Files and inputs

- `IBD.R` — IBD inference, pairwise fractions, network summaries, and figures.
- `run_IBD.slurm` — VSC SLURM launcher for `IBD.R` (3 CPU cores, 72 hours).
- Default genotype inputs: `/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/New/77/77.19298.ped` and `.map`.
- Default metadata: `/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/SNP_matrix_with_metadata.csv`, with `Sample` or `id`, `travel`, and optionally `community` columns. Sample identifiers must match the filtered PED by FID/IID, IID, or FID. Every retained sample needs valid travel metadata.

The script applies MAF ≥0.01, isolate missingness ≤0.40, and SNP missingness ≤0.60 through `getGenotypes()`. It assigns MOI=2 to all isolates, as the supplied scripts did. That is a **modeling assumption**, not an MOI estimate. Set `IBD_MOI_MODE=ped` to use curated 1/2 values from PED column 5 instead. By default, MAP genetic positions are recalculated as `bp / 13700` cM, matching the supplied scripts. Set `IBD_MAP_MODE=input` to use existing cM values in MAP column 3.

IBD segments are called with at least 450 SNPs, at least 700,000 bp, and genotyping error 0.001. The pairwise genome fraction is the sum of detected segment lengths for a sample pair divided by the total span from the first to last retained SNP on each chromosome. Networks join pairs whose fraction meets or exceeds the threshold; all filtered samples, including isolated ones, remain as vertices. This follows the fraction and threshold definition in [isoRelate's `getIBDpclusters()` implementation](https://rdrr.io/github/bahlolab/isoRelate/src/R/get_ibd_p_clusters.R).

## Run on VSC

Place `IBD.R` and `run_IBD.slurm` in `/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/`, then submit:

```bash
sbatch /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/run_IBD.slurm
```

The R bundle must provide `isoRelate` and `igraph` (and isoRelate's plotting dependencies). To use another dataset or directory, export `IBD_INPUT_DIR`, `IBD_METADATA`, `IBD_OUT_DIR`, or `IBD_PREFIX` before submitting. Other supported overrides are `IBD_EXPECTED_SAMPLES` (default 56; 0 disables the check), `IBD_CORES`, `IBD_ISOLATE_MAX_MISSING`, `IBD_MOI_MODE`, `IBD_MAP_MODE`, and `IBD_BP_PER_CM`. SLURM normally exports these variables to the job.

## Outputs

By default, files go to `/scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/results/`:

| File or folder | Contents |
| --- | --- |
| `77.19298_segments.csv` / `.rds` | Detected IBD segments. |
| `77.19298_pairwise_ibd_fraction.csv` and `_pairwise_ibd_fraction_matrix.csv` | Pairwise genome fractions used for network thresholds and the heatmap. |
| `77.19298_snp_ibd_proportion.csv` | **Per-SNP proportion of sample pairs IBD**, a different quantity from the pairwise genome fraction. |
| `thr05/`, `thr20/`, `thr40/`, `thr60/`, `thr80/` | Edge, node, component, and travel-mixing tables; RDS/GraphML networks; network PDFs. |
| `77.19298_threshold_network_summary.csv` and `_threshold_travel_mixing.csv` | Cross-threshold network and travel summaries. |
| `77.19298_threshold_node_metrics.csv` and `_threshold_component_summary.csv` | Sample and component details across thresholds. |
| `77.19298_pairwise_ibd_heatmap.pdf` and `_threshold_overview.pdf` | Pairwise sharing and threshold comparison figures. |
| `77.19298_run_settings.txt` and `_session_info.txt` | Recorded settings and R package session. |

The supplied short IBD script used isolate missingness 0.52 and genotyping error 0.01, whereas the longer threshold script used 0.40 and 0.001. **This version uses the longer threshold script's values.** Confirm them against the actual analysis log and manuscript before describing a new run as a reproduction of published results. The network layout is illustrative; pairwise fractions are in the tables.
