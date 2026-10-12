    # IBD sample heatmap

    This script creates the sample-level pairwise identity-by-descent (IBD) heatmap used to visualize IBD sharing between *Plasmodium vivax* isolates. It is a downstream plotting script: it does not re-estimate IBD segments, but uses the pairwise IBD proportions produced by the IBD analysis.

    In the manuscript workflow, this script belongs after the `isoRelate` IBD estimation described in the Methods section. It supports the IBD Results section by showing how pairwise IBD sharing is distributed across samples, communities, and traveler-associated infections.

    ## Script

    ```text
    ibd_sample_heatmap.R
    ```

    ## Input files

    The script uses the annotated pairwise IBD dataset and the cleaned sample metadata:

    ```text
    /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/03_pairwise_ibd_dataset/77.19298_pairwise_ibd_annotated.csv
    /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/02_metadata_clean/metadata_clean.csv
    ```

    The pairwise IBD file must contain:

    ```text
    sample1, sample2, ibd_prop
    ```

    The metadata file must contain:

    ```text
    sample_id, community, travel_status
    ```

    `ibd_prop` is the pairwise proportion of the genome inferred to be IBD between two isolates. Values range from `0` to `1`, with larger values indicating greater IBD sharing.

    ## Method

    The script first reads the annotated pairwise IBD file and the cleaned metadata file. The pairwise file contains one row per sample pair, while the metadata file provides the community and travel-status information used to annotate the heatmap.

    The heatmap matrix is then built by converting the pairwise IBD table into a symmetric sample-by-sample matrix. Each pairwise value is added in both directions, so that `sample1` versus `sample2` and `sample2` versus `sample1` have the same IBD value. The diagonal is set to `1`, representing each sample compared with itself.

    Samples are ordered by community and then by travel status. The community order used in the script is:

    ```text
    Libertad, Gamitanacocha, UrcoMirano
    ```

    Within each community, samples are ordered as:

    ```text
    Non-traveler, Traveler
    ```

    This ordering makes the heatmap directly interpretable for the manuscript results: blocks along the diagonal show within-community IBD sharing, off-diagonal blocks show between-community IBD sharing, and the travel-status symbols show whether traveler-associated infections cluster separately or are embedded among other infections from the same communities.

    The plot uses a continuous IBD color scale from `0` to `1`. The legend includes `0.50` as a reference value, matching the IBD heatmap interpretation used in the manuscript. Community boundaries are shown with separator lines, and community labels are added to the top and left sides of the heatmap.

    Travel status is shown with symbols on the top and left margins of the heatmap:

    ```text
    circle   = Non-traveler
    triangle = Traveler
    ```

    The script also writes the exact plotting dataset used to create the heatmap. This file is useful for checking the sample order, IBD values, and metadata annotations used in the final figure.

    ## Outputs

    Outputs are written to:

    ```text
    /scratch/antwerpen/208/vsc20843/WGS/results/gatk/pv/IBD/ibd_pairwise_relatedness_analysis/12_final_selected_figures/
    ```

    The main outputs are:

    ```text
    final_sample_level_ibd_heatmap_metadata_ordered_plot_data.csv
    final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.pdf
    final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.png
    final_sample_level_ibd_heatmap_metadata_ordered_annotated_corrected.tiff
    ```

    The PDF is used for vector-quality inspection, while the PNG and TIFF files are publication-ready raster versions of the same heatmap.

    ## Recommended repository location

    Keep this script and README in a subfolder of the IBD analysis:

    ```text
    population_genetics/IBD/ibd_sample_heatmap/
    ```

    Recommended folder contents:

    ```text
    population_genetics/IBD/ibd_sample_heatmap/ibd_sample_heatmap.R
    population_genetics/IBD/ibd_sample_heatmap/README.md
    ```
