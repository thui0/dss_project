# Epithelial-specific loss of Smad4 in the Colon Enhances the Wound Healing Response and Alleviates the Fibrotic Response in an Acute DSS Mouse Model

This repository contains the code used for the RNA-seq analysis and figure generation in R, as described in the corresponding publication.

## Repository Structure
`data/`: Includes RNA-seq counts data and the data used for boxplots and heatmaps
`code/`: Contains R scripts and Rmd Notebooks for data processing, differential gene expression analysis, and visualization.
`code/figures/`: R scripts to reproduce the figures in the manuscript.
`figures/`: Contains the outputs generated from `code/figures/` scripts.
`results/`: Outputs from analyses, including normalized expression matrices, DESeq2 results, and enrichment analyses.

## Requirements
- The analysis was performed using R (version 4.3.1) with the following key R packages:
  - `DESeq2`
  - `clusterProfiler`
  - `ggplot2`
  - `ComplexHeatmap`
For a more complete list of package dependencies, use `renv::dependencies()` in the project directory.
- GSEA was performed using the desktop application (version 4.3.2).

Note: R may run out of memory when running the `code/03_ora.R` script. To avoid this, it is recommended to run the script manually in smaller chunks.

## Data Availability
Raw RNA-seq data is not included in this repository. The raw data can be downloaded from GEO: [GSE252864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252864).
