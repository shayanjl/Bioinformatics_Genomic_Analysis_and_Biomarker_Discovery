# Gastric Cancer Biomarker Discovery Pipeline

A comprehensive pipeline for the discovery and validation of gastric cancer‑specific circulating miRNA biomarkers, integrating transcriptomic data from public repositories with functional enrichment and machine learning approaches.

## Table of Contents

* [Project Overview](#project-overview)
* [Repository Structure](#repository-structure)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Data Sources](#data-sources)
* [Usage](#usage)

  * [1. Microarray Analysis](#1-microarray-analysis)
  * [2. TCGA-STAD Data Analysis](#2-tcga-stad-data-analysis)
  * [3. TCGA-ESCA Data Analysis](#3-tcga-esca-data-analysis)
  * [4. Functional Enrichment & Visualization](#4-functional-enrichment--visualization)
* [Results](#results)
* [Manuscript](#manuscript)
* [Contributing](#contributing)
* [License](#license)

## Project Overview

Late diagnosis of gastric cancer (GC) leads to poor survival rates. This project develops a non‑invasive miRNA biomarker panel by:

1. Performing differential expression analysis on gastric cancer tissue (TCGA‑STAD) and blood miRNA (GEO) datasets.
2. Identifying GC‑specific pathways shared between tissue‑derived DEGs and circulating miRNA targets.
3. Selecting mechanistically relevant miRNAs and refining them via machine learning feature selection.
4. Validating predictive performance on external datasets.

The pipeline is implemented in R and produces publication‑quality figures and tables.

## Repository Structure

```
├── data/                          # Raw and processed data files
│   ├── raw/                       # Downloaded GEO and TCGA data
│   └── processed/                 # Filtered & normalized expression matrices
├── scripts/                       # Analysis scripts
│   ├── MicroarrayAnalysis_ShayanJL.R   # GEO miRNA microarray analysis
│   ├── TCGA_Data_Analysis_Using_TCGAbiolinks_and_DEA_ShayanJL.R  # TCGA‑STAD analysis
│   ├── Download_to_DEGs_by_TCGAbiolinks_COMPLETE_FILE.R           # TCGA‑ESCA analysis
│   └── 103.R                      # Functional enrichment barplot by gene ratio
├── results/                       # Generated figures and tables
│   ├── box.plot.pdf
│   ├── Cor.heatmap.pdf
│   ├── PCA.pdf
│   ├── VolcanoPlot.pdf
│   ├── Heatmap.DEGs.pdf
│   └── 20-103_by_gene_ratio_categorized_padj.jpg
├── manuscript/                    # Draft manuscript and supplementary files
│   └── After_second_revision.docx
└── README.md                      # This file
```

## Prerequisites

* **R** (version ≥ 4.0)
* **Rtools** (on Windows)
* Internet access (to download TCGA/GEO data)

### R Packages

The following R packages are required (each script installs missing packages automatically):

* **Data retrieval & preprocessing**: `GEOquery`, `BiocManager`, `TCGAbiolinks`, `SummarizedExperiment`, `EDASeq`, `edgeR`, `limma`
* **Visualization**: `ggplot2`, `pheatmap`, `gplots`, `EnhancedVolcano`, `forcats`, `scales`
* **Data manipulation**: `dplyr`, `readxl`

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/<username>/momeni-mentor2.git
   cd momeni-mentor2
   ```
2. Open R or RStudio and set the working directory to the project root.
3. The first time you run each script, required packages will be installed automatically.

## Data Sources

* **TCGA‑STAD** (Stomach Adenocarcinoma) transcriptome profiling via `TCGAbiolinks`.
* **TCGA‑ESCA** (Esophageal Carcinoma) transcriptome profiling via `TCGAbiolinks`.
* **GEO** miRNA microarray dataset `GSE106817` for gastric cancer serum miRNAs.
* **External validation** blood miRNA dataset `GSE164174`.

## Usage

### 1. Microarray Analysis

Performs QC, normalization, differential expression, and visualization on GEO miRNA data (`GSE106817`).

```bash
# In R or terminal:
Rscript scripts/MicroarrayAnalysis_ShayanJL.R
```

**Outputs**:

* Boxplots (`results/box.plot.pdf`)
* Correlation heatmap (`results/Cor.heatmap.pdf`)
* PCA (`results/PCA.pdf`)
* Volcano plot (`results/VolcanoPlot.pdf`)
* Heatmap of top 50 DEGs (`results/Heatmap.DEGs.pdf`)
* Lists of up‐ and downregulated miRNAs in `results/cancer.up.genes.txt` and `cancer.down.genes.txt`.

### 2. TCGA‑STAD Data Analysis

Downloads TCGA‑STAD RNA‑Seq data, preprocesses, and runs differential expression using `edgeR` and `limma`.

```bash
Rscript scripts/TCGA_Data_Analysis_Using_TCGAbiolinks_and_DEA_ShayanJL.R
```

**Outputs**:

* Raw and preprocessed RData files
* CSV tables of DEGs (`dataDEGs*.csv`)

### 3. TCGA‑ESCA Data Analysis

Analogous pipeline for TCGA‑ESCA, including normalization, filtering, DEA, plotting, clinical data extraction, and WGCNA preprocessing.

```bash
Rscript scripts/Download_to_DEGs_by_TCGAbiolinks_COMPLETE_FILE.R
```

**Outputs**:

* DEGs CSVs for multiple thresholds
* Volcano and enrichment plots
* Prepared clinical tables for downstream analysis

### 4. Functional Enrichment & Visualization

Generates a barplot of functional enrichment results by gene ratio for top pathways (from external Excel file).

```bash
Rscript scripts/103.R
```

**Inputs**:

* `data/processed/103_top20_smallest_padj.xlsx` (place the enrichment results Excel here)

**Outputs**:

* Barplot image: `results/20-103_by_gene_ratio_categorized_padj.jpg`

## Results

All figures, heatmaps, ROC curves, and DEG tables generated by the pipeline are stored in the `results/` directory. Key outputs include:

* **DEG heatmaps & volcano plots** for GEO and TCGA datasets.
* **Functional enrichment barplots** highlighting gene‐ratio categories.
* **ROC curves** and performance metrics for selected miRNA biomarker panels (see manuscript Figure 3).

## Manuscript

The draft manuscript (**After\_second\_revision.docx**) contains:

* Introduction, Materials & Methods, Results, Discussion
* Supplementary tables (S1–S9) and figures (Venn diagrams, ROC curves).

Please refer to this file for detailed study findings and interpretation.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for:

* Bug fixes or improvements to scripts
* Additional dataset integration
* Upstream/downstream analysis modules (e.g., WGCNA, survival analysis)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

*For questions or feedback, contact Shayan Jalali.*
