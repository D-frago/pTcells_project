# pTcells: Cell-specific regulatory feature extraction from public scRNA-seq data

## Table of contents
* [Biological framework](#biological-framework)
* [Pipeline overview](#pipeline-overview)
* [Installation and execution](#installation-and-execution)
* [Tutorial](#tutorial)
* [Bibliography](#bibliography)

---

## Biological framework
**pTcells** allows the systematic analysis of scRNA-seq datasets obtained from public repositories, optimized specifically for the [CellxGene](https://cellxgene.cziscience.com/) database. It is designed to generate high-quality training data consisting of DNA sequences from regulatory elements, which can be used for the fine-tuning of foundational DNA Language Models like **EVO2**.

This pipeline uses Differential Expression Analysis (DEA) tools and the Ensembl database to extract regulatory elements from genes that are overexpressed in T cells. By contrasting T cells against another cell type of interest (e.g., B cells), we can identify genes specific to T cells. These genes are often associated with unique regulatory elements that drive T cell-specific gene expression. By analyzing the genomic regions flanking or inside these genes, the pipeline identifies and extracts potential regulatory elements (such as promoters or enhancers) that dictate T cell identity.

---

## Pipeline overview

![Pipeline Overview](Pipeline%20Overview.png)

The pipeline is split into three core functional layers:

### 1. Job Orchestration (`ptcells_pipeline.sh`)
* **SLURM Integration:** Encompasses both the raw data extraction parameters and environment variable settings. It manages computational resource allocation, handling initialization and argument passing to sequentially execute the downstream modules without manual intervention.

### 2. Data Extraction (`data_prep.py`)
* **Input:** Accepts scRNA-seq data in the `.h5ad` format (e.g., from CellxGene).
* **Pre-processing:** Filters out samples of non-interest cell types and merges cell types into groups of interest.
* **Processing:** Parses the AnnData object structure (`obs`, `var`, and `.X` matrix) to isolate and export the un-normalized, raw count data (`raw.h5ad`), ensuring clean downstream statistical analysis.

### 3. Analysis & quality control (`scRNAseq_QC.R`)
* **Pre-processing & Quality Control:** Filters out low-quality data based on the percentage of mitochondrial genes, normalizes the data, detects doublets, removes batch effects, and performs dimensionality reduction.

### 4. Downstream Processing (`scRNAseq_downstream.R`)
* **Differential Expression Analysis (DEA):** Employs **DESeq2** to find statistically significant, cell-specific differentially expressed genes, which are further cross-validated using Fast Gene Set Enrichment Analysis (**FGSEA**).
* **Sequence Extraction (BioMART):** Uses Ensembl's BioMART API to scan genomic regions associated with the target genes, locate annotated regulatory features, and extract the precise nucleotide sequences.

### Expected outputs
Upon successful completion of the workflow, the pipeline yields two primary categories of results:
* **Visualizations & Diagnostics:** UMAP & QC plots for data inspection, alongside Volcano plots to evaluate the statistical significance of the DEA results.
* **Target Genomic Sequences:** The final extracted regulatory DNA sequences can delivered in standard `.fa` (FASTA) format for a more compressed result or in an `.xlsx` spreadsheet, along with more relevant metrics about the analysis.

---

## Installation and execution

### Prerequisites

### Installation
pTcells is a pipeline developed using versions of Python/3.11.3 and R/4.5.1. To install pTcells in your local computer execute the following commands:
```bash
git clone [https://github.com/mariaartlle/pTcells_repo/David/pTcells.git]([https://github.com/mariaartlle/pTcells_repo/David/pTcells.git)
cd pTcells
```

## Tutorial
