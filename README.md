# pTcells: Cell-specific regulatory feature extraction from public scRNA-seq data

## Table of contents
* [Biological framework](#biological-framework)
* [Pipeline overview](#pipeline-overview)
* [Installation](#installation)
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

### 1. Job Orchestration (`ptcells.sh`)
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
Upon successful completion of the workflow, the pipeline yields two primary categories of results within an Outputs/ folder:
* **Visualizations & Diagnostics:** UMAP & QC plots for data inspection, alongside Volcano plots to evaluate the statistical significance of the DEA results.
* **Target Genomic Sequences:** The final extracted regulatory DNA sequences can delivered in standard `.fa` (FASTA) format for a more compressed result or in an `.xlsx` spreadsheet, along with more relevant metrics about the analysis.

---

## Installation
pTcells is a pipeline developed using versions of Python/3.11.3 and R/4.5.1. To install pTcells in your local computer execute the following commands:
```bash
git clone [https://github.com/D-frago/pTcells.git]
cd pTcells
```
To ensure the portability of pTcells, both a python environment and an R library are provided with all the packages needed to execute the pipeline without running into dependencies errors.

## Tutorial

### Prerequisites
This pipeline is optimized for Python 3.11.3 and R 4.5.1; alternative versions may trigger dependency issues. Furthermore, given the significant memory overhead required to process large-scale scRNA-seq datasets (>$10^6$ cells), executing this pipeline on a High-Performance Computing (HPC) cluster is strongly advised.

### Configuration Variables

The pipeline script (`pTcells_pipeline.sh`) contains a dedicated **USER CONFIG** section. Below are the variables you can modify to adapt the pipeline to different datasets or biological questions.

#### 1. Dataset & File Paths
| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `DATASET` | `"Edgar_2025"` | The name identifier for your current run. Used for naming output files. |
| `WORK_DIR` | `"/scratch/.../${DATASET}"` | The absolute path to the directory where input files sit and results will be saved. |
| `INPUT_FILE` | `"${DATASET}.h5ad"` | Name of the primary input AnnData file containing processed/filtered data. |
| `RAW_FILE` | `"${DATASET}_raw_counts.h5ad"` | Name of the file containing raw, unnormalized counts (required for differential expression). |
| `PROM_FILE` | `"/scratch/.../List of...xlsx"` | Absolute path to the Excel sheet detailing known tissue-specific promoters. |

#### 2. Biological Metadata & Comparisons
| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `REF_CT` | `"Tcell"` | The **target/reference** cell type you are investigating. |
| `CONTRAST_CTS` | `("Bcell" "hepatocyte" "Kupffer_cells")` | A Bash array of cell types you want to compare *individually* against your `REF_CT`. |
| `CT_COLUMN` | `"cell_type"` | The column name in the `.obs` metadata slot of your AnnData object that contains cell labels. |
| `BATCH_COLUMN` | `"donor_id"` | The column name in metadata representing batch effects (e.g., patient ID, sample pool, or sequencing lane). |
| `REG_FEATURES`| `"Promoter,Enhancer"` | Comma-separated genomic features you wish to extract and analyze for your target genes. |

#### 3. Genome Annotation (Ensembl)
| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `ENSEMBL_GEN` | `"/scratch/.../gene_annot.csv"` | Path to the local reference CSV mapping gene IDs to names. |
| `ENSEMBL_REG` | `"/scratch/.../feature_annot.csv"`| Path to local regulatory feature annotations (promoters/enhancers). |
| `ENSEMBL_VERS`| `"114,115,113"` | Comma-separated priority list of Ensembl database versions to look up. |

#### 4. Quality Control & Cell Filtering
| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `NMADS_QC` | `"3"` | Number of Median Absolute Deviations (MADs) allowed for adaptive filtering thresholds. |
| `ONLY_NORMAL` | `"TRUE"` | If set to `TRUE`, the pipeline filters out tumor/disease-state cells if that metadata exists. |
| `FILTER_CTS` | `"T cell, B cell, Hepa..."` | Comma-separated keywords. Only cells whose annotation has these strings in them will pass the initial python script. |
| `MERGE_MAP` | `'{"Tcell": [ "T cell",  "naive T cell", ...]}` | A JSON-formatted mapping dictionary used to aggregate highly specifically labeled cells into clean, high-level groups (e.g., merging "naive CD4+ T cell" into "Tcell"). |

#### 5. Downstream Analysis Tuning
| Variable | Default Value | Description |
| :--- | :--- | :--- |
| `TOP_PERCENT` | `"5"` | The top percentage cutoff of Differentially Expressed (DE) genes to extract regulatory features for. |
| `UPSTREAM_BP` | `"2000"` | Distance (in base pairs) upstream of the Transcription Start Site (TSS) to consider when searching for regulatory features. |
| `DOWNSTREAM_BP`| `"100"` | Distance (in base pairs) downstream of the TSS to consider when searching for regulatory features. |


### Basic execution example

To run the pipeline on your High-Performance Computing (HPC) cluster, follow these steps. 

#### 1. Prepare Your Data
Ensure your dataset files are placed inside your directory of choice and match the structure expected by the script. For the default example `Edgar_2025`:

```bash
# Check that your input files are in your dataset folder
ls /scratch/.../scRNA-seq_datasets/Edgar_2025/
# Expected output: Edgar_2025.h5ad
```

#### 2. Adapt variables to your dataset
Before submitting the job, make sure that the variables from the **USER CONFIG** are adapted to your dataset. The most essential variables include:
* DATASET: Name of the dataset you're analysisng
* WORK_DIR: Directory where the input file is and where the outputs will be generated.
* REF_CT: Choose the cell type of interest.
* CONTRAST_CTS: Choose the cell type(s) that you wan't to contrast the reference cell type against.

It is also heavily recommended to adapt:
* FILTER_CTS: Filter cells that are not of interest for a clearer UMAP. It also helps reduce the dataset's size which can make the pipeline faster.
* MERGE_MAP: Really useful in case you want to merge two cell types into a single group.
_Note: if the no merge map is included, the ref and contrast cell types must be annotated exactly as in the input dataset's metadata._ 

#### Run 
Navigate to your repository directory where pTcells.sh is located and submit the script to the scheduler using sbatch:

```bash
cd /.../pTcells_repo
sbatch pTcells.sh
```


## Bibliography

CZ CELLxGENE Discover: A single-cell data platform for scalable exploration, analysis and modeling of aggregated data CZI Single-Cell Biology, et al. bioRxiv 2023.10.30; doi: https://doi.org/10.1101/2023.10.30.563174

CELLxGENE: a performant, scalable exploration platform for high dimensional sparse matrices. CZI Single-Cell Biology, et al. bioRxiv 2021.04.05; doi: https://doi.org/10.1101/2021.04.05.438318

Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, Waldron L, Pagès H, Smith ML, Huber W, Morgan M, Gottardo R, Hicks SC. Orchestrating single-cell analysis with Bioconductor. Nature Methods, 2020. doi: 10.1038/s41592-019-0654-x

Deconinck L, Zappia L, Cannoodt R, Morgan M, scverse core, Virshup I, Sang-aram C, Bredikhin D, Seurinck R, Saeys Y (2025). “anndataR improves interoperability between R and Python in single-cell transcriptomics.” bioRxiv, 2025.08.18.669052. doi:10.1101/2025.08.18.669052.

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). Fast gene set enrichment analysis. bioRxiv. https://doi.org/10.1101/060012

Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191. doi:10.1038/nprot.2009.97.
