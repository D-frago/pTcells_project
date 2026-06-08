# 1. GLOBAL CONFIGURATION (SLURM Strict Mode) ----------------------------------
options(timeout = 1200)
options(rhdf5.bit64conversion = "double")
# cat(paste(.libPaths()))
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- gsub("^--", "", args[i])
    val <- args[i + 1]
    out[[key]] <- val
  }
  out
}

p <- parse_args(args)


# Extract variables
dataset        <- p$dataset
work_directory <- p$work_dir
raw_file       <- p$raw_file
ensembl_gene_path <- p$ensembl_gen
ensembl_reg_path <- p$ensembl_reg
ref_ct         <- p$ref_ct
contrast_ct    <- p$contrast_ct
ct_column      <- p$ct_column
batch_column <- p$batch_column

ensembl_vers <- as.integer(strsplit(p$ensembl_vers, ",")[[1]])
reg_features <- strsplit(as.character(p$reg_features), ",")[[1]]

# QC & Filtering
nmads_qc       <- as.numeric(p$nmads_qc)
only_normal    <- as.logical(p$only_normal)



# Logging for SLURM .out file
cat("--- Configuration Loaded ---\n")
cat(paste0("Working Dir: ", work_directory, "\n"))
cat(paste0("Input File:  ", raw_file, "\n"))
cat("---------------------------\n")



# 2. SETUP AND LIBRARIES -------------------------------------------------------
suppressPackageStartupMessages({
  library(BiocManager)
  library(BiocParallel)
})

# Use all available cores for installation
bp <- MulticoreParam(workers = parallel::detectCores())
register(bp)

# List of required packages
pkgs <- c(
  # --- Bioconductor core infrastructure ---
  "BiocGenerics",  "S4Vectors",  "IRanges",  "GenomeInfoDb",  "SummarizedExperiment",

  # --- Single-cell core ---
  "SingleCellExperiment",  "scuttle",  "scater",  "scran",  "DropletUtils", "batchelor",

  # --- Bulk RNA-seq ---
  "DESeq2",

  # --- Annotation & enrichment ---
  "reactome.db",  "msigdbr",  "biomaRt",  "fgsea",  "GEOquery",

  # --- zellkonverter & anndataR ---
  "zellkonverter", 

  # --- CRAN packages (safe to install after Bioc) ---
  "ggplot2",  "ggrepel",  "RColorBrewer",  "hexbin",  "tidyr",  "dplyr",  "factoextra",  "cowplot",  "readxl",  "R.utils",  "umap",  "car",  "gridExtra",  "data.table",  "pheatmap",  "httr"
)


# Install missing packages
missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, Ncpus = parallel::detectCores())
} else {
  message("All packages already installed.")
}

cat("--- STEP 1: LOADING LIBRARIES ---\n")
suppressPackageStartupMessages({
  # clea; 
  library(DESeq2); library(vsn)
  library(ggplot2); library(ggrepel); library(RColorBrewer)
  library(hexbin); library(tidyr); library(dplyr)
  library(factoextra); library(BiocParallel); library(cowplot)
  library(readxl); library(R.utils); library(umap)
  library(car); library(gridExtra); library(fgsea)
  library(data.table); library(msigdbr); library(reactome.db)
  library(biomaRt); library(httr); library(GEOquery)
  library(scuttle); library(scater); library(scran)
  library(scDblFinder); library(DropletUtils); library(SingleCellExperiment)
  library(batchelor); library(edgeR); library(pheatmap); library(zellkonverter)
})


# 3. CUSTOM FUNCTIONS ----------------------------------------------------------
cat("--- STEP 2: DEFINING HELPER FUNCTIONS ---\n")

safe_useEnsembl <- function(versions, biomart = "genes", dataset = "hsapiens_gene_ensembl") {
  for (v in versions) {
    message("Trying Ensembl version: ", v)
    out <- try(useEnsembl(biomart = biomart, dataset = dataset, version = v), silent = TRUE)
    if (!inherits(out, "try-error")) {
      message("Success with version ", v); return(out)
    }
  }
  stop("All versions failed")
}

# 4. DATA LOADING & PRE-PROCESSING ---------------------------------------------
cat("--- STEP 3: DATA LOADING AND PRE-PROCESSING ---\n")
setwd(work_directory)

cat("Reading H5AD file...\n")
sce <- readH5AD(raw_file)
assay(sce, "counts") <- assay(sce, "X")
print(sce)

# Handle Gene Annotations
if (!is.null(ensembl_gene_path)) {
  cat("Loading gene annotations from CSV...\n")
  gene_annot <- read.csv(ensembl_gene_path, row.names = 1, check.names = FALSE)
  
} else {
  cat("Connecting to gene BioMART...\n")
  ensembl <- safe_useEnsembl(versions = ensembl_vers)
  gene_annot <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), mart = ensembl)
  
  cat("Saving gene annotations to disk as CSV...\n")
  write.csv(gene_annot, file = file.path("/scratch/lab_mguell/dfragoso/gene_annot.csv"), row.names = TRUE)
}




# Removing non-normal samples from diseased individuals
if (only_normal && "disease" %in% colnames(colData(sce))) {
  cat("Filtering for normal samples...\n")
  sce <- sce[, sce$disease == "normal"]
}



rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- gene_annot$external_gene_name[match(rowData(sce)$ENSEMBL, gene_annot$ensembl_gene_id)]
rowData(sce)$SYMBOL[rowData(sce)$SYMBOL == "" | is.na(rowData(sce)$SYMBOL)] <- rowData(sce)$ENSEMBL [rowData(sce)$SYMBOL == "" | is.na(rowData(sce)$SYMBOL)]
rownames(sce) <- rowData(sce)$SYMBOL



# 5. QUALITY CONTROL -----------------------------------------------------------
cat("--- STEP 4: QUALITY CONTROL ---\n")


is.mito <- which(gene_annot$chromosome_name[match(rowData(sce)$ENSEMBL, gene_annot$ensembl_gene_id)] == "MT")
df_qc <- perCellQCMetrics(sce, subsets = list(Mito = is.mito))
colData(sce) <- cbind(colData(sce), df_qc)

reasons <- perCellQCFilters(df_qc, sub.fields = "subsets_Mito_percent", nmads = nmads_qc)
sce$discard <- reasons$discard

MitoPercent <- plotColData(sce, y = "subsets_Mito_percent", colour_by = "discard")
MitoPercent_sum <- plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "discard")
sce <- sce[, !sce$discard]

# 6. NORMALIZATION & DIM RED ---------------------------------------------------
cat("--- STEP 5: NORMALIZATION AND DIMENSION REDUCTION ---\n")
sizeFactors(sce) <- librarySizeFactors(sce)
sce <- logNormCounts(sce)

dec.sce <- modelGeneVar(sce)
hvg.sce.var <- getTopHVGs(dec.sce, n = 1000)
sce <- runPCA(sce, subset_row = hvg.sce.var)
sce <- runUMAP(sce, dimred = "PCA")

cat("Detecting doublets...\n")
dbl.dens <- computeDoubletDensity(sce, subset.row = hvg.sce.var, dims = ncol(reducedDim(sce)))
sce$DoubletScore <- dbl.dens
sce$doublet <- doubletThresholding(data.frame(score = dbl.dens), method = "griffiths", returnType = "call")
DoubletScore <- plotColData(sce, y = "DoubletScore", colour_by = "doublet")
DoubletDetected_sum <- plotColData(sce, "detected", "sum", colour_by = "doublet")
sce <- sce[, sce$doublet != "doublet"]

# 7. BATCH CORRECTION & DESEQ2 -------------------------------------------------
cat("--- STEP 6: BATCH CORRECTION AND DIMRED ---\n")
cat("Selecting hvgs...\n")
dec <- modelGeneVar(sce, block = sce[[batch_column]])
chosen.hvgs <- dec$bio > 0
sce <- runPCA(sce, subset_row = chosen.hvgs, ntop = 1000, ncomponents = 50)


set.seed(10102)
cat("")
merged <- correctExperiments(sce, batch = sce[[batch_column]], subset.row = chosen.hvgs, PARAM = FastMnnParam(d = 50))
merged <- runUMAP(merged, dimred = "corrected")
UMAP <- plotUMAP(sce, colour_by = ct_column) + 
                 scale_fill_brewer(palette = "Spectral")
UMAP_corrected <- plotUMAP(merged, colour_by = ct_column) + 
                  scale_fill_brewer(palette = "Spectral")
PCA_corrected <- plotPCA(merged, dimred = "PCA", colour_by = ct_column)



# Saving results
dir.create("Outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/Figures", recursive = TRUE, showWarnings = FALSE)


cat("Saving processed sce object...\n")
# Save sce object in outputs folder in working outputs
saveRDS(merged, file = file.path("Outputs", paste0(dataset, "_processed_sce.rds")))

cat("Saving figures...\n \n")
ggsave(paste0("Outputs/Figures/", dataset, "_MitoPercent.png"), MitoPercent)
ggsave(paste0("Outputs/Figures/", dataset, "_MitoPercent_sum.png"), MitoPercent_sum)
ggsave(paste0("Outputs/Figures/", dataset, "_DoubletScore.png"), DoubletScore)
ggsave(paste0("Outputs/Figures/", dataset, "_DoubletDetected_sum.png"), DoubletDetected_sum)
ggsave(paste0("Outputs/Figures/", dataset, "_UMAP.png"), UMAP)
ggsave(paste0("Outputs/Figures/", dataset, "_UMAP_corrected.png"), UMAP_corrected)
ggsave(paste0("Outputs/Figures/", dataset, "_PCA_corrected.png"), PCA_corrected)

# Create extra results folder
# setwd("..")
# setwd("..")
# dir.create("Results/Figures/UMAP", recursive = TRUE, showWarnings = FALSE)
# ggsave(paste0("Results/Figures/UMAP/", dataset, "_UMAP_corrected.png"), UMAP_corrected)
