# 1. GLOBAL CONFIGURATION (SLURM Strict Mode) ----------------------------------
cat(paste(.libPaths()))
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
promoter_file  <- p$promoter_file

reg_features <- strsplit(as.character(p$reg_features), ",")[[1]]

# Analysis Params
n_top_percent  <- as.numeric(p$top_percent)
upstream_bp    <- as.numeric(p$upstream_bp)
downstream_bp  <- as.numeric(p$downstream_bp)


# Logging for SLURM .out file
cat("--- Configuration Loaded ---\n")
cat(paste0("Comparison:  ", ref_ct, " vs ", contrast_ct, "\n"))
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

  # --- zellkonverter + basilisk ---
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

get_results <- function(dds, column, ct1, ct2, alpha = 0.05, important_genes = NULL, dotsize = 1, textsize = 3, return = FALSE, only_significant = TRUE) {
  res <- results(dds, contrast = c(column, ct1, ct2), alpha = alpha)
  res <- res[!is.na(res$padj),]
  if (only_significant) res <- res[res$padj < 0.05, ]
  if (return) return(res)
  
  res$is_promoter <- if (!is.null(important_genes)) rownames(res) %in% important_genes else FALSE
  
  p <- ggplot(transform(as.data.frame(res), gene = rownames(res), 
                        type = ifelse(padj < 0.05 & log2FoldChange > 2, "Up", 
                                      ifelse(padj < 0.05 & log2FoldChange < -2, "Down", "NS"))),
              aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = type), size = dotsize)
  
  if (!is.null(important_genes)) {
    p <- p + geom_point(data = subset(transform(as.data.frame(res), gene = rownames(res), is_promoter = rownames(res) %in% important_genes), is_promoter), color = "red", size = dotsize)
  }
  
  p + geom_text_repel(data = function(d) subset(d, is_promoter), aes(label = gene), color = "red", size = textsize, max.overlaps = Inf) +
    geom_text_repel(data = function(d) {top10 = d[order(-d$log2FoldChange), ][1:10, ]; subset(top10, !is_promoter)}, aes(label = gene), color = "darkcyan", size = textsize, max.overlaps = Inf) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") + 
    scale_color_manual(values = c("Down" = "darkgoldenrod2", "Up" = "darkcyan", "NS" = "gray20")) + 
    labs(x = "log2 Fold Change", y = "-log10(p-adj)", title = paste(ct1, "vs", ct2)) +
    theme(legend.position = "none", axis.title = element_text(size = 14), plot.title = element_text(size = 16, face = "bold"))
}

arrange_res <- function(res, gene_annot, crop_df = TRUE, top = 5, order_by = 'log2FoldChange', only_significant = FALSE){
  res <- res[order(res[[order_by]], decreasing = TRUE), ]
  if (only_significant) res <- res[res$padj < 0.05,]
  res_df <- data.frame(Gene_symbol = rownames(res), log2FoldChange = res$log2FoldChange, padj = res$padj, baseMean = res$baseMean)
  if (!is.null(gene_annot)) res_df <- res_df %>% left_join(gene_annot %>% dplyr::select(ensembl_gene_id, external_gene_name, description), by = c("Gene_symbol" = "external_gene_name"))
  if (crop_df) res_df <- head(res_df, nrow(res_df)*top/100)
  return(res_df)
}

top_percent_genes <- function(df, crop_df = TRUE, top = 5, order_by = "log2FoldChange", only_pathways = FALSE, only_ensembl = FALSE, bottom = FALSE) {
  if (bottom) df <- df[order(df[[order_by]], decreasing = FALSE), ] else df <- df[order(df[[order_by]], decreasing = TRUE), ]
  df <- df[df$padj < 0.05,]; if (only_ensembl) df <- df[df$ensembl_gene_id != "",]
  if (only_pathways) df <- df[!is.na(df$pathways),]; if (crop_df) df <- head(df, nrow(df)*top/100)
  return(df)
}

make_fgsea_stats <- function(df, gene_col = "Gene_symbol", logfc_col = "log2FoldChange") {
  stats <- df |> dplyr::select(all_of(c(gene_col, logfc_col))) |> dplyr::filter(!is.na(.data[[gene_col]]), !is.na(.data[[logfc_col]])) |>
    dplyr::distinct(.data[[gene_col]], .keep_all = TRUE) |> tibble::deframe()
  return(sort(stats, decreasing = TRUE))
}

top_pathways_plot <- function(res_df, stats){
  topPathwaysUp <- res_df[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- res_df[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(reactome[topPathways], stats, res_df, gseaParam=0.5)
}

pathways_into_df <- function(gene_df, res_fgsea){
  res_fgsea_tbl <- as_tibble(res_fgsea) 
  gene_to_pathway <- res_fgsea_tbl %>% dplyr::select(pathway, leadingEdge) %>% tidyr::unnest(cols = leadingEdge) %>% dplyr::rename(gene = leadingEdge)
  gene_to_pathway_collapsed <- gene_to_pathway %>% group_by(gene) %>% summarise(pathways = paste(pathway, collapse = ", "))
  return(gene_df %>% left_join(gene_to_pathway_collapsed, by = c("Gene_symbol" = "gene")))
}

add_mean_expression <- function(df, dds_obj, group_var, ref_celltype, contrast_celltype) {
  norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
  if (is.null(colnames(norm_counts))) colnames(norm_counts) <- rownames(colData(dds_obj))
  meta <- as.data.frame(colData(dds_obj))
  ref_samples <- rownames(meta[meta[[group_var]] == ref_celltype, , drop = FALSE])
  contrast_samples <- rownames(meta[meta[[group_var]] == contrast_celltype, , drop = FALSE])
  mean_ref <- rowMeans(norm_counts[, ref_samples, drop = FALSE])
  mean_contrast <- rowMeans(norm_counts[, contrast_samples, drop = FALSE])
  mean_df <- data.frame(Gene_symbol = rownames(norm_counts), meanExp_ref = mean_ref, meanExp_contrast = mean_contrast, diff_meanExp = mean_ref - mean_contrast, check.names = FALSE)
  names(mean_df)[2:3] <- c(paste0("meanExp_", ref_celltype), paste0("meanExp_", contrast_celltype))
  return(left_join(df, mean_df, by = "Gene_symbol"))
}



get_feature_coords <- function(gene_df, id_column = "Gene_symbol", upstream_bp = 1000, downstream_bp = 100, 
                               gene_annot_df = gene_annot, reg_annot_df = feature_annot, 
                               features = c("Promoter", "Enhancer")) {
  match_col <- if (id_column == "Gene_symbol") "external_gene_name" else "ensembl_gene_id"
  results_list <- list()
  
  for (i in seq_len(nrow(gene_df))) {
    gene_id <- gene_df[[id_column]][i]
    message("Processing: ", gene_id)
    gene_info <- gene_annot_df[gene_annot_df[[match_col]] == gene_id, ]
    
    if (nrow(gene_info) == 0) {message("Gene not found locally: ", gene_id)
      next
    }

    gene_info <- gene_info[1, ]
    chr <- gene_info$chromosome_name
    strand <- gene_info$strand

    if (is.na(strand) || is.na(gene_info$start_position) || is.na(gene_info$end_position)) {
      message("Skipping ", gene_id, ": Incomplete coordinate or strand data.")
      next
    }

    if (strand == 1) { r_start <- gene_info$start_position - upstream_bp; r_end <- gene_info$start_position + downstream_bp } else { r_start <- gene_info$end_position - downstream_bp; r_end <- gene_info$end_position + upstream_bp }
    
    message("Looking for regulatory features in chr", chr, ":", r_start, "-", r_end)

    reg_elements <- reg_annot_df[
      reg_annot_df$chromosome_name == chr & 
      reg_annot_df$chromosome_start <= r_end & 
      reg_annot_df$chromosome_end >= r_start &
      reg_annot_df$feature_type_name %in% features, 
    ]
    
    if (nrow(reg_elements) == 0) next

    reg_elements$length <- reg_elements$chromosome_end - reg_elements$chromosome_start + 1
    
    results_list[[gene_id]] <- cbind(
      gene_df[rep(i, nrow(reg_elements)), , drop = FALSE], 
      feature_type = reg_elements$feature_type_name, 
      regulatory_id = reg_elements$regulatory_stable_id, 
      chromosome_name = reg_elements$chromosome_name, 
      element_start = reg_elements$chromosome_start, 
      element_end = reg_elements$chromosome_end, 
      element_length = reg_elements$length
    )
  }
  
  return(do.call(rbind, results_list))
}

get_feature_fasta <- function(df) {
  all_fastas <- c(); for (i in 1:nrow(df)) {
    id <- df$regulatory_id[i]; gene_id <- df$Gene_symbol[i]
    region <- paste0(df$chromosome_name[i], ":", df$element_start[i], "..", df$element_end[i])
    url <- paste0("https://rest.ensembl.org/sequence/region/human/", region)
    response <- GET(url, content_type("text/x-fasta"))
    if (status_code(response) == 200) {
      parts <- strsplit(content(response, as = "text", encoding = "UTF-8"), "\n")[[1]]
      new_header <- sub(">", paste0(">", id, " | gene:", gene_id, " | "), parts[1])
      all_fastas <- c(all_fastas, paste0(new_header, "\n", paste(parts[-1], collapse = ""), "\n"))
    }
    Sys.sleep(0.1)
  }; return(paste(all_fastas, collapse = ""))
}

get_feature_seqs <- function(df) {
  df$Sequence <- NA_character_; for (i in 1:nrow(df)) {
    region <- paste0(df$chromosome_name[i], ":", df$element_start[i], "..", df$element_end[i])
    url <- paste0("https://rest.ensembl.org/sequence/region/human/", region)
    response <- GET(url, content_type("text/x-fasta"))
    if (status_code(response) == 200) {
      parts <- strsplit(content(response, as = "text", encoding = "UTF-8"), "\n")[[1]]
      df$Sequence[i] <- paste(parts[-1], collapse = "")
    }
    Sys.sleep(0.1)
  }; return(df)
}


# 4. DATA LOADING & PRE-PROCESSING ---------------------------------------------
cat("--- STEP 3: DATA LOADING AND PRE-PROCESSING ---\n")
setwd(work_directory)

cat("Reading processed sce...\n")
sce <- readRDS(file.path("Outputs", paste0(dataset, "_processed_sce.rds")))
summed <- aggregateAcrossCells(sce, id = colData(sce)[, c(ct_column, batch_column)])
rownames(summed) <- rowData(summed)$SYMBOL

# Skipping if celltype not found
if (!contrast_ct %in% sce[[ct_column]]) {
    message(
        "\n[WARNING] ", contrast_ct, " not found in dataset ", dataset, ". Skipping to next dataset...\n"
    )
    quit(save = "no", status = 0)
}

# 1. Handle Gene Annotations
if (!is.null(ensembl_gene_path)) {
  cat("Loading existing gene annotations from CSV...\n")
  gene_annot <- read.csv(ensembl_gene_path, row.names = 1, check.names = FALSE)
  
} else {
  cat("Connecting to gene BioMART...\n")
  ensembl <- safe_useEnsembl(versions = ensembl_vers)
  gene_annot <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"), mart = ensembl)
  
  cat("Saving gene annotations to disk as CSV...\n")
  write.csv(gene_annot, file = file.path("/scratch/lab_mguell/dfragoso/gene_annot.csv"), row.names = TRUE)
}


# 2. Handle Regulatory Annotations
if (!is.null(ensembl_reg_path)) {
  cat("Loading existing regulatory annotations from CSV...\n")
  feature_annot <- read.csv(ensembl_reg_path, row.names = 1, check.names = FALSE)
  
} else {
  cat("Connecting to regulatory BioMART...\n")
  ensembl_reg <- safe_useEnsembl(versions = ensembl_vers, biomart = "regulation", dataset = "hsapiens_regulatory_feature")

  chrs <- c(1:22, "X", "Y", "MT")
  feature_list <- list()

  cat("Starting feature download...\n")

  for (chr in chrs) {
    cat(paste0("--- Fetching Chromosome: ", chr, " --- \n"))
    
    # Try the download for this specific chromosome
    dat <- tryCatch({
      getBM(
        attributes = c('chromosome_name', 'chromosome_start', 'chromosome_end', 
                      'feature_type_name', 'regulatory_stable_id'), 
        filters = 'chromosome_name',
        values = chr,
        mart = ensembl_reg
      )
    }, error = function(e) {
      message(paste0("Failed to download Chromosome ", chr, ": ", e$message))
      return(NULL)
    })
    
    if (!is.null(dat)) {
      feature_list[[chr]] <- dat
      cat(paste0("Successfully retrieved ", nrow(dat), " features for Chr ", chr, "\n"))
    }
  }

  # Combine all pieces into one table
  cat("Merging all chromosomes...\n")
  feature_annot <- do.call(rbind, feature_list)

  # Check if we actually got data
  if (nrow(feature_annot) > 0) {
    cat(paste0("Total features retrieved: ", nrow(feature_annot), "\n"))
    cat("Saving to CSV...\n")
    write.csv(feature_annot, file = "/scratch/lab_mguell/dfragoso/feature_annot.csv", row.names = TRUE)
  } else {
    stop("No data was retrieved. Check your BioMart connection.")
  }
}


if (!is.null(promoter_file)) {
  cat("Loading specific genes list...\n")
  specific_prom <- read_xlsx(promoter_file)
  specific_prom <- specific_prom[!is.na(specific_prom$`PROMOTER ID`),]
  promoter_genes <- specific_prom$`GENE SYMBOL`
}


cat("--- STEP 4: STARTING DESeq2 ---\n")
dds <- DESeqDataSet(summed, design = as.formula(paste("~", ct_column)))
dds[[ct_column]] <- relevel(dds[[ct_column]], ref = ref_ct)
new_ids <- make.unique(paste(colData(dds)[[batch_column]], colData(dds)[[ct_column]]))
rownames(colData(dds)) <- new_ids; colnames(dds) <- new_ids

dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)


cat("--- STEP 5: EXTRACTING RESULTS AND PATHWAY ANALYSIS ---\n")
cat("Getting results...\n")
T_vs_cct <- get_results(dds, ct_column, ref_ct, contrast_ct, return = TRUE, only_significant = FALSE)
T_vs_cct_volcano <- get_results(dds, ct_column, ref_ct, contrast_ct, only_significant = FALSE, important_genes = promoter_genes)

T_vs_cct_full_df <- T_vs_cct_df <- arrange_res(T_vs_cct, gene_annot, crop_df = FALSE, only_significant = FALSE)
T_vs_cct_df <- arrange_res(T_vs_cct, gene_annot, crop_df = FALSE, only_significant = TRUE)
T_vs_cct_stats <- make_fgsea_stats(T_vs_cct_df)

cat("Performing FGSEA...\n")
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>% split(x = .$gene_symbol, f = .$gs_name)

resfgsea <- fgsea(pathways = reactome, stats = T_vs_cct_stats)
top_pathways_plt <- top_pathways_plot(resfgsea, T_vs_cct_stats)

cat("Arranging dataframes...\n")
T_vs_cct_df <- pathways_into_df(T_vs_cct_df, resfgsea)
T_vs_cct_df <- add_mean_expression(T_vs_cct_df, dds, ct_column, ref_ct, contrast_ct)

T_vs_cct_full_df <- pathways_into_df(T_vs_cct_full_df, resfgsea)
T_vs_cct_full_df <- add_mean_expression(T_vs_cct_full_df, dds, ct_column, ref_ct, contrast_ct)

# 9. TOP/BOTTOM PROMOTERS ------------------------------------------------------
cat("--- STEP 6: FETCHING PROMOTER SEQUENCES ---\n")
T_vs_cct_top5 <- top_percent_genes(T_vs_cct_df, top = n_top_percent, only_ensembl = TRUE)
T_vs_cct_bot5 <- top_percent_genes(T_vs_cct_df, top = n_top_percent, only_ensembl = TRUE, bottom = TRUE)

cat("Exctracting top coordinates...\n")
prom_coords_top5 <- get_feature_coords(T_vs_cct_top5, upstream_bp = upstream_bp, downstream_bp = downstream_bp)
cat("Exctracting bot coordinates...\n")
prom_coords_bot5 <- get_feature_coords(T_vs_cct_bot5, upstream_bp = upstream_bp, downstream_bp = downstream_bp)

cat("Exctracting fasta seqs...\n")
top_seqs_fasta <- get_feature_fasta(prom_coords_top5)
bot_seqs_fasta <- get_feature_fasta(prom_coords_bot5)

cat("Exctracting seqs in df...\n")
top_promoters_df <- get_feature_seqs(prom_coords_top5)
bot_promoters_df <- get_feature_seqs(prom_coords_bot5)

# 10. EXPORT --------------------------------------------------------------------
cat("--- STEP 7: EXPORTING RESULTS ---\n \n")
dir.create("Outputs/DataFrames", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/Sequences", recursive = TRUE, showWarnings = FALSE)

writexl::write_xlsx(T_vs_cct_df, paste0("Outputs/DataFrames/", dataset, "_df_", ref_ct, "_vs_", contrast_ct, ".xlsx"))
writexl::write_xlsx(T_vs_cct_full_df, paste0("Outputs/DataFrames/", dataset, "_full_df_", ref_ct, "_vs_", contrast_ct, ".xlsx"))
writexl::write_xlsx(T_vs_cct_top5, paste0("Outputs/DataFrames/", dataset, "_top_gen_", ref_ct, "_vs_", contrast_ct, ".xlsx"))
writexl::write_xlsx(T_vs_cct_bot5, paste0("Outputs/DataFrames/", dataset, "_bot_gen_", ref_ct, "_vs_", contrast_ct, ".xlsx"))
writexl::write_xlsx(top_promoters_df, paste0("Outputs/DataFrames/", dataset, "_top_features_", ref_ct, "_vs_", contrast_ct, ".xlsx"))
writexl::write_xlsx(bot_promoters_df, paste0("Outputs/DataFrames/", dataset, "_bot_features_", ref_ct, "_vs_", contrast_ct, ".xlsx"))

ggsave(paste0("Outputs/Figures/", dataset, "_", ref_ct, "_vs_", contrast_ct, "_volcano.png"), T_vs_cct_volcano)
ggsave(paste0("Outputs/Figures/", dataset, "_", ref_ct, "_vs_", contrast_ct, "_top_pathways.png"), top_pathways_plt)

writeLines(top_seqs_fasta, paste0("Outputs/Sequences/", dataset, "_", ref_ct, "_vs_", contrast_ct, "_top_promoters.fasta"))
writeLines(bot_seqs_fasta, paste0("Outputs/Sequences/", dataset, "_", ref_ct, "_vs_", contrast_ct, "_bot_promoters.fasta"))

# Saving in a general folder
setwd("..")
setwd("..")

dir.create("Results/DataFrames/top_features", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/DataFrames/bot_features", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Figures/Volcano", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Sequences", recursive = TRUE, showWarnings = FALSE)

ggsave(paste0("Results/Figures/Volcano/", dataset, "_", ref_ct, "_vs_", contrast_ct, "_volcano.png"), T_vs_cct_volcano)

write.csv(top_promoters_df, 
          file = paste0("Results/DataFrames/top_features/", dataset, "_top_features_", ref_ct, "_vs_", contrast_ct, ".csv"), 
          row.names = FALSE)
write.csv(bot_promoters_df, 
          file = paste0("Results/DataFrames/bot_features/", dataset, "_bot_features_", ref_ct, "_vs_", contrast_ct, ".csv"), 
          row.names = FALSE)