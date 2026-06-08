#!/bin/bash -l

#SBATCH --job-name=multi-pTcells
#SBATCH --output=multi-pTcells.out
#SBATCH --error=multi-pTcells.err
#SBATCH --time=24:00:00
#SBATCH --partition=haswell
#SBATCH --mem=512G

set -euo pipefail


########################
# USER CONFIG
########################

BASE_DIR="/scratch/lab_mguell/dfragoso/scRNA-seq_datasets"

# <<< CHOOSE DATASETS HERE >>>
DATASETS=("Xu_2023_blood" "Xu_2023_liver" "Edgar_2025" "Edgar_2025_2" "MacParland_2018" "Suo_2022" "Trebo_2025" "Wells_2025" "Fachrul_2026" "Kock_2025" "Gong_2025")

# Promoter file 
PROM_FILE="/scratch/lab_mguell/dfragoso/RNA-seq_datasets/List of tissue-specific promoters.xlsx"

# Biology / metadata
REF_CT="Tcell"          # <-- set to your reference cell type
CONTRAST_CTS=("Bcell" "hepatocyte" "Kupffer_cells" "myeloid")    # <-- set to your contrast cell types
CT_COLUMN="cell_type"      # <-- column in metadata with cell types
BATCH_COLUMN="donor_id"     # <-- column in metadata with patient IDs
REG_FEATURES="Promoter,Enhancer"

# Ensembl connection
ENSEMBL_GEN="/scratch/lab_mguell/dfragoso/gene_annot.csv"
ENSEMBL_REG="/scratch/lab_mguell/dfragoso/feature_annot.csv"
ENSEMBL_VERS="114,115,113"

# QC & filtering
NMADS_QC="3"           # <-- higher = more strict 
ONLY_NORMAL="TRUE"
FILTER_CTS="T cell, B cell, Hepa, liver, Kupffer, monocyte, macroph, dendri, neutr, granu, myeloid" # <-- key words you want the remaining cells to have
MERGE_MAP='{
  "Tcell": [
    "T cell", "alpha-beta T cell", "naive T cell", "CD4+ T cell", "CD8+ T cell", "memory regulatory T cell",
    "naive regulatory T cell", "naive thymus-derived CD4-positive, alpha-beta T cell",
    "T lymphocyte", "CD4 T", "CD8 T", "regulatory T cell", "CD4-positive, alpha-beta T cell", 
    "CD4-positive alpha-beta memory T cell", "central memory CD4-positive alpha-beta T cell", 
    "effector memory CD4-positive alpha-beta T cell", "naive thymus-derived CD4-positive alpha-beta T cell", 
    "CD4-positive, alpha-beta cytotoxic T cell", "central memory CD4-positive, alpha-beta T cell",
	  "effector memory CD4-positive, alpha-beta T cell", "CD4-positive, alpha-beta memory T cell",
	  "CD8-positive, alpha-beta T cell", "central memory CD8-positive, alpha-beta T cell",
    "CD8-positive alpha-beta memory T cell", "central memory CD8-positive alpha-beta T cell", 
    "effector memory CD8-positive alpha-beta T cell", "naive thymus-derived CD8-positive alpha-beta T cell", 
    "CD8-positive, alpha-beta cytotoxic T cell", "effector memory CD8-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta memory T cell", "CD8-positive, alpha-beta regulatory T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell"
  ],

  "Bcell": [
    "B cell", "B lymphocyte", "naive B cell", "memory B cell", "fraction A pre-pro B cell", 
    "pro-B cell", "B-2 B cell", "B-1 B cell", "unswitched memory B cell", 
    "germinal center B cell", "immature B cell", "mature B cell", "transitional stage B cell",
	"class switched memory B cell"
  ],

  "hepatocyte": [
    "Hepatocyte", "hepatocytes", "mature hepatocyte", "midzonal region hepatocyte", 
    "centrilobular region hepatocyte", "periportal region hepatocyte"
  ],

  "Kupffer_cells": [
    "Kup", "Kupffer", "Kupffer cell", "liver macrophage"
  ],

  "other_Tcells": [
    "gamma-delta T cell", "gdT", "gd T cell", "invariant natural killer T", "iNKT", 
    "NKT cell", "mature gamma-delta T cell", "mucosal invariant T cell", "MAIT",
    "effector memory CD8-positive alpha-beta T cell terminally differentiated",
    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
    "mucosal invariant T cell", "memory CCR4-positive regulatory T cell",
    "effector memory CD4-positive, alpha-beta T cell, terminally differentiated",
    "CD8-alpha-alpha-positive, alpha-beta intraepithelial T cell"
  ],

  "other_liver_cells": [
    "LSEC", "liver sinusoidal endothelial cell", "hepatic stellate cell", "HSC", 
    "cholangiocyte", "biliary epithelial cell", "endothelial cell of pericentral hepatic sinusoid", 
    "endothelial cell of periportal hepatic sinusoid", "endothelial cell of hepatic sinusoid"
  ],

  "myeloid": [
    "monocyte", "macrophage", "dendritic cell", "plasmacytoid dendritic cell", "DC", 
    "classical monocyte", "non-classical monocyte", "neutrophil", "granulocyte", 
    "myeloid cell", "inflammatory macrophage", "cycling myeloid cell", 
    "conventional dendritic cell", "dendritic cell, human", "Axl+ dendritic cell, human",
    "CD14-low, CD16-positive monocyte", "CD14-positive, CD16-negative classical monocyte",
    "common myeloid progenitor", "intermediate monocyte", "mast cell", "Axl+ dendritic cell",
    "granulocyte monocyte progenitor cell", "pre-conventional dendritic cell"
  ]
}'


# Analysis params
TOP_PERCENT="5"     # <-- percent of DE genes whose reg features will be exctracted
UPSTREAM_BP="2000"  # <-- upstream region from TSS
DOWNSTREAM_BP="100" # <-- downstream region from TSS



# Redirect logs for each dataset
exec > "${BASE_DIR}/multi-pTcells.out"
exec 2> "${BASE_DIR}/multi-pTcells.err"

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job started on $(hostname)"

########################
# ENVIRONMENT
########################

ml purge
ml modulepath/haswell
ml R/4.5.1-gfbf-2023a
ml Python/3.11.3-GCCcore-12.3.0

source /homes/users/dfragoso/scratch/scRNA-seq_datasets/.pTcells/bin/activate


################################
# LOOP THROUGH SELECTED DATASETS
################################

for DATASET in "${DATASETS[@]}"; do
    WORK_DIR="${BASE_DIR}/${DATASET}"
    INPUT_FILE="${DATASET}.h5ad"
    RAW_FILE="${DATASET}_raw_counts.h5ad"

    echo ""
    echo ""
    echo "############################################################"
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Processing dataset: $DATASET"
    echo "############################################################"
	  echo ""

    #######################
    # STEP 1: Preprocessing
    #######################
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running data_prep.py"

    python data_prep.py \
	  --work_dir "$WORK_DIR" \
	  --input_file "$INPUT_FILE" \
	  --raw_file "$RAW_FILE" \
	  --ct_column "$CT_COLUMN" \
	  --filter_cts "$FILTER_CTS" \
	  --merge_map "$MERGE_MAP"

    ########################
    # STEP 2: Analysis & QC
    ########################

    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running scRNAseq_QC.R"

	Rscript scRNAseq_QC.R \
	  --dataset "$DATASET" \
	  --work_dir "$WORK_DIR" \
	  --raw_file "$RAW_FILE" \
	  --ct_column "$CT_COLUMN" \
  	--batch_column "$BATCH_COLUMN" \
	  --ensembl_gen "$ENSEMBL_GEN" \
	  --ensembl_reg "$ENSEMBL_REG" \
	  --ensembl_vers $ENSEMBL_VERS \
	  --reg_features "$REG_FEATURES" \
	  --nmads_qc $NMADS_QC \
	  --only_normal $ONLY_NORMAL 

	########################
    # STEP 3: DEA & seq extraction
    ########################

    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running scRNAseq_downstream.R"

	for CELLTYPE in "${CONTRAST_CTS[@]}"; do
		CONTRAST_CT="${CELLTYPE}"

		echo "-----------------------------------------------------------------"
		echo "[$(date +"%Y-%m-%d %H:%M:%S")] Processing $REF_CT vs $CONTRAST_CT"
		echo "-----------------------------------------------------------------"

		Rscript scRNAseq_downstream.R \
		  --dataset "$DATASET" \
		  --work_dir "$WORK_DIR" \
		  --raw_file "$RAW_FILE" \
		  --ct_column "$CT_COLUMN" \
		  --batch_column "$BATCH_COLUMN" \
		  --ensembl_gen "$ENSEMBL_GEN" \
		  --ensembl_reg "$ENSEMBL_REG" \
	  	--ensembl_vers $ENSEMBL_VERS \
		  --ref_ct "$REF_CT" \
		  --contrast_ct "$CONTRAST_CT" \
		  --reg_features "$REG_FEATURES" \
		  --promoter_file "$PROM_FILE" \
		  --top_percent $TOP_PERCENT \
		  --upstream_bp $UPSTREAM_BP \
		  --downstream_bp $DOWNSTREAM_BP


	done
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Finished dataset: $DATASET"
done

echo "[$(date +"%Y-%m-%d %H:%M:%S")] All selected datasets completed"

##############################
# Merge all analyised datasets
##############################
echo ""
echo "---------------------------------------------------"
echo "[$(date +"%Y-%m-%d %H:%M:%S")] Running merge_dfs.py"
echo "---------------------------------------------------"
echo ""

python merge_dfs.py \
    --base_dir "$BASE_DIR" \
    --datasets "${DATASETS[@]}" \
    --ref_ct "$REF_CT" \
    --contrast_cts "${CONTRAST_CTS[@]}"

echo "[$(date +"%Y-%m-%d %H:%M:%S")] pTcells pipeline completed successfully"
