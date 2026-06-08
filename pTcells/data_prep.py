import argparse
import os
import scanpy as sc
import anndata as ad
import re
import json


# 1. Correct Argument Parsing
parser = argparse.ArgumentParser()
parser.add_argument("--work_dir", type=str, required=True)
parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--raw_file", type=str, required=True)
parser.add_argument("--ct_column", type=str, required=True)
# Use action='store_true'. If the flag is present, it's True. If absent, it's False.
parser.add_argument("--filter_cts", type=str, default="", help="Comma-separated cell types to keep") 
parser.add_argument("--merge_map", type=str, default="{}")
args = parser.parse_args()

work_dir = args.work_dir
input_file = args.input_file
filter_cts = args.filter_cts
raw_file = args.raw_file
ct_column = args.ct_column

print("Reading .h5ad file...", flush=True)
adata = sc.read_h5ad(os.path.join(work_dir, input_file))

# Filter cts of interest
if args.filter_cts:
    print("Filtering for specific lineages...", flush=True)
    
    target_col = ct_column 
    
    # 1. Split the comma-separated string into a list and strip whitespace
    cts_list = [ct.strip() for ct in args.filter_cts.split(",") if ct.strip()]
    
    # 2. Join them with a pipe '|' to create a regex OR pattern
    regex_pattern = "|".join(cts_list).replace(" ", " ?")
    
    # 3. Apply the dynamic filter
    mask = adata.obs[target_col].str.contains(regex_pattern, case=False, na=False)
    adata = adata[mask].copy()

    remaining = sorted(adata.obs[target_col].unique())
    print(f"Remaining cell types: {', '.join(remaining)}", flush=True)
 

# Merge cts
if args.merge_map:
    print("Merging cell type groups...", flush=True)
    
    # 1. Parse the Slurm JSON string
    slurm_dict = json.loads(args.merge_map)

    # 2. Invert the dict and force all incoming Slurm keys to lowercase
    mapping_dict = {}
    for group, names_list in slurm_dict.items():
        # If it's a list, iterate; if it's a string (backwards compatibility), split
        items = names_list if isinstance(names_list, list) else names_list.split(",")
        for name in items:
            mapping_dict[name.strip().lower()] = group
    
    # # 2. Invert the dict and force all incoming Slurm keys to lowercase
    # mapping_dict = {}
    # for group, old_names in slurm_dict.items():
    #     clean_names = old_names.replace("\n", "")
    #     for name in clean_names.split(","):
    #         mapping_dict[name.strip().lower()] = group

    # 3. Get the original cell types as strings
    old_series = adata.obs[ct_column].astype(str)
    
    # 4. Map using a lowercase version of the data
    merged_series = old_series.str.lower().map(mapping_dict)
    
    # Replace NaNs (unmapped cell types) with "other"
    merged_series = merged_series.fillna("other")
    
    adata.obs[ct_column] = merged_series.astype("category")

    # 5. Print Summary
    print("Summary of merged cell types:", flush=True)
    # Include "other" in the summary printout
    all_groups = list(slurm_dict.keys()) + ["other"]
    
    for group in all_groups:
        mapped_old_types = old_series[adata.obs[ct_column] == group].unique()
        if len(mapped_old_types) > 0:
            print(f"  {group} <- {', '.join(sorted(mapped_old_types))}", flush=True)


# Extract raw counts
if adata.raw is not None:

	print("Exctracting raw counts...")
	
	adata_raw = ad.AnnData(
		X = adata.raw.X,
		var = adata.raw.var,
		obs = adata.obs
	)
	raw_counts = adata_raw.X

	print("Raw matrix shape:", raw_counts.shape)
	print("Raw matrix type:", type(raw_counts))
	print("Raw matrix dtype:", raw_counts.dtype)

else:
	print("raw counts not found")
	adata_raw = None

from scipy.sparse import csc_matrix



# Convert X to CSC format (standard for R)
print("Converting matrix to CSC...")
adata_raw.X = csc_matrix(adata_raw.X)

# Clean up 'uns' (Unstructured metadata often crashes R readers)
adata_raw.uns = {} 

# Save output
if adata_raw is not None:
	out_path = os.path.join(work_dir,  raw_file)
	adata_raw.write_h5ad(out_path)

	
	print(f"Saved raw counts to: {out_path}\n")