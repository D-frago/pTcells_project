import pandas as pd
import os
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Merge, filter, and deduplicate scRNA-seq feature tables.")
    parser.add_argument("--base_dir", type=str, required=True)
    parser.add_argument("--datasets", nargs='+', required=True)
    parser.add_argument("--ref_ct", type=str, required=True)
    parser.add_argument("--contrast_cts", nargs='+', required=True)
    return parser.parse_args()

def merge_files(args, feature_type):
    dfs = []
    for dataset in args.datasets:
        for contrast_ct in args.contrast_cts:
            contrast_label = f"{args.ref_ct}_vs_{contrast_ct}"
            
            file_path = os.path.join(
                args.base_dir, 
                dataset, 
                "Outputs/DataFrames/", 
                f"{dataset}_{feature_type}_features_{contrast_label}.xlsx"
            )
            
            if os.path.exists(file_path):
                # Reading the Excel file
                temp_df = pd.read_excel(file_path)
                # Adding metadata columns
                temp_df['dataset'] = dataset
                temp_df['contrast'] = contrast_label
                dfs.append(temp_df)
            else:
                print(f"File not found: {file_path}")

    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

def finalize_dataframe(df, name):
    """Subsets specific columns and ensures unique ensembl_gene_id."""
    if df.empty:
        return df
    
    cols_to_keep = ["Gene_symbol", "ensembl_gene_id", "Sequence"]
    # Filter only for columns that actually exist in the dataframe
    available_cols = [c for c in cols_to_keep if c in df.columns]
    
    # 1. Select Columns | 2. Drop duplicates based on ensembl_gene_id
    new_df = df[available_cols].copy()
    new_df = new_df.drop_duplicates(subset=['ensembl_gene_id'], keep='first')
    
    print(f"Deduplication ({name}): Kept {len(new_df)} unique IDs.")
    return new_df

if __name__ == "__main__":
    args = get_args()

    # --- 1. Load the two Full DataFrames ---
    df_full_top = merge_files(args, "top")
    df_full_bot = merge_files(args, "bot")

    # --- 2. Cross-Filter Logic ---
    # Remove rows from the TOP dataframe if the ID appears anywhere in the BOTTOM dataframe
    if not df_full_top.empty and not df_full_bot.empty:
        bot_ids_blacklist = df_full_bot['ensembl_gene_id'].unique()
        df_full_top = df_full_top[~df_full_top['ensembl_gene_id'].isin(bot_ids_blacklist)]
        print(f"Scrubbed IDs found in 'bot' from 'top' dataframe.")

    # --- 3. Create the Cleaned/Unique DataFrames ---
    df_unique_top = finalize_dataframe(df_full_top, "Top Features Cleaned")
    df_unique_bot = finalize_dataframe(df_full_bot, "Bottom Features Cleaned")

    # --- 4. Exporting All Four Files ---
    
    # Save the full versions (all columns, all occurrences after cross-filtering)
    df_full_top.to_csv("all_top_features.csv", index=False)
    df_full_bot.to_csv("all_bot_features.csv", index=False)
    
    # Save the cleaned versions (selected columns, unique IDs)
    df_unique_top.to_csv("all_top_features_UNIQUE.csv", index=False)
    df_unique_bot.to_csv("all_bot_features_UNIQUE.csv", index=False)
    
    print("\nProcess finished. Created:")
    print("- all_top_features.csv / all_bot_features.csv")
    print("- all_top_features_UNIQUE.csv / all_bot_features_UNIQUE.csv")