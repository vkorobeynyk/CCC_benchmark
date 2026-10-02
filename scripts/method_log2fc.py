import liana as li
import scanpy as sc
import pandas as pd
import numpy as np
from liana.method import logfc

# =============================================================================
# 1. Load Input Data
# =============================================================================
# Paths are provided via Snakemake wildcards and parameters
counts_path = snakemake.input["inflated_normalized_counts"]
metadata_path = snakemake.input["metadata_processed"]
db_path = snakemake.params["LR_database"]

# Load normalized counts and prepare AnnData object
counts_df = pd.read_csv(counts_path, sep="\t").T
adata = sc.AnnData(counts_df)
adata.raw = adata.copy()

# Integrate cell type metadata
metadata = pd.read_csv(metadata_path, sep="\t")
adata.obs["Celltype"] = metadata["Celltype"].values

# Load the curated Ligand-Receptor database using space separator
lr_database = pd.read_csv(db_path, sep=" ")

# =============================================================================
# 2. Run Log2FC (via LIANA)
# =============================================================================
results = logfc(
    adata, 
    groupby='Celltype', 
    resource=lr_database, 
    expr_prop=0, 
    min_cells=0, 
    inplace=False, 
    verbose=True
)

# =============================================================================
# 3. Post-processing and NaN/Infinite Value Handling
# =============================================================================
# Replace infinite values with NaN and remove them to ensure valid statistics
results.replace([np.inf, -np.inf], np.nan, inplace=True)
results.dropna(inplace=True)

# Filter for the specific Sender-to-Receiver interactions used in the benchmark
results = results[(results["source"] == "Sender") & (results["target"] == "Receiver")]

# Format the final output dataframe
output_df = pd.DataFrame({
    "ligand_receptor": results["ligand_complex"] + "_" + results["receptor_complex"],
    "statistics": results["lr_logfc"]
})

# Define significance (Log2FC > 1, representing a ~2x fold increase)
output_df["significant"] = output_df["statistics"] > 1

# Sort results: higher Log2FC values (stronger signals) appear first
output_df = output_df.sort_values("statistics", ascending=False)

# =============================================================================
# 4. Save Results
# =============================================================================
output_path = snakemake.output["significant_interactions"]
output_df.to_csv(output_path, sep="\t", index=False)
