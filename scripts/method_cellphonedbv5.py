import liana as li
import scanpy as sc
import pandas as pd
from liana.method import cellphonedb

# =============================================================================
# 1. Load Input Data
# =============================================================================
# Paths are provided by the Snakemake workflow wildcards
counts_path = snakemake.input["inflated_normalized_counts"]
metadata_path = snakemake.input["metadata_processed"]
db_path = snakemake.params["LR_database"]

# Load normalized counts and create AnnData object
counts_df = pd.read_csv(counts_path, sep="\t").T
adata = sc.AnnData(counts_df)
adata.raw = adata.copy()

# Add cell type metadata
metadata = pd.read_csv(metadata_path, sep="\t")
adata.obs["Celltype"] = metadata["Celltype"].values

# Load the curated Ligand-Receptor database
lr_database = pd.read_csv(db_path, sep=" ")

# =============================================================================
# 2. Run CellPhoneDB (via LIANA)
# =============================================================================
results = cellphonedb(
  adata, 
  groupby='Celltype', 
  resource=lr_database, 
  expr_prop=0, 
  min_cells=0, 
  inplace=False, 
  verbose=True
)

# =============================================================================
# 3. Post-processing and Formatting
# =============================================================================
# Filter for specific Sender-Receiver interactions defined in the simulation
results = results[(results["source"] == "Sender") & (results["target"] == "Receiver")]

# Format output dataframe
output_df = pd.DataFrame({
  "ligand_receptor": results["ligand_complex"] + "_" + results["receptor_complex"],
  "statistics": results["cellphone_pvals"]
})

# Define significance (p-value < 0.05)
output_df["significant"] = output_df["statistics"] < 0.05

# Sort results: lower p-values (more significant) appear first
output_df = output_df.sort_values("statistics", ascending=True)

# =============================================================================
# 4. Save Results
# =============================================================================
output_path = snakemake.output["significant_interactions"]
output_df.to_csv(output_path, sep="\t", index=False)
