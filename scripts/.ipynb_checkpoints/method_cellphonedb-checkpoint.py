import liana as li
import scanpy as sc
import pandas as pd
# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean, rank_aggregate, scseqcomm

#############
### INPUT ###
#############
inflated_normalized_counts_path = snakemake.input["inflated_normalized_counts"]
cellmetadata_path = snakemake.input["metadata_processed"]

adata = sc.AnnData(pd.read_csv(inflated_normalized_counts_path, sep="\t").T)
adata.raw = adata.copy()
cm = pd.read_csv(cellmetadata_path, sep="\t")
adata.obs["Celltype"] = cm["Celltype"]

##############
### OUTPUT ###
##############
significant_interactions_path = snakemake.output["significant_interactions"]

##############
### Params ###
##############
LR_database_path = snakemake.params["LR_database"]
LR_database = pd.read_csv(LR_database_path, sep=" ")

# run cellphonedb
cpdb_results = cellphonedb(adata,
            groupby='Celltype', 
            resource=LR_database,
            expr_prop=0,
            min_cells = 0,
            inplace = False,
            verbose=True, key_added='cpdb_res')

LRdata_df = pd.DataFrame({"ligand_receptor" : cpdb_results["ligand"] + "_" + cpdb_results["receptor"], "pval": cpdb_results["cellphone_pvals"]})
LRdata_df["significant"] = LRdata_df["pval"] < 0.05
LRdata_df = LRdata_df.rename({"pval":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics")

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
