import liana as li
import scanpy as sc
import pandas as pd
# import all individual methods
from liana.method import scseqcomm

#############
### INPUT ###
#############
inflated_normalized_counts_path = snakemake.input["inflated_normalized_counts"]
cellmetadata_path = snakemake.input["metadata_processed"]

adata = sc.AnnData(pd.read_csv(inflated_normalized_counts_path, sep="\t").T)
adata.raw = adata.copy()
cm = pd.read_csv(cellmetadata_path, sep="\t")
adata.obs["Celltype"] = cm["Celltype"].values

##############
### OUTPUT ###
##############
significant_interactions_path = snakemake.output["significant_interactions"]

##############
### Params ###
##############
LR_database_path = snakemake.params["LR_database"]
LR_database = pd.read_csv(LR_database_path, sep=" ")

# run 
x = scseqcomm(adata,
            groupby='Celltype', 
            resource=LR_database,
            expr_prop=0,
            min_cells = 0,
            inplace = False,
            verbose=True)
            
# select only CT1-CT2 interaction
x = x[(x["source"] == "CT1") & (x["target"] == "CT2")]
LRdata_df = pd.DataFrame({"ligand_receptor" : x["ligand"] + "_" + x["receptor"], "inter_score": x["inter_score"]})
LRdata_df["significant"] = LRdata_df["inter_score"] > 0.8
LRdata_df = LRdata_df.rename({"inter_score":"statistics"},axis=1)
LRdata_df = LRdata_df.sort_values("statistics", ascending=False) # higher vaue, the more significant the interaction

LRdata_df.to_csv(significant_interactions_path, sep = "\t")
