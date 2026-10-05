library(Seurat)
library(magrittr)
library(dplyr)
library(tibble)

# -------------------------- #
# 10x download from https://zenodo.org/records/13893852/files/seurat_Chromium_All.rds?download=1 file = "seurat_Chromium_All.rds"
# -------------------------- #

SO  = readRDS("/home/vkorob/Downloads/seurat_Chromium_All.rds")

# remove 75% of cells 
set.seed(1)
Idents(SO) = "Celltype" 
sampled_cells = sample(Cells(SO), size = floor(ncol(SO) * 0.20))
SO = subset(x = SO, cells = sampled_cells)

# replace astrocyte with "Sender" and microglia with "Receiver"
SO$Celltype %<>% gsub("Astrocyte", "Sender",. ) %>% gsub("Endothelial", "Receiver",. )

meta = SO@meta.data %>% 
  select(Celltype) %>% 
  rownames_to_column("cell_ID")

counts = GetAssayData(SO, "RNA" , "counts")

# Filter to decrease size
counts = counts[rowSums(counts != 0) > 10,]
counts = counts[rowSums(counts) > 200,]

write.table(meta,"data/10x/raw_metadata.tsv") # save metadata
write.table(counts,"data/10x/raw_counts.tsv") # save counts

# -------------------------- #
# SMARTseq2 download from https://github.com/vkorobeynyk/scRNA-seq_analysis_jessberger_lab/tree/master/data file = "SO_Gli_Ascl_comb_int.Robj"
# -------------------------- #
load("/home/vkorob/Documents/data_scripts/RNA_seq/GliAscl1_jessberger/SO_Gli_Ascl_comb_int.Robj")

meta = pbmc.combined@meta.data
meta$cluster %<>% gsub("ndNSC", "Sender",. ) %>% gsub("dNSC", "Receiver",. )

meta = data.frame(cell_ID = rownames(meta), Celltype = meta$cluster)

write.table(meta,"data/SMARTseq2/raw_metadata.tsv") # save metadata
write.table(GetAssayData(pbmc.combined, "RNA", "counts"),"data/SMARTseq2/raw_counts.tsv") # save counts

# -------------------------- #
# 10x Immune (GSE195848)
# -------------------------- #
SO = readRDS("GSE195848_Seurat_object.RDS")

SOb1 = subset(SO, batch == "Batch1")

set.seed(1)
Idents(SOb1) = "cell_type" 
sampled_cells = sample(Cells(SOb1), size = floor(ncol(SOb1) * 0.50))
SO_subset = subset(x = SOb1, cells = sampled_cells)

# replace astrocyte with "Sender" and microglia with "Receiver"
SO_subset$cell_type %<>% gsub("Oligodendrocyte precursor cells", "Sender",. ) %>% gsub("Dendritic cells", "Receiver",. )

meta = SO_subset@meta.data %>% 
  select(cell_type) %>% 
  tibble::rownames_to_column("cell_ID")
colnames(meta)[2] = "Celltype" 

counts = GetAssayData(SO_subset, "RNA" , "counts")
rownames(counts) = toupper(rownames(counts))
# Filter to decrease size
counts = counts[rowSums(counts != 0) > 10,]
counts = counts[rowSums(counts) > 200,]

write.table(meta,"data/10x_immune_R1/raw_metadata.tsv") # save metadata
write.table(counts,"data/10x_immune_R1/raw_counts.tsv") # save counts

# -------------------------- #
# BD_wholeTranscriptome liver https://cellxgene.cziscience.com/collections/be679cb1-35f0-46c9-9a2d-30691862a54a
# -------------------------- #

'
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("b425976f-9d73-4388-95dd-e7cd0f8caca0.h5ad")
# Subsample to 50% (0.5 fraction)
sc.pp.subsample(adata, fraction=0.5, random_state=42)

adata_P1 = adata[adata.obs["donor_id"] == "P2"]

# remove genes expressed in less than 10 cells
sc.pp.filter_genes(adata_P1, min_cells=10)
# remove genes with less than 200 total counts across all cells
sc.pp.filter_genes(adata_P1, min_counts=200)

# get list of genes
ensemblgenes_to_keep = adata_P1.var_names
genesymbol_to_keep = adata_P1.var["gene_symbols"]
# subset
raw_subset_df_P1 = adata_P1.raw.to_adata()[:, ensemblgenes_to_keep].to_df()
# change genes from ensembl to symbols
raw_subset_df_P1.columns = genesymbol_to_keep
raw_subset_df_P1 = raw_subset_df_P1.transpose()

meta = pd.DataFrame(adata_P1.obs[["cell_type"]])
meta["cell_ID"] = meta.index
meta = meta.rename(columns={"cell_type":"Celltype"})

mapping = {
    "mononuclear phagocyte": "Sender",
    "natural killer cell": "Receiver"
}

# Apply to the column
meta["celltype"] = meta["celltype"].replace(mapping)

# save 
meta.to_csv("/home/vkorob/Documents/git/CCC_benchmark/data/BD_wholeTranscriptome_P1/raw_metadata.tsv", sep="")
raw_subset_df_P1.to_csv("/home/vkorob/Documents/git/CCC_benchmark/data/BD_wholeTranscriptome_P1/raw_counts.tsv", sep="")

adata_P2 = adata[adata.obs["donor_id"] == "P2"]

# remove genes expressed in less than 10 cells
sc.pp.filter_genes(adata_P2, min_cells=10)
# remove genes with less than 200 total counts across all cells
sc.pp.filter_genes(adata_P2, min_counts=200)

# get list of genes
ensemblgenes_to_keep = adata_P2.var_names
genesymbol_to_keep = adata_P2.var["gene_symbols"]
# subset
raw_subset_df_P2 = adata_P2.raw.to_adata()[:, ensemblgenes_to_keep].to_df()
# change genes from ensembl to symbols
raw_subset_df_P2.columns = genesymbol_to_keep
raw_subset_df_P2 = raw_subset_df_P2.transpose()

meta = pd.DataFrame(adata_P2.obs[["cell_type"]])
meta["cell_ID"] = meta.index
meta = meta.rename(columns={"cell_type":"Celltype"})

mapping = {
    "mononuclear phagocyte": "Sender",
    "natural killer cell": "Receiver"
}

# Apply to the column
meta["celltype"] = meta["celltype"].replace(mapping)

# save 
meta.to_csv("/home/vkorob/Documents/git/CCC_benchmark/data/BD_wholeTranscriptome_P2/raw_metadata.tsv", sep="")
raw_subset_df_P2.to_csv("/home/vkorob/Documents/git/CCC_benchmark/data/BD_wholeTranscriptome_P2/raw_counts.tsv", sep="")

## here one has to load in R using counts = read.table("data/BD_wholeTranscriptome_P1/raw_counts.tsv", row.names = 1, header = T) and the save again
# I also transformed cell names so that its not just numbers
'

# -------------------------- #
# VASAseq is in-house dataset
# -------------------------- #
