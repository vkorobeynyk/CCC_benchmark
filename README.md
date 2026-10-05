# CCC_benchmark

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

A Snakemake pipeline for benchmarking **cell-cell communication (CCC) inference methods** on single-cell RNA-seq data with a **semi-simulation** framework.

For each dataset we simulate a set of ligand-receptor (LR) interactions between a *Sender* and a *Receiver* cell type. The pipeline then runs each CCC method on the modified data and scores how well it recovers that ground truth. Two parameters set how strong the simulated signal is:

- **PCE** (*Percentage of Cells Expressing*): the fraction of cells in a cell type that express the inflated ligand/receptor gene. This controls **detection** (sparsity).
- **FC** (*Fold Change*): a multiplier on the gene's estimated negative-binomial mean. This controls **expression magnitude**.

PCE and FC are decoupled. PCE-selected cells get a "+1" detection floor, and FC-scaled signal is added on top of it.

The data folder used in the benchmark, along with the output folder can be found in ................
---

## Contact

Vladyslav Korobeynyk, HIFO / DMLS, University of Zurich (Jessberger lab / Mark D. Robinson lab)
Questions and bug reports: please open a [GitHub issue](https://github.com/vkorobeynyk/CCC_benchmark/issues).

## Generative AI statement
Generative AI was used throughout this entire benchmark to make code nicer to read and more efficient. The entire logic of the benchmark was created by myself and I assume responsability of the content within this repo.

## Table of contents

- [Pipeline overview](#pipeline-overview)
- [Repository structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input data](#input-data)
- [Configuration](#configuration)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Benchmarked methods](#benchmarked-methods)
- [Metrics](#metrics)
- [Adding a new method](#adding-a-new-method)
- [Adding a new dataset](#adding-a-new-dataset)
- [License](#license)

---

## Pipeline overview

```
data/{dataset}/raw_counts.tsv + raw_metadata.tsv
        │
        ▼  1. run_processing                      (edgeR.sif)
data/processed/{dataset}/  counts, metadata, genemetadata.RDS
        │      filtering, edgeR dispersion + per-cell-type means,
        │      selection of LR pairs to simulate
        ▼  2. run_semiSimulation_inflateCounts    (edgeR.sif)   × FC × PCE
output/{dataset}_semiSimulation_NB/  inflated counts, realized PCE / FC
        │
        ▼  3. run_normalization                   (edgeR.sif)   scater::logNormCounts
output/{dataset}_semiSimulation_NB/  inflated normalized counts
        │
        ▼  4. run_methods                         (lianaPlus.sif) × method
output/{dataset}/{method}/significant_interactions_FC_{FC}_PercCellsExpressing_{PCE}.csv
        │
        ▼  5. run_metric_f1score_rankingLRgenes   (edgeR.sif)
output/final_scores.RDS
        │
        ▼  all → scripts/visualization.R          (edgeR.sif)
output/results/
```

Every dataset is run over the full grid of `FC × PCE` values and every method listed in `config.yaml`.

---

## Repository structure

```
CCC_benchmark/
├── Snakefile                  # workflow definition
├── config.yaml                # datasets, methods, FC/PCE grid, LR database
├── sbatch_submit              # example SLURM submission script
├── fix_indentation.sh         # replaces tabs with spaces in Snakefile/config.yaml
├── scripts/
│   ├── processing_dataset.R                    # step 1
│   ├── processing_semisimulation_nbsampling.R  # step 2
│   ├── processing_normalization.R              # step 3
│   ├── method_<name>.py                        # step 4, one script per CCC method
│   ├── metric_f1score_rankingLRgenes.R         # step 5
│   ├── visualization.R                         # final plots
│   └── helper_functions.R                      # semi-simulation + plotting helpers
├── sing_container/
│   ├── README.md              # how to build the containers
│   ├── get_def_files.sh
│   └── defs/
│       ├── edgeR.def          # R / Bioconductor environment
│       └── lianaPlus.def      # Python / LIANA+ environment
├── data/                      # not tracked: see "Input data"
│   ├── LR_database.tsv
│   ├── InfoAbout_LR_database.txt
│   └── download_Datasets.R
└── output/                    # not tracked: created by the pipeline
```

`data/`, `output/`, `log/` and `sing_container/*.sif` are in `.gitignore`.

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [Apptainer](https://apptainer.org/) or [SingularityCE](https://sylabs.io/singularity/). Every step runs inside a container, so you don't need a local R or Python installation.

The final `all` rule calls `singularity exec` directly, so a `singularity` command must be on your `PATH`.

---

## Installation

```bash
git clone git@github.com:vkorobeynyk/CCC_benchmark.git
cd CCC_benchmark
```

Build the two containers (see [`sing_container/README.md`](sing_container/README.md) for details):

---

## Input data

The datasets are not distributed with this repository. Each dataset needs its own folder under `data/`:

```
data/{dataset}/raw_counts.tsv      # genes × cells, raw counts; gene names as row names
data/{dataset}/raw_metadata.tsv    # one row per cell
```

`raw_metadata.tsv` must contain:

- `cell_ID`: matches the column names of `raw_counts.tsv`, in the same order
- `Celltype`: the two cell types between which interactions are simulated must be labelled **`Sender`** and **`Receiver`**

### Datasets used in the benchmark

`data/download_Datasets.R` documents how each dataset was obtained, subsampled and relabelled:

| Dataset | Technology | Source | Sender → Receiver |
|---|---|---|---|
| `10x` | 10x Chromium | [Zenodo 13893852](https://zenodo.org/records/13893852) | Astrocyte → Endothelial |
| `10x_immune_R1` | 10x Chromium | [GEO GSE195848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195848) (Batch1) | OPC → Dendritic cells |
| `SMARTseq2` | Smart-seq2 | [scRNA-seq_analysis_jessberger_lab](https://github.com/vkorobeynyk/scRNA-seq_analysis_jessberger_lab/tree/master/data) | ndNSC → dNSC |
| `BD_wholeTranscriptome_P1` | BD Rhapsody | [CELLxGENE collection](https://cellxgene.cziscience.com/collections/be679cb1-35f0-46c9-9a2d-30691862a54a) | Mononuclear phagocyte → NK cell |
| `VASAseq` | VASA-seq | in-house (not public) | — |

### Ligand-receptor database

`data/LR_database.tsv` is a curated subset of CellPhoneDB in LIANA format. It keeps only interactions also present in `CellChatDB.human`, and removes self-interactions and genes containing hyphens. The rows are shuffled because the simulation picks LR pairs in order. Complexes use the `L_R1_R2` format, and we add signal to every subunit of a complex. See `data/InfoAbout_LR_database.txt` for the full curation code.

---

## Configuration

All run parameters are set in `config.yaml`:

Every combination of `dataset × FC × PCE × method` becomes one job.

---

## Running the pipeline

Do a dry run first:

```bash
snakemake -n
```

Run locally:

```bash
# Snakemake with 4 cores
snakemake --cores 4 --use-singularity
```

### On a SLURM cluster

`sbatch_submit` is an example submission script. Adjust the module names, the bind path and the profile to your cluster:

```bash
mkdir -p log
sbatch sbatch_submit
```

Per-rule memory requirements are set in the `resources:` field of each rule in the `Snakefile`. The simulation step is the most memory-intensive (≈15 GB).

> **Tip:** if Snakemake complains about mixed tabs and spaces, run `./fix_indentation.sh Snakefile` (the SLURM script does this automatically).

---

## Outputs

| Path | Content |
|---|---|
| `data/processed/{dataset}/genemetadata.RDS` | edgeR means/dispersions, offsets, simulated LR pairs (ground truth), dataset summary statistics |
| `output/{dataset}_semiSimulation_NB/` | inflated raw and normalized counts, realized PCE per gene, realized FC after simulation |
| `output/{dataset}/{method}/significant_interactions_*.csv` | per-method results: `ligand_receptor`, `statistics`, `significant` (Sender → Receiver only, ranked from strongest to weakest) |
| `output/final_scores.RDS` | nested list `[[dataset]][[method]][[run]]` of metric tables |
| `output/results/` | diagnostic and benchmark figures |

---

## Benchmarked methods

All methods run through [LIANA+](https://github.com/saezlab/liana-py) using the same LR database, with `expr_prop = 0` and `min_cells = 0` so that methods don't apply their own pre-filtering.

| Config name | Method (LIANA+) | Score used | Significant if |
|---|---|---|---|
| `cellphonedbv5` | CellPhoneDB v5 | `cellphone_pvals` | < 0.05 |
| `cellchat` | CellChat | `cellchat_pvals` | < 0.05 |
| `geometric_mean` | Geometric mean | `gmean_pvals` | < 0.05 |
| `connectome` | Connectome | `scaled_weight` | > 0.5 |
| `natmi_specificity` | NATMI | `spec_weight` | > 0.1 |
| `singlecellsignalR` | SingleCellSignalR | `lrscore` | > 0.5 |
| `scseqcomm` | scSeqComm | `inter_score` | > 0.8 |
| `log2fc` | log2FC | `lr_logfc` | > 1 |

---

## Metrics

`scripts/metric_f1score_rankingLRgenes.R` compares each method's output with the simulated ground truth:

- **TP / FP / FN, precision, recall, F1 score**, computed on interactions flagged as `significant`
- **% significant among simulated**: the share of planted interactions that the method called significant
- **Detection rate**: the fraction of planted interactions that appear anywhere in the method's output
- **Conditional normalized rank**: the mean normalized rank (0 = top, 1 = bottom) of planted interactions, computed **only over those detected**
- **Penalized normalized rank**: the same, but undetected interactions are given the worst rank

Ranking quality is reported as two separate numbers (detection rate and conditional rank) rather than one blended score.

---

## Adding a new method

1. Create `scripts/method_<name>.py`. It reads the Snakemake inputs `inflated_normalized_counts` and `metadata_processed` and the param `LR_database`.
2. Keep only `Sender → Receiver` results and write a tab-separated file to `snakemake.output["significant_interactions"]` with the columns:
   - `ligand_receptor`: `<ligand_complex>_<receptor_complex>`
   - `statistics`: the method's score
   - `significant`: `True` / `False`

   Sort rows from strongest to weakest interaction.
3. Add `<name>` to `methods` in `config.yaml`.

If the method needs packages that aren't in `lianaPlus.sif`, add them to `defs/lianaPlus.def` and rebuild, or give the method its own container in the `Snakefile`.

## Adding a new dataset

1. Create `data/<name>/raw_counts.tsv` and `data/<name>/raw_metadata.tsv` as described in [Input data](#input-data), labelling the two cell types of interest `Sender` and `Receiver`.
2. Add `<name>` to `datasets` in `config.yaml`.

---

## License

Copyright (C) 2026 Vladyslav Korobeynyk

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file for the full text of the GNU General Public License v3.0.

The CCC methods, LR databases and datasets used by the pipeline are distributed under their own licenses; please refer to the respective sources.
