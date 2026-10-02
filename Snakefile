import os

# Load configuration parameters
configfile: "config.yaml"

# =============================================================================
# TARGET RULE
# =============================================================================
# This rule defines the final desired outputs of the workflow.
# It also executes the final visualization script.
rule all:
    input:
        final_scores = "output/final_scores.RDS"
    output:
        visualization = directory("output/results")
    shell:
        """
        singularity exec --no-home sing_container/edgeR.sif \
        Rscript scripts/visualization.R \
            --config.yaml config.yaml \
            --path_output_dir output \
            --path_results_dir "output/results"
        """

# =============================================================================
# STEP 1: DATA PROCESSING
# =============================================================================
# Cleans raw counts and metadata, and estimates dispersion using edgeR.
rule run_processing:
    threads: 1
    resources:
        mem_mb = 10000
    input:
        raw_counts = "data/{dataset}/raw_counts.tsv",
        raw_metadata = "data/{dataset}/raw_metadata.tsv"
    output:
        processed_counts = "data/processed/{dataset}/counts_{dataset}_processed.tsv",
        metadata_processed = "data/processed/{dataset}/metadata_{dataset}_processed.tsv",
        genemetadata = "data/processed/{dataset}/genemetadata.RDS"
    params:
        nLR_per_CTCTcomb = config["nLR_per_CTCTcomb"],
        LR_database = config["LR_database"]
    container:
        "sing_container/edgeR.sif"
    script:
        "scripts/processing_dataset.R"

# =============================================================================
# STEP 2: SEMI-SIMULATION (NB SAMPLING)
# =============================================================================
# Inflates gene expression counts based on specific Fold Change (FC) 
# and Percentage of Cells Expressing (PCE) parameters.
rule run_semiSimulation_inflateCounts:
    threads: 1
    resources:
        mem_mb = 15000
    input:
        processed_counts = "data/processed/{dataset}/counts_{dataset}_processed.tsv",
        metadata_processed = "data/processed/{dataset}/metadata_{dataset}_processed.tsv",
        genemetadata = "data/processed/{dataset}/genemetadata.RDS"
    output:
        inflated_counts = "output/{dataset}_semiSimulation_NB/sc_inflated_counts_FC_{FC}_PercCellsExpressing_{PCE}.tsv",
        metadata_processed = "output/{dataset}_semiSimulation_NB/sc_metadata_FC_{FC}_PercCellsExpressing_{PCE}.tsv",
        PCE_perGene = "output/{dataset}_semiSimulation_NB/PCE_perGene_FC_{FC}_PercCellsExpressing_{PCE}.RDS",
        FC_after_simulation = "output/{dataset}_semiSimulation_NB/realFC_aftersimulation_FC_{FC}_PercCellsExpressing_{PCE}.RDS"
    params:
        nLR_per_CTCTcomb = config["nLR_per_CTCTcomb"],
        LR_database = config["LR_database"]
    container:
        "sing_container/edgeR.sif"
    script:
        "scripts/processing_semisimulation_nbsampling.R"

# =============================================================================
# STEP 3: NORMALIZATION
# =============================================================================
# Normalizes the inflated counts for downstream CCC method compatibility.
rule run_normalization:
    threads: 1
    resources:
        mem_mb = 10000
    input:
        inflated_counts = "output/{dataset}_semiSimulation_NB/sc_inflated_counts_FC_{FC}_PercCellsExpressing_{PCE}.tsv"
    output:
        inflated_normalized_counts = "output/{dataset}_semiSimulation_NB/sc_inflated_normalized_counts_FC_{FC}_PercCellsExpressing_{PCE}.tsv"
    container:
        "sing_container/edgeR.sif"
    script:
        "scripts/processing_normalization.R"

# =============================================================================
# STEP 4: CCC METHOD EXECUTION
# =============================================================================
# Iterates through every method defined in the config to find significant interactions.
rule run_methods:
    threads: 1
    resources:
        mem_mb = 10000
    input:
        inflated_normalized_counts = "output/{dataset}_semiSimulation_NB/sc_inflated_normalized_counts_FC_{FC}_PercCellsExpressing_{PCE}.tsv",
        metadata_processed = "output/{dataset}_semiSimulation_NB/sc_metadata_FC_{FC}_PercCellsExpressing_{PCE}.tsv"
    output:
        significant_interactions = "output/{dataset}/{method}/significant_interactions_FC_{FC}_PercCellsExpressing_{PCE}.csv"
    params:
        LR_database = config["LR_database"]
    container:
        "sing_container/lianaPlus.sif"
    script:
        "scripts/method_{wildcards.method}.py"

# =============================================================================
# STEP 5: METRIC CALCULATION
# =============================================================================
# Calculates F1 scores and ranking performance across all datasets and methods.
rule run_metric_f1score_rankingLRgenes:
    threads: 1
    resources:
        mem_mb = 1000
    input:
        significant_interactions = expand(
            "output/{dataset}/{method}/significant_interactions_FC_{FC}_PercCellsExpressing_{PCE}.csv",
            method = config["methods"], 
            FC = config["semiSimulation"]["FC"], 
            dataset = config["datasets"], 
            PCE = config["semiSimulation"]["PCE"]
        ),
        genemetadata = expand(
            "data/processed/{dataset}/genemetadata.RDS",
            dataset = config["datasets"]
        )
    output:
        final_scores = "output/final_scores.RDS"
    container:
        "sing_container/edgeR.sif"
    script:
        "scripts/metric_f1score_rankingLRgenes.R"
