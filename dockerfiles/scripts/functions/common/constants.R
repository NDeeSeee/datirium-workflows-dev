#!/usr/bin/env Rscript

# --- Common Constants ---
# This file contains constants used across multiple DESeq analysis scripts.

# Column names and identifiers
READ_COL <- "reads"                  # Pattern to identify read count columns
RPKM_COL <- "RPKM"                   # Pattern to identify RPKM columns
RPKM_TREATED_ALIAS <- "TREATED_RPKM"  # Alias for treated RPKM values
RPKM_UNTREATED_ALIAS <- "UNTREATED_RPKM" # Alias for untreated RPKM values
INTERSECT_BY <- "GeneId"             # Column to use for matching between datasets

# Default thresholds
DEFAULT_FDR_CUTOFF <- 0.1            # Default FDR threshold
DEFAULT_LFC_THRESHOLD <- 0.585       # Default log2 fold change threshold (~1.5 fold)
DEFAULT_RPKM_CUTOFF <- 1.0           # Default RPKM expression cutoff

# Scaling type options for clustering
SCALING_TYPES <- c("none", "row", "column", "both")
DEFAULT_SCALING <- "row"             # Default scaling type

# Distance metrics for clustering
DISTANCE_METRICS <- c("euclidean", "maximum", "manhattan", "canberra", 
                     "binary", "pearson", "abspearson", "correlation", 
                     "abscorrelation", "spearman", "kendall")
DEFAULT_DISTANCE <- "euclidean"      # Default distance metric

# Clustering methods
CLUSTERING_METHODS <- c("none", "row", "column", "both")
DEFAULT_CLUSTERING <- "none"         # Default clustering method

# Batch correction methods
BATCH_CORRECTION_METHODS <- c("none", "model", "combat")
DEFAULT_BATCH_CORRECTION <- "none"   # Default batch correction method

# File paths and prefixes
DEFAULT_OUTPUT_PREFIX <- "deseq2_results" # Default prefix for output files

# Analysis types
ANALYSIS_TYPES <- c("standard", "LRT")
DEFAULT_ANALYSIS_TYPE <- "standard"  # Default analysis type

# Export options
DEFAULT_GCT_EXPORT <- TRUE           # Whether to export GCT files by default
DEFAULT_CLS_EXPORT <- TRUE           # Whether to export CLS files by default
DEFAULT_CSV_EXPORT <- TRUE           # Whether to export CSV files by default

# Visualization options
DEFAULT_WIDTH <- 8                   # Default plot width in inches
DEFAULT_HEIGHT <- 6                  # Default plot height in inches
DEFAULT_DPI <- 300                   # Default plot resolution

# Other parameters
DEFAULT_DIGITS <- 4                  # Default number of digits for rounding
DEFAULT_THREADS <- 1                 # Default number of threads to use
DEFAULT_SEED <- 123                  # Default random seed for reproducibility

# DiffBind constants (required for ATAC-seq analysis)
# These must match the values defined in DiffBind package
DBA_CONDITION <- 4         # DiffBind condition constant
DBA_DESEQ2 <- "DESeq2"     # DiffBind DESeq2 method constant (string, not numeric!)
DBA_SCORE_READS <- 3       # DiffBind score type: raw reads (CORRECTED)
DBA_SCORE_RPKM <- 1        # DiffBind score type: RPKM (value is 1, not 2)
DBA_SCORE_TMM_MINUS_FULL <- 6  # DiffBind score type: TMM normalized

# Export this to environment when sourced
.common_constants_loaded <- TRUE 