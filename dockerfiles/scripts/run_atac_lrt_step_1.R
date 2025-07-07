#!/usr/bin/env Rscript

# Main runner script for ATAC-seq LRT Step 1 Analysis

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
})

# Source the workflow functions
source("/usr/local/bin/functions/atac_lrt_step_1/workflow.R")

# Run the main workflow - CRITICAL: No parameters since function calls get_args() internally
tryCatch({
  initialize_environment()
  main_workflow()
}, error = function(e) {
  cat("ERROR: ATAC-seq LRT Step 1 analysis failed\n")
  cat("Error message:", conditionMessage(e), "\n")
  traceback()
  quit(status = 1)
})

message("ATAC-seq LRT Step 1 analysis completed successfully!")
