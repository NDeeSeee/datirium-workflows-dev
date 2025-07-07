#!/usr/bin/env Rscript
#
# Main entry point for ATAC-seq LRT Step 2 Analysis
#

# Source required function files
source("/usr/local/bin/functions/common/error_handling.R")
source("/usr/local/bin/functions/common/output_utils.R")
source("/usr/local/bin/functions/atac_lrt_step_2/workflow.R")

# Set up error handling
options(error = function(e) handle_error(e, "ATAC-seq LRT Step 2"))

# Run the workflow
tryCatch({
    # Initialize the environment and load required packages
    initialize_environment()
    
    # Execute the main workflow with memory management
    results <- main_with_memory_management()
    
    message("ATAC-seq LRT Step 2 analysis completed successfully.")
}, error = function(e) {
    handle_error(e, "ATAC-seq LRT Step 2")
}) 