#!/usr/bin/env Rscript
#
# Main entry point for DESeq2 LRT Step 1 Analysis
#

# Source required function files
source("/usr/local/bin/functions/common/error_handling.R")
source("/usr/local/bin/functions/common/output_utils.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/workflow.R")

# Set up error handling
options(error = function(e) handle_error(e, "DESeq2 LRT Step 1"))

# Run the workflow
tryCatch({
    # Initialize the environment and load required packages
    initialize_environment()
    
    # Execute the main workflow
    run_deseq2_lrt_workflow()
    
    message("DESeq2 LRT Step 1 analysis completed successfully.")
}, error = function(e) {
    handle_error(e, "DESeq2 LRT Step 1")
})
