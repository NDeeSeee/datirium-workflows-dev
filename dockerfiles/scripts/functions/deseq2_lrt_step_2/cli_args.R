#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq2 LRT Step 2
#

# Try to source helper functions with robust path resolution
if (exists("source_cli_helpers")) {
  source_cli_helpers()
} else {
  # Fallback to basic sourcing
  tryCatch({
    possible_paths <- c(
      file.path(dirname(getwd()), "common", "cli_helpers.R"),
      "../common/cli_helpers.R",
      "/usr/local/bin/functions/common/cli_helpers.R"
    )
    sourced <- FALSE
    for (path in possible_paths) {
      if (file.exists(path)) {
        source(path)
        sourced <- TRUE
        break
      }
    }
    if (!sourced) {
      message("CLI helpers not available, using manual parsing")
    }
  }, error = function(e) {
    message("CLI helpers not available, using manual parsing")
  })
}

#' Define and parse command line arguments
#'
#' @return Parsed arguments list
#' @export
get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 analysis using contrasts from previous LRT step")
  parser$add_argument(
    "--dsq_obj_data",
    help = "RDS file containing contrasts and expression data from step 1",
    required = TRUE,
    type = "character"
  )
  parser$add_argument(
    "--contrast_df",
    required = TRUE,
    help = "TSV file containing contrasts data",
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = "Batch correction method to use: 'none', 'combatseq', or 'model'",
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  parser$add_argument(
    "--contrast_indices",
    help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "FDR cutoff for significance filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Use lfcthreshold as the null hypothesis value in the results function call.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison: 'both', 'up', or 'down'.",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts.",
      "Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Type of scaling for expression data: 'minmax' or 'zscore'.",
      "Default: zscore"
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = paste(
      "Distance metric for HOPACH row clustering.",
      "Default: cosangle"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "cosangle"
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering.",
      "Default: euclid"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "euclid"
  )
  parser$add_argument(
    "--k",
    help = "Number of levels for Hopach clustering (1-15). Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help = "Maximum number of clusters at each level (2-9). Default: 5.",
    type = "integer",
    default = 5
  )
  parser$add_argument(
    "--output",
    help = "Output prefix. Default: deseq-lrt-step-2",
    type = "character",
    default = "deseq-lrt-step-2"
  )
  parser$add_argument(
    "--threads",
    help = "Number of threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run in test mode (first 500 rows only)",
    action = "store_true",
    default = FALSE
  )
  
  # Parse arguments with better error handling
  args <- tryCatch({
    parser$parse_args(commandArgs(trailingOnly = TRUE))
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Use helper functions if available, otherwise fallback to manual parsing
    if (exists("cli_helpers") && is.environment(cli_helpers)) {
      message("Using CLI helper functions for manual parsing")
      
      # Parse using helpers
      manual_args <- list()
      
      # Required single-value arguments
      required_args <- c("dsq_obj_data", "contrast_df", "contrast_indices")
      for (arg in required_args) {
        manual_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg)
      }
      
      # Optional single-value arguments with defaults
      optional_args <- list(
        batchcorrection = "none",
        regulation = "both",
        cluster = "none",
        scaling_type = "zscore",
        rowdist = "cosangle",
        columndist = "euclid",
        output = "deseq-lrt-step-2"
      )
      
      for (arg in names(optional_args)) {
        manual_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg, optional_args[[arg]])
      }
      
      # Numeric arguments
      numeric_args <- list(fdr = "double", lfcthreshold = "double", k = "integer", kmax = "integer", threads = "integer")
      numeric_defaults <- list(fdr = 0.1, lfcthreshold = 0.59, k = 3, kmax = 5, threads = 1)
      numeric_values <- cli_helpers$parse_numeric_args(all_args, numeric_args, numeric_defaults)
      manual_args <- c(manual_args, numeric_values)
      
      # Boolean flags
      boolean_flags <- c("use_lfc_thresh", "test_mode")
      boolean_values <- cli_helpers$parse_boolean_flags(all_args, boolean_flags)
      manual_args <- c(manual_args, boolean_values)
      
    } else {
      # Fallback to original manual parsing logic
      message("CLI helpers not available, using original manual parsing")
      manual_args <- list()
      
      # Original manual parsing code as fallback...
      required_args <- c("dsq_obj_data", "contrast_df", "contrast_indices")
      for (req_arg in required_args) {
        req_flag <- paste0("--", req_arg)
        arg_idx <- which(all_args == req_flag)
        if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
          manual_args[[req_arg]] <- all_args[arg_idx[1] + 1]
        }
      }
      
      # Simple flag processing for other arguments
      i <- 1
      while (i <= length(all_args)) {
        current_arg <- all_args[i]
        if (grepl("^--", current_arg)) {
          arg_name <- sub("^--", "", current_arg)
          if (i < length(all_args) && !grepl("^--", all_args[i + 1])) {
            manual_args[[arg_name]] <- all_args[i + 1]
            i <- i + 2
          } else {
            manual_args[[arg_name]] <- TRUE
            i <- i + 1
          }
        } else {
          i <- i + 1
        }
      }
    }
    
    message("Manually parsed arguments using helpers")
    return(manual_args)
  })
  
  # Validate required arguments
  required_args <- c("dsq_obj_data", "contrast_df", "contrast_indices")
  missing_args <- required_args[!required_args %in% names(args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Process contrast indices
  if (!is.null(args$contrast_indices)) {
    args$contrast_indices <- as.numeric(unlist(strsplit(args$contrast_indices, ",")))
  }
  
  # Convert boolean string values if needed
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      if (is.character(args[[arg_name]])) {
        args[[arg_name]] <- toupper(args[[arg_name]]) %in% c("TRUE", "T", "YES", "Y", "1")
      }
    }
  }
  
  # Convert numeric values
  for (arg_name in c("fdr", "lfcthreshold", "k", "kmax", "threads")) {
    if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
      if (grepl("^[0-9.]+$", args[[arg_name]])) {
        args[[arg_name]] <- as.numeric(args[[arg_name]])
      }
    }
  }
  
  return(args)
} 