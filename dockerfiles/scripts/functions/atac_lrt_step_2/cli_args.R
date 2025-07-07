#!/usr/bin/env Rscript
#
# Command-line argument handling for ATAC-seq LRT Step 2
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
  # Try to load argparse with error handling
  tryCatch({
    if (!require("argparse", quietly = TRUE)) {
      stop("argparse package not available")
    }
    
    parser <- argparse::ArgumentParser(description = "Run ATAC-seq analysis using contrasts from previous LRT step")
    parser$add_argument(
      "--atac_obj_data",
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
      help = "Batch correction method to use: 'none', 'combatseq', 'limmaremovebatcheffect', or 'model'",
      type = "character",
      choices = c("none", "combatseq", "limmaremovebatcheffect", "model"),
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
        "Log2 fold change threshold for determining significant differential accessibility.",
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
        "Direction of differential accessibility comparison: 'both', 'up', or 'down'.",
        "Default: both"
      ),
      type = "character",
      choices = c("both", "up", "down"),
      default = "both"
    )
    parser$add_argument(
      "--cluster",
      help = paste(
        "Hopach clustering method to be run on normalized accessibility counts.",
        "Default: none"
      ),
      type = "character",
      choices = c("row", "column", "both", "none"),
      default = "none"
    )
    parser$add_argument(
      "--scaling_type",
      help = paste(
        "Type of scaling for accessibility data: 'minmax' or 'zscore'.",
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
      help = "Output prefix. Default: atac-lrt-step-2",
      type = "character",
      default = "atac-lrt-step-2"
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
    
    # Parse arguments
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return(args)
    
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    return(parse_manually_atac_step2())
  })
}

#' Manual argument parsing fallback for ATAC LRT Step 2
#'
#' @return Parsed arguments list
parse_manually_atac_step2 <- function() {
  all_args <- commandArgs(trailingOnly = TRUE)
  args <- list()
  
  # Set defaults
  args$batchcorrection <- "none"
  args$fdr <- 0.1
  args$lfcthreshold <- 0.59
  args$use_lfc_thresh <- FALSE
  args$regulation <- "both"
  args$cluster <- "none"
  args$scaling_type <- "zscore"
  args$rowdist <- "cosangle"
  args$columndist <- "euclid"
  args$k <- 3
  args$kmax <- 5
  args$output <- "atac-lrt-step-2"
  args$threads <- 1
  args$test_mode <- FALSE
  
  # Parse arguments
  i <- 1
  while (i <= length(all_args)) {
    current_arg <- all_args[i]
    
    if (current_arg == "--atac_obj_data" && i < length(all_args)) {
      args$atac_obj_data <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--contrast_df" && i < length(all_args)) {
      args$contrast_df <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--contrast_indices" && i < length(all_args)) {
      args$contrast_indices <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--batchcorrection" && i < length(all_args)) {
      args$batchcorrection <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--fdr" && i < length(all_args)) {
      args$fdr <- as.numeric(all_args[i + 1])
      i <- i + 2
    } else if (current_arg == "--lfcthreshold" && i < length(all_args)) {
      args$lfcthreshold <- as.numeric(all_args[i + 1])
      i <- i + 2
    } else if (current_arg == "--regulation" && i < length(all_args)) {
      args$regulation <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--cluster" && i < length(all_args)) {
      args$cluster <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--scaling_type" && i < length(all_args)) {
      args$scaling_type <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--rowdist" && i < length(all_args)) {
      args$rowdist <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--columndist" && i < length(all_args)) {
      args$columndist <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--k" && i < length(all_args)) {
      args$k <- as.integer(all_args[i + 1])
      i <- i + 2
    } else if (current_arg == "--kmax" && i < length(all_args)) {
      args$kmax <- as.integer(all_args[i + 1])
      i <- i + 2
    } else if (current_arg == "--output" && i < length(all_args)) {
      args$output <- all_args[i + 1]
      i <- i + 2
    } else if (current_arg == "--threads" && i < length(all_args)) {
      args$threads <- as.integer(all_args[i + 1])
      i <- i + 2
    } else if (current_arg == "--use_lfc_thresh") {
      if (i < length(all_args) && !startsWith(all_args[i + 1], "--")) {
        args$use_lfc_thresh <- toupper(all_args[i + 1]) %in% c("TRUE", "T", "YES", "Y", "1")
        i <- i + 2
      } else {
        args$use_lfc_thresh <- TRUE
        i <- i + 1
      }
    } else if (current_arg == "--test_mode") {
      if (i < length(all_args) && !startsWith(all_args[i + 1], "--")) {
        args$test_mode <- toupper(all_args[i + 1]) %in% c("TRUE", "T", "YES", "Y", "1")
        i <- i + 2
      } else {
        args$test_mode <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  # Validate required arguments
  required_args <- c("atac_obj_data", "contrast_df", "contrast_indices")
  missing_args <- required_args[!required_args %in% names(args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Process contrast indices
  if (!is.null(args$contrast_indices)) {
    args$contrast_indices <- as.numeric(unlist(strsplit(args$contrast_indices, ",")))
  }
  
  message("Successfully parsed arguments manually for ATAC LRT Step 2")
  return(args)
}