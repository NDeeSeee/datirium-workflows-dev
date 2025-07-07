#!/usr/bin/env Rscript
#
# Command-line argument handling for ATAC-seq Pairwise Analysis
#

# Load required libraries
# library(argparse) - Made optional, using manual parsing fallback

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

#' Define and parse command line arguments for ATAC-seq pairwise analysis
#' Uses common ATAC argument base functions to reduce duplication
#'
#' @return Parsed arguments list
#' @export
get_args <- function() {
  # Use manual parsing as fallback if ArgumentParser is not available
  if (!requireNamespace("argparse", quietly = TRUE)) {
    args <- parse_args_manual_atac()
    # Map to expected format (same as argparse version)
    args$condition1 <- args$treated_name
    args$condition2 <- args$untreated_name
    args$condition_column <- "condition"
    args$input_files <- c(args$untreated_files, args$treated_files)
    args$name <- c(args$untreated_sample_names, args$treated_sample_names)
    return(args)
  }
  
  # Try using argparse first, with fallback to manual parsing
  tryCatch({
    # Simple argparse approach without helper functions
    parser <- argparse::ArgumentParser(description = "ATAC-seq pairwise differential accessibility analysis")
    
    # Basic required arguments for pairwise comparison
    parser$add_argument("-t", "--treated_files", nargs = "+", required = TRUE, help = "Treated/condition1 files")
    parser$add_argument("-u", "--untreated_files", nargs = "+", required = TRUE, help = "Untreated/condition2 files")
    parser$add_argument("-tn", "--treated_name", required = TRUE, help = "Treated group name")
    parser$add_argument("-un", "--untreated_name", required = TRUE, help = "Untreated group name")
    parser$add_argument("-ta", "--treated_sample_names", nargs = "+", required = TRUE, help = "Treated sample names")
    parser$add_argument("-ua", "--untreated_sample_names", nargs = "+", required = TRUE, help = "Untreated sample names")
    parser$add_argument("-o", "--output", default = "atac-pairwise", help = "Output prefix")
    
    # Analysis parameters
    parser$add_argument("--fdr", type = "double", default = 0.05, help = "FDR threshold")
    parser$add_argument("--lfcthreshold", type = "double", default = 0.5, help = "Log fold change threshold")
    parser$add_argument("--use_lfc_thresh", action = "store_true", default = FALSE, help = "Use LFC threshold")
    parser$add_argument("--regulation", default = "both", choices = c("both", "up", "down"), help = "Regulation direction")
    parser$add_argument("--batchcorrection", default = "none", help = "Batch correction method")
    parser$add_argument("--cluster", default = "row", help = "Clustering method")
    parser$add_argument("--scaling_type", default = "zscore", help = "Scaling type")
    parser$add_argument("--k", type = "integer", default = 3, help = "Hopach k parameter")
    parser$add_argument("--kmax", type = "integer", default = 5, help = "Hopach kmax parameter")
    parser$add_argument("-p", "--threads", type = "integer", default = 1, help = "Number of threads")
    parser$add_argument("--test_mode", action = "store_true", default = FALSE, help = "Test mode")
    
    # Parse arguments
    args <- parser$parse_args()
    
    # Map to expected format
    args$condition1 <- args$treated_name
    args$condition2 <- args$untreated_name
    args$condition_column <- "condition"
    args$input_files <- c(args$untreated_files, args$treated_files)
    args$name <- c(args$untreated_sample_names, args$treated_sample_names)
    
  }, error = function(e) {
    message("Warning: argparse failed, using manual argument parsing")
    message(paste("Error:", e$message))
    
    # Fallback to manual parsing
    all_args <- commandArgs(trailingOnly = TRUE)
    args <- parse_manually_pairwise(all_args)
    
    # Map to expected format (same as argparse version)
    args$condition1 <- args$treated_name
    args$condition2 <- args$untreated_name
    args$condition_column <- "condition"
    args$input_files <- c(args$untreated_files, args$treated_files)
    args$name <- c(args$untreated_sample_names, args$treated_sample_names)
  })
  
  # Add pairwise-specific required arguments to validation
  pairwise_required <- c("condition1", "condition2")
  missing_args <- pairwise_required[!pairwise_required %in% names(args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required pairwise arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Set default values for pairwise-specific arguments
  if (is.null(args$condition_column)) {
    args$condition_column <- "condition"
  }
  
  return(args)
}

#' Manual parsing fallback for ATAC pairwise arguments
#' @param all_args Command line arguments vector
#' @return Parsed arguments list
parse_manually_pairwise <- function(all_args) {
  message("Using manual parsing for ATAC pairwise arguments")
  
  # Initialize with defaults
  args <- list(
    condition_column = "condition",
    fdr = 0.05,
    lfcthreshold = 0.5,
    use_lfc_thresh = FALSE,
    regulation = "both",
    batchcorrection = "none",
    cluster = "row",
    scaling_type = "zscore",
    k = 3,
    kmax = 5,
    output = "atac-pairwise",
    threads = 1,
    test_mode = FALSE
  )
  
  # Parse arguments manually
  i <- 1
  while (i <= length(all_args)) {
    arg <- all_args[i]
    
    if (arg == "-t" && i + 1 <= length(all_args)) {
      # Parse treated files (multiple values)
      i <- i + 1
      treated_files <- c()
      while (i <= length(all_args) && !startsWith(all_args[i], "-")) {
        treated_files <- c(treated_files, all_args[i])
        i <- i + 1
      }
      args$treated_files <- treated_files
      next
    } else if (arg == "-u" && i + 1 <= length(all_args)) {
      # Parse untreated files (multiple values)
      i <- i + 1
      untreated_files <- c()
      while (i <= length(all_args) && !startsWith(all_args[i], "-")) {
        untreated_files <- c(untreated_files, all_args[i])
        i <- i + 1
      }
      args$untreated_files <- untreated_files
      next
    } else if (arg == "-tn" && i + 1 <= length(all_args)) {
      args$treated_name <- all_args[i + 1]
      # Map to condition1 (treatment)
      args$condition1 <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "-un" && i + 1 <= length(all_args)) {
      args$untreated_name <- all_args[i + 1]
      # Map to condition2 (reference)
      args$condition2 <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "-ta" && i + 1 <= length(all_args)) {
      # Parse treated sample names (multiple values)
      i <- i + 1
      treated_sample_names <- c()
      while (i <= length(all_args) && !startsWith(all_args[i], "-")) {
        treated_sample_names <- c(treated_sample_names, all_args[i])
        i <- i + 1
      }
      args$treated_sample_names <- treated_sample_names
      next
    } else if (arg == "-ua" && i + 1 <= length(all_args)) {
      # Parse untreated sample names (multiple values)
      i <- i + 1
      untreated_sample_names <- c()
      while (i <= length(all_args) && !startsWith(all_args[i], "-")) {
        untreated_sample_names <- c(untreated_sample_names, all_args[i])
        i <- i + 1
      }
      args$untreated_sample_names <- untreated_sample_names
      next
    } else if (arg == "-o" && i + 1 <= length(all_args)) {
      args$output <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "--fdr" && i + 1 <= length(all_args)) {
      args$fdr <- as.numeric(all_args[i + 1])
      i <- i + 2
    } else if (arg == "--lfcthreshold" && i + 1 <= length(all_args)) {
      args$lfcthreshold <- as.numeric(all_args[i + 1])
      i <- i + 2
    } else if (arg == "--use_lfc_thresh") {
      args$use_lfc_thresh <- TRUE
      i <- i + 1
    } else if (arg == "--test_mode") {
      args$test_mode <- TRUE
      i <- i + 1
    } else if (arg == "--regulation" && i + 1 <= length(all_args)) {
      args$regulation <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "--scaling_type" && i + 1 <= length(all_args)) {
      args$scaling_type <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "--cluster" && i + 1 <= length(all_args)) {
      args$cluster <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "--batchcorrection" && i + 1 <= length(all_args)) {
      args$batchcorrection <- all_args[i + 1]
      i <- i + 2
    } else if (arg == "--k" && i + 1 <= length(all_args)) {
      args$k <- as.integer(all_args[i + 1])
      i <- i + 2
    } else if (arg == "--kmax" && i + 1 <= length(all_args)) {
      args$kmax <- as.integer(all_args[i + 1])
      i <- i + 2
    } else if (arg == "-p" && i + 1 <= length(all_args)) {
      args$threads <- as.integer(all_args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  # Create metadata information based on the parsed arguments
  # This mimics what the original pairwise workflow expects
  args$meta <- NULL  # Will be created from sample information
  args$input_files <- c(args$untreated_files, args$treated_files)
  args$name <- c(args$untreated_sample_names, args$treated_sample_names)
  
  message(paste("Parsed condition1 (treatment):", args$condition1))
  message(paste("Parsed condition2 (reference):", args$condition2))
  message(paste("Parsed", length(args$input_files), "input files"))
  message(paste("Parsed", length(args$name), "sample names"))
  
  return(args)
}

#' Validate command line arguments for ATAC-seq pairwise analysis
#' Uses common validation functions to reduce duplication
#'
#' @param args Parsed arguments
#' @return Validated args
#' @export
validate_args <- function(args) {
  log_message("Validating ATAC-seq pairwise arguments")
  
  # Basic validation without helper functions
  errors <- c()
  
  # Check required arguments
  if (is.null(args$condition1) || args$condition1 == "") {
    errors <- c(errors, "condition1 (treated_name) is required")
  }
  
  if (is.null(args$condition2) || args$condition2 == "") {
    errors <- c(errors, "condition2 (untreated_name) is required")
  }
  
  # Pairwise-specific validation
  if (!is.null(args$condition1) && !is.null(args$condition2) && args$condition1 == args$condition2) {
    errors <- c(errors, "Condition1 and condition2 cannot be the same")
  }
  
  # Validate numeric parameters
  if (!is.null(args$fdr) && (args$fdr <= 0 || args$fdr >= 1)) {
    errors <- c(errors, "FDR must be between 0 and 1")
  }
  
  if (!is.null(args$lfcthreshold) && args$lfcthreshold < 0) {
    errors <- c(errors, "Log fold change threshold must be non-negative")
  }
  
  if (!is.null(args$k) && (args$k < 1 || args$k > 15)) {
    errors <- c(errors, "k must be between 1 and 15")
  }
  
  if (!is.null(args$kmax) && (args$kmax < 2 || args$kmax > 9)) {
    errors <- c(errors, "kmax must be between 2 and 9")
  }
  
  # Check file arguments
  if (is.null(args$input_files) || length(args$input_files) == 0) {
    errors <- c(errors, "No input files provided")
  }
  
  if (is.null(args$name) || length(args$name) == 0) {
    errors <- c(errors, "No sample names provided")
  }
  
  if (!is.null(args$input_files) && !is.null(args$name) && 
      length(args$input_files) != length(args$name)) {
    errors <- c(errors, "Number of input files must match number of sample names")
  }
  
  # Convert boolean values if they're strings
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      if (is.character(args[[arg_name]])) {
        args[[arg_name]] <- toupper(args[[arg_name]]) %in% c("TRUE", "T", "YES", "Y", "1")
      }
    }
  }
  
  if (length(errors) > 0) {
    stop("Argument validation failed:\n", paste(errors, collapse = "\n"))
  }
  
  log_message("ATAC-seq pairwise arguments validated successfully")
  return(args)
}

#' Print command line arguments for debugging
#' Uses common print function with pairwise-specific additions
#'
#' @param args Parsed arguments
#' @export
print_args <- function(args) {
  # Use common print function
  print_atac_args(args, "ATAC-seq Pairwise")
  
  # Add pairwise-specific information
  log_message("  Pairwise-specific parameters:")
  log_message(paste("    Condition column:", args$condition_column))
  log_message(paste("    Condition 1 (treatment):", args$condition1))
  log_message(paste("    Condition 2 (reference):", args$condition2))
}

#' Manual argument parsing function as fallback when argparse is not available
#'
#' @return Parsed argument list
parse_args_manual_atac <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize result with defaults
  result <- list(
    treated_files = NULL,
    untreated_files = NULL,
    treated_name = "treated",
    untreated_name = "untreated",
    treated_sample_names = NULL,
    untreated_sample_names = NULL,
    batchcorrection = "none",
    fdr = 0.05,
    lfcthreshold = 0.5,
    use_lfc_thresh = FALSE,
    regulation = "both",
    scaling_type = "zscore",
    cluster = "row",
    rowdist = "cosangle",
    columndist = "euclid",
    k = 3,
    kmax = 5,
    output = "atac_pairwise",
    threads = 1,
    test_mode = FALSE
  )
  
  # Parse arguments manually
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    
    if (arg %in% c("-t", "--treated_files")) {
      result$treated_files <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$treated_files <- c(result$treated_files, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-u", "--untreated_files")) {
      result$untreated_files <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$untreated_files <- c(result$untreated_files, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-ta", "--treated_sample_names")) {
      result$treated_sample_names <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$treated_sample_names <- c(result$treated_sample_names, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-ua", "--untreated_sample_names")) {
      result$untreated_sample_names <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$untreated_sample_names <- c(result$untreated_sample_names, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-tn", "--treated_name")) {
      result$treated_name <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-un", "--untreated_name")) {
      result$untreated_name <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-o", "--output")) {
      result$output <- args[i + 1]
      i <- i + 2
    } else if (arg == "--fdr") {
      result$fdr <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg == "--lfcthreshold") {
      result$lfcthreshold <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (arg %in% c("-p", "--threads")) {
      result$threads <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--test_mode") {
      result$test_mode <- TRUE
      i <- i + 1
    } else if (arg == "--use_lfc_thresh") {
      result$use_lfc_thresh <- TRUE
      i <- i + 1
    } else if (arg == "--batchcorrection") {
      result$batchcorrection <- args[i + 1]
      i <- i + 2
    } else if (arg == "--regulation") {
      result$regulation <- args[i + 1]
      i <- i + 2
    } else if (arg == "--scaling_type") {
      result$scaling_type <- args[i + 1]
      i <- i + 2
    } else if (arg == "--cluster") {
      result$cluster <- args[i + 1]
      i <- i + 2
    } else {
      # Skip unknown arguments
      i <- i + 1
    }
  }
  
  return(result)
}