#!/usr/bin/env Rscript
#
# Base CLI argument definitions for ATAC-seq workflows
# This module provides common argument parsing functions to reduce duplication
#

#' Get base ATAC-seq arguments common to all workflows
#'
#' @param parser ArgumentParser object to add arguments to
#' @return Modified parser with base arguments
#' @export
add_base_atac_args <- function(parser) {
  # Analysis parameters
  parser$add_argument(
    "--fdr",
    help = "FDR cutoff for significance filtering. Default: 0.1",
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = "Log2 fold change threshold for significance. Default: 0.59 (1.5 fold)",
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = "Use lfcthreshold as null hypothesis in results function",
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--regulation",
    help = "Direction of differential accessibility: 'both', 'up', or 'down'",
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  
  # Batch correction
  parser$add_argument(
    "--batchcorrection",
    help = "Batch correction method: 'none', 'combatseq', 'limmaremovebatcheffect', 'model'",
    type = "character",
    choices = c("none", "combatseq", "limmaremovebatcheffect", "model"),
    default = "none"
  )
  
  # Clustering options
  parser$add_argument(
    "--cluster",
    help = "Hopach clustering method to be run on normalized accessibility counts",
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = "Type of scaling for accessibility data: 'minmax' or 'zscore'",
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = "Distance metric for HOPACH row clustering",
    type = "character",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    default = "cosangle"
  )
  parser$add_argument(
    "--columndist",
    help = "Distance metric for HOPACH column clustering",
    type = "character",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
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
  
  # Output options
  parser$add_argument(
    "--output",
    help = "Output prefix. Default varies by workflow type",
    type = "character",
    default = NULL
  )
  parser$add_argument(
    "--threads",
    help = "Number of threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run in test mode (first 500 peaks only)",
    action = "store_true",
    default = FALSE
  )
  
  return(parser)
}

#' Get DiffBind-specific arguments
#'
#' @param parser ArgumentParser object to add arguments to
#' @return Modified parser with DiffBind arguments
#' @export
add_diffbind_args <- function(parser) {
  # Input files
  parser$add_argument(
    "--input_files",
    help = "Peak files for all samples (space-separated)",
    nargs = "+",
    required = TRUE,
    type = "character"
  )
  parser$add_argument(
    "--bamfiles",
    help = "BAM files for all samples (space-separated)",
    nargs = "+",
    required = TRUE,
    type = "character"
  )
  parser$add_argument(
    "--name",
    help = "Sample names (space-separated)",
    nargs = "+",
    required = TRUE,
    type = "character"
  )
  parser$add_argument(
    "--meta",
    help = "Metadata file with sample information",
    required = TRUE,
    type = "character"
  )
  
  # DiffBind parameters
  parser$add_argument(
    "--peakformat",
    help = "Peak file format",
    type = "character",
    choices = c("bed", "narrow", "macs"),
    default = "narrow"
  )
  parser$add_argument(
    "--peakcaller",
    help = "Peak caller used",
    type = "character",
    choices = c("macs", "bed", "narrow"),
    default = "macs"
  )
  parser$add_argument(
    "--scorecol",
    help = "Column to use for peak scores",
    type = "integer",
    default = 5
  )
  parser$add_argument(
    "--minoverlap",
    help = "Minimum overlap for consensus peaks",
    type = "integer",
    default = 2
  )
  
  return(parser)
}

#' Parse arguments with fallback handling
#'
#' @param parser Configured ArgumentParser object
#' @param cli_helpers CLI helper environment (optional)
#' @return Parsed arguments list
#' @export
parse_args_with_fallback <- function(parser, cli_helpers = NULL) {
  args <- tryCatch({
    parser$parse_args(commandArgs(trailingOnly = TRUE))
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    all_args <- commandArgs(trailingOnly = TRUE)
    
    if (!is.null(cli_helpers) && exists("cli_helpers") && is.environment(cli_helpers)) {
      message("Using CLI helper functions for manual parsing")
      return(parse_with_helpers(all_args, cli_helpers))
    } else {
      message("CLI helpers not available, using basic manual parsing")
      return(parse_manually(all_args))
    }
  })
  
  return(args)
}

#' Parse arguments using CLI helper functions
#'
#' @param all_args Command line arguments vector
#' @param cli_helpers CLI helper environment
#' @return Parsed arguments list
parse_with_helpers <- function(all_args, cli_helpers) {
  manual_args <- list()
  
  # Required single-value arguments (common pattern)
  required_args <- c("meta")
  for (arg in required_args) {
    manual_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg)
  }
  
  # Multi-value arguments (common pattern)
  array_args <- c("input_files", "bamfiles", "name")
  for (arg in array_args) {
    if (paste0("--", arg) %in% all_args) {
      manual_args[[arg]] <- cli_helpers$parse_multi_value_arg(all_args, arg)
    }
  }
  
  # Optional single-value arguments with defaults
  optional_args <- list(
    peakformat = "narrow",
    peakcaller = "macs",
    regulation = "both",
    batchcorrection = "none",
    cluster = "none",
    scaling_type = "zscore",
    rowdist = "cosangle",
    columndist = "euclid"
  )
  
  for (arg in names(optional_args)) {
    manual_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg, optional_args[[arg]])
  }
  
  # Numeric arguments
  numeric_args <- list(
    scorecol = "integer", minoverlap = "integer", 
    fdr = "double", lfcthreshold = "double", 
    k = "integer", kmax = "integer", threads = "integer"
  )
  numeric_defaults <- list(
    scorecol = 5, minoverlap = 2, fdr = 0.1, lfcthreshold = 0.59, 
    k = 3, kmax = 5, threads = 1
  )
  numeric_values <- cli_helpers$parse_numeric_args(all_args, numeric_args, numeric_defaults)
  manual_args <- c(manual_args, numeric_values)
  
  # Boolean flags
  boolean_flags <- c("use_lfc_thresh", "test_mode")
  boolean_values <- cli_helpers$parse_boolean_flags(all_args, boolean_flags)
  manual_args <- c(manual_args, boolean_values)
  
  return(manual_args)
}

#' Basic manual argument parsing fallback
#'
#' @param all_args Command line arguments vector
#' @return Parsed arguments list
parse_manually <- function(all_args) {
  manual_args <- list()
  
  # Simple flag processing
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
  
  return(manual_args)
}

#' Validate base ATAC arguments
#'
#' @param args Parsed arguments
#' @return Validated args
#' @export
validate_base_atac_args <- function(args) {
  log_message("Validating base ATAC-seq arguments")
  
  # Validate numeric parameters
  if (!is.null(args$fdr) && (args$fdr <= 0 || args$fdr >= 1)) {
    stop("FDR must be between 0 and 1")
  }
  
  if (!is.null(args$lfcthreshold) && args$lfcthreshold < 0) {
    stop("Log fold change threshold must be non-negative")
  }
  
  if (!is.null(args$minoverlap) && args$minoverlap < 1) {
    stop("Minimum overlap must be at least 1")
  }
  
  if (!is.null(args$k) && (args$k < 1 || args$k > 15)) {
    stop("k must be between 1 and 15")
  }
  
  if (!is.null(args$kmax) && (args$kmax < 2 || args$kmax > 9)) {
    stop("kmax must be between 2 and 9")
  }
  
  # Convert boolean string values if needed
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      if (is.character(args[[arg_name]])) {
        args[[arg_name]] <- toupper(args[[arg_name]]) %in% c("TRUE", "T", "YES", "Y", "1")
      }
    }
  }
  
  # Convert numeric values if they're strings
  numeric_args <- c("scorecol", "minoverlap", "fdr", "lfcthreshold", "k", "kmax", "threads")
  for (arg_name in numeric_args) {
    if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
      if (grepl("^[0-9.]+$", args[[arg_name]])) {
        args[[arg_name]] <- as.numeric(args[[arg_name]])
      }
    }
  }
  
  log_message("Base ATAC-seq arguments validated successfully")
  return(args)
}

#' Validate DiffBind arguments
#'
#' @param args Parsed arguments
#' @return Validated args
#' @export
validate_diffbind_args <- function(args) {
  log_message("Validating DiffBind arguments")
  
  # Check required arguments exist
  required_args <- c("input_files", "bamfiles", "name", "meta")
  missing_args <- required_args[!required_args %in% names(args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required DiffBind arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Validate input file counts
  if (length(args$input_files) != length(args$bamfiles)) {
    stop("Number of peak files must match number of BAM files")
  }
  
  if (length(args$input_files) != length(args$name)) {
    stop("Number of input files must match number of sample names")
  }
  
  # Check file existence
  for (file in args$input_files) {
    if (!file.exists(file)) {
      stop(paste("Peak file not found:", file))
    }
  }
  
  for (file in args$bamfiles) {
    if (!file.exists(file)) {
      stop(paste("BAM file not found:", file))
    }
  }
  
  if (!file.exists(args$meta)) {
    stop(paste("Metadata file not found:", args$meta))
  }
  
  log_message("DiffBind arguments validated successfully")
  return(args)
}

#' Print ATAC arguments for debugging
#'
#' @param args Parsed arguments
#' @param workflow_type Type of ATAC workflow
#' @export
print_atac_args <- function(args, workflow_type = "ATAC-seq") {
  log_message(paste(workflow_type, "Analysis Arguments:"))
  
  if (!is.null(args$input_files)) {
    log_message(paste("  Input files:", length(args$input_files), "peak files"))
  }
  if (!is.null(args$bamfiles)) {
    log_message(paste("  BAM files:", length(args$bamfiles), "files"))
  }
  if (!is.null(args$name)) {
    log_message(paste("  Sample names:", paste(args$name, collapse=", ")))
  }
  if (!is.null(args$meta)) {
    log_message(paste("  Metadata file:", args$meta))
  }
  
  log_message(paste("  Peak format:", args$peakformat %||% "default"))
  log_message(paste("  Peak caller:", args$peakcaller %||% "default"))
  log_message(paste("  Score column:", args$scorecol %||% "default"))
  log_message(paste("  Minimum overlap:", args$minoverlap %||% "default"))
  log_message(paste("  FDR threshold:", args$fdr %||% "default"))
  log_message(paste("  LFC threshold:", args$lfcthreshold %||% "default"))
  log_message(paste("  Regulation:", args$regulation %||% "default"))
  log_message(paste("  Batch correction:", args$batchcorrection %||% "default"))
  log_message(paste("  Clustering:", args$cluster %||% "default"))
  log_message(paste("  Output prefix:", args$output %||% "default"))
  log_message(paste("  Test mode:", args$test_mode %||% "default"))
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x