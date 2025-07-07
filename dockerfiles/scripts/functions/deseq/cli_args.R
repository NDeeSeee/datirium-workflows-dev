#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for parsing and validating command-line arguments
# for the main DESeq analysis workflow.
#
# Version: 0.1.0

# Try to source helper functions (optional, with fallback)
tryCatch({
  source_path <- file.path(dirname(getwd()), "common", "cli_helpers.R")
  if (file.exists(source_path)) {
    source(source_path)
  } else {
    # Try Docker path
    docker_path <- "/usr/local/bin/functions/common/cli_helpers.R"
    if (file.exists(docker_path)) {
      source(docker_path)
    }
  }
}, error = function(e) {
  # Helpers not available, continue with manual parsing
  message("CLI helpers not available, using manual parsing")
})

#' Assert and validate command line arguments
#'
#' @param args The parsed arguments from ArgumentParser
#' @return Modified args with validated and processed values
assert_args <- function(parsed_arguments) {
  log_message("Checking input parameters")
  
  # Process aliases if not provided
  if (is.null(parsed_arguments$untreated_sample_names) | is.null(parsed_arguments$treated_sample_names)) {
    log_message("--untreated_sample_names or --treated_sample_names were not set, using default values based on expression file names")
    
    parsed_arguments$untreated_sample_names <- character(0)
    for (i in 1:length(parsed_arguments$untreated_files)) {
      parsed_arguments$untreated_sample_names <- append(parsed_arguments$untreated_sample_names, head(unlist(
        strsplit(basename(parsed_arguments$untreated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
    
    parsed_arguments$treated_sample_names <- character(0)
    for (i in 1:length(parsed_arguments$treated_files)) {
      parsed_arguments$treated_sample_names <- append(parsed_arguments$treated_sample_names, head(unlist(
        strsplit(basename(parsed_arguments$treated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
  } else {
    # Verify correct number of aliases; if mismatched, regenerate from filenames
    if ((length(parsed_arguments$untreated_sample_names) != length(parsed_arguments$untreated_files)) |
        (length(parsed_arguments$treated_sample_names) != length(parsed_arguments$treated_files))) {
      log_warning("Length mismatch between sample alias vectors and file lists; regenerating aliases from file basenames")
      parsed_arguments$untreated_sample_names <- sapply(parsed_arguments$untreated_files, function(x) sub("\\..*$", "", basename(x)))
      parsed_arguments$treated_sample_names   <- sapply(parsed_arguments$treated_files,   function(x) sub("\\..*$", "", basename(x)))
    }
  }

  # Check for minimum file requirements
  if (length(parsed_arguments$treated_files) == 1 || length(parsed_arguments$untreated_files) == 1) {
    log_warning("Only one file in a group. DESeq2 requires at least two replicates for accurate analysis.")
    parsed_arguments$batch_file <- NULL # reset batch_file to NULL. We don't need it for DESeq even if it was provided
  }

  # Process batch file if provided
  if (!is.null(parsed_arguments$batch_file)) {
    batch_metadata <- with_error_handling({
      read.table(
        parsed_arguments$batch_file,
        sep = get_file_type(parsed_arguments$batch_file),
        row.names = 1,
        col.names = c("name", "batch"),
        header = FALSE,
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(batch_metadata)) {
      log_error("Failed to read batch metadata file")
      parsed_arguments$batch_file <- NULL
      return(parsed_arguments)
    }
    
    log_message("Loaded batch metadata")
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    
    if (all(is.element(c(parsed_arguments$untreated_sample_names, parsed_arguments$treated_sample_names), rownames(batch_metadata)))) {
      parsed_arguments$batch_file <- batch_metadata # dataframe
    } else {
      log_warning("Missing values in batch metadata file. Skipping multi-factor analysis")
      log_debug(paste("Expected:", paste(c(parsed_arguments$untreated_sample_names, parsed_arguments$treated_sample_names), collapse=", ")))
      log_debug(paste("Found:", paste(rownames(batch_metadata), collapse=", ")))
      parsed_arguments$batch_file <- NULL
    }
  }

  # Convert boolean string values if they came as strings
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(parsed_arguments[[arg_name]])) {
      parsed_arguments[[arg_name]] <- convert_to_boolean(parsed_arguments[[arg_name]], FALSE)
    }
  }

  # Map argument names for compatibility with workflow
  parsed_arguments$treated <- parsed_arguments$treated_files
  parsed_arguments$untreated <- parsed_arguments$untreated_files
  parsed_arguments$talias <- parsed_arguments$treated_sample_names  
  parsed_arguments$ualias <- parsed_arguments$untreated_sample_names
  parsed_arguments$tname <- parsed_arguments$treated_name
  parsed_arguments$uname <- parsed_arguments$untreated_name
  parsed_arguments$output <- parsed_arguments$output_prefix
  parsed_arguments$batchfile <- parsed_arguments$batch_file

  return(parsed_arguments)
}

#' Parse command line arguments for DESeq analysis
#'
#' @return Parsed and validated argument list
get_args <- function() {
  # Use manual parsing as fallback if ArgumentParser is not available
  if (!requireNamespace("argparse", quietly = TRUE)) {
    parsed_args <- parse_args_manual()
    # Apply same validation and mapping as argparse version
    args <- assert_args(parsed_args)
    return(args)
  }
  
  parser <- ArgumentParser(description = "Run DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2)")
  
  # Input file parameters
  parser$add_argument(
    "-u", "--untreated_files",
    help = "Untreated (condition 1) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-t", "--treated_files",
    help = "Treated (condition 2) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-ua", "--untreated_sample_names",
    help = "Unique aliases for untreated (condition 1) expression files. Default: basenames of -u without extensions",
    type = "character",
    nargs = "*"
  )
  parser$add_argument(
    "-ta", "--treated_sample_names",
    help = "Unique aliases for treated (condition 2) expression files. Default: basenames of -t without extensions",
    type = "character",
    nargs = "*"
  )
  
  # Condition naming parameters
  parser$add_argument(
    "-un", "--untreated_name",
    help = "Name for untreated (condition 1), use only letters and numbers",
    type = "character",
    default = "untreated"
  )
  parser$add_argument(
    "-tn", "--treated_name",
    help = "Name for treated (condition 2), use only letters and numbers",
    type = "character",
    default = "treated"
  )
  
  # Batch correction parameters
  parser$add_argument(
    "-bf", "--batch_file",
    help = paste(
      "Metadata file for multi-factor analysis. Headerless TSV/CSV file.",
      "First column - names from --untreated_sample_names and --treated_sample_names, second column - batch group name.",
      "Default: None"
    ),
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.",
      "- 'model' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  
  # Statistical and filtering parameters
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--rpkm_cutoff",
    help = paste(
      "RPKM cutoff for filtering genes. Genes with RPKM values below this threshold will be excluded from the analysis.",
      "Default: NULL (no filtering)"
    ),
    type = "integer",
    default = NULL
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold for greaterAbs and |β| < lfcThreshold for lessAbs, with p-values being two-tailed or maximum of the upper and lower tests, respectively); ",
      "'up' for upregulated genes (β > lfcThreshold in condition2 compared to condition1); ",
      "'down' for downregulated genes (β < -lfcThreshold in condition2 compared to condition1). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  
  # Clustering parameters
  parser$add_argument(
    "--cluster_method",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Specifies the type of scaling to be applied to the expression data.",
      "- 'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].",
      "- 'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.",
      "- Default: zscore"
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--row_distance",
    help = paste(
      "Distance metric for HOPACH row clustering. Ignored if --cluster_method is not",
      "provided. Default: cosangle"
    ),
    type = "character",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--column_distance",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster_method is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--k",
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5.",
    type = "integer",
    default = 5
  )
  
  # Testing parameters
  parser$add_argument(
    "--test_mode",
    help = "Enable test mode for faster processing with reduced data. Default: FALSE",
    action = "store_true",
    default = FALSE
  )
  
  # Output parameters
  parser$add_argument(
    "-o", "--output_prefix",
    help = "Output prefix. Default: deseq",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-d", "--digits",
    help = "Precision, number of digits to print. Default: 3",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "-p", "--threads",
    help = "Number of threads to use for parallel processing. Default: 1",
    type = "integer",
    default = 1
  )
  
  # Parse arguments with better error handling
  parsed_args <- tryCatch({
    parser$parse_args()
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Use helper functions if available, otherwise fallback to manual parsing
    if (exists("cli_helpers") && is.environment(cli_helpers)) {
      message("Using CLI helper functions for manual parsing")
      
      # Parse using helpers
      parsed_args <- list()
      
      # Multi-value arguments  
      parsed_args$untreated_files <- cli_helpers$parse_multi_value_args(all_args, "untreated_files")
      parsed_args$treated_files <- cli_helpers$parse_multi_value_args(all_args, "treated_files")
      parsed_args$untreated_sample_names <- cli_helpers$parse_multi_value_args(all_args, "untreated_sample_names")
      parsed_args$treated_sample_names <- cli_helpers$parse_multi_value_args(all_args, "treated_sample_names")
      
      # Try short flags too
      if (length(parsed_args$untreated_files) == 0) {
        parsed_args$untreated_files <- cli_helpers$parse_multi_value_args(all_args, "u")
      }
      if (length(parsed_args$treated_files) == 0) {
        parsed_args$treated_files <- cli_helpers$parse_multi_value_args(all_args, "t")
      }
      
      # Optional single-value arguments with defaults
      optional_args <- list(
        untreated_name = "untreated",
        treated_name = "treated",
        batchcorrection = "none",
        regulation = "both", 
        cluster_method = "none",
        scaling_type = "zscore",
        row_distance = "cosangle",
        column_distance = "euclid",
        output_prefix = "./deseq"
      )
      
      for (arg in names(optional_args)) {
        parsed_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg, optional_args[[arg]])
      }
      
      # Try short flags for some args
      if (is.null(parsed_args$batch_file)) {
        parsed_args$batch_file <- cli_helpers$parse_single_value_arg(all_args, "bf")
      }
      if (parsed_args$output_prefix == "./deseq") {
        parsed_args$output_prefix <- cli_helpers$parse_single_value_arg(all_args, "o", "./deseq") 
      }
      
      # Numeric arguments
      numeric_args <- list(fdr = "double", lfcthreshold = "double", k = "integer", kmax = "integer", threads = "integer", digits = "integer", rpkm_cutoff = "integer")
      numeric_defaults <- list(fdr = 0.1, lfcthreshold = 0.59, k = 3, kmax = 5, threads = 1, digits = 3, rpkm_cutoff = NULL)
      numeric_values <- cli_helpers$parse_numeric_args(all_args, numeric_args, numeric_defaults)
      parsed_args <- c(parsed_args, numeric_values)
      
      # Try short flag for threads
      if (parsed_args$threads == 1) {
        threads_val <- cli_helpers$parse_single_value_arg(all_args, "p", "1")
        parsed_args$threads <- as.integer(threads_val)
      }
      
      # Boolean flags
      boolean_flags <- c("use_lfc_thresh", "test_mode")
      boolean_values <- cli_helpers$parse_boolean_flags(all_args, boolean_flags)
      parsed_args <- c(parsed_args, boolean_values)
      
    } else {
      # Fallback to complete manual parsing
      message("CLI helpers not available, using complete manual parsing")
      parsed_args <- list()
      
      # Parse multi-value arguments manually
      parse_multi_args <- function(all_args, flag_names) {
        result <- character(0)
        for (flag in flag_names) {
          flag_idx <- which(all_args == flag)
          if (length(flag_idx) > 0) {
            start_idx <- flag_idx[1] + 1
            end_idx <- start_idx
            while (end_idx <= length(all_args) && !startsWith(all_args[end_idx], "--") && !startsWith(all_args[end_idx], "-")) {
              end_idx <- end_idx + 1
            }
            if (start_idx < end_idx) {
              result <- c(result, all_args[start_idx:(end_idx - 1)])
            }
          }
        }
        return(result)
      }
      
      # Parse single-value arguments manually
      parse_single_arg <- function(all_args, flag_names, default_val = NULL) {
        for (flag in flag_names) {
          flag_idx <- which(all_args == flag)
          if (length(flag_idx) > 0 && flag_idx[1] < length(all_args) && !startsWith(all_args[flag_idx[1] + 1], "--")) {
            return(all_args[flag_idx[1] + 1])
          }
        }
        return(default_val)
      }
      
      # Parse boolean flags manually
      parse_bool_arg <- function(all_args, flag_names) {
        for (flag in flag_names) {
          flag_idx <- which(all_args == flag)
          if (length(flag_idx) > 0) {
            if (flag_idx[1] < length(all_args) && !startsWith(all_args[flag_idx[1] + 1], "--")) {
              val <- all_args[flag_idx[1] + 1]
              return(toupper(val) == "TRUE")
            } else {
              return(TRUE)
            }
          }
        }
        return(FALSE)
      }
      
      # Extract all required and optional arguments
      parsed_args$untreated_files <- parse_multi_args(all_args, c("-u", "--untreated_files"))
      parsed_args$treated_files <- parse_multi_args(all_args, c("-t", "--treated_files"))
      parsed_args$untreated_sample_names <- parse_multi_args(all_args, c("-ua", "--untreated_sample_names"))
      parsed_args$treated_sample_names <- parse_multi_args(all_args, c("-ta", "--treated_sample_names"))
      
      # Single-value arguments with defaults
      parsed_args$untreated_name <- parse_single_arg(all_args, c("-un", "--untreated_name"), "untreated")
      parsed_args$treated_name <- parse_single_arg(all_args, c("-tn", "--treated_name"), "treated")
      parsed_args$batch_file <- parse_single_arg(all_args, c("-bf", "--batch_file"))
      parsed_args$batchcorrection <- parse_single_arg(all_args, c("--batchcorrection"), "none")
      parsed_args$regulation <- parse_single_arg(all_args, c("--regulation"), "both")
      parsed_args$cluster_method <- parse_single_arg(all_args, c("--cluster_method"), "none")
      parsed_args$scaling_type <- parse_single_arg(all_args, c("--scaling_type"), "zscore")
      parsed_args$row_distance <- parse_single_arg(all_args, c("--row_distance"), "cosangle")
      parsed_args$column_distance <- parse_single_arg(all_args, c("--column_distance"), "euclid")
      parsed_args$output_prefix <- parse_single_arg(all_args, c("-o", "--output_prefix"), "./deseq")
      
      # Numeric arguments
      parsed_args$fdr <- as.numeric(parse_single_arg(all_args, c("--fdr"), "0.1"))
      parsed_args$lfcthreshold <- as.numeric(parse_single_arg(all_args, c("--lfcthreshold"), "0.59"))
      parsed_args$k <- as.integer(parse_single_arg(all_args, c("--k"), "3"))
      parsed_args$kmax <- as.integer(parse_single_arg(all_args, c("--kmax"), "5"))
      parsed_args$threads <- as.integer(parse_single_arg(all_args, c("-p", "--threads"), "1"))
      parsed_args$digits <- as.integer(parse_single_arg(all_args, c("-d", "--digits"), "3"))
      rpkm_val <- parse_single_arg(all_args, c("--rpkm_cutoff"))
      parsed_args$rpkm_cutoff <- if (is.null(rpkm_val)) NULL else as.integer(rpkm_val)
      
      # Boolean flags
      parsed_args$use_lfc_thresh <- parse_bool_arg(all_args, c("--use_lfc_thresh"))
      parsed_args$test_mode <- parse_bool_arg(all_args, c("--test_mode"))
      
      message("Directly extracted required argument: untreated_files = ", paste(parsed_args$untreated_files, collapse=", "))
      message("Directly extracted required argument: treated_files = ", paste(parsed_args$treated_files, collapse=", "))
      if (length(parsed_args$untreated_sample_names) > 0) {
        message("Extracted optional array argument: untreated_sample_names = ", paste(parsed_args$untreated_sample_names, collapse=", "))
      }
      if (length(parsed_args$treated_sample_names) > 0) {
        message("Extracted optional array argument: treated_sample_names = ", paste(parsed_args$treated_sample_names, collapse=", "))
      }
    }
    
    message("Manually parsed arguments using helpers")
    return(parsed_args)
  })
  
  # Debug: Check what parsed_args is before calling assert_args
  message("DEBUG: Before assert_args - parsed_args class = ", class(parsed_args))
  message("DEBUG: Before assert_args - parsed_args type = ", typeof(parsed_args))
  message("DEBUG: Before assert_args - is.function(parsed_args) = ", is.function(parsed_args))
  message("DEBUG: Before assert_args - is.list(parsed_args) = ", is.list(parsed_args))
  
  # Validate arguments and set defaults
  args <- assert_args(parsed_args)
  
  # Convert numeric values
  for (arg_name in c("fdr", "lfcthreshold", "threads")) {
    if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
      if (grepl("^[0-9.]+$", args[[arg_name]])) {
        args[[arg_name]] <- as.numeric(args[[arg_name]])
      }
    }
  }
  
  return(args)
}

#' Manual argument parsing function as fallback when argparse is not available
#'
#' @return Parsed argument list
parse_args_manual <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize result with defaults
  result <- list(
    untreated_files = NULL,
    treated_files = NULL, 
    untreated_sample_names = NULL,
    treated_sample_names = NULL,
    untreated_name = "untreated",
    treated_name = "treated",
    batch_file = NULL,
    batchcorrection = "none",
    fdr = 0.1,
    lfcthreshold = 0.59,
    use_lfc_thresh = FALSE,
    regulation = "both",
    scaling_type = "zscore",
    cluster_method = "none",
    row_distance = "cosangle",
    column_distance = "euclid",
    k = 3,
    kmax = 5,
    rpkm_cutoff = NULL,
    output_prefix = "deseq",
    threads = 1,
    digits = 3,
    test_mode = FALSE
  )
  
  # Parse arguments manually
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    
    if (arg %in% c("-u", "--untreated_files")) {
      result$untreated_files <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$untreated_files <- c(result$untreated_files, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-t", "--treated_files")) {
      result$treated_files <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$treated_files <- c(result$treated_files, args[j])
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
    } else if (arg %in% c("-ta", "--treated_sample_names")) {
      result$treated_sample_names <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[j], "-")) {
        result$treated_sample_names <- c(result$treated_sample_names, args[j])
        j <- j + 1
      }
      i <- j
    } else if (arg %in% c("-un", "--untreated_name")) {
      result$untreated_name <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-tn", "--treated_name")) {
      result$treated_name <- args[i + 1]
      i <- i + 2
    } else if (arg %in% c("-o", "--output_prefix")) {
      result$output_prefix <- args[i + 1]
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
    } else {
      # Skip unknown arguments
      i <- i + 1
    }
  }
  
  return(result)
} 