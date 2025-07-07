#!/usr/bin/env Rscript

# --- Command line argument parsing functions ---

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

# Function to parse command line arguments
get_args <- function() {
  # Get raw command line args for backup
  raw_args <- commandArgs(trailingOnly = TRUE)
  
  parser <- argparse::ArgumentParser(
    description = "Run DESeq2 analysis with Likelihood Ratio Test (LRT)",
    formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
  )
  
  # Input data arguments
  parser$add_argument(
    "--input",
    type = "character",
    required = TRUE,
    nargs = "+",
    help = "List of input files with expression data (CSV or TSV format)"
  )
  parser$add_argument(
    "--name",
    type = "character",
    required = TRUE,
    nargs = "+",
    help = "Names for input files (in the same order as input files)"
  )
  parser$add_argument(
    "--meta",
    type = "character",
    required = TRUE,
    help = "Metadata file in CSV or TSV format"
  )
  
  # Analysis parameters
  parser$add_argument(
    "--design",
    type = "character",
    required = TRUE,
    help = "Design formula for DESeq2 (e.g., '~condition+batch')"
  )
  parser$add_argument(
    "--reduced",
    type = "character",
    required = TRUE,
    help = "Reduced design formula for LRT (e.g., '~batch')"
  )
  
  # CWL-aligned parameters
  parser$add_argument(
    "--batchcorrection",
    type = "character",
    default = "none",
    choices = c("none", "combatseq", "model"),
    help = "Batch correction method: 'none', 'combatseq', or 'model'"
  )
  
  parser$add_argument(
    "--scaling_type",
    type = "character",
    default = "zscore",
    choices = c("minmax", "zscore"),
    help = "Scaling type for expression data: 'minmax' or 'zscore'"
  )
  
  parser$add_argument(
    "--fdr",
    type = "double",
    default = 0.1,
    help = "FDR threshold for significance"
  )
  
  parser$add_argument(
    "--lfcthreshold",
    type = "double",
    default = 0.59,
    help = "Log2 fold change threshold for determining significant differential expression"
  )
  
  parser$add_argument(
    "--use_lfc_thresh",
    action = "store_true",
    default = FALSE,
    help = "Use lfcthreshold as the null hypothesis value in the results function call"
  )
  
  parser$add_argument(
    "--rpkm_cutoff",
    type = "integer",
    default = NULL,
    help = "Integer cutoff for filtering rows in the expression data"
  )
  
  # Using the names directly as in CWL
  parser$add_argument(
    "--cluster",
    type = "character",
    default = "none",
    choices = c("row", "column", "both", "none"),
    help = "Hopach clustering method to be run on normalized read counts"
  )
  
  parser$add_argument(
    "--rowdist",
    type = "character",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH row clustering"
  )
  
  parser$add_argument(
    "--columndist",
    type = "character",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH column clustering"
  )
  
  parser$add_argument(
    "--k",
    type = "integer",
    default = 3,
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15"
  )
  
  parser$add_argument(
    "--kmax",
    type = "integer",
    default = 5,
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9"
  )
  
  # Output arguments
  parser$add_argument(
    "--output",
    type = "character",
    default = "./deseq_lrt_step_1",
    help = "Output prefix for generated files"
  )
  
  parser$add_argument(
    "--threads",
    type = "integer",
    default = 1,
    help = "Number of threads to use for parallel processing"
  )
  
  parser$add_argument(
    "--lrt_only_mode",
    action = "store_true",
    default = FALSE,
    help = "Run LRT only, no contrasts"
  )
  
  parser$add_argument(
    "--test_mode",
    action = "store_true",
    default = FALSE,
    help = "Run for test, only first 500 rows"
  )
  
  # Parse arguments safely with error handling
  tryCatch({
    args <- parser$parse_args()
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Use helper functions if available, otherwise fallback to manual parsing
    if (exists("cli_helpers") && is.environment(cli_helpers)) {
      message("Using CLI helper functions for manual parsing")
      
      # Parse using helpers
      args <- list()
      
      # Multi-value arguments
      args$input <- cli_helpers$parse_multi_value_args(all_args, "input")
      args$name <- cli_helpers$parse_multi_value_args(all_args, "name")
      
      # Required single-value arguments
      args$meta <- cli_helpers$parse_single_value_arg(all_args, "meta")
      args$design <- cli_helpers$parse_single_value_arg(all_args, "design")
      args$reduced <- cli_helpers$parse_single_value_arg(all_args, "reduced")
      
      # Optional single-value arguments with defaults
      optional_args <- list(
        batchcorrection = "none",
        scaling_type = "zscore",
        cluster = "none",
        rowdist = "cosangle",
        columndist = "euclid",
        output = "./deseq_lrt_step_1"
      )
      
      for (arg in names(optional_args)) {
        args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg, optional_args[[arg]])
      }
      
      # Numeric arguments
      numeric_args <- list(fdr = "double", lfcthreshold = "double", k = "integer", kmax = "integer", threads = "integer", rpkm_cutoff = "integer")
      numeric_defaults <- list(fdr = 0.1, lfcthreshold = 0.59, k = 3, kmax = 5, threads = 1, rpkm_cutoff = NULL)
      numeric_values <- cli_helpers$parse_numeric_args(all_args, numeric_args, numeric_defaults)
      args <- c(args, numeric_values)
      
      # Boolean flags
      boolean_flags <- c("use_lfc_thresh", "lrt_only_mode", "test_mode")
      boolean_values <- cli_helpers$parse_boolean_flags(all_args, boolean_flags)
      args <- c(args, boolean_values)
      
    } else {
      # Fallback to original defaults
      message("CLI helpers not available, using minimal manual parsing")
      args <- list(
        input = character(0),
        name = character(0),
        meta = NULL,
        design = NULL,
        reduced = NULL,
        batchcorrection = "none",
        scaling_type = "zscore",
        fdr = 0.1,
        lfcthreshold = 0.59,
        use_lfc_thresh = FALSE,
        rpkm_cutoff = NULL,
        cluster = "none",
        rowdist = "cosangle",
        columndist = "euclid",
        k = 3,
        kmax = 5,
        output = "./deseq_lrt_step_1",
        threads = 1,
        lrt_only_mode = FALSE,
        test_mode = FALSE
      )
      
      # Simple fallback parsing for required args only
      for (req_arg in c("meta", "design", "reduced")) {
        req_flag <- paste0("--", req_arg)
        arg_idx <- which(all_args == req_flag)
        if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
          args[[req_arg]] <- all_args[arg_idx[1] + 1]
        }
      }
    }
    
    message("Manually parsed arguments using helpers")
  })
  
  # Validate arguments
  args <- assert_args(args)
  
  # Return the parsed and validated arguments
  return(args)
}

# Function to validate command line arguments
assert_args <- function(args) {
  # Check for required arguments
  required_args <- c("meta", "design", "reduced")
  missing_required <- required_args[!required_args %in% names(args)]
  
  if (length(missing_required) > 0) {
    # Before failing, check if any required arguments may have been parsed 
    # correctly but are not accessible due to variable scope issues
    # Try to extract these from the original command arguments if possible
    all_args <- commandArgs(trailingOnly = TRUE)
    
    for (missing_arg in missing_required) {
      arg_flag <- paste0("--", missing_arg)
      arg_idx <- which(all_args == arg_flag)
      if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
        args[[missing_arg]] <- all_args[arg_idx[1] + 1]
        message(paste("Recovered missing required argument:", missing_arg, "=", args[[missing_arg]]))
        # Remove from missing list
        missing_required <- setdiff(missing_required, missing_arg)
      }
    }
    
    # If we still have missing arguments, abort
    if (length(missing_required) > 0) {
      message(paste("Error: Missing required arguments:", paste(missing_required, collapse=", ")))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }
  
  # Create empty arrays for input/name if they don't exist
  if (is.null(args$input)) args$input <- character(0)
  if (is.null(args$name)) args$name <- character(0)
  
  # Validate input and name parameters
  if (length(args$input) == 0 || length(args$name) == 0) {
    message("Error: --input and --name parameters are required")
    message("Input files found: ", length(args$input))
    message("Sample names found: ", length(args$name))
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  if (length(args$input) != length(args$name)) {
    message("Error: --input and --name have different number of values")
    message(paste("Number of input files:", length(args$input)))
    message(paste("Number of sample names:", length(args$name)))
    message("Input files:")
    print(args$input)
    message("Sample names:")
    print(args$name)
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate design formula
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(args$design)
    },
    error = function(e) {
      message(paste0("Error: failed to load --design ", args$design, " as formula"))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  )
  
  # Validate reduced formula
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(args$reduced)
    },
    error = function(e) {
      message(paste0("Error: failed to load --reduced ", args$reduced, " as formula"))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  )
  
  # Check if metadata file exists
  if (!file.exists(args$meta)) {
    message(paste("Error: Metadata file does not exist:", args$meta))
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate batch correction
  if (!is.null(args$batchcorrection) && args$batchcorrection != "none") {
    # Try to read the metadata file
    meta_data <- tryCatch(
      {
        # Try to read the file
        if (endsWith(args$meta, ".csv")) {
          read.csv(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        } else if (endsWith(args$meta, ".tsv")) {
          read.delim(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        } else {
          # Default to CSV
          read.csv(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        }
      },
      error = function(e) {
        message(paste0("Error: failed to read metadata file: ", e$message))
        quit(save = "no", status = 1, runLast = FALSE)
      }
    )
    
    if (!"batch" %in% colnames(meta_data)) {
      message("Warning: batch correction requested but 'batch' column not found in metadata file")
      message("Proceeding without batch correction")
      args$batchcorrection <- "none"
    } else if (!is.numeric(meta_data$batch)) {
      message("Warning: 'batch' column in metadata file is not numeric")
      message("Converting batch to numeric values")
      # No need to quit here, we'll handle conversion elsewhere
    }
  }
  
  # Validate rpkm_cutoff (if provided)
  if (!is.null(args$rpkm_cutoff)) {
    if (!is.numeric(args$rpkm_cutoff) || args$rpkm_cutoff < 0) {
      message("Warning: --rpkm_cutoff must be a non-negative integer, using default NULL")
      args$rpkm_cutoff <- NULL
    }
  }
  
  # Validate k parameter (1-15)
  if (!is.null(args$k) && (args$k < 1 || args$k > 15)) {
    message("Warning: --k must be between 1 and 15, using default value 3")
    args$k <- 3
  }
  
  # Validate kmax parameter (2-9)
  if (!is.null(args$kmax) && (args$kmax < 2 || args$kmax > 9)) {
    message("Warning: --kmax must be between 2 and 9, using default value 5")
    args$kmax <- 5
  }
  
  # Validate fdr parameter
  if (!is.null(args$fdr) && (args$fdr < 0 || args$fdr > 1)) {
    message("Warning: --fdr must be between 0 and 1, using default value 0.1")
    args$fdr <- 0.1
  }

  return(args)
}

# This namespace contains exported parameter handling functions
params <- new.env()

#' Get command line arguments with standard parsing
#' @export
params$get_cli_args <- function() {
  return(get_args())
}

params$assert_args <- assert_args

# Export the params namespace
assign("params", params, envir = .GlobalEnv)

# Command line argument parsing for DESeq2 LRT Step 1
#
# This module handles all command-line argument parsing and validation

#' Parse command line arguments using base R functions
#' @return List of parsed arguments
#' @export
parse_cli_args <- function() {
  cat("Parsing command line arguments...\n")
  
  # Get raw command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize default values
  result <- list(
    input = NULL,
    name = NULL,
    meta = NULL,
    design = NULL,
    reduced = NULL,
    batchcorrection = "none",
    scaling_type = "zscore",
    fdr = 0.1,
    lfcthreshold = 0.59,
    use_lfc_thresh = FALSE,
    rpkm_cutoff = NULL,
    cluster = "none",
    rowdist = "cosangle",
    columndist = "euclid",
    k = 3,
    kmax = 5,
    output = "./deseq_lrt_step_1",
    threads = 1,
    lrt_only_mode = FALSE,
    test_mode = FALSE
  )
  
  # Manual parsing of command line arguments
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    
    # Handle flags with multiple values (--input, --name)
    if (arg == "--input" || arg == "--name") {
      param_name <- substring(arg, 3)  # Remove leading --
      values <- c()
      j <- i + 1
      
      # Collect all values until next flag
      while (j <= length(args) && !startsWith(args[j], "--")) {
        values <- c(values, args[j])
        j <- j + 1
      }
      
      # Store values and update index
      result[[param_name]] <- values
      i <- j
    }
    # Handle boolean flags
    else if (arg == "--use_lfc_thresh" || arg == "--lrt_only_mode" || arg == "--test_mode") {
      param_name <- substring(arg, 3)  # Remove leading --
      result[[param_name]] <- TRUE
      i <- i + 1
    }
    # Handle flags with single values
    else if (startsWith(arg, "--") && i < length(args) && !startsWith(args[i+1], "--")) {
      param_name <- substring(arg, 3)  # Remove leading --
      value <- args[i+1]
      
      # Convert numeric values
      if (param_name %in% c("fdr", "lfcthreshold")) {
        value <- as.numeric(value)
      }
      else if (param_name %in% c("rpkm_cutoff", "k", "kmax", "threads")) {
        value <- as.integer(value)
      }
      
      result[[param_name]] <- value
      i <- i + 2
    }
    else {
      # Skip unknown arguments
      cat("WARNING: Skipping unknown argument:", arg, "\n")
      i <- i + 1
    }
  }
  
  # Validate required parameters
  missing_params <- c()
  for (param in c("input", "name", "meta", "design", "reduced")) {
    if (is.null(result[[param]])) {
      missing_params <- c(missing_params, param)
    }
  }
  
  if (length(missing_params) > 0) {
    stop(paste("Missing required parameters:", paste(missing_params, collapse=", ")))
  }
  
  # Ensure input files and sample names match in length
  if (length(result$input) != length(result$name)) {
    cat(paste("WARNING: Number of input files (", length(result$input), 
              ") does not match number of sample names (", length(result$name), ")\n", sep=""))
  }
  
  cat("Command line arguments parsed successfully\n")
  
  # Print the parsed arguments
  cat("Parsed parameters:\n")
  for (name in names(result)) {
    if (name %in% c("input", "name")) {
      cat(paste("  ", name, ":", paste(result[[name]][1:min(3, length(result[[name]]))], collapse=", "), 
                if(length(result[[name]]) > 3) "..." else "", "\n", sep=""))
    } else {
      cat(paste("  ", name, ":", result[[name]], "\n", sep=""))
    }
  }
  
  return(result)
} 