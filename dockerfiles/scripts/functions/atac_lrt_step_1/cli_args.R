#!/usr/bin/env Rscript

# --- Command line argument parsing functions for ATAC-seq LRT Step 1 ---

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
    description = "Run ATAC-seq analysis with Likelihood Ratio Test (LRT) using DiffBind and DESeq2",
    formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
  )
  
  # Input data arguments - adapted for ATAC-seq
  parser$add_argument(
    "--input_files",
    type = "character",
    required = TRUE,
    nargs = "+",
    help = "List of input files with peak data (CSV or TSV format)"
  )
  parser$add_argument(
    "--name",
    type = "character",
    required = TRUE,
    nargs = "+",
    help = "Names for input files (in the same order as input files)"
  )
  parser$add_argument(
    "--bamfiles",
    type = "character",
    required = TRUE,
    nargs = "+",
    help = "List of BAM files corresponding to peak files"
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
  
  # DiffBind-specific parameters
  parser$add_argument(
    "--minoverlap",
    type = "integer",
    default = 2,
    help = "Minimum overlap for consensus peaks (default: 2)"
  )
  parser$add_argument(
    "--peakformat",
    type = "character",
    default = "macs",
    choices = c("csv", "bed", "narrow", "macs"),
    help = "Peak file format: csv, bed, narrow, macs (default: macs for .xls files)"
  )
  parser$add_argument(
    "--peakcaller",
    type = "character",
    default = "macs",
    choices = c("macs", "spp", "bed", "bayes", "peakseq"),
    help = "Peak caller used: macs, spp, bed, bayes, peakseq (default: macs)"
  )
  parser$add_argument(
    "--scorecol",
    type = "integer",
    default = 6,
    help = "Column number containing peak scores (default: 6)"
  )
  
  # CWL-aligned parameters
  parser$add_argument(
    "--batchcorrection",
    type = "character",
    default = "none",
    choices = c("none", "combatseq", "limmaremovebatcheffect"),
    help = "Batch correction method (default: none)"
  )
  
  parser$add_argument(
    "--scaling_type",
    type = "character",
    default = "zscore",
    choices = c("minmax", "zscore"),
    help = "Scaling type for accessibility data: 'minmax' or 'zscore'"
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
    help = "Log2 fold change threshold for determining significant differential accessibility"
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
    help = "Integer cutoff for filtering rows in the accessibility data"
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
    default = "./atac_lrt_step_1",
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
    help = "Run for test, only first 500 peaks"
  )
  
  # Add interaction design flag
  parser$add_argument(
    "--with_interaction",
    action = "store_true",
    default = FALSE,
    help = "Use interaction design (e.g., ~Condition*Tissue)"
  )
  
  # Enhance batch correction options
  parser$add_argument(
    "--cluster_method",
    type = "character",
    default = "row",
    choices = c("none", "row", "column", "both"),
    help = "Clustering method for heatmaps (default: row)"
  )
  
  parser$add_argument(
    "--use_limma_correction",
    action = "store_true",
    default = FALSE,
    help = "Apply limma batch correction with interaction design"
  )
  
  # Add score type for DiffBind counting
  parser$add_argument(
    "--score_type",
    type = "character",
    default = "DBA_SCORE_RPKM",
    choices = c("DBA_SCORE_READS", "DBA_SCORE_RPKM", "DBA_SCORE_TMM_MINUS_FULL"),
    help = "Score type for DiffBind counting (default: DBA_SCORE_RPKM)"
  )
  
  # Parse arguments safely with error handling
  args <- tryCatch({
    parser$parse_args()
  }, error = function(e) {
    message("Warning: Argument parsing error. Attempting to handle arguments manually.")
    
    all_args <- commandArgs(trailingOnly = TRUE)
    
    # Use helper functions if available, otherwise fallback to manual parsing
    if (exists("cli_helpers") && is.environment(cli_helpers)) {
      message("Using CLI helper functions for manual parsing")
      
      # Parse using helpers
      parsed_args <- list()
      
      # Multi-value arguments specific to ATAC
      parsed_args$input_files <- cli_helpers$parse_multi_value_args(all_args, "input_files")
      parsed_args$name <- cli_helpers$parse_multi_value_args(all_args, "name")
      parsed_args$bamfiles <- cli_helpers$parse_multi_value_args(all_args, "bamfiles")
      
      # Required single-value arguments
      parsed_args$meta <- cli_helpers$parse_single_value_arg(all_args, "meta")
      parsed_args$design <- cli_helpers$parse_single_value_arg(all_args, "design")
      parsed_args$reduced <- cli_helpers$parse_single_value_arg(all_args, "reduced")
      
      # Optional single-value arguments with defaults
      optional_args <- list(
        peakformat = "macs",
        peakcaller = "macs",
        batchcorrection = "none",
        scaling_type = "zscore",
        cluster = "none",
        rowdist = "cosangle",
        columndist = "euclid",
        score_type = "DBA_SCORE_RPKM",
        output = "./atac_lrt_step_1"
      )
      
      for (arg in names(optional_args)) {
        parsed_args[[arg]] <- cli_helpers$parse_single_value_arg(all_args, arg, optional_args[[arg]])
      }
      
      # Numeric arguments
      numeric_args <- list(
        minoverlap = "integer", scorecol = "integer", fdr = "double", 
        lfcthreshold = "double", rpkm_cutoff = "integer", k = "integer", 
        kmax = "integer", threads = "integer"
      )
      numeric_defaults <- list(
        minoverlap = 2, scorecol = 6, fdr = 0.1, lfcthreshold = 0.59,
        rpkm_cutoff = NULL, k = 3, kmax = 5, threads = 1
      )
      numeric_values <- cli_helpers$parse_numeric_args(all_args, numeric_args, numeric_defaults)
      parsed_args <- c(parsed_args, numeric_values)
      
      # Boolean flags
      boolean_flags <- c("use_lfc_thresh", "lrt_only_mode", "test_mode", "with_interaction", "use_limma_correction")
      boolean_values <- cli_helpers$parse_boolean_flags(all_args, boolean_flags)
      parsed_args <- c(parsed_args, boolean_values)
      
    } else {
      # Fallback to original defaults
      message("CLI helpers not available, using minimal manual parsing")
      parsed_args <- list(
        input_files = character(0),
        name = character(0),
        bamfiles = character(0),
        meta = NULL,
        design = NULL,
        reduced = NULL,
        minoverlap = 2,
        peakformat = "macs",
        peakcaller = "macs",
        scorecol = 6,
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
        output = "./atac_lrt_step_1",
        threads = 1,
        lrt_only_mode = FALSE,
        test_mode = FALSE,
        with_interaction = FALSE,
        use_limma_correction = FALSE,
        score_type = "DBA_SCORE_RPKM"
      )
      
      # Simple fallback parsing for required args only
      for (req_arg in c("meta", "design", "reduced")) {
        req_flag <- paste0("--", req_arg)
        arg_idx <- which(all_args == req_flag)
        if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
          parsed_args[[req_arg]] <- all_args[arg_idx[1] + 1]
        }
      }
      
      # Parse input_files manually (CRITICAL FIX)
      input_idx <- which(all_args == "--input_files")
      if (length(input_idx) > 0) {
        start_idx <- input_idx[1] + 1
        end_idx <- start_idx
        while (end_idx <= length(all_args) && !startsWith(all_args[end_idx], "--")) {
          end_idx <- end_idx + 1
        }
        parsed_args$input_files <- all_args[start_idx:(end_idx - 1)]
      }
      
      # Find --name arguments
      name_idx <- which(all_args == "--name")
      if (length(name_idx) > 0) {
        start_idx <- name_idx[1] + 1
        end_idx <- start_idx
        while (end_idx <= length(all_args) && !startsWith(all_args[end_idx], "--")) {
          end_idx <- end_idx + 1
        }
        parsed_args$name <- all_args[start_idx:(end_idx - 1)]
      }
      
      # Find --bamfiles arguments
      bam_idx <- which(all_args == "--bamfiles")
      if (length(bam_idx) > 0) {
        start_idx <- bam_idx[1] + 1
        end_idx <- start_idx
        while (end_idx <= length(all_args) && !startsWith(all_args[end_idx], "--")) {
          end_idx <- end_idx + 1
        }
        parsed_args$bamfiles <- all_args[start_idx:(end_idx - 1)]
      }
      
      # Handle boolean flags (both --flag and --flag TRUE/FALSE formats)
      boolean_flags <- c("use_lfc_thresh", "lrt_only_mode", "test_mode", "with_interaction", 
                        "use_limma_correction")
      for (flag in boolean_flags) {
        flag_name <- paste0("--", flag)
        flag_idx <- which(all_args == flag_name)
        if (length(flag_idx) > 0) {
          # Check if there's a value after the flag
          if (flag_idx[1] < length(all_args) && !startsWith(all_args[flag_idx[1] + 1], "--")) {
            # Has a value (TRUE/FALSE)
            val <- all_args[flag_idx[1] + 1]
            parsed_args[[flag]] <- toupper(val) == "TRUE"
          } else {
            # Just the flag present (store_true action)
            parsed_args[[flag]] <- TRUE
          }
        }
      }
    }
    
    message("Manually parsed arguments using helpers")
    
    # Convert args to match expected format and return
    return(as.list(parsed_args))
  })
  
  return(args)
}

# Function to validate arguments
validate_args <- function(args) {
  message("DEBUG: Starting validate_args")
  message("DEBUG: args class:", class(args))
  message("DEBUG: args names:", paste(names(args), collapse = ", "))
  
  errors <- character(0)
  
  message("DEBUG: Checking required arguments...")
  
  # Check required arguments
  if (length(args$input_files) == 0) {
    errors <- c(errors, "Input peak files are required")
  }
  
  if (length(args$name) == 0) {
    errors <- c(errors, "Sample names are required")
  }
  
  if (length(args$bamfiles) == 0) {
    errors <- c(errors, "BAM files are required")
  }
  
  if (is.null(args$meta)) {
    errors <- c(errors, "Metadata file is required")
  }
  
  if (is.null(args$design)) {
    errors <- c(errors, "Design formula is required")
  }
  
  if (is.null(args$reduced)) {
    errors <- c(errors, "Reduced formula is required")
  }
  
  message("DEBUG: Checking array lengths...")
  
  # Check that input, name, and bamfiles have the same length
  if (length(args$input_files) != length(args$name)) {
    errors <- c(errors, "Number of input files must match number of names")
  }
  
  if (length(args$input_files) != length(args$bamfiles)) {
    errors <- c(errors, "Number of input files must match number of BAM files")
  }
  
  message("DEBUG: Checking file existence...")
  
  # Check file existence
  for (file in args$input_files) {
    if (!file.exists(file)) {
      errors <- c(errors, paste("Peak file does not exist:", file))
    }
  }
  
  for (file in args$bamfiles) {
    if (!file.exists(file)) {
      errors <- c(errors, paste("BAM file does not exist:", file))
    }
  }
  
  if (!file.exists(args$meta)) {
    errors <- c(errors, paste("Metadata file does not exist:", args$meta))
  }
  
  message("DEBUG: Checking parameter ranges...")
  
  # Check parameter ranges
  if (args$k < 1 || args$k > 15) {
    errors <- c(errors, "k must be between 1 and 15")
  }
  
  if (args$kmax < 2 || args$kmax > 9) {
    errors <- c(errors, "kmax must be between 2 and 9")
  }
  
  if (args$fdr < 0 || args$fdr > 1) {
    errors <- c(errors, "FDR must be between 0 and 1")
  }
  
  if (args$minoverlap < 1) {
    errors <- c(errors, "minoverlap must be at least 1")
  }
  
  if (args$scorecol < 1) {
    errors <- c(errors, "scorecol must be at least 1")
  }
  
  message("DEBUG: validate_args completed successfully")
  
  if (length(errors) > 0) {
    stop("Argument validation failed:\n", paste(errors, collapse = "\n"))
  }
  
  return(TRUE)
}

# Function to print arguments summary
print_args <- function(args) {
  cat("ATAC-seq LRT Step 1 Analysis Parameters:\n")
  cat("=====================================\n")
  cat("Peak files:", length(args$input_files), "files\n")
  cat("BAM files:", length(args$bamfiles), "files\n")
  cat("Sample names:", paste(args$name, collapse = ", "), "\n")
  cat("Metadata file:", args$meta, "\n")
  cat("Design formula:", args$design, "\n")
  cat("Reduced formula:", args$reduced, "\n")
  cat("Peak format:", args$peakformat, "\n")
  cat("Peak caller:", args$peakcaller, "\n")
  cat("Score column:", args$scorecol, "\n")
  cat("Minimum overlap:", args$minoverlap, "\n")
  cat("Batch correction:", args$batchcorrection, "\n")
  cat("FDR threshold:", args$fdr, "\n")
  cat("LFC threshold:", args$lfcthreshold, "\n")
  cat("Use LFC threshold:", args$use_lfc_thresh, "\n")
  cat("Clustering method:", args$cluster, "\n")
  cat("Row distance:", args$rowdist, "\n")
  cat("Column distance:", args$columndist, "\n")
  cat("Hopach k:", args$k, "\n")
  cat("Hopach kmax:", args$kmax, "\n")
  cat("Output prefix:", args$output, "\n")
  cat("Threads:", args$threads, "\n")
  cat("LRT only mode:", args$lrt_only_mode, "\n")
  cat("Test mode:", args$test_mode, "\n")
  cat("With interaction:", args$with_interaction, "\n")
  cat("Cluster method:", args$cluster, "\n")
  cat("Use limma correction:", args$use_limma_correction, "\n")
  cat("Score type:", args$score_type, "\n")
  cat("=====================================\n")
} 