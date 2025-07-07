#!/usr/bin/env Rscript
#
# Base workflow functions for ATAC-seq workflows
# This module provides common workflow patterns and environment setup
#

#' Initialize ATAC-seq environment with common setup
#'
#' @param workflow_type Type of ATAC workflow for logging
#' @param specific_modules Vector of workflow-specific module paths to source
#' @export
initialize_atac_environment <- function(workflow_type = "ATAC-seq", specific_modules = NULL) {
  # Display startup message
  message(paste("Starting", workflow_type, "Analysis"))
  message("Working directory:", getwd())
  
  # Print command line arguments for debugging purposes
  raw_args <- commandArgs(trailingOnly = TRUE)
  message("Command line arguments received:")
  message(paste(raw_args, collapse = " "))
  
  # First, make sure we have the utilities module
  if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
    message("Loading utilities from Docker path: /usr/local/bin/functions/common/utilities.R")
    source("/usr/local/bin/functions/common/utilities.R")
  } else if (file.exists("functions/common/utilities.R")) {
    message("Loading utilities from relative path: functions/common/utilities.R")
    source("functions/common/utilities.R")
  } else {
    # Try one more location
    script_dir <- tryCatch({
      dirname(sys.frame(1)$ofile)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(script_dir)) {
      potential_path <- file.path(script_dir, "../common/utilities.R")
      if (file.exists(potential_path)) {
        message(paste("Loading utilities from script relative path:", potential_path))
        source(potential_path)
      } else {
        stop("Could not find utilities.R file")
      }
    } else {
      stop("Could not find utilities.R file")
    }
  }
  
  # Now we have access to source_with_fallback and other utilities
  # Source common functions
  source_with_fallback("functions/common/constants.R", "/usr/local/bin/functions/common/constants.R")
  source_with_fallback("functions/common/output_utils.R", "/usr/local/bin/functions/common/output_utils.R")
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
  source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
  source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")

  # Source ATAC common functions
  source_with_fallback("functions/common/atac_common/cli_args_base.R", "/usr/local/bin/functions/common/atac_common/cli_args_base.R")
  source_with_fallback("functions/common/atac_common/data_processing_base.R", "/usr/local/bin/functions/common/atac_common/data_processing_base.R")
  source_with_fallback("functions/common/atac_common/diffbind_utils.R", "/usr/local/bin/functions/common/atac_common/diffbind_utils.R")
  
  # Source workflow-specific modules if provided
  if (!is.null(specific_modules)) {
    for (module in specific_modules) {
      # Try both relative and Docker paths
      relative_path <- paste0("functions/", module)
      docker_path <- paste0("/usr/local/bin/functions/", module)
      source_with_fallback(relative_path, docker_path)
    }
  }
  
  # Load required libraries with clear error messages
  load_atac_libraries()

  # Configure R options
  configure_r_options()
  
  # Configure plot theme
  configure_plot_theme()
  
  log_message(paste("Environment initialized for", workflow_type, "analysis"))
}

#' Load required libraries for ATAC-seq analysis
#'
#' @export
load_atac_libraries <- function() {
  tryCatch({
    message("Loading required libraries...")
    suppressPackageStartupMessages({
      required_packages <- c(
        "DiffBind",
        "DESeq2",
        "BiocParallel",
        "data.table",
        "ggplot2",
        "plotly",
        "limma",
        "hopach",
        "stringr",
        "GenomicRanges",
        "rtracklayer",
        "Rsamtools"
      )
      
      for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          stop(paste("Required package not found:", pkg))
        }
        message(paste("Loading package:", pkg))
        library(pkg, character.only = TRUE)
      }
    })
  }, error = function(e) {
    stop(paste("Error loading libraries:", e$message))
  })
}

#' Main workflow wrapper with memory management for ATAC-seq
#'
#' @param workflow_function Main workflow function to execute
#' @return Results from the workflow function
#' @export
main_atac_workflow_wrapper <- function(workflow_function) {
  # Report initial memory usage
  report_memory_usage("Initial")
  
  # Get command-line arguments
  args <- get_args()
  
  # Run the main workflow
  results <- tryCatch({
    workflow_function(args)
  }, error = function(e) {
    log_message(paste("Error in ATAC-seq workflow execution:", e$message), "ERROR")
    # Force garbage collection
    gc(verbose = FALSE)
    # Re-raise the error
    stop(e)
  })
  
  # Final memory usage report
  report_memory_usage("Final")
  
  return(results)
}

#' Run base DiffBind analysis pipeline common to all ATAC workflows
#'
#' @param sample_sheet DiffBind sample sheet
#' @param args Command-line arguments
#' @param test_mode_handling Whether to include special test mode handling
#' @return List containing DBA object and consensus peaks
#' @export
run_base_diffbind_analysis <- function(sample_sheet, args, test_mode_handling = TRUE) {
  log_message("Running base DiffBind analysis pipeline...")
  
  # Pre-validation: Check file accessibility before calling DiffBind
  if (test_mode_handling) {
    validate_files_for_diffbind(sample_sheet, args$test_mode)
  }
  
  # Create DBA object with comprehensive error handling
  log_message("Creating DBA object...")
  log_message("DiffBind parameters:")
  log_message(paste("  - peakFormat:", args$peakformat))
  log_message(paste("  - peakCaller:", args$peakcaller))
  log_message(paste("  - scoreCol:", args$scorecol))
  
  dba_obj <- tryCatch({
    dba(
      sampleSheet = sample_sheet,
      peakFormat = args$peakformat,
      peakCaller = args$peakcaller,
      scoreCol = args$scorecol
    )
  }, error = function(e) {
    handle_dba_creation_error(e, sample_sheet, args)
  })
  
  # Handle test mode with special DBA objects
  if (test_mode_handling && args$test_mode && !inherits(dba_obj, "DBA")) {
    return(handle_test_mode_diffbind(dba_obj, args))
  }
  
  # Create consensus peaks
  log_message("Creating consensus peaks...")
  dba_consensus <- dba.peakset(dba_obj, consensus = DBA_CONDITION, minOverlap = args$minoverlap)
  dba_consensus <- dba(dba_consensus, mask = dba_consensus$masks$Consensus, minOverlap = 1)
  
  # Get consensus peaks
  consensus_peaks <- dba.peakset(dba_consensus, bRetrieve = TRUE, minOverlap = 1)
  
  # Count reads in consensus peaks
  log_message("Counting reads in consensus peaks...")
  dba_obj <- count_reads_with_error_handling(dba_obj, consensus_peaks, args, test_mode_handling)
  
  # Apply test mode filtering if requested
  if (args$test_mode && test_mode_handling) {
    dba_obj <- apply_test_mode_filtering(dba_obj, args)
  }
  
  return(list(dba_obj = dba_obj, consensus_peaks = consensus_peaks))
}

#' Validate files for DiffBind with test mode considerations
#'
#' @param sample_sheet DiffBind sample sheet
#' @param test_mode Whether running in test mode
validate_files_for_diffbind <- function(sample_sheet, test_mode = FALSE) {
  message("Pre-validating files before DiffBind analysis...")
  
  for (i in 1:nrow(sample_sheet)) {
    # Check BAM file
    bam_file <- sample_sheet$bamReads[i]
    if (!file.exists(bam_file)) {
      stop(sprintf("BAM file not found: %s (sample: %s)", bam_file, sample_sheet$SampleID[i]))
    }
    
    # Check peak file
    peak_file <- sample_sheet$Peaks[i]
    if (!file.exists(peak_file)) {
      stop(sprintf("Peak file not found: %s (sample: %s)", peak_file, sample_sheet$SampleID[i]))
    }
    
    # Check file sizes (warn about very small files)
    bam_size <- file.info(bam_file)$size
    peak_size <- file.info(peak_file)$size
    
    if (test_mode) {
      if (bam_size < 1000) {
        message(sprintf("Test mode: BAM file %s is very small (%d bytes) - may be a dummy file", 
                        basename(bam_file), bam_size))
      }
      
      if (peak_size < 100) {
        message(sprintf("Test mode: Peak file %s is very small (%d bytes) - may be minimal", 
                        basename(peak_file), peak_size))
      }
    } else {
      if (bam_size < 1000) {
        warning(sprintf("BAM file %s is very small (%d bytes) - may be a dummy file", 
                        basename(bam_file), bam_size))
      }
      
      if (peak_size < 100) {
        warning(sprintf("Peak file %s is very small (%d bytes) - may be empty", 
                        basename(peak_file), peak_size))
      }
    }
  }
}

#' Handle DBA object creation errors
#'
#' @param e Error object
#' @param sample_sheet DiffBind sample sheet
#' @param args Command-line arguments
#' @return DBA object or mock object for test mode
handle_dba_creation_error <- function(e, sample_sheet, args) {
  message("ERROR in dba() function:")
  message("Error message: ", e$message)
  
  # For test mode with dummy files, create a minimal mock DBA object
  if (args$test_mode) {
    message("Test mode: creating mock DBA object due to file issues...")
    
    # Create a minimal mock DBA object structure
    mock_dba <- list(
      samples = sample_sheet,
      class = c("DBA"),
      config = list(fragmentSize = 0, th = 0.05)
    )
    
    return(mock_dba)
  }
  
  stop("Failed to create DBA object. See error details above.")
}

#' Handle test mode DiffBind analysis
#'
#' @param mock_dba Mock DBA object
#' @param args Command-line arguments
#' @return Mock results for testing
handle_test_mode_diffbind <- function(mock_dba, args) {
  message("Test mode: bypassing DiffBind analysis due to dummy files")
  message("Creating mock results for testing...")
  
  # Create minimal mock results
  mock_results <- list(
    dba_obj = mock_dba,
    consensus_peaks = GenomicRanges::GRanges(
      seqnames = rep("chr1", 100),
      ranges = IRanges::IRanges(start = seq(1000, 100000, 1000), width = 500)
    )
  )
  
  return(mock_results)
}

#' Count reads with error handling for test mode
#'
#' @param dba_obj DBA object
#' @param consensus_peaks Consensus peaks
#' @param args Command-line arguments
#' @param test_mode_handling Whether to handle test mode specially
#' @return Updated DBA object with counts
count_reads_with_error_handling <- function(dba_obj, consensus_peaks, args, test_mode_handling) {
  if (test_mode_handling && args$test_mode) {
    message("Test mode: checking BAM file validity...")
    
    # Check if BAM files are dummy/invalid
    bam_files <- dba_obj$samples$bamReads
    bam_sizes <- sapply(bam_files, function(f) file.info(f)$size)
    
    if (any(bam_sizes < 1000)) {
      message("Test mode: dummy BAM files detected, using minimal peak counting...")
      
      # Create a minimal DBA object for testing without full BAM processing
      tryCatch({
        dba_obj <- dba.count(dba_obj, peaks = consensus_peaks, minOverlap = 1)
      }, error = function(e) {
        message("BAM processing failed as expected with dummy files. Creating mock count matrix for testing...")
        
        # Create mock count data
        n_peaks <- length(consensus_peaks)
        n_samples <- nrow(dba_obj$samples)
        
        mock_counts <- matrix(
          sample(1:100, n_peaks * n_samples, replace = TRUE),
          nrow = n_peaks,
          ncol = n_samples
        )
        colnames(mock_counts) <- dba_obj$samples$SampleID
        
        # Add the mock data to the DBA object structure
        dba_obj$binding <- mock_counts
        dba_obj$peaks <- as.data.frame(consensus_peaks)[1:n_peaks, ]
        
        message(paste("Created mock binding matrix with", n_peaks, "peaks and", n_samples, "samples"))
        return(dba_obj)
      })
    } else {
      # Real BAM files in test mode
      dba_obj <- dba.count(dba_obj, peaks = consensus_peaks, minOverlap = 1)
    }
  } else {
    # Production mode - full BAM processing
    dba_obj <- dba.count(dba_obj, peaks = consensus_peaks, minOverlap = 1)
  }
  
  return(dba_obj)
}

#' Apply test mode filtering to reduce dataset size
#'
#' @param dba_obj DBA object
#' @param args Command-line arguments
#' @return Filtered DBA object
apply_test_mode_filtering <- function(dba_obj, args) {
  message("Test mode: reducing to first 500 peaks for faster processing...")
  
  # Get binding matrix
  binding_matrix <- dba_obj$binding
  n_peaks <- min(500, nrow(binding_matrix))
  
  # Filter the binding matrix
  dba_obj$binding <- binding_matrix[1:n_peaks, , drop = FALSE]
  
  # Filter peaks data if it exists
  if (!is.null(dba_obj$peaks)) {
    dba_obj$peaks <- dba_obj$peaks[1:n_peaks, , drop = FALSE]
  }
  
  message(paste("Reduced dataset to", n_peaks, "peaks for test mode"))
  return(dba_obj)
}

#' Export ATAC-seq results with standard formatting
#'
#' @param results_df Results data frame
#' @param output_prefix Output file prefix
#' @param workflow_type Type of workflow for file naming
#' @export
export_atac_results <- function(results_df, output_prefix, workflow_type = "atac") {
  log_message("Exporting ATAC-seq results...")
  
  # Export main results table
  results_file <- paste0(output_prefix, "_", workflow_type, "_results.tsv")
  write_atac_results(results_df, results_file)
  log_message(paste("Exported ATAC-seq results table to", results_file))
  
  return(results_file)
}

#' Create standard ATAC-seq visualizations
#'
#' @param dds DESeq2 object
#' @param results DESeq2 results (optional)
#' @param args Command-line arguments
#' @param title_prefix Prefix for plot titles
#' @export
create_atac_visualizations <- function(dds, results = NULL, args, title_prefix = "ATAC-seq") {
  log_message("Creating ATAC-seq visualizations...")
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Create MDS plot
  col_data <- as.data.frame(colData(dds))
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  
  if (length(factor_names) > 0) {
    color_by <- factor_names[1]
  } else {
    color_by <- NULL
  }
  
  generate_mds_plot_html(
    norm_counts,
    paste0(args$output, "_mds_plot.html"),
    metadata = col_data,
    color_by = color_by,
    title = paste(title_prefix, "MDS Plot")
  )
  log_message("Created MDS plot")
  
  # Create additional plots if results are provided
  if (!is.null(results)) {
    # MA plot
    ma_plot <- function() {
      DESeq2::plotMA(results, main = paste(title_prefix, "MA Plot"))
    }
    save_plot(ma_plot, args$output, "ma_plot")
    log_message("Created MA plot")
    
    # Volcano plot
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      results_df <- as.data.frame(results)
      results_df$significance <- ifelse(
        results_df$padj < args$fdr & abs(results_df$log2FoldChange) > args$lfcthreshold,
        "Significant", 
        "Not Significant"
      )
      
      volcano_plot <- ggplot2::ggplot(
        results_df, 
        ggplot2::aes(
          x = log2FoldChange, 
          y = -log10(pvalue),
          color = significance
        )
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
        ggplot2::labs(
          title = paste(title_prefix, "Volcano Plot"),
          x = "log2 Fold Change",
          y = "-log10(p-value)"
        ) +
        ggplot2::theme_minimal()
      
      save_plot(volcano_plot, args$output, "volcano_plot")
      log_message("Created volcano plot")
    }
  }
}