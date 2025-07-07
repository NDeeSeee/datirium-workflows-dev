#!/usr/bin/env Rscript
#
# Workflow functions for ATAC-seq Pairwise Analysis
#

#' Main workflow function with memory management
#' 
#' Wraps the run_workflow function with memory tracking
#' 
#' @param args Command line arguments (optional, will be parsed if not provided)
#' @return Results from the analysis workflow
#' @export
main_with_memory_management <- function(args = NULL) {
  # Report initial memory usage
  report_memory_usage("Initial")
  
  # Get command-line arguments if not provided
  if (is.null(args)) {
    args <- get_args()
  }
  
  # Run the main workflow
  results <- tryCatch({
    run_workflow(args)
  }, error = function(e) {
    log_message(paste("Error in ATAC-seq pairwise workflow execution:", e$message), "ERROR")
    # Force garbage collection
    gc(verbose = FALSE)
    # Re-raise the error
    stop(e)
  })
  
  # Final memory usage report
  report_memory_usage("Final")
  
  return(results)
}

#' Initialize the environment for ATAC-seq pairwise analysis
#'
#' Loads required libraries, sources dependency files, and configures environment
initialize_environment <- function() {
  # Display startup message
  message("Starting ATAC-seq Pairwise Analysis")
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

  # Source ATAC-seq pairwise specific functions
  source_with_fallback("functions/atac_pairwise/cli_args.R", "/usr/local/bin/functions/atac_pairwise/cli_args.R")
  source_with_fallback("functions/atac_pairwise/data_processing.R", "/usr/local/bin/functions/atac_pairwise/data_processing.R")
  source_with_fallback("functions/atac_pairwise/diffbind_analysis.R", "/usr/local/bin/functions/atac_pairwise/diffbind_analysis.R")
  
  # Load required libraries with clear error messages
  tryCatch({
    message("Loading required libraries...")
    suppressPackageStartupMessages({
      required_packages <- c(
        # "argparse", # Optional - using manual parsing fallback
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

  # Configure R options
  configure_r_options()
  
  # Configure plot theme
  configure_plot_theme()
  
  log_message("Environment initialized for ATAC-seq pairwise analysis")
}

#' Main workflow function for ATAC-seq pairwise analysis
#'
#' @param args Command-line arguments
#' @return Final results list
#' @export
run_workflow <- function(args) {
  log_message("Starting ATAC-seq pairwise workflow")
  
  # Parse and validate arguments
  validate_args(args)
  print_atac_args(args, "ATAC-seq Pairwise")
  
  # Load and validate metadata for pairwise comparison
  metadata_df <- load_and_validate_pairwise_metadata(args)
  
  # Load and validate ATAC-seq data
  sample_sheet <- load_and_validate_pairwise_atac_data(args, metadata_df)
  
  # Run DiffBind analysis
  diffbind_results <- run_pairwise_diffbind_analysis(sample_sheet, args)
  
  # Handle test mode differently
  if (args$test_mode && !is.null(diffbind_results$method) && diffbind_results$method == "Mock (Test Mode)") {
    log_message("Test mode: generating simplified mock results and summary...")
    
    # Generate a simple mock summary directly
    generate_test_mode_summary(args)
    
    # Return early for test mode
    log_message("Test mode: ATAC-seq pairwise analysis completed with mock results", "SUCCESS")
    return(list(test_mode = TRUE, method = "Mock"))
  }
  
  dba_obj <- diffbind_results$dba_obj
  consensus_peaks <- diffbind_results$consensus_peaks
  
  # Run DESeq2 pairwise analysis
  deseq2_results <- run_pairwise_deseq2_analysis(dba_obj, args)
  dds <- deseq2_results$dds
  pairwise_results <- deseq2_results$pairwise_results
  pairwise_peakset <- deseq2_results$pairwise_peakset
  counts_mat <- deseq2_results$counts_mat
  contrast_name <- deseq2_results$contrast_name
  
  # Get normalized counts from DESeq2
  normCounts <- counts(dds, normalized = TRUE)
  log_message(paste("Extracted normalized counts with dimensions:", nrow(normCounts), "x", ncol(normCounts)))
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    log_message(paste("Applying", args$batchcorrection, "batch correction to normalized ATAC-seq counts"))
    normCounts <- apply_pairwise_batch_correction(
      count_data = normCounts,
      metadata_df = metadata_df,
      batch_method = args$batchcorrection
    )
    log_message("Batch correction completed for ATAC-seq pairwise data")
  } else {
    log_message("No batch correction applied to normalized ATAC-seq counts")
  }
  
  # Export significant peaks
  export_pairwise_significant_peaks(pairwise_peakset, pairwise_results, args)
  
  # Create expression data frame for results
  expDataDf <- data.frame(
    RefseqId = rownames(pairwise_results),
    GeneId = rownames(pairwise_results),
    chromosome = seqnames(pairwise_peakset),
    start = start(pairwise_peakset),
    end = end(pairwise_peakset),
    peak_width = width(pairwise_peakset),
    baseMean = pairwise_results$baseMean,
    stringsAsFactors = FALSE
  )
  rownames(expDataDf) <- rownames(pairwise_results)
  
  # Add pairwise comparison results
  expDataDf[[paste0(contrast_name, "_LFC")]] <- pairwise_results$log2FoldChange
  expDataDf[[paste0(contrast_name, "_pvalue")]] <- pairwise_results$pvalue
  expDataDf[[paste0(contrast_name, "_FDR")]] <- pairwise_results$padj
  
  # Prepare final results
  final_results <- list(
    dds = dds,
    normCounts = normCounts,
    expDataDf = expDataDf,
    pairwise_results = pairwise_results,
    contrast_name = contrast_name,
    metadata = metadata_df
  )
  
  # Save results
  output_file <- paste0(args$output, "_pairwise_results.rds")
  log_message(paste("Saving ATAC-seq pairwise results to", output_file))
  saveRDS(final_results, file = output_file)
  
  # Export reports
  log_message("Exporting ATAC-seq pairwise reports and visualizations", "STEP")
  
  # Export pairwise results table
  results_file <- paste0(args$output, "_pairwise_results.tsv")
  write_pairwise_atac_results(expDataDf, results_file, args)
  write_pairwise_atac_results(expDataDf, paste0(args$output, "_report.tsv"), args)
  log_message(paste("Exported ATAC-seq pairwise results table to", results_file), "INFO")
  
  # Export normalized counts in GCT format
  write_gct_file(normCounts, paste0(args$output, "_normalized_counts.gct"))
  write_gct_file(normCounts, paste0(args$output, "_counts_all.gct"))
  log_message("Exported normalized ATAC-seq counts to GCT format", "INFO")
  
  # Generate phenotypes CLS file for GSEA/CWL
  condition_vector <- colData(dds)[[args$condition_column]]
  write_cls_file(condition_vector, paste0(args$output, "_phenotypes.cls"))
  
  # Generate summary markdown
  log_message("Generating analysis summary", "INFO")
  summary_file <- generate_deseq_summary(
    pairwise_results,
    get_output_filename(args$output, "summary", "md"),
    title = "ATAC-seq Pairwise Analysis Summary",
    parameters = list(
      "Condition 1" = args$uname,
      "Condition 2" = args$tname,
      "Batch correction" = args$batchcorrection,
      "Analysis Method" = "ATAC-seq Pairwise"
    )
  )
  log_message(paste("Generated analysis summary:", summary_file), "INFO")
  
  # Create visualizations
  create_pairwise_visualizations(dds, pairwise_results, args)
  
  # Perform clustering if requested
  if (args$cluster != "none") {
    log_message("Performing clustering analysis...")
    perform_clustering(normCounts, pairwise_results, args)
  }
  
  # Generate summary statistics
  summary_stats <- generate_pairwise_summary(pairwise_results, args)
  
  # Verify output files
  verify_pairwise_outputs(args$output, fail_on_missing = FALSE)
  
  log_message("ATAC-seq pairwise workflow completed successfully", "SUCCESS")
  
  return(final_results)
}

#' Generate test mode summary for ATAC-seq pairwise analysis
#'
#' @param args Command-line arguments
#' @export
generate_test_mode_summary <- function(args) {
  log_message("Generating test mode summary for ATAC-seq pairwise analysis...")
  
  # Create simple mock summary statistics
  mock_total_peaks <- 1000
  mock_significant_peaks <- 100
  mock_upregulated <- 60
  mock_downregulated <- 40
  
  # Get condition names
  condition1 <- if (!is.null(args$treated_name)) args$treated_name else "Condition1"
  condition2 <- if (!is.null(args$untreated_name)) args$untreated_name else "Condition2"
  
  # Create markdown summary content
  md_content <- c(
    "# ATAC-seq Pairwise Analysis Summary (Test Mode)",
    "",
    paste("**Comparison:**", condition1, "vs", condition2),
    paste("**Total peaks analyzed:**", mock_total_peaks),
    paste("**Significant peaks (FDR <", args$fdr, "):**", mock_significant_peaks),
    paste("**Upregulated in", condition1, ":**", mock_upregulated),
    paste("**Downregulated in", condition1, ":**", mock_downregulated),
    paste("**Percentage significant:**", round(100 * mock_significant_peaks / mock_total_peaks, 2), "%"),
    "",
    "## Analysis Parameters",
    paste("- FDR threshold:", args$fdr),
    paste("- Log2 fold change threshold:", args$lfcthreshold),
    paste("- Analysis mode: Test Mode (Mock Results)"),
    "",
    "## Note",
    "This analysis was run in test mode with mock data for validation purposes.",
    "",
    paste("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  
  # Write the summary file that CWL expects
  summary_file <- paste0(args$output, "_summary.md")
  writeLines(md_content, summary_file)
  log_message(paste("Generated test mode summary:", summary_file), "INFO")
  
  # Also create a simple normalized counts file for CWL
  mock_counts <- matrix(
    rpois(1000 * 4, lambda = 50),
    nrow = 1000,
    ncol = 4
  )
  colnames(mock_counts) <- paste0("sample", 1:4)
  rownames(mock_counts) <- paste0("peak_", 1:1000)
  
  # Write GCT files for CWL outputs
  write_gct_file(mock_counts, paste0(args$output, "_counts_all.gct"))
  write_gct_file(mock_counts, paste0(args$output, "_counts_filtered.gct"))
  log_message("Generated test mode GCT files", "INFO")
  
  # Create mock report.tsv file
  mock_report <- data.frame(
    RefseqId = paste0("peak_", 1:100),
    GeneId = paste0("peak_", 1:100),
    Chrom = paste0("chr", sample(1:22, 100, replace = TRUE)),
    TxStart = sample(1000000:200000000, 100),
    TxEnd = sample(1000000:200000000, 100),
    Strand = sample(c("+", "-"), 100, replace = TRUE),
    TotalReads = sample(100:1000, 100),
    BaseMean = runif(100, 10, 1000),
    log2FoldChange = rnorm(100, 0, 2),
    lfcSE = runif(100, 0.1, 1),
    stat = rnorm(100, 0, 2),
    pvalue = runif(100, 0, 1),
    padj = runif(100, 0, 1)
  )
  write.table(mock_report, paste0(args$output, "_report.tsv"), sep = "\t", 
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  log_message("Generated test mode report file", "INFO")
  
  # Create mock phenotypes.cls file
  phenotypes_content <- c(
    paste(4, 2, 1),
    paste("#", condition1, condition2),
    paste(c(rep(condition2, 2), rep(condition1, 2)), collapse = " ")
  )
  writeLines(phenotypes_content, paste0(args$output, "_phenotypes.cls"))
  log_message("Generated test mode phenotypes file", "INFO")
  
  # Create simple mock plots (empty PNG files for now)
  mock_plot_files <- c(
    paste0(args$output, "_ma_plot.png"),
    paste0(args$output, "_expression_heatmap.png"), 
    paste0(args$output, "_pca_plot.png"),
    paste0(args$output, "_ma_plot.pdf"),
    paste0(args$output, "_expression_heatmap.pdf"),
    paste0(args$output, "_pca_plot.pdf")
  )
  
  for (plot_file in mock_plot_files) {
    # Create minimal PNG/PDF headers (for file detection)
    if (grepl("\\.png$", plot_file)) {
      # Minimal PNG header
      writeBin(as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A)), plot_file)
    } else {
      # Minimal PDF header
      writeLines("%PDF-1.4", plot_file)
    }
  }
  log_message("Generated test mode plot files", "INFO")
  
  # Create MDS plot HTML file
  mds_html <- c(
    "<!DOCTYPE html>",
    "<html><head><title>Test Mode MDS Plot</title></head>",
    "<body><h1>ATAC-seq Pairwise MDS Plot (Test Mode)</h1>",
    "<p>This is a mock MDS plot generated in test mode.</p>",
    "</body></html>"
  )
  writeLines(mds_html, paste0(args$output, "_mds_plot.html"))
  log_message("Generated test mode MDS plot", "INFO")
  
  # Print summary to console
  cat("=== ATAC-seq Pairwise Analysis Summary (Test Mode) ===\n")
  cat("Comparison:", condition1, "vs", condition2, "\n")
  cat("Total peaks analyzed:", mock_total_peaks, "\n")
  cat("Significant peaks (FDR <", args$fdr, "):", mock_significant_peaks, "\n")
  cat("Upregulated in", condition1, ":", mock_upregulated, "\n")
  cat("Downregulated in", condition1, ":", mock_downregulated, "\n")
  cat("Percentage significant:", round(100 * mock_significant_peaks / mock_total_peaks, 2), "%\n")
  cat("Analysis mode: Test Mode (Mock Results)\n")
  cat("==================================================\n")
}

#' Perform clustering analysis on pairwise results
#'
#' @param norm_counts Normalized count matrix
#' @param results DESeq2 results
#' @param args Command-line arguments
#' @export
perform_clustering <- function(norm_counts, results, args) {
  log_message("Performing clustering analysis for pairwise comparison...")
  
  # Filter for significant peaks for clustering
  significant_mask <- !is.na(results$padj) & results$padj < args$fdr
  
  if (sum(significant_mask) < 10) {
    log_warning("Too few significant peaks for clustering analysis")
    return()
  }
  
  sig_counts <- norm_counts[significant_mask, ]
  
  # Apply scaling
  if (args$scaling_type == "zscore") {
    scaled_counts <- t(scale(t(sig_counts)))
  } else if (args$scaling_type == "minmax") {
    scaled_counts <- apply(sig_counts, 1, function(x) (x - min(x)) / (max(x) - min(x)))
    scaled_counts <- t(scaled_counts)
  } else {
    scaled_counts <- sig_counts
  }
  
  # Perform clustering using hopach
  if (requireNamespace("hopach", quietly = TRUE)) {
    log_message("Running HOPACH clustering...")
    
    # Run clustering based on specified method
    if (args$cluster %in% c("row", "both")) {
      row_clustering <- hopach::hopach(scaled_counts, 
                                      dmat = hopach::distancematrix(scaled_counts, d = args$rowdist),
                                      K = args$k, kmax = args$kmax)
      
      # Save row clustering results
      clustering_file <- paste0(args$output, "_row_clustering.rds")
      saveRDS(row_clustering, clustering_file)
      log_message(paste("Row clustering results saved to", clustering_file))
    }
    
    if (args$cluster %in% c("column", "both")) {
      col_clustering <- hopach::hopach(t(scaled_counts),
                                      dmat = hopach::distancematrix(t(scaled_counts), d = args$columndist),
                                      K = args$k, kmax = args$kmax)
      
      # Save column clustering results  
      clustering_file <- paste0(args$output, "_column_clustering.rds")
      saveRDS(col_clustering, clustering_file)
      log_message(paste("Column clustering results saved to", clustering_file))
    }
  } else {
    log_warning("hopach package not available - skipping clustering analysis")
  }
}