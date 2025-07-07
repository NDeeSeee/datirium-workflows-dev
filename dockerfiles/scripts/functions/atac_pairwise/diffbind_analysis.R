#!/usr/bin/env Rscript
#
# DiffBind analysis functions for ATAC-seq Pairwise Analysis
#

#' Run DiffBind analysis pipeline for pairwise comparison
#'
#' @param sample_sheet DiffBind sample sheet
#' @param args Command-line arguments
#' @return List containing DBA object and analysis results
#' @export
run_pairwise_diffbind_analysis <- function(sample_sheet, args) {
  log_message("Running DiffBind analysis pipeline for pairwise comparison...")
  
  # Set default values for DiffBind parameters if not provided
  peakformat <- if (is.null(args$peakformat)) "csv" else args$peakformat
  peakcaller <- if (is.null(args$peakcaller)) "macs" else args$peakcaller
  scorecol <- if (is.null(args$scorecol)) 6 else args$scorecol
  minoverlap <- if (is.null(args$minoverlap)) 2 else args$minoverlap
  
  # Auto-detect MACS2 format and adjust parameters
  peak_files <- c(args$untreated_files, args$treated_files)
  if (any(grepl("\\.xls$", peak_files))) {
    log_message("Detected MACS2 .xls peak files - adjusting format parameters...")
    log_message("  - Auto-corrected peakFormat to: macs (for .xls files)")
    log_message("  - Auto-corrected scoreCol to: 6 (pileup column)")
    peakformat <- "macs"  # Use MACS format for .xls files
    scorecol <- 6  # Column 6 is pileup in MACS2 output
  }
  
  log_message("DiffBind parameters:")
  log_message(paste("  - peakFormat:", peakformat))
  log_message(paste("  - peakCaller:", peakcaller))
  log_message(paste("  - scoreCol:", scorecol))
  
  # Test mode bypass for DiffBind analysis - check FIRST before any DiffBind operations
  if (args$test_mode) {
    log_message("Test mode: bypassing DiffBind analysis and generating mock results")
    return(generate_mock_diffbind_results(sample_sheet, args))
  }
  
  # Create DBA object
  log_message("Creating DBA object...")
  dba_obj <- dba(
    sampleSheet = sample_sheet,
    peakFormat = peakformat,
    peakCaller = peakcaller,
    scoreCol = scorecol
  )
  
  # Create consensus peaks
  log_message("Creating consensus peaks...")
  dba_consensus <- dba.peakset(dba_obj, consensus = DBA_CONDITION, minOverlap = minoverlap)
  dba_consensus <- dba(dba_consensus, mask = dba_consensus$masks$Consensus, minOverlap = 1)
  
  # Get consensus peaks
  consensus_peaks <- dba.peakset(dba_consensus, bRetrieve = TRUE, minOverlap = 1)
  
  # Count reads in consensus peaks
  log_message("Counting reads in consensus peaks...")
  dba_obj <- dba.count(dba_obj, peaks = consensus_peaks, minOverlap = 1)
  
  # Apply test mode filtering if requested
  if (args$test_mode) {
    log_message("Test mode: reducing to first 500 peaks for faster processing...")
    # Get binding matrix
    binding_matrix <- dba_obj$binding
    n_peaks <- min(500, nrow(binding_matrix))
    
    # Create subset of peaks
    test_peaks <- consensus_peaks[1:n_peaks]
    
    # Recreate DBA object with subset
    dba_obj <- dba.count(dba_obj, peaks = test_peaks, minOverlap = 1)
  }
  
  return(list(dba_obj = dba_obj, consensus_peaks = consensus_peaks))
}

#' Run DESeq2 pairwise analysis on DiffBind results
#'
#' @param dba_obj DBA object
#' @param args Command-line arguments
#' @return List containing analysis results
#' @export
run_pairwise_deseq2_analysis <- function(dba_obj, args) {
  log_message("Running DiffBind differential analysis for pairwise comparison...")
  
  # Simple pairwise design
  design_formula <- paste0("~", args$condition_column)
  
  # Run DiffBind analysis
  dba_analyzed <- dba.analyze(
    dba_obj, 
    method = DBA_DESEQ2, 
    design = design_formula
  )
  
  # Save the DBA object
  saveRDS(dba_analyzed, file.path(args$output, paste0(args$output, "_dba_pairwise_analyzed.rds")))
  
  # Extract the DESeq2 object for custom analysis
  log_message("Extracting DESeq2 object for pairwise analysis...")
  dds <- dba_analyzed$DESeq2$DEdata
  
  # Run standard DESeq2 analysis (Wald test for pairwise comparison)
  log_message("Running DESeq2 Wald test for pairwise comparison...")
  dds <- DESeq(dds, test = "Wald")
  
  # Get pairwise results
  contrast_name <- paste(args$condition1, "vs", args$condition2, sep = "_")
  pairwise_results <- results(
    dds, 
    contrast = c(args$condition_column, args$condition1, args$condition2),
    alpha = args$fdr
  )
  
  # Extract count matrix and peak information
  log_message("Extracting count matrix and genomic coordinates...")
  
  # Get the peakset with proper genomic coordinates
  pairwise_peakset <- dba.peakset(dba_analyzed, bRetrieve = TRUE)
  
  # Extract count matrix
  counts_mat <- as.matrix(mcols(pairwise_peakset))
  rownames(counts_mat) <- paste0(seqnames(pairwise_peakset), ":", 
                                 start(pairwise_peakset), "-", 
                                 end(pairwise_peakset))
  colnames(counts_mat) <- dba_analyzed$samples$SampleID
  
  return(list(
    dba_analyzed = dba_analyzed,
    dds = dds,
    pairwise_results = pairwise_results,
    pairwise_peakset = pairwise_peakset,
    counts_mat = counts_mat,
    contrast_name = contrast_name
  ))
}

#' Export significant peaks for pairwise comparison
#'
#' @param peakset Genomic ranges object with peaks
#' @param results DESeq2 results object
#' @param args Command-line arguments
#' @export
export_pairwise_significant_peaks <- function(peakset, results, args) {
  log_message("Exporting significant peaks for pairwise comparison...")
  
  # Filter for significant peaks
  significant_mask <- !is.na(results$padj) & results$padj < args$fdr & 
                     abs(results$log2FoldChange) > args$lfcthreshold
  
  if (sum(significant_mask) == 0) {
    log_warning("No significant peaks found for the given thresholds")
    return()
  }
  
  # Get significant peaks
  sig_peaks <- peakset[significant_mask]
  sig_results <- results[significant_mask, ]
  
  # Create results data frame with genomic coordinates
  peaks_df <- data.frame(
    chromosome = seqnames(sig_peaks),
    start = start(sig_peaks),
    end = end(sig_peaks),
    peak_id = paste0(seqnames(sig_peaks), ":", start(sig_peaks), "-", end(sig_peaks)),
    baseMean = sig_results$baseMean,
    log2FoldChange = sig_results$log2FoldChange,
    lfcSE = sig_results$lfcSE,
    stat = sig_results$stat,
    pvalue = sig_results$pvalue,
    padj = sig_results$padj,
    stringsAsFactors = FALSE
  )
  
  # Sort by adjusted p-value
  peaks_df <- peaks_df[order(peaks_df$padj), ]
  
  # Export to file
  output_file <- paste0(args$output, "_significant_peaks.tsv")
  write.table(peaks_df, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  log_message(paste("Exported", nrow(peaks_df), "significant peaks to", output_file))
  
  # Also export as BED file for visualization
  bed_file <- paste0(args$output, "_significant_peaks.bed")
  bed_df <- data.frame(
    chrom = peaks_df$chromosome,
    chromStart = peaks_df$start - 1,  # BED is 0-based
    chromEnd = peaks_df$end,
    name = peaks_df$peak_id,
    score = round(-log10(peaks_df$padj), 2),
    strand = ".",
    stringsAsFactors = FALSE
  )
  
  write.table(bed_df, file = bed_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  log_message(paste("Exported significant peaks in BED format to", bed_file))
}

#' Create visualization plots for pairwise comparison
#'
#' @param dds DESeq2 object
#' @param results DESeq2 results
#' @param args Command-line arguments
#' @export
create_pairwise_visualizations <- function(dds, results, args) {
  log_message("Creating visualizations for pairwise comparison...")
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Create MA plot
  ma_plot <- function() {
    DESeq2::plotMA(results, main = paste("ATAC-seq MA Plot:", args$condition1, "vs", args$condition2))
  }
  ma_files <- save_plot(ma_plot, args$output, "pairwise_ma_plot")
  log_message("Created pairwise MA plot")
  
  # Create volcano plot
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
        title = paste("ATAC-seq Volcano Plot:", args$condition1, "vs", args$condition2),
        x = "log2 Fold Change",
        y = "-log10(p-value)"
      ) +
      ggplot2::theme_minimal()
    
    volcano_files <- save_plot(volcano_plot, args$output, "pairwise_volcano_plot")
    log_message("Created pairwise volcano plot")
  }
  
  # Create MDS plot
  generate_atac_mds_plot(
    norm_counts,
    "atac_pairwise_mds_plot.html",
    metadata = as.data.frame(colData(dds)),
    color_by = args$condition_column,
    title = paste("ATAC-seq Pairwise MDS Plot:", args$condition1, "vs", args$condition2)
  )
  log_message("Created pairwise MDS plot")
}

#' Generate summary statistics for pairwise comparison
#'
#' @param results DESeq2 results
#' @param args Command-line arguments
#' @return Summary statistics list
#' @export
generate_pairwise_summary <- function(results, args) {
  log_message("Generating summary statistics for pairwise comparison...")
  
  # Calculate basic statistics
  total_peaks <- nrow(results)
  significant_peaks <- sum(!is.na(results$padj) & results$padj < args$fdr)
  upregulated <- sum(!is.na(results$padj) & results$padj < args$fdr & results$log2FoldChange > args$lfcthreshold)
  downregulated <- sum(!is.na(results$padj) & results$padj < args$fdr & results$log2FoldChange < -args$lfcthreshold)
  
  # Create summary
  summary_stats <- list(
    total_peaks = total_peaks,
    significant_peaks = significant_peaks,
    upregulated_peaks = upregulated,
    downregulated_peaks = downregulated,
    fdr_threshold = args$fdr,
    lfc_threshold = args$lfcthreshold,
    condition1 = ifelse(is.null(args$condition1), args$untreated_name, args$condition1),
    condition2 = ifelse(is.null(args$condition2), args$treated_name, args$condition2),
    comparison = paste(ifelse(is.null(args$condition1), args$untreated_name, args$condition1), "vs", ifelse(is.null(args$condition2), args$treated_name, args$condition2))
  )
  
  # Save summary
  saveRDS(summary_stats, file.path(args$output, paste0(args$output, "_pairwise_summary.rds")))
  
  # Create markdown summary file (required by CWL)
  md_content <- c(
    "# ATAC-seq Pairwise Analysis Summary",
    "",
    paste("**Comparison:**", ifelse(is.null(args$condition1), args$untreated_name, args$condition1), "vs", ifelse(is.null(args$condition2), args$treated_name, args$condition2)),
    paste("**Total peaks analyzed:**", total_peaks),
    paste("**Significant peaks (FDR <", args$fdr, "):**", significant_peaks),
    paste("**Upregulated in", ifelse(is.null(args$condition1), args$untreated_name, args$condition1), ":**", upregulated),
    paste("**Downregulated in", ifelse(is.null(args$condition1), args$untreated_name, args$condition1), ":**", downregulated),
    paste("**Percentage significant:**", round(100 * significant_peaks / total_peaks, 2), "%"),
    "",
    "## Analysis Parameters",
    paste("- FDR threshold:", args$fdr),
    paste("- Log2 fold change threshold:", args$lfcthreshold),
    "",
    paste("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  
  # Write markdown summary (required by CWL glob pattern)
  summary_file <- paste0(args$output, "_summary.md")
  writeLines(md_content, summary_file)
  
  # Print summary
  cat("=== ATAC-seq Pairwise Analysis Summary ===\n")
  cat("Comparison:", args$condition1, "vs", args$condition2, "\n")
  cat("Total peaks analyzed:", total_peaks, "\n")
  cat("Significant peaks (FDR <", args$fdr, "):", significant_peaks, "\n")
  cat("Upregulated in", args$condition1, ":", upregulated, "\n")
  cat("Downregulated in", args$condition1, ":", downregulated, "\n")
  cat("Percentage significant:", round(100 * significant_peaks / total_peaks, 2), "%\n")
  cat("=========================================\n")
  
  return(summary_stats)
}

#' Generate mock DiffBind results for test mode
#' @param sample_sheet Sample sheet data frame
#' @param args Command line arguments
#' @return Mock DiffBind results list
generate_mock_diffbind_results <- function(sample_sheet, args) {
  log_message("Generating mock DiffBind results for test mode")
  
  # Create mock peak data based on sample sheet
  n_peaks <- 1000  # Mock number of peaks
  
  # Generate mock differential binding results
  mock_results <- data.frame(
    seqnames = paste0("chr", sample(1:22, n_peaks, replace = TRUE)),
    start = sample(1000000:200000000, n_peaks),
    end = sample(1000000:200000000, n_peaks),
    width = sample(200:2000, n_peaks),
    strand = "*",
    Conc = runif(n_peaks, 0, 10),
    Conc_condition1 = runif(n_peaks, 0, 10),
    Conc_condition2 = runif(n_peaks, 0, 10),
    Fold = rnorm(n_peaks, 0, 2),
    p.value = runif(n_peaks, 0, 1),
    FDR = runif(n_peaks, 0, 1)
  )
  
  # Make some peaks "significant"
  sig_indices <- sample(1:n_peaks, min(100, n_peaks * 0.1))
  mock_results$FDR[sig_indices] <- runif(length(sig_indices), 0, 0.05)
  
  # Create mock binding matrix
  n_samples <- nrow(sample_sheet)
  mock_binding_matrix <- matrix(
    rpois(n_peaks * n_samples, lambda = 50),
    nrow = n_peaks,
    ncol = n_samples
  )
  colnames(mock_binding_matrix) <- sample_sheet$SampleID
  
  # Return results in expected format
  list(
    results = mock_results,
    binding_matrix = mock_binding_matrix,
    sample_sheet = sample_sheet,
    method = "Mock (Test Mode)"
  )
}