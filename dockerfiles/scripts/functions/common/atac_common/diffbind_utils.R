#!/usr/bin/env Rscript
#
# DiffBind utility functions for ATAC-seq workflows
# This module provides common DiffBind-specific functionality
#

#' Run DiffBind differential analysis with DESeq2 backend
#'
#' @param dba_obj DBA object with count data
#' @param design_formula Design formula for analysis
#' @param test_type Type of statistical test ("Wald" or "LRT")
#' @param reduced_formula Reduced formula for LRT (optional)
#' @param output_prefix Output prefix for saving results
#' @return List containing analysis results
#' @export
run_diffbind_deseq2_analysis <- function(dba_obj, design_formula, test_type = "Wald", 
                                        reduced_formula = NULL, output_prefix = NULL) {
  log_message("Running DiffBind differential analysis with DESeq2...")
  
  # Run DiffBind analysis
  dba_analyzed <- dba.analyze(
    dba_obj, 
    method = DBA_DESEQ2, 
    design = design_formula
  )
  
  # Save the DBA object if output prefix provided
  if (!is.null(output_prefix)) {
    saveRDS(dba_analyzed, paste0(output_prefix, "_dba_analyzed.rds"))
    log_message(paste("Saved DBA analyzed object to", paste0(output_prefix, "_dba_analyzed.rds")))
  }
  
  # Extract the DESeq2 object for custom analysis
  log_message("Extracting DESeq2 object for custom analysis...")
  dds <- dba_analyzed$DESeq2$DEdata
  
  # Run DESeq2 analysis based on test type
  if (test_type == "LRT") {
    if (is.null(reduced_formula)) {
      stop("Reduced formula required for LRT test")
    }
    log_message("Running DESeq2 LRT test...")
    dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
  } else {
    log_message("Running DESeq2 Wald test...")
    dds <- DESeq(dds, test = "Wald")
  }
  
  # Extract peakset with proper genomic coordinates
  log_message("Extracting genomic coordinates...")
  peakset <- dba.peakset(dba_analyzed, bRetrieve = TRUE)
  
  # Extract count matrix
  counts_mat <- as.matrix(mcols(peakset))
  rownames(counts_mat) <- paste0(seqnames(peakset), ":", 
                                 start(peakset), "-", 
                                 end(peakset))
  colnames(counts_mat) <- dba_analyzed$samples$SampleID
  
  return(list(
    dba_analyzed = dba_analyzed,
    dds = dds,
    peakset = peakset,
    counts_mat = counts_mat
  ))
}

#' Get contrast results from DESeq2 object
#'
#' @param dds DESeq2 object
#' @param contrast_vector Contrast vector (e.g., c("condition", "treatment", "control"))
#' @param contrast_name Name for the contrast (optional)
#' @param alpha Significance threshold
#' @param lfcThreshold Log fold change threshold
#' @param altHypothesis Alternative hypothesis for testing
#' @return DESeq2 results object
#' @export
get_deseq2_contrast_results <- function(dds, contrast_vector = NULL, contrast_name = NULL, 
                                       alpha = 0.1, lfcThreshold = 0, altHypothesis = "greaterAbs") {
  log_message("Extracting contrast results from DESeq2 object...")
  
  if (!is.null(contrast_vector)) {
    log_message(paste("Getting results for contrast:", paste(contrast_vector, collapse = " vs ")))
    results_obj <- results(
      dds,
      contrast = contrast_vector,
      alpha = alpha,
      lfcThreshold = lfcThreshold,
      independentFiltering = TRUE,
      altHypothesis = altHypothesis
    )
  } else if (!is.null(contrast_name)) {
    log_message(paste("Getting results for contrast name:", contrast_name))
    results_obj <- results(
      dds,
      name = contrast_name,
      alpha = alpha,
      lfcThreshold = lfcThreshold,
      independentFiltering = TRUE,
      altHypothesis = altHypothesis
    )
  } else {
    log_message("Getting results for last coefficient")
    results_obj <- results(
      dds,
      alpha = alpha,
      lfcThreshold = lfcThreshold,
      independentFiltering = TRUE,
      altHypothesis = altHypothesis
    )
  }
  
  log_message("DESeq2 contrast results obtained successfully")
  return(results_obj)
}

#' Export significant peaks from ATAC-seq analysis
#'
#' @param peakset GenomicRanges object with peak coordinates
#' @param results DESeq2 results object
#' @param output_prefix Output file prefix
#' @param fdr_threshold FDR significance threshold
#' @param lfc_threshold Log fold change threshold
#' @return Number of significant peaks exported
#' @export
export_significant_atac_peaks <- function(peakset, results, output_prefix, 
                                         fdr_threshold = 0.1, lfc_threshold = 0) {
  log_message("Exporting significant peaks...")
  
  # Filter for significant peaks
  significant_mask <- !is.na(results$padj) & results$padj < fdr_threshold & 
                     abs(results$log2FoldChange) > lfc_threshold
  
  n_significant <- sum(significant_mask)
  
  if (n_significant == 0) {
    log_warning("No significant peaks found for the given thresholds")
    return(0)
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
  
  # Export to TSV file
  tsv_file <- paste0(output_prefix, "_significant_peaks.tsv")
  write.table(peaks_df, file = tsv_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  log_message(paste("Exported", nrow(peaks_df), "significant peaks to", tsv_file))
  
  # Also export as BED file for visualization
  bed_file <- paste0(output_prefix, "_significant_peaks.bed")
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
  
  return(n_significant)
}

#' Apply score type transformations to DBA object
#'
#' @param dba_obj DBA object
#' @param score_type Score type to apply (e.g., "DBA_SCORE_RPKM")
#' @return Updated DBA object
#' @export
apply_dba_score_type <- function(dba_obj, score_type) {
  if (is.null(score_type) || score_type == "none") {
    log_message("No score type transformation applied")
    return(dba_obj)
  }
  
  log_message(paste("Applying", score_type, "scoring for count normalization..."))
  
  # Apply the specified score type
  if (score_type == "DBA_SCORE_RPKM") {
    dba_obj <- dba.count(dba_obj, peaks = NULL, score = DBA_SCORE_RPKM)
  } else if (score_type == "DBA_SCORE_RPKM_FOLD") {
    dba_obj <- dba.count(dba_obj, peaks = NULL, score = DBA_SCORE_RPKM_FOLD)
  } else if (score_type == "DBA_SCORE_TMM_MINUS_FULL") {
    dba_obj <- dba.count(dba_obj, peaks = NULL, score = DBA_SCORE_TMM_MINUS_FULL)
  } else {
    log_warning(paste("Unknown score type:", score_type, "- using default"))
  }
  
  return(dba_obj)
}

#' Generate summary statistics for ATAC-seq analysis
#'
#' @param results DESeq2 results object
#' @param fdr_threshold FDR significance threshold
#' @param lfc_threshold Log fold change threshold
#' @param comparison_name Name of the comparison (optional)
#' @return List of summary statistics
#' @export
generate_atac_summary_stats <- function(results, fdr_threshold = 0.1, lfc_threshold = 0, 
                                       comparison_name = NULL) {
  log_message("Generating ATAC-seq summary statistics...")
  
  # Calculate basic statistics
  total_peaks <- nrow(results)
  
  # Significant peaks
  significant_mask <- !is.na(results$padj) & results$padj < fdr_threshold
  significant_peaks <- sum(significant_mask)
  
  # Directional changes
  upregulated <- sum(significant_mask & results$log2FoldChange > lfc_threshold)
  downregulated <- sum(significant_mask & results$log2FoldChange < -lfc_threshold)
  
  # P-value distribution
  valid_pvals <- !is.na(results$pvalue)
  mean_pvalue <- mean(results$pvalue[valid_pvals])
  
  # Create summary
  summary_stats <- list(
    total_peaks = total_peaks,
    significant_peaks = significant_peaks,
    upregulated_peaks = upregulated,
    downregulated_peaks = downregulated,
    percentage_significant = round(100 * significant_peaks / total_peaks, 2),
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold,
    mean_pvalue = round(mean_pvalue, 4),
    comparison = comparison_name
  )
  
  # Print summary
  cat("=== ATAC-seq Analysis Summary ===\n")
  if (!is.null(comparison_name)) {
    cat("Comparison:", comparison_name, "\n")
  }
  cat("Total peaks analyzed:", total_peaks, "\n")
  cat("Significant peaks (FDR <", fdr_threshold, "):", significant_peaks, "\n")
  cat("Upregulated peaks:", upregulated, "\n")
  cat("Downregulated peaks:", downregulated, "\n")
  cat("Percentage significant:", summary_stats$percentage_significant, "%\n")
  cat("Mean p-value:", summary_stats$mean_pvalue, "\n")
  cat("================================\n")
  
  return(summary_stats)
}

#' Create DBA plot for sample correlation
#'
#' @param dba_obj DBA object
#' @param output_prefix Output file prefix
#' @param plot_type Type of plot ("heatmap" or "PCA")
#' @export
create_dba_plot <- function(dba_obj, output_prefix, plot_type = "heatmap") {
  log_message(paste("Creating DBA", plot_type, "plot..."))
  
  tryCatch({
    if (plot_type == "heatmap") {
      # Create correlation heatmap
      plot_file <- paste0(output_prefix, "_correlation_heatmap.pdf")
      pdf(plot_file, width = 8, height = 6)
      dba.plotHeatmap(dba_obj)
      dev.off()
      log_message(paste("Created correlation heatmap:", plot_file))
      
    } else if (plot_type == "PCA") {
      # Create PCA plot
      plot_file <- paste0(output_prefix, "_pca_plot.pdf")
      pdf(plot_file, width = 8, height = 6)
      dba.plotPCA(dba_obj)
      dev.off()
      log_message(paste("Created PCA plot:", plot_file))
    }
  }, error = function(e) {
    log_warning(paste("Failed to create DBA plot:", e$message))
  })
}

#' Extract peak annotation information if available
#'
#' @param peakset GenomicRanges object
#' @param annotation_file Optional annotation file path
#' @return Data frame with peak annotations
#' @export
annotate_atac_peaks <- function(peakset, annotation_file = NULL) {
  log_message("Annotating ATAC-seq peaks...")
  
  # Create basic annotation data frame
  annotation_df <- data.frame(
    peak_id = paste0(seqnames(peakset), ":", start(peakset), "-", end(peakset)),
    chromosome = seqnames(peakset),
    start = start(peakset),
    end = end(peakset),
    width = width(peakset),
    stringsAsFactors = FALSE
  )
  
  # Add external annotation if provided
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    log_message(paste("Loading external annotation from", annotation_file))
    # Implementation would depend on annotation file format
    # This is a placeholder for future enhancement
  }
  
  # Add basic genomic features if possible
  # This could be enhanced with promoter/enhancer annotation
  
  return(annotation_df)
}