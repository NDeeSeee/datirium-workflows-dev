#
# EdgeR Analysis Functions for DESeq Workflow
# Used for single replicate differential expression analysis
#

#' Run EdgeR-based differential expression analysis for single replicates
#'
#' @param count_data Count matrix with genes as rows and samples as columns
#' @param col_data Sample metadata data frame
#' @param condition_names Named vector with condition names
#' @param args Command line arguments list
#' @return List containing results and normalized counts
#' @export
run_edger_analysis <- function(count_data, col_data, condition_names, args) {
  log_message("Starting EdgeR analysis for single replicate comparison")
  
  # Load EdgeR if not already loaded
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("EdgeR package is required for single replicate analysis")
  }
  
  tryCatch({
    # Create DGEList object
    log_message("Creating DGEList object")
    dge <- edgeR::DGEList(counts = count_data, group = col_data$condition)
    
    # Calculate normalization factors
    log_message("Calculating normalization factors")
    dge <- edgeR::calcNormFactors(dge)
    
    # Since we have single replicates, use a fixed dispersion
    # This is a common approach for single replicate comparisons
    log_message("Setting fixed dispersion for single replicate analysis")
    bcv <- 0.4  # Biological coefficient of variation (40% is conservative)
    
    # Perform exact test for differential expression
    log_message("Performing EdgeR exact test")
    et <- edgeR::exactTest(dge, dispersion = bcv^2)
    
    # Get results table
    log_message("Extracting results table")
    results_table <- edgeR::topTags(et, n = Inf, sort.by = "PValue")$table
    
    # Add gene names if available
    if (!is.null(rownames(count_data))) {
      results_table$gene_id <- rownames(results_table)
    }
    
    # Calculate normalized counts (CPM - counts per million)
    log_message("Calculating normalized counts")
    norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
    
    # Rename columns to match expected format
    colnames(results_table)[colnames(results_table) == "logFC"] <- "log2FoldChange"
    colnames(results_table)[colnames(results_table) == "PValue"] <- "pvalue"
    colnames(results_table)[colnames(results_table) == "FDR"] <- "padj"
    
    # Add baseMean column (average of normalized counts)
    results_table$baseMean <- rowMeans(norm_counts)
    
    # Filter by FDR and log fold change if specified
    if (!is.null(args$fdr)) {
      log_message(paste("Filtering results by FDR <", args$fdr))
      significant_genes <- sum(results_table$padj < args$fdr, na.rm = TRUE)
      log_message(paste("Found", significant_genes, "genes with FDR <", args$fdr))
    }
    
    if (!is.null(args$lfcthreshold) && args$use_lfc_thresh) {
      log_message(paste("Filtering by log2 fold change threshold:", args$lfcthreshold))
      lfc_genes <- sum(abs(results_table$log2FoldChange) > args$lfcthreshold, na.rm = TRUE)
      log_message(paste("Found", lfc_genes, "genes with |log2FC| >", args$lfcthreshold))
    }
    
    log_message("EdgeR analysis completed successfully")
    
    # Return results in format compatible with DESeq2 workflow
    return(list(
      res = results_table,
      norm_counts = norm_counts,
      dge = dge,
      method = "EdgeR"
    ))
    
  }, error = function(e) {
    log_message(paste("Error in EdgeR analysis:", e$message))
    stop(paste("EdgeR analysis failed:", e$message))
  })
}