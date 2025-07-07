#!/usr/bin/env Rscript

# --- ATAC-seq specific analysis functions ---

#' Export LRT results to file
#'
#' @param lrt_results DESeq2 results object from LRT test
#' @param consensus_peaks GRanges object with consensus peaks
#' @param args Command line arguments
export_lrt_results <- function(lrt_results, consensus_peaks, args) {
  message("Exporting LRT results...")
  
  # Convert results to data frame
  results_df <- as.data.frame(lrt_results)
  results_df$peak_id <- rownames(results_df)
  
  # Add genomic coordinates
  peak_coords <- as.data.frame(consensus_peaks)
  results_df$chr <- peak_coords$seqnames
  results_df$start <- peak_coords$start
  results_df$end <- peak_coords$end
  results_df$width <- peak_coords$width
  
  # Reorder columns for better readability
  col_order <- c("peak_id", "chr", "start", "end", "width", 
                 "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  results_df <- results_df[, col_order[col_order %in% colnames(results_df)]]
  
  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$pvalue), ]
  
  # Add significance column
  results_df$significant <- !is.na(results_df$padj) & results_df$padj < args$fdr
  
  # Export to file
  output_file <- paste0(args$output, "_lrt_results.tsv")
  write.table(results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("LRT results exported to: %s", output_file))
  message(sprintf("Found %d significant peaks (FDR < %g)", 
                  sum(results_df$significant, na.rm = TRUE), args$fdr))
  
  return(results_df)
}

#' Export DESeq2 object
#'
#' @param dds_lrt DESeq2 object
#' @param args Command line arguments
export_deseq_object <- function(dds_lrt, args) {
  message("Exporting DESeq2 object...")
  
  output_file <- paste0(args$output, "_deseq_object.rds")
  saveRDS(dds_lrt, output_file)
  
  message(sprintf("DESeq2 object exported to: %s", output_file))
}

#' Export normalized counts
#'
#' @param dds_lrt DESeq2 object
#' @param args Command line arguments
export_normalized_counts <- function(dds_lrt, args) {
  message("Exporting normalized counts...")
  
  # Get normalized counts
  norm_counts <- counts(dds_lrt, normalized = TRUE)
  
  # Add peak coordinates as row names or separate columns
  output_file <- paste0(args$output, "_normalized_counts.tsv")
  write.table(norm_counts, output_file, sep = "\t", quote = FALSE)
  
  message(sprintf("Normalized counts exported to: %s", output_file))
  
  # Also export in GCT format for GSEA compatibility
  export_gct_format(norm_counts, args)
}

#' Export counts in GCT format
#'
#' @param norm_counts Normalized counts matrix
#' @param args Command line arguments
export_gct_format <- function(norm_counts, args) {
  message("Exporting GCT format...")
  
  # Create GCT header
  gct_header <- c("#1.2", paste(nrow(norm_counts), ncol(norm_counts), sep = "\t"))
  
  # Create data frame with NAME and Description columns
  gct_data <- data.frame(
    NAME = rownames(norm_counts),
    Description = paste("Peak", rownames(norm_counts)),
    norm_counts,
    stringsAsFactors = FALSE
  )
  
  # Export GCT file
  output_file <- paste0(args$output, "_normalized_counts.gct")
  
  # Write header
  writeLines(gct_header, output_file)
  
  # Append data
  write.table(gct_data, output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, append = TRUE)
  
  message(sprintf("GCT format exported to: %s", output_file))
}

#' Generate visualizations
#'
#' @param dds_lrt DESeq2 object
#' @param lrt_results LRT results
#' @param args Command line arguments
generate_visualizations <- function(dds_lrt, lrt_results, args) {
  message("Generating visualizations...")
  
  # Generate PCA plot
  generate_pca_plot(dds_lrt, args)
  
  # Generate MA plot
  generate_ma_plot(lrt_results, args)
  
  # Generate volcano plot
  generate_volcano_plot(lrt_results, args)
  
  # Generate heatmap of top peaks
  generate_heatmap(dds_lrt, lrt_results, args)
}

#' Generate PCA plot
#'
#' @param dds_lrt DESeq2 object
#' @param args Command line arguments
generate_pca_plot <- function(dds_lrt, args) {
  message("Generating PCA plot...")
  
  # Apply variance stabilizing transformation
  vsd <- vst(dds_lrt, blind = FALSE)
  
  # Generate PCA plot
  pca_plot <- plotPCA(vsd, intgroup = names(colData(dds_lrt))[1:2], returnData = FALSE)
  
  # Save plot
  output_file <- paste0(args$output, "_pca_plot.png")
  ggsave(output_file, pca_plot, width = 8, height = 6, dpi = 300)
  
  message(sprintf("PCA plot saved to: %s", output_file))
}

#' Generate MA plot
#'
#' @param lrt_results LRT results
#' @param args Command line arguments
generate_ma_plot <- function(lrt_results, args) {
  message("Generating MA plot...")
  
  # Create MA plot
  ma_plot <- ggplot(as.data.frame(lrt_results), aes(x = baseMean, y = log2FoldChange)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_point(data = subset(as.data.frame(lrt_results), padj < args$fdr), 
               color = "red", alpha = 0.7, size = 0.8) +
    scale_x_log10() +
    labs(title = "MA Plot - ATAC-seq LRT Results",
         x = "Mean of Normalized Counts (log10)",
         y = "Log2 Fold Change") +
    theme_minimal()
  
  # Save plot
  output_file <- paste0(args$output, "_ma_plot.png")
  ggsave(output_file, ma_plot, width = 8, height = 6, dpi = 300)
  
  message(sprintf("MA plot saved to: %s", output_file))
}

#' Generate volcano plot
#'
#' @param lrt_results LRT results
#' @param args Command line arguments
generate_volcano_plot <- function(lrt_results, args) {
  message("Generating volcano plot...")
  
  # Prepare data
  results_df <- as.data.frame(lrt_results)
  results_df$neg_log10_pval <- -log10(results_df$pvalue)
  results_df$significant <- !is.na(results_df$padj) & results_df$padj < args$fdr
  
  # Create volcano plot
  volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = neg_log10_pval)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_point(data = subset(results_df, significant), 
               color = "red", alpha = 0.7, size = 0.8) +
    labs(title = "Volcano Plot - ATAC-seq LRT Results",
         x = "Log2 Fold Change",
         y = "-Log10(p-value)") +
    theme_minimal()
  
  # Save plot
  output_file <- paste0(args$output, "_volcano_plot.png")
  ggsave(output_file, volcano_plot, width = 8, height = 6, dpi = 300)
  
  message(sprintf("Volcano plot saved to: %s", output_file))
}

#' Generate heatmap of top peaks
#'
#' @param dds_lrt DESeq2 object
#' @param lrt_results LRT results
#' @param args Command line arguments
generate_heatmap <- function(dds_lrt, lrt_results, args) {
  message("Generating heatmap...")
  
  # Get top peaks (by p-value)
  top_peaks <- head(order(lrt_results$pvalue), 50)
  
  # Get normalized counts for top peaks
  norm_counts <- counts(dds_lrt, normalized = TRUE)
  top_counts <- norm_counts[top_peaks, ]
  
  # Log transform and scale
  log_counts <- log2(top_counts + 1)
  scaled_counts <- t(scale(t(log_counts)))
  
  # Create heatmap
  heatmap_file <- paste0(args$output, "_heatmap.png")
  
  png(heatmap_file, width = 10, height = 8, units = "in", res = 300)
  heatmap(scaled_counts, 
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Top 50 Variable Peaks",
          cexRow = 0.8, cexCol = 0.8)
  dev.off()
  
  message(sprintf("Heatmap saved to: %s", heatmap_file))
}

#' Perform clustering analysis
#'
#' @param dds_lrt DESeq2 object
#' @param lrt_results LRT results
#' @param args Command line arguments
perform_clustering <- function(dds_lrt, lrt_results, args) {
  message("Performing clustering analysis...")
  
  # Get significant peaks
  sig_peaks <- which(!is.na(lrt_results$padj) & lrt_results$padj < args$fdr)
  
  if (length(sig_peaks) < 10) {
    warning("Too few significant peaks for clustering analysis")
    return(NULL)
  }
  
  # Get normalized counts for significant peaks
  norm_counts <- counts(dds_lrt, normalized = TRUE)
  sig_counts <- norm_counts[sig_peaks, ]
  
  # Log transform and scale
  log_counts <- log2(sig_counts + 1)
  
  if (args$cluster %in% c("row", "both")) {
    perform_row_clustering(log_counts, args)
  }
  
  if (args$cluster %in% c("column", "both")) {
    perform_column_clustering(log_counts, args)
  }
}

#' Perform row clustering (peaks)
#'
#' @param log_counts Log-transformed counts matrix
#' @param args Command line arguments
perform_row_clustering <- function(log_counts, args) {
  message("Performing row clustering...")
  
  # Calculate distance matrix
  dist_matrix <- dist(log_counts, method = args$rowdist)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix)
  
  # Save clustering results
  cluster_file <- paste0(args$output, "_row_clustering.rds")
  saveRDS(hc, cluster_file)
  
  # Create dendrogram plot
  dendro_file <- paste0(args$output, "_row_dendrogram.png")
  png(dendro_file, width = 12, height = 8, units = "in", res = 300)
  plot(hc, main = "Peak Clustering Dendrogram", xlab = "", sub = "")
  dev.off()
  
  message(sprintf("Row clustering results saved to: %s", cluster_file))
  message(sprintf("Row dendrogram saved to: %s", dendro_file))
}

#' Perform column clustering (samples)
#'
#' @param log_counts Log-transformed counts matrix
#' @param args Command line arguments
perform_column_clustering <- function(log_counts, args) {
  message("Performing column clustering...")
  
  # Calculate distance matrix
  dist_matrix <- dist(t(log_counts), method = args$columndist)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix)
  
  # Save clustering results
  cluster_file <- paste0(args$output, "_column_clustering.rds")
  saveRDS(hc, cluster_file)
  
  # Create dendrogram plot
  dendro_file <- paste0(args$output, "_column_dendrogram.png")
  png(dendro_file, width = 10, height = 6, units = "in", res = 300)
  plot(hc, main = "Sample Clustering Dendrogram", xlab = "", sub = "")
  dev.off()
  
  message(sprintf("Column clustering results saved to: %s", cluster_file))
  message(sprintf("Column dendrogram saved to: %s", dendro_file))
}

#' Export contrasts table
#'
#' @param contrasts_list List of contrast results
#' @param args Command line arguments
export_contrasts_table <- function(contrasts_list, args) {
  message("Exporting contrasts table...")
  
  if (is.null(contrasts_list) || length(contrasts_list) == 0) {
    message("No contrasts to export")
    return(NULL)
  }
  
  # Combine all contrasts into a single table
  all_contrasts <- data.frame()
  
  for (contrast_name in names(contrasts_list)) {
    contrast_data <- as.data.frame(contrasts_list[[contrast_name]])
    contrast_data$contrast <- contrast_name
    contrast_data$peak_id <- rownames(contrast_data)
    
    all_contrasts <- rbind(all_contrasts, contrast_data)
  }
  
  # Export combined table
  output_file <- paste0(args$output, "_contrasts_table.tsv")
  write.table(all_contrasts, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("Contrasts table exported to: %s", output_file))
  
  return(all_contrasts)
}

# Apply advanced limma batch correction with interaction design
apply_limma_batch_correction <- function(counts_mat, dba_analyzed, lrt_results, args) {
  message("Applying advanced limma batch correction with interaction design...")
  
  # Get significant peaks
  significant_peaks <- which(!is.na(lrt_results$padj) & lrt_results$padj < args$fdr)
  
  if (length(significant_peaks) == 0) {
    warning("No significant peaks found for batch correction")
    return(NULL)
  }
  
  message(sprintf("Applying batch correction to %d significant peaks", length(significant_peaks)))
  
  # Subset to significant peaks
  sub_counts <- counts_mat[significant_peaks, , drop = FALSE]
  
  # Create column data for design matrix
  coldata <- dba_analyzed$samples
  rownames(coldata) <- coldata$SampleID
  
  # Create batch factor
  batch <- as.factor(coldata$Replicate)
  
  # Create interaction design matrix matching the example
  if (args$with_interaction) {
    design <- model.matrix(~ Condition + Tissue + Condition:Tissue, coldata)
  } else {
    design <- model.matrix(~ Condition + Tissue, coldata)
  }
  
  # Apply limma batch correction
  log2_sub_counts <- log2(sub_counts + 1)
  limma_corrected_data <- removeBatchEffect(log2_sub_counts, batch = batch, design = design)
  
  # Scale the corrected data
  scale_min_max <- function(x, min_range = -2, max_range = 2) {
    min_val <- min(x)
    max_val <- max(x)
    scaled_x <- (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
    return(scaled_x)
  }
  
  scaled_sub_counts_limma <- t(apply(limma_corrected_data, 1, FUN = function(x) scale_min_max(x)))
  
  # Perform hierarchical clustering on corrected data
  d_limma <- dist(as.matrix(scaled_sub_counts_limma))
  hc_limma <- hclust(d_limma)
  
  # Save results
  saveRDS(hc_limma, file.path(args$output, paste0(args$alias, "_hc_clustering_limma_corrected.rds")))
  saveRDS(scaled_sub_counts_limma, file.path(args$output, paste0(args$alias, "_limma_corrected_scaled_counts.rds")))
  
  return(list(
    limma_corrected_data = limma_corrected_data,
    scaled_sub_counts_limma = scaled_sub_counts_limma,
    hc_limma = hc_limma,
    significant_peaks = significant_peaks
  ))
}

# Export significant peaks with genomic coordinates (matching example approach)
export_significant_peaks_genomic <- function(interaction_peakset, lrt_results, args) {
  message("Exporting significant peaks with genomic coordinates...")
  
  # Get significant peaks
  significant_peaks <- which(!is.na(lrt_results$padj) & lrt_results$padj < args$fdr)
  
  if (length(significant_peaks) > 0) {
    # Export the genomic ranges for significant peaks
    # This preserves proper chromosome ordering (avoiding Chr10 -> Chr2 sorting issues)
    significant_genomic_peaks <- interaction_peakset[significant_peaks]
    
    saveRDS(significant_genomic_peaks, 
            file.path(args$output, paste0(args$alias, "_significant_peaks_genomic.rds")))
    
    # Also export as a data frame for easier inspection
    peaks_df <- data.frame(
      chr = seqnames(significant_genomic_peaks),
      start = start(significant_genomic_peaks),
      end = end(significant_genomic_peaks),
      log2FoldChange = lrt_results$log2FoldChange[significant_peaks],
      pvalue = lrt_results$pvalue[significant_peaks],
      padj = lrt_results$padj[significant_peaks]
    )
    
    write.csv(peaks_df, 
              file.path(args$output, paste0(args$alias, "_significant_peaks.csv")),
              row.names = FALSE)
    
    message(sprintf("Exported %d significant peaks", length(significant_peaks)))
  } else {
    message("No significant peaks found")
  }
} 