#!/usr/bin/env Rscript
#
# Data processing functions for ATAC-seq LRT Step 2
#

#' Validate metadata consistency with step 1
#'
#' @param col_data Column data from DESeq2 object
#' @param batchcorrection Batch correction method
#' @param design_formula Design formula
#' @return Validated col_data
#' @export
validate_metadata <- function(col_data, batchcorrection, design_formula) {
  log_message("Validating metadata consistency with step 1")
  
  # Check for required columns based on design formula
  formula_vars <- all.vars(design_formula)
  missing_vars <- formula_vars[!formula_vars %in% colnames(col_data)]
  
  if (length(missing_vars) > 0) {
    stop("Missing variables in metadata: ", paste(missing_vars, collapse = ", "))
  }
  
  # Check batch correction requirements
  if (batchcorrection != "none" && !"batch" %in% colnames(col_data)) {
    warning("Batch correction requested but no 'batch' column found in metadata")
  }
  
  # Ensure factor columns are properly formatted
  for (var in formula_vars) {
    if (!is.factor(col_data[[var]])) {
      log_message(paste("Converting", var, "to factor"))
      col_data[[var]] <- factor(col_data[[var]])
    }
  }
  
  log_message("Metadata validation completed")
  return(col_data)
}

#' Verify consistency between LRT steps
#'
#' @param col_data Column data from DESeq2 object
#' @param design_formula Design formula
#' @param args Command-line arguments
#' @return NULL (validation function)
#' @export
verify_step_consistency <- function(col_data, design_formula, args) {
  log_message("Verifying consistency between LRT steps")
  
  # Check sample numbers
  log_message(paste("Number of samples in step 2:", nrow(col_data)))
  
  # Check design formula variables
  design_vars <- all.vars(design_formula)
  log_message(paste("Design variables:", paste(design_vars, collapse = ", ")))
  
  # Verify factor levels
  for (var in design_vars) {
    if (is.factor(col_data[[var]])) {
      levels_count <- length(levels(col_data[[var]]))
      log_message(paste("Factor", var, "has", levels_count, "levels:", 
                       paste(levels(col_data[[var]]), collapse = ", ")))
    }
  }
  
  log_message("Step consistency verification completed")
}

#' Apply batch correction to count data
#'
#' @param count_data Count matrix (genes x samples)
#' @param metadata_df Sample metadata
#' @param batch_method Batch correction method
#' @param normalized Whether counts are already normalized
#' @return Batch-corrected count matrix
#' @export
apply_batch_correction <- function(count_data, metadata_df, batch_method, normalized = TRUE) {
  log_message(paste("Applying", batch_method, "batch correction"))
  
  if (batch_method == "none") {
    log_message("No batch correction applied")
    return(count_data)
  }
  
  # Check for batch column
  if (!"batch" %in% colnames(metadata_df)) {
    log_warning("No 'batch' column found in metadata - skipping batch correction")
    return(count_data)
  }
  
  batch_factor <- factor(metadata_df$batch)
  
  if (batch_method == "limmaremovebatcheffect") {
    log_message("Applying limma removeBatchEffect")
    
    if (!requireNamespace("limma", quietly = TRUE)) {
      log_warning("limma package not available - skipping batch correction")
      return(count_data)
    }
    
    # For ATAC-seq data, we typically work with log2-transformed counts
    log_counts <- log2(count_data + 1)
    
    # Apply batch correction
    corrected_log_counts <- limma::removeBatchEffect(log_counts, batch = batch_factor)
    
    # Convert back to count scale
    corrected_counts <- 2^corrected_log_counts - 1
    corrected_counts[corrected_counts < 0] <- 0
    
    return(corrected_counts)
    
  } else if (batch_method == "combatseq") {
    log_message("Applying ComBat-Seq batch correction")
    
    if (!requireNamespace("sva", quietly = TRUE)) {
      log_warning("sva package not available - skipping batch correction")
      return(count_data)
    }
    
    # ComBat-Seq expects integer counts
    integer_counts <- round(count_data)
    
    # Apply ComBat-Seq
    corrected_counts <- sva::ComBat_seq(integer_counts, batch = batch_factor)
    
    return(corrected_counts)
    
  } else if (batch_method == "model") {
    log_message("Batch correction handled in model design")
    # For model-based correction, no preprocessing needed
    return(count_data)
  }
  
  log_warning(paste("Unknown batch correction method:", batch_method))
  return(count_data)
}

#' Rebuild DESeq dataset with updated reference levels
#'
#' @param dds_base Original DESeq dataset
#' @param factor_ref_list Named list of factors and their reference levels
#' @return Rebuilt DESeq dataset
#' @export
rebuild_dds <- function(dds_base, factor_ref_list) {
  log_message("Rebuilding DESeq dataset with updated reference levels")
  
  if (length(factor_ref_list) == 0) {
    log_message("No reference levels to update")
    return(dds_base)
  }
  
  # Get column data
  col_data <- colData(dds_base)
  
  # Update reference levels
  for (factor_name in names(factor_ref_list)) {
    ref_level <- factor_ref_list[[factor_name]]
    
    if (factor_name %in% colnames(col_data)) {
      current_levels <- levels(col_data[[factor_name]])
      
      if (ref_level %in% current_levels) {
        log_message(paste("Setting reference level for", factor_name, "to", ref_level))
        col_data[[factor_name]] <- relevel(col_data[[factor_name]], ref = ref_level)
      } else {
        log_warning(paste("Reference level", ref_level, "not found in factor", factor_name,
                         "Available levels:", paste(current_levels, collapse = ", ")))
      }
    } else {
      log_warning(paste("Factor", factor_name, "not found in column data"))
    }
  }
  
  # Update column data in DESeq dataset
  colData(dds_base) <- col_data
  
  # Rebuild the DESeq dataset to ensure proper design matrix
  dds_rebuilt <- DESeqDataSet(dds_base, design = design(dds_base))
  
  log_message("DESeq dataset rebuilt successfully")
  return(dds_rebuilt)
}

#' Export ATAC-seq results in standard format
#'
#' @param results_df Results data frame
#' @param output_file Output file path
#' @return NULL
#' @export
write_atac_results <- function(results_df, output_file) {
  log_message(paste("Writing ATAC-seq results to", output_file))
  
  # Ensure required columns exist
  if (!"RefseqId" %in% colnames(results_df)) {
    results_df$RefseqId <- rownames(results_df)
  }
  
  if (!"GeneId" %in% colnames(results_df)) {
    results_df$GeneId <- rownames(results_df)
  }
  
  # Write results
  write.table(results_df, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  log_message(paste("ATAC-seq results written successfully to", output_file))
}

#' Create interactive MDS plot for ATAC-seq data
#'
#' @param count_matrix Count matrix (peaks x samples)
#' @param output_file Output HTML file path
#' @param metadata Sample metadata
#' @param color_by Column to color points by
#' @param title Plot title
#' @return Output file path
#' @export
generate_atac_mds_plot <- function(count_matrix, output_file, metadata = NULL, 
                                   color_by = NULL, title = "ATAC-seq MDS Plot") {
  log_message("Generating ATAC-seq MDS plot")
  
  # Use the common visualization function
  return(generate_mds_plot_html(
    count_matrix,
    output_file,
    metadata = metadata,
    color_by = color_by,
    title = title
  ))
}

#' Verify ATAC-seq outputs
#'
#' @param output_prefix Output prefix
#' @param workflow_type Type of workflow
#' @param fail_on_missing Whether to fail if files are missing
#' @return NULL
#' @export
verify_outputs <- function(output_prefix, workflow_type, fail_on_missing = FALSE) {
  log_message("Verifying ATAC-seq output files")
  
  # Expected output files for ATAC-seq LRT step 2
  expected_files <- c(
    paste0(output_prefix, "_atac_results_table.tsv"),
    paste0(output_prefix, "_normalized_counts.gct"),
    "mds_plot.html"
  )
  
  missing_files <- c()
  
  for (file in expected_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
      log_warning(paste("Expected output file not found:", file))
    } else {
      log_message(paste("Found expected output file:", file))
    }
  }
  
  if (length(missing_files) > 0 && fail_on_missing) {
    stop(paste("Missing required output files:", paste(missing_files, collapse = ", ")))
  }
  
  log_message("Output verification completed")
}