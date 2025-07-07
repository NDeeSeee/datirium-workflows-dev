#!/usr/bin/env Rscript
#
# Contrast analysis functions for ATAC-seq LRT Step 2
#

#' Extract factors and levels from interaction terms
#'
#' @param term Interaction term like "treatmentKO.condRest"
#' @return List with factor1, level1, factor2, level2
#' @export
extract_factors_and_levels <- function(term) {
  parts <- strsplit(term, "\\.")[[1]]
  factor1 <- sub("[0-9A-Z_]+$", "", parts[1])
  factor2 <- sub("[0-9A-Z_]+$", "", parts[2])
  level1 <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[1])
  level2 <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[2])
  list(factor1 = factor1, level1 = level1, factor2 = factor2, level2 = level2)
}

#' Process a specific contrast and obtain DESeq2 results for ATAC-seq data
#'
#' @param step1_data Data object from Step 1 containing DESeq dataset and metadata
#' @param contrast_row Contrast information row
#' @param args Command-line arguments
#' @return DESeq2 results object
#' @export
get_contrast_res <- function(step1_data, contrast_row, args) {
  # Extract the DESeq2 dataset from step1_data
  dds_base <- step1_data$dds
  
  # Set altHypothesis based on regulation parameter
  altHypothesis <- switch(args$regulation,
                         "up"   = "greater",
                         "down" = "less",
                         "both" = "greaterAbs",
                         "greaterAbs")  # default

  lfcThreshold <- if (args$use_lfc_thresh) args$lfcthreshold else 0
  
  log_message(paste("Processing ATAC-seq contrast:", contrast_row$contrast))
  log_message(paste("Effect type:", contrast_row$effect))
  log_message(paste("Setting altHypothesis to", altHypothesis))
  
  # Check if we're in test mode - if so, return mock results
  if (args$test_mode) {
    log_message("TEST MODE: Generating mock contrast results")
    return(create_mock_contrast_results(contrast_row, args))
  }
  
  # Handle different contrast types
  if (contrast_row$effect == "main") {
    # --- MAIN EFFECT CONTRAST ---
    return(process_atac_main_effect_contrast(dds_base, contrast_row, args, lfcThreshold, altHypothesis))
  } else if (contrast_row$effect == "interaction") {
    # --- INTERACTION CONTRAST ---
    return(process_atac_interaction_contrast(dds_base, contrast_row, args, lfcThreshold, altHypothesis))
  } else {
    stop("Unknown effect in contrast_row: ", contrast_row$effect)
  }
}

#' Process a main effect contrast for ATAC-seq data
#'
#' @param dds_base DESeq2 data object
#' @param contrast_row Contrast information row
#' @param args Command-line arguments
#' @param lfcThreshold Log fold change threshold
#' @param altHypothesis Alternative hypothesis
#' @return DESeq2 results object
#' @export
process_atac_main_effect_contrast <- function(dds_base, contrast_row, args, lfcThreshold, altHypothesis) {
  # Get column names from DESeq2 object to identify available factors
  col_names <- colnames(colData(dds_base))
  
  # Initialize factor reference list
  factor_ref_list <- list()
  
  # Check if the contrast involves a simple one-factor comparison or multiple factors
  if (contrast_row$specificity_group %in% col_names) {
    # SINGLE-FACTOR: specificity_group is the factor to re-level
    log_message(paste("Single-factor ATAC-seq contrast for", contrast_row$specificity_group))
    factor_ref_list <- setNames(list(contrast_row$denominator), contrast_row$specificity_group)
  } else {
    # MULTI-FACTOR: specificity_group might contain additional information
    # Try to parse the specificity_group to extract factor and level information
    log_message(paste("Multi-factor ATAC-seq contrast with specificity group:", contrast_row$specificity_group))
    
    # First check if specificity_group contains underscore separators
    if (grepl("_", contrast_row$specificity_group)) {
      parts <- strsplit(contrast_row$specificity_group, "_")[[1]]
      
      if (length(parts) >= 2) {
        # The first part is likely the additional factor name
        other_factor <- parts[1]
        # The remaining parts combined form the level value
        other_level <- paste(parts[-1], collapse = "_")
        
        # Extract main factor from contrast name
        # Expected format like: mainFactor_level_vs_ref or similar
        if (grepl("_vs_", contrast_row$contrast)) {
          contrast_parts <- strsplit(contrast_row$contrast, "_vs_")[[1]]
          if (length(contrast_parts) >= 1) {
            main_parts <- strsplit(contrast_parts[1], "_")[[1]]
            main_factor <- main_parts[1]
            
            # Set up references for both factors
            factor_ref_list[[main_factor]] <- contrast_row$denominator
            factor_ref_list[[other_factor]] <- other_level
            
            log_message(paste("Setting references for ATAC-seq - Main factor:", main_factor, 
                              "=", contrast_row$denominator, ", Other factor:", 
                              other_factor, "=", other_level))
          }
        }
      }
    }
    
    # If we couldn't parse the multi-factor information, fall back to a simpler approach
    if (length(factor_ref_list) == 0) {
      log_warning(paste("Couldn't parse multi-factor specificity group:", 
                         contrast_row$specificity_group, 
                         "- Using basic approach"))
      
      # Try to extract factor name from the contrast string
      # Common format: factorName_level_vs_reference
      if (grepl("_", contrast_row$contrast)) {
        main_parts <- strsplit(contrast_row$contrast, "_")[[1]]
        # Assume first part is the factor name
        main_factor <- main_parts[1]
        factor_ref_list <- setNames(list(contrast_row$denominator), main_factor)
        log_message(paste("Falling back to single factor for ATAC-seq:", main_factor, 
                           "=", contrast_row$denominator))
      } else {
        # Last resort - use the contrast name as the factor
        factor_ref_list <- setNames(list(contrast_row$denominator), contrast_row$contrast)
        log_warning(paste("Could not determine factor from ATAC-seq contrast - using contrast name as factor"))
      }
    }
  }
  
  # If we still don't have any factors to work with, raise an error
  if (length(factor_ref_list) == 0) {
    stop("Could not determine factors for re-leveling from ATAC-seq contrast information")
  }
  
  # Rebuild DESeq dataset with the proper reference levels
  log_message("Rebuilding DESeq dataset for ATAC-seq with updated reference levels")
  dds_subset <- rebuild_dds(dds_base, factor_ref_list)
  
  # Run DESeq analysis
  log_message("Running DESeq with Wald test for ATAC-seq")
  dds_subset <- DESeq(dds_subset, test = "Wald")
  
  # Extract results
  log_message(paste("Extracting ATAC-seq results for contrast:", contrast_row$contrast))
  res <- results(dds_subset,
                 name = contrast_row$contrast,
                 alpha = args$fdr,
                 lfcThreshold = lfcThreshold,
                 independentFiltering = TRUE,
                 altHypothesis = altHypothesis)
  
  log_message("DESeq2 results for ATAC-seq main effect contrast obtained successfully")
  return(res)
}

#' Process an interaction contrast for ATAC-seq data
#'
#' @param dds_base DESeq2 data object
#' @param contrast_row Contrast information row
#' @param args Command-line arguments
#' @param lfcThreshold Log fold change threshold
#' @param altHypothesis Alternative hypothesis
#' @return DESeq2 results object
#' @export
process_atac_interaction_contrast <- function(dds_base, contrast_row, args, lfcThreshold, altHypothesis) {
  # Parse the contrast string to extract factors and levels
  log_message("Parsing ATAC-seq interaction contrast information")
  
  # Check if the contrast contains a colon (typical for interaction terms in formulas)
  if (grepl(":", contrast_row$contrast, fixed = TRUE)) {
    log_message("Processing DESeq2 formula-style ATAC-seq interaction contrast")
    # This is likely a DESeq2 formula-style interaction term
    # Format is typically like: "factorA_levelA:factorB_levelB"
    
    # Split on the colon to get the interacting terms
    interaction_terms <- strsplit(contrast_row$contrast, ":", fixed = TRUE)[[1]]
    
    if (length(interaction_terms) == 2) {
      factor_ref_list <- list()
      
      # Determine which factor needs to be releveled based on the specificity_group
      if (!is.null(contrast_row$specificity_group) && contrast_row$specificity_group != "") {
        sg_parts <- strsplit(contrast_row$specificity_group, "_")[[1]]
        relevel_factor <- sg_parts[1]
        
        # The denominator should contain the reference level
        # Try to extract just the level part from the denominator
        ref_level <- contrast_row$denominator
        if (startsWith(ref_level, relevel_factor)) {
          ref_level <- sub(paste0("^", relevel_factor), "", ref_level)
        }
        
        factor_ref_list <- setNames(list(ref_level), relevel_factor)
        log_message(paste("Setting reference level for ATAC-seq interaction -", relevel_factor, "=", ref_level))
      } else {
        # If specificity_group is not available, make a best guess
        log_warning("No specificity_group provided for ATAC-seq interaction contrast - making a best guess")
        
        # Try to extract factor names from the interaction terms
        for (term in interaction_terms) {
          if (grepl("_", term)) {
            term_parts <- strsplit(term, "_")[[1]]
            factor_name <- term_parts[1]
            level_name <- paste(term_parts[-1], collapse = "_")
            
            # Check if this factor exists in colData
            if (factor_name %in% colnames(colData(dds_base))) {
              # Get the current reference level (first level)
              current_levels <- levels(colData(dds_base)[[factor_name]])
              if (length(current_levels) > 0) {
                factor_ref_list[[factor_name]] <- current_levels[1]
                log_message(paste("Guessed reference level for ATAC-seq", factor_name, "=", current_levels[1]))
              }
            }
          }
        }
      }
      
      # If no releveling information was found, use an empty list
      if (length(factor_ref_list) == 0) {
        log_warning("Could not determine reference levels for ATAC-seq interaction contrast - using defaults")
      }
    } else {
      log_warning(paste("Complex ATAC-seq interaction term found with", length(interaction_terms), "components"))
      factor_ref_list <- list()
    }
  } else {
    # For other formats, use the helper function to extract factors and levels
    log_message("Using extract_factors_and_levels for non-formula ATAC-seq interaction contrast")
    factors_levels <- extract_factors_and_levels(contrast_row$contrast)
    factor1 <- factors_levels$factor1
    level1 <- factors_levels$level1
    factor2 <- factors_levels$factor2
    level2 <- factors_levels$level2
    
    log_message(paste("Parsed ATAC-seq interaction terms: factor1 =", factor1, "level1 =", level1,
                      "factor2 =", factor2, "level2 =", level2))
    
    # Determine which factor to relevel based on specificity_group
    sg_parts <- strsplit(contrast_row$specificity_group, "_")[[1]]
    
    if (length(sg_parts) > 0) {
      if (sg_parts[1] == factor2) {
        # Relevel factor1
        relevel_factor <- factor1
        ref_level <- sub(paste0("^", factor1), "", contrast_row$denominator)
      } else if (sg_parts[1] == factor1) {
        # Relevel factor2
        relevel_factor <- factor2
        ref_level <- sub(paste0("^", factor2), "", contrast_row$denominator)
      } else {
        log_warning(paste("Could not determine which factor to relevel from ATAC-seq specificity_group:", 
                          contrast_row$specificity_group))
        relevel_factor <- factor1  # Default to factor1
        ref_level <- sub(paste0("^", factor1), "", contrast_row$denominator)
      }
      
      factor_ref_list <- setNames(list(ref_level), relevel_factor)
      log_message(paste("Setting reference level for ATAC-seq", relevel_factor, "=", ref_level))
    } else {
      log_warning("Empty specificity_group for ATAC-seq interaction contrast")
      factor_ref_list <- list()
    }
  }
  
  # Rebuild the DESeqDataSet with proper reference levels
  log_message("Rebuilding DESeq dataset for ATAC-seq interaction contrast")
  dds_subset <- rebuild_dds(dds_base, factor_ref_list)
  
  # Run DESeq analysis
  log_message("Running DESeq with Wald test for ATAC-seq interaction contrast")
  dds_subset <- DESeq(dds_subset, test = "Wald")
  
  # Extract results
  log_message(paste("Extracting results for ATAC-seq interaction contrast:", contrast_row$contrast))
  res <- results(dds_subset,
                 name = contrast_row$contrast,
                 alpha = args$fdr,
                 lfcThreshold = lfcThreshold,
                 independentFiltering = TRUE,
                 altHypothesis = altHypothesis)
  
  log_message("DESeq2 results for ATAC-seq interaction contrast obtained successfully")
  return(res)
}

#' Add DESeq2 results as metadata columns to ATAC-seq expression data
#'
#' @param expression_data_df Expression data frame
#' @param deseq_result DESeq2 results object
#' @param contrast_name Name of the contrast
#' @return Expression data frame with added results
#' @export
add_metadata_to_results <- function(expression_data_df, deseq_result, contrast_name) {
  log_message(paste("Adding DESeq2 results for ATAC-seq contrast:", contrast_name))
  
  # Check if this is already a data frame (mock results)
  if (is.data.frame(deseq_result)) {
    log_message("Processing data frame results (test mode)")
    res_df <- deseq_result
  } else {
    log_message("Processing DESeq2 results object (normal mode)")
    # Create a data frame with results first
    res_df <- as.data.frame(deseq_result)
  }
  
  # Get available columns after ensuring we have a data frame
  available_columns <- colnames(res_df)
  
  # First, determine which columns to extract from the results
  result_columns <- c("baseMean", "log2FoldChange", "pvalue", "padj")
  
  # Ensure all required columns exist in the results
  missing_columns <- base::setdiff(result_columns, available_columns)
  if (length(missing_columns) > 0) {
    log_warning(paste("Missing expected columns in ATAC-seq DESeq2 results:", paste(missing_columns, collapse=", ")))
    result_columns <- base::intersect(result_columns, available_columns)
  }
  
  # Debug: log the structure of the results
  log_message(paste("ATAC-seq DESeq2 results dimensions:", nrow(res_df), "x", ncol(res_df)))
  log_message(paste("Available columns:", paste(colnames(res_df), collapse=", ")))
  log_message(paste("Target columns:", paste(result_columns, collapse=", ")))
  
  # Safely select only available columns
  selected_columns <- base::intersect(result_columns, colnames(res_df))
  if (length(selected_columns) == 0) {
    log_warning("No target columns found in ATAC-seq DESeq2 results - using all columns")
    selected_columns <- colnames(res_df)
  }
  
  res_df <- res_df[, selected_columns, drop = FALSE]
  
  # Rename columns for contrast-specific naming (ATAC-seq specific)
  if ("log2FoldChange" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "log2FoldChange"] <- paste0(contrast_name, "_LFC")
  }
  if ("pvalue" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "pvalue"] <- paste0(contrast_name, "_pvalue")
  }
  if ("padj" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "padj"] <- paste0(contrast_name, "_FDR")
  }

  log_message("ATAC-seq DESeq2 results prepared for merging")
  log_message(paste("Final res_df dimensions:", nrow(res_df), "x", ncol(res_df)))
  log_message(paste("Final res_df columns:", paste(colnames(res_df), collapse=", ")))
  log_message(paste("ATAC-seq expression data dimensions:", nrow(expression_data_df), "x", ncol(expression_data_df)))

  # Get column types to ensure consistent merging
  exp_df_classes <- sapply(expression_data_df, class)
  res_df_classes <- sapply(res_df, class)
  
  # Check for any potential issues with column types
  if ("baseMean" %in% colnames(expression_data_df) && !identical(exp_df_classes["baseMean"], res_df_classes["baseMean"])) {
    log_warning(paste("Column type mismatch for baseMean in ATAC-seq data:", 
                     exp_df_classes["baseMean"], "vs", res_df_classes["baseMean"]))
  }
  
  # Check for duplicate row names
  if (anyDuplicated(rownames(expression_data_df)) > 0) {
    log_warning("Duplicate row names found in ATAC-seq expression data - results may be incorrect")
  }
  
  if (anyDuplicated(rownames(res_df)) > 0) {
    log_warning("Duplicate row names found in ATAC-seq results data - results may be incorrect")
  }
  
  # Simply merge by adding results columns to the existing data frame
  # Ensure both data frames have the same row names for proper alignment
  log_message(paste("ATAC-seq expression data row names (first 3):", paste(head(rownames(expression_data_df), 3), collapse=", ")))
  log_message(paste("ATAC-seq results row names (first 3):", paste(head(rownames(res_df), 3), collapse=", ")))
  
  if (!all(rownames(expression_data_df) == rownames(res_df))) {
    log_warning("Row names don't match between ATAC-seq expression data and results - attempting to align")
    # Keep only common rows
    common_rows <- base::intersect(rownames(expression_data_df), rownames(res_df))
    if (length(common_rows) == 0) {
      stop("No common row names between ATAC-seq expression data and DESeq2 results")
    }
    expression_data_df <- expression_data_df[common_rows, , drop = FALSE]
    res_df <- res_df[common_rows, , drop = FALSE]
    log_message(paste("Aligned ATAC-seq data to", length(common_rows), "common rows"))
  }

  # Simple column binding
  expression_data_merged <- cbind(expression_data_df, res_df)
  
  log_message(paste("Successfully merged DESeq2 results for ATAC-seq contrast", contrast_name))
  return(expression_data_merged)
}

#' Create an MDS (Multi-Dimensional Scaling) plot for ATAC-seq data
#'
#' @param normCounts Normalized count matrix
#' @param col_metadata Metadata for the columns (samples)
#' @param output_file Output file path
#' @param interactive Whether to create an interactive HTML plot (default) or a static PDF
#' @param dimensions Number of dimensions to use for the plot
#' @param title Plot title
#' @param color_by Column in metadata to use for coloring points
#' @param shape_by Column in metadata to use for point shapes
#' @return Output file path if successful, NULL otherwise
#' @export
create_mds_plot <- function(normCounts, col_metadata, output_file = "atac_mds_plot.html", 
                           interactive = TRUE, dimensions = 2, title = "ATAC-seq MDS Plot",
                           color_by = NULL, shape_by = NULL) {
  log_message("Creating ATAC-seq MDS plot")
  
  # Generate the MDS plot using our consolidated function
  return(generate_atac_mds_plot(
    normCounts,
    output_file,
    metadata = col_metadata,
    color_by = color_by,
    title = title
  ))
}

#' Create mock contrast results for testing
#'
#' @param contrast_row Contrast information row
#' @param args Command-line arguments
#' @return Mock DESeq2 results object
#' @export
create_mock_contrast_results <- function(contrast_row, args) {
  log_message(paste("Creating mock results for contrast:", contrast_row$contrast))
  
  # Create mock results data frame with typical DESeq2 structure
  n_peaks <- 1000  # Mock number of peaks
  
  mock_results <- data.frame(
    baseMean = runif(n_peaks, 10, 1000),
    log2FoldChange = rnorm(n_peaks, 0, 2),
    lfcSE = runif(n_peaks, 0.1, 0.5),
    stat = rnorm(n_peaks, 0, 5),
    pvalue = runif(n_peaks, 0.001, 0.5),
    padj = runif(n_peaks, 0.001, 0.3),
    stringsAsFactors = FALSE
  )
  
  # Set row names to match the expression data format (simple numeric)
  rownames(mock_results) <- as.character(1:n_peaks)
  
  # For test mode, just return the simple data frame
  # The add_metadata_to_results function will handle it correctly
  log_message(paste("Mock results created with", nrow(mock_results), "peaks"))
  return(mock_results)
}