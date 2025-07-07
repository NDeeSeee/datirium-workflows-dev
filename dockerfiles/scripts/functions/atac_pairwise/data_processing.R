#!/usr/bin/env Rscript
#
# Data processing functions for ATAC-seq Pairwise Analysis
#

#' Load and validate metadata for pairwise comparison
#'
#' @param args Command-line arguments
#' @return Validated metadata data frame
#' @export
load_and_validate_pairwise_metadata <- function(args) {
  log_message("Loading metadata for ATAC-seq pairwise analysis...")
  
  # For pairwise analysis, create metadata from sample information
  if (is.null(args$meta)) {
    log_message("Creating metadata from sample information for pairwise analysis")
    
    # Create metadata data frame from sample names and conditions
    all_sample_names <- args$name
    all_conditions <- c(
      rep(args$condition2, length(args$untreated_sample_names)),  # Reference condition first
      rep(args$condition1, length(args$treated_sample_names))     # Treatment condition second
    )
    
    metadata_df <- data.frame(
      condition = all_conditions,
      stringsAsFactors = FALSE,
      row.names = all_sample_names
    )
    
    # Set condition column name
    args$condition_column <- "condition"
    
  } else {
    # Load metadata from file if provided
    log_message("Loading metadata from file...")
    
    # Get the file delimiter
    delimiter <- check_file_delimiter(args$meta)
    
    # Load metadata
    metadata_df <- read.table(
      args$meta,
      sep = delimiter,
      header = TRUE,
      stringsAsFactors = FALSE,
      row.names = 1
    )
    
    # Clean metadata column and row names
    colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
    rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  }
  
  log_message(paste("Loaded metadata for", nrow(metadata_df), "samples with", 
                   ncol(metadata_df), "covariates"))
  
  # Validate condition column exists
  if (!args$condition_column %in% colnames(metadata_df)) {
    stop(paste("Condition column", args$condition_column, "not found in metadata"))
  }
  
  # Check that both conditions exist in the data
  available_conditions <- unique(metadata_df[[args$condition_column]])
  
  if (!args$condition1 %in% available_conditions) {
    stop(paste("Condition1", args$condition1, "not found in", args$condition_column, 
               "Available conditions:", paste(available_conditions, collapse=", ")))
  }
  
  if (!args$condition2 %in% available_conditions) {
    stop(paste("Condition2", args$condition2, "not found in", args$condition_column,
               "Available conditions:", paste(available_conditions, collapse=", ")))
  }
  
  # Filter metadata to only include samples from the two conditions of interest
  condition_mask <- metadata_df[[args$condition_column]] %in% c(args$condition1, args$condition2)
  metadata_df <- metadata_df[condition_mask, , drop = FALSE]
  
  log_message(paste("Filtered to", nrow(metadata_df), "samples for pairwise comparison"))
  log_message(paste("Condition1 (", args$condition1, "):", 
                   sum(metadata_df[[args$condition_column]] == args$condition1), "samples"))
  log_message(paste("Condition2 (", args$condition2, "):", 
                   sum(metadata_df[[args$condition_column]] == args$condition2), "samples"))
  
  # Ensure condition column is a factor with condition2 as reference
  metadata_df[[args$condition_column]] <- factor(
    metadata_df[[args$condition_column]], 
    levels = c(args$condition2, args$condition1)
  )
  
  # Add design formula for simple pairwise comparison
  design_formula <- as.formula(paste("~", args$condition_column))
  attr(metadata_df, "design_formula") <- design_formula
  
  return(metadata_df)
}

#' Load and validate ATAC-seq data for pairwise comparison
#'
#' @param args Command-line arguments
#' @param metadata_df Validated metadata
#' @return DiffBind sample sheet
#' @export
load_and_validate_pairwise_atac_data <- function(args, metadata_df) {
  log_message("Loading ATAC-seq data for pairwise analysis...")
  
  # Clean sample names for consistency
  clean_names <- clean_sample_names(args$name)
  
  # Filter files and names to only include samples in the two conditions
  sample_mask <- clean_names %in% rownames(metadata_df)
  
  if (!any(sample_mask)) {
    stop("No sample names match metadata row names")
  }
  
  # Filter to only relevant samples
  filtered_input_files <- args$input_files[sample_mask]
  filtered_names <- clean_names[sample_mask]
  
  # Derive BAM paths next to each peak file (CWL secondaryFiles staging)
  log_message("Deriving BAM paths from peak file locations (no --bamfiles supplied)")
  filtered_bamfiles <- vapply(seq_along(filtered_input_files), function(i) {
    peak_path <- filtered_input_files[i]
    peak_dir  <- dirname(peak_path)
    peak_base <- basename(peak_path)
    # Handle both *_peaks.csv and generic names
    bam_base  <- sub("_peaks\\.csv$", ".bam", peak_base)
    file.path(peak_dir, bam_base)
  }, character(1))
  
  log_message(paste("Filtered to", length(filtered_names), "samples matching metadata"))
  
  # Create DiffBind sample sheet for pairwise comparison
  sample_sheet <- create_pairwise_diffbind_sample_sheet(
    filtered_input_files, filtered_bamfiles, filtered_names, metadata_df, args
  )
  
  log_message(paste("Created DiffBind sample sheet for", nrow(sample_sheet), "samples"))
  
  return(sample_sheet)
}

#' Create DiffBind sample sheet for pairwise comparison
#'
#' @param input_files Peak files
#' @param bamfiles BAM files
#' @param clean_names Sample names
#' @param metadata_df Metadata
#' @param args Command-line arguments
#' @return DiffBind sample sheet
#' @export
create_pairwise_diffbind_sample_sheet <- function(input_files, bamfiles, clean_names, metadata_df, args) {
  # Set default values for missing arguments
  peakcaller <- if (is.null(args$peakcaller)) "macs" else args$peakcaller
  
  # Create the basic sample sheet structure required by DiffBind
  sample_sheet <- data.frame(
    SampleID = clean_names,
    Tissue = rep("N", length(clean_names)),  # Default tissue type
    Factor = rep("ATAC", length(clean_names)),  # Factor being analyzed
    Condition = metadata_df[clean_names, args$condition_column],
    Treatment = metadata_df[clean_names, args$condition_column],
    Replicate = seq_along(clean_names),  # Sequential replicate numbers
    bamReads = bamfiles,
    Peaks = input_files,
    PeakCaller = rep(peakcaller, length(clean_names)),
    stringsAsFactors = FALSE
  )
  
  # Add all metadata columns to the sample sheet
  for (col_name in colnames(metadata_df)) {
    if (col_name != args$condition_column) {  # Avoid duplication
      sample_sheet[[col_name]] <- metadata_df[clean_names, col_name]
    }
  }
  
  # Remove any rows with missing metadata
  complete_rows <- complete.cases(sample_sheet)
  if (!all(complete_rows)) {
    warning(paste("Removing", sum(!complete_rows), "samples with missing metadata"))
    sample_sheet <- sample_sheet[complete_rows, , drop = FALSE]
  }
  
  # Check that we still have samples from both conditions
  remaining_conditions <- unique(sample_sheet$Condition)
  if (length(remaining_conditions) < 2) {
    stop("After filtering, samples from both conditions are not available")
  }
  
  # Validate file accessibility
  if (!is.null(args$test_mode) && args$test_mode) {
    log_message("Test mode enabled – skipping BAM existence validation")
  } else {
    # Validate file accessibility
    for (i in 1:nrow(sample_sheet)) {
      # Only check BAM files if they're real files (not dummy paths)
      if (file.exists(dirname(sample_sheet$bamReads[i])) && !grepl("^[^/]*\\.bam$", sample_sheet$bamReads[i])) {
        if (!file.exists(sample_sheet$bamReads[i])) {
          stop(paste("BAM file not found:", sample_sheet$bamReads[i]))
        }
      }
      if (!file.exists(sample_sheet$Peaks[i])) {
        stop(paste("Peak file not found:", sample_sheet$Peaks[i]))
      }
    }
  }
  
  return(sample_sheet)
}

#' Apply batch correction for pairwise ATAC-seq data
#'
#' @param count_data Count matrix
#' @param metadata_df Sample metadata
#' @param batch_method Batch correction method
#' @return Batch-corrected count matrix
#' @export
apply_pairwise_batch_correction <- function(count_data, metadata_df, batch_method) {
  log_message(paste("Applying", batch_method, "batch correction for pairwise analysis"))
  
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
    log_message("Applying limma removeBatchEffect for pairwise comparison")
    
    if (!requireNamespace("limma", quietly = TRUE)) {
      log_warning("limma package not available - skipping batch correction")
      return(count_data)
    }
    
    # For ATAC-seq data, work with log2-transformed counts
    log_counts <- log2(count_data + 1)
    
    # Apply batch correction
    corrected_log_counts <- limma::removeBatchEffect(log_counts, batch = batch_factor)
    
    # Convert back to count scale
    corrected_counts <- 2^corrected_log_counts - 1
    corrected_counts[corrected_counts < 0] <- 0
    
    return(corrected_counts)
    
  } else if (batch_method == "combatseq") {
    log_message("Applying ComBat-Seq batch correction for pairwise comparison")
    
    if (!requireNamespace("sva", quietly = TRUE)) {
      log_warning("sva package not available - skipping batch correction")
      return(count_data)
    }
    
    # ComBat-Seq expects integer counts
    integer_counts <- round(count_data)
    
    # Apply ComBat-Seq
    corrected_counts <- sva::ComBat_seq(integer_counts, batch = batch_factor)
    
    return(corrected_counts)
  }
  
  log_warning(paste("Unknown batch correction method:", batch_method))
  return(count_data)
}

#' Export ATAC-seq pairwise results
#'
#' @param results_df Results data frame
#' @param output_file Output file path
#' @param args Command-line arguments
#' @export
write_pairwise_atac_results <- function(results_df, output_file, args) {
  log_message(paste("Writing ATAC-seq pairwise results to", output_file))
  
  # Ensure required columns exist
  if (!"RefseqId" %in% colnames(results_df)) {
    results_df$RefseqId <- rownames(results_df)
  }
  
  if (!"GeneId" %in% colnames(results_df)) {
    results_df$GeneId <- rownames(results_df)
  }
  
  # Add comparison information as a comment
  comparison_info <- paste0("# ATAC-seq Pairwise Comparison: ", args$condition1, " vs ", args$condition2)
  
  # Write header comment and results
  writeLines(comparison_info, output_file)
  write.table(results_df, file = output_file, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE, append = TRUE)
  
  log_message(paste("ATAC-seq pairwise results written successfully to", output_file))
}

#' Verify pairwise ATAC-seq outputs
#'
#' @param output_prefix Output prefix
#' @param fail_on_missing Whether to fail if files are missing
#' @export
verify_pairwise_outputs <- function(output_prefix, fail_on_missing = FALSE) {
  log_message("Verifying ATAC-seq pairwise output files")
  
  # Expected output files for ATAC-seq pairwise analysis
  expected_files <- c(
    paste0(output_prefix, "_pairwise_results.tsv"),
    paste0(output_prefix, "_normalized_counts.gct"),
    "atac_pairwise_mds_plot.html"
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
  
  log_message("ATAC-seq pairwise output verification completed")
}