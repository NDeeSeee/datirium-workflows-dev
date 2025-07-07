#!/usr/bin/env Rscript
#
# Base data processing functions for ATAC-seq workflows
# This module provides common data loading and validation functions
#

#' Load and validate ATAC-seq metadata with common patterns
#'
#' @param meta_file Path to metadata file
#' @param required_columns Vector of required column names
#' @param clean_names Whether to clean sample names (default TRUE)
#' @return Validated metadata data frame
#' @export
load_atac_metadata <- function(meta_file, required_columns = NULL, clean_names = TRUE) {
  log_message("Loading ATAC-seq metadata...")
  
  # Get the file delimiter
  delimiter <- check_file_delimiter(meta_file)
  
  # Load metadata
  metadata_df <- read.table(
    meta_file,
    sep = delimiter,
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  # Clean metadata column and row names if requested
  if (clean_names) {
    colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
    rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  }
  
  log_message(paste("Loaded metadata for", nrow(metadata_df), "samples with", 
                   ncol(metadata_df), "covariates"))
  
  # Validate required columns exist
  if (!is.null(required_columns)) {
    missing_cols <- required_columns[!required_columns %in% colnames(metadata_df)]
    if (length(missing_cols) > 0) {
      stop(paste("Required columns missing from metadata:", paste(missing_cols, collapse=", ")))
    }
  }
  
  return(metadata_df)
}

#' Validate ATAC-seq input files
#'
#' @param input_files Vector of peak file paths
#' @param bamfiles Vector of BAM file paths
#' @param sample_names Vector of sample names
#' @param check_sizes Whether to check file sizes (default TRUE)
#' @return NULL (validation function)
#' @export
validate_atac_files <- function(input_files, bamfiles, sample_names, check_sizes = TRUE) {
  log_message("Validating ATAC-seq input files...")
  
  # Check file counts match
  if (length(input_files) != length(bamfiles)) {
    stop("Number of peak files must match number of BAM files")
  }
  
  if (length(input_files) != length(sample_names)) {
    stop("Number of input files must match number of sample names")
  }
  
  # Check peak file existence
  for (i in seq_along(input_files)) {
    file <- input_files[i]
    if (!file.exists(file)) {
      stop(paste("Peak file not found:", file, "(sample:", sample_names[i], ")"))
    }
    
    if (check_sizes) {
      file_size <- file.info(file)$size
      if (file_size < 100) {
        warning(paste("Peak file", basename(file), "is very small (", file_size, "bytes) - may be empty"))
      }
    }
  }
  
  # Check BAM file existence
  for (i in seq_along(bamfiles)) {
    file <- bamfiles[i]
    if (!file.exists(file)) {
      stop(paste("BAM file not found:", file, "(sample:", sample_names[i], ")"))
    }
    
    if (check_sizes) {
      file_size <- file.info(file)$size
      if (file_size < 1000) {
        warning(paste("BAM file", basename(file), "is very small (", file_size, "bytes) - may be a dummy file"))
      }
    }
  }
  
  log_message("ATAC-seq input files validated successfully")
}

#' Create DiffBind sample sheet with common structure
#'
#' @param input_files Vector of peak file paths
#' @param bamfiles Vector of BAM file paths
#' @param sample_names Vector of sample names
#' @param metadata_df Metadata data frame
#' @param peakcaller Peak caller name
#' @param condition_column Column name for condition (optional)
#' @return DiffBind sample sheet data frame
#' @export
create_atac_sample_sheet <- function(input_files, bamfiles, sample_names, metadata_df, 
                                    peakcaller = "macs", condition_column = NULL) {
  log_message("Creating DiffBind sample sheet...")
  
  # Filter to only include samples in metadata
  sample_mask <- sample_names %in% rownames(metadata_df)
  
  if (!any(sample_mask)) {
    stop("No sample names match metadata row names")
  }
  
  if (!all(sample_mask)) {
    log_warning(paste("Filtering out", sum(!sample_mask), "samples not found in metadata"))
    input_files <- input_files[sample_mask]
    bamfiles <- bamfiles[sample_mask]
    sample_names <- sample_names[sample_mask]
  }
  
  # Create the basic sample sheet structure required by DiffBind
  sample_sheet <- data.frame(
    SampleID = sample_names,
    Tissue = "N",  # Default tissue type
    Factor = "ATAC",  # Factor being analyzed
    Condition = if (!is.null(condition_column)) metadata_df[sample_names, condition_column] else "unknown",
    Treatment = if (!is.null(condition_column)) metadata_df[sample_names, condition_column] else "unknown",
    Replicate = seq_along(sample_names),  # Sequential replicate numbers
    bamReads = bamfiles,
    Peaks = input_files,
    PeakCaller = peakcaller,
    stringsAsFactors = FALSE
  )
  
  # Add all metadata columns to the sample sheet
  for (col_name in colnames(metadata_df)) {
    if (!col_name %in% colnames(sample_sheet)) {
      sample_sheet[[col_name]] <- metadata_df[sample_names, col_name]
    }
  }
  
  # Remove any rows with missing critical metadata
  complete_rows <- complete.cases(sample_sheet[, c("SampleID", "bamReads", "Peaks")])
  if (!all(complete_rows)) {
    warning(paste("Removing", sum(!complete_rows), "samples with missing critical data"))
    sample_sheet <- sample_sheet[complete_rows, , drop = FALSE]
  }
  
  # Validate final sample sheet
  if (nrow(sample_sheet) == 0) {
    stop("No valid samples remain after filtering")
  }
  
  log_message(paste("Created DiffBind sample sheet for", nrow(sample_sheet), "samples"))
  
  return(sample_sheet)
}

#' Apply batch correction for ATAC-seq data with common methods
#'
#' @param count_data Count matrix (peaks x samples)
#' @param metadata_df Sample metadata
#' @param batch_method Batch correction method
#' @param batch_column Column name for batch information (default "batch")
#' @return Batch-corrected count matrix
#' @export
apply_atac_batch_correction <- function(count_data, metadata_df, batch_method, batch_column = "batch") {
  log_message(paste("Applying", batch_method, "batch correction for ATAC-seq data"))
  
  if (batch_method == "none") {
    log_message("No batch correction applied")
    return(count_data)
  }
  
  # Check for batch column
  if (!batch_column %in% colnames(metadata_df)) {
    log_warning(paste("No", batch_column, "column found in metadata - skipping batch correction"))
    return(count_data)
  }
  
  batch_factor <- factor(metadata_df[[batch_column]])
  
  # Check if we have multiple batches
  if (length(levels(batch_factor)) < 2) {
    log_warning("Only one batch found - skipping batch correction")
    return(count_data)
  }
  
  if (batch_method == "limmaremovebatcheffect") {
    log_message("Applying limma removeBatchEffect")
    
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
    
    log_message("limma batch correction completed")
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
    
    log_message("ComBat-Seq batch correction completed")
    return(corrected_counts)
    
  } else if (batch_method == "model") {
    log_message("Batch correction will be handled in model design")
    # For model-based correction, no preprocessing needed
    return(count_data)
  }
  
  log_warning(paste("Unknown batch correction method:", batch_method))
  return(count_data)
}

#' Filter metadata to specific conditions for pairwise analysis
#'
#' @param metadata_df Metadata data frame
#' @param condition_column Column name for conditions
#' @param condition1 First condition name
#' @param condition2 Second condition name
#' @param reference_condition Which condition should be the reference (default: condition2)
#' @return Filtered and processed metadata
#' @export
filter_metadata_for_pairwise <- function(metadata_df, condition_column, condition1, condition2, 
                                        reference_condition = condition2) {
  log_message("Filtering metadata for pairwise comparison...")
  
  # Validate condition column exists
  if (!condition_column %in% colnames(metadata_df)) {
    stop(paste("Condition column", condition_column, "not found in metadata"))
  }
  
  # Check that both conditions exist in the data
  available_conditions <- unique(metadata_df[[condition_column]])
  
  if (!condition1 %in% available_conditions) {
    stop(paste("Condition1", condition1, "not found in", condition_column, 
               "Available conditions:", paste(available_conditions, collapse=", ")))
  }
  
  if (!condition2 %in% available_conditions) {
    stop(paste("Condition2", condition2, "not found in", condition_column,
               "Available conditions:", paste(available_conditions, collapse=", ")))
  }
  
  # Filter metadata to only include samples from the two conditions of interest
  condition_mask <- metadata_df[[condition_column]] %in% c(condition1, condition2)
  filtered_metadata <- metadata_df[condition_mask, , drop = FALSE]
  
  log_message(paste("Filtered to", nrow(filtered_metadata), "samples for pairwise comparison"))
  log_message(paste("Condition1 (", condition1, "):", 
                   sum(filtered_metadata[[condition_column]] == condition1), "samples"))
  log_message(paste("Condition2 (", condition2, "):", 
                   sum(filtered_metadata[[condition_column]] == condition2), "samples"))
  
  # Ensure condition column is a factor with proper reference level
  filtered_metadata[[condition_column]] <- factor(
    filtered_metadata[[condition_column]], 
    levels = c(reference_condition, setdiff(c(condition1, condition2), reference_condition))
  )
  
  return(filtered_metadata)
}

#' Create expression data frame with genomic coordinates
#'
#' @param peak_ranges GenomicRanges object with peak coordinates
#' @param deseq_results DESeq2 results object (optional)
#' @return Expression data frame with coordinates
#' @export
create_atac_expression_df <- function(peak_ranges, deseq_results = NULL) {
  log_message("Creating ATAC-seq expression data frame...")
  
  # Create base data frame with genomic coordinates
  expr_df <- data.frame(
    RefseqId = names(peak_ranges) %||% paste0("peak_", seq_along(peak_ranges)),
    GeneId = names(peak_ranges) %||% paste0("peak_", seq_along(peak_ranges)),
    chromosome = seqnames(peak_ranges),
    start = start(peak_ranges),
    end = end(peak_ranges),
    peak_width = width(peak_ranges),
    stringsAsFactors = FALSE
  )
  
  # Add peak coordinate as row names
  peak_ids <- paste0(seqnames(peak_ranges), ":", start(peak_ranges), "-", end(peak_ranges))
  rownames(expr_df) <- peak_ids
  
  # Add DESeq2 results if provided
  if (!is.null(deseq_results)) {
    # Ensure row names match
    if (all(rownames(deseq_results) == peak_ids)) {
      expr_df$baseMean <- deseq_results$baseMean
      expr_df$log2FoldChange <- deseq_results$log2FoldChange
      expr_df$lfcSE <- deseq_results$lfcSE
      expr_df$stat <- deseq_results$stat
      expr_df$pvalue <- deseq_results$pvalue
      expr_df$padj <- deseq_results$padj
    } else {
      log_warning("Row names don't match between peak ranges and DESeq2 results")
    }
  }
  
  log_message(paste("Created expression data frame with", nrow(expr_df), "peaks"))
  return(expr_df)
}

#' Write ATAC-seq results with standard format
#'
#' @param results_df Results data frame
#' @param output_file Output file path
#' @param comparison_info Optional comparison information to add as header
#' @export
write_atac_results <- function(results_df, output_file, comparison_info = NULL) {
  log_message(paste("Writing ATAC-seq results to", output_file))
  
  # Ensure required columns exist
  if (!"RefseqId" %in% colnames(results_df)) {
    results_df$RefseqId <- rownames(results_df)
  }
  
  if (!"GeneId" %in% colnames(results_df)) {
    results_df$GeneId <- rownames(results_df)
  }
  
  # Write header comment if provided
  if (!is.null(comparison_info)) {
    writeLines(paste0("# ", comparison_info), output_file)
    write.table(results_df, file = output_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = TRUE, append = TRUE)
  } else {
    write.table(results_df, file = output_file, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = TRUE)
  }
  
  log_message(paste("ATAC-seq results written successfully to", output_file))
}

#' Verify standard ATAC-seq output files
#'
#' @param output_prefix Output prefix
#' @param expected_files Vector of expected file suffixes
#' @param fail_on_missing Whether to fail if files are missing
#' @export
verify_atac_outputs <- function(output_prefix, expected_files = NULL, fail_on_missing = FALSE) {
  log_message("Verifying ATAC-seq output files")
  
  # Default expected files if not provided
  if (is.null(expected_files)) {
    expected_files <- c(
      "_results.tsv",
      "_normalized_counts.gct",
      "_mds_plot.html"
    )
  }
  
  # Add prefix to file names
  full_expected_files <- paste0(output_prefix, expected_files)
  
  missing_files <- c()
  
  for (file in full_expected_files) {
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
  
  log_message("ATAC-seq output verification completed")
  return(missing_files)
}