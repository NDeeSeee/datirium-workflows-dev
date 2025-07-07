#!/usr/bin/env Rscript

# --- Data processing functions for ATAC-seq LRT Step 1 ---

#' Load and process peak files for ATAC-seq analysis
#'
#' @param args Command line arguments containing input files and parameters
#' @return List containing processed peak data
load_peak_data <- function(args) {
  message("Loading peak files...")
  
  peak_files <- args$input_files
  sample_names <- clean_sample_names(args$name)
  
  # Validate file existence
  missing_files <- peak_files[!file.exists(peak_files)]
  if (length(missing_files) > 0) {
    stop(paste("Peak files not found:", paste(missing_files, collapse = ", ")))
  }
  
  # Load and process each peak file
  peak_data_list <- list()
  
  for (i in seq_along(peak_files)) {
    file_path <- peak_files[i]
    sample_name <- sample_names[i]
    
    message(sprintf("Loading peak file %d/%d: %s", i, length(peak_files), basename(file_path)))
    
    # Determine file format and load accordingly
    peak_data <- load_single_peak_file(file_path, args$peakformat, args$scorecol)
    
    # Add sample information
    peak_data$SampleID <- sample_name
    
    peak_data_list[[sample_name]] <- peak_data
  }
  
  message(sprintf("Successfully loaded %d peak files", length(peak_data_list)))
  return(peak_data_list)
}

#' Load a single peak file
#'
#' @param file_path Path to the peak file
#' @param peak_format Format of the peak file (csv, bed, narrow, macs)
#' @param score_col Column number containing peak scores
#' @return Data frame with standardized peak information
load_single_peak_file <- function(file_path, peak_format = "csv", score_col = 6) {
  
  if (peak_format == "csv") {
    # Handle CSV format with potential headers and different separators
    delimiter <- check_file_delimiter(file_path)
    
    # Try to read with header first
    peak_data <- tryCatch({
      read.table(file_path, sep = delimiter, header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
      # If header reading fails, try without header
      read.table(file_path, sep = delimiter, header = FALSE, stringsAsFactors = FALSE)
    })
    
    # Standardize column names based on expected format
    peak_data <- standardize_peak_columns(peak_data, peak_format)
    
  } else if (peak_format == "bed") {
    # Standard BED format
    peak_data <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(peak_data)[1:3] <- c("chr", "start", "end")
    
    # Add score column if it exists
    if (ncol(peak_data) >= score_col) {
      colnames(peak_data)[score_col] <- "score"
    } else {
      peak_data$score <- 1  # Default score
    }
    
  } else if (peak_format == "narrow") {
    # narrowPeak format (ENCODE standard)
    peak_data <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    narrow_cols <- c("chr", "start", "end", "name", "score", "strand", 
                     "signalValue", "pValue", "qValue", "peak")
    
    # Assign column names up to the number of columns available
    n_cols <- min(length(narrow_cols), ncol(peak_data))
    colnames(peak_data)[1:n_cols] <- narrow_cols[1:n_cols]
    
  } else if (peak_format == "macs") {
    # MACS format - usually tab-delimited with specific columns
    peak_data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    peak_data <- standardize_peak_columns(peak_data, peak_format)
  }
  
  # Ensure essential columns exist
  essential_cols <- c("chr", "start", "end")
  missing_cols <- essential_cols[!essential_cols %in% colnames(peak_data)]
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing essential columns in peak file %s: %s", 
                 file_path, paste(missing_cols, collapse = ", ")))
  }
  
  # Convert coordinates to numeric
  peak_data$start <- as.numeric(peak_data$start)
  peak_data$end <- as.numeric(peak_data$end)
  
  # Add peak width if not present
  if (!"width" %in% colnames(peak_data)) {
    peak_data$width <- peak_data$end - peak_data$start
  }
  
  # Ensure score column exists
  if (!"score" %in% colnames(peak_data)) {
    if (score_col <= ncol(peak_data)) {
      colnames(peak_data)[score_col] <- "score"
    } else {
      peak_data$score <- 1  # Default score
    }
  }
  
  # Convert score to numeric
  peak_data$score <- as.numeric(peak_data$score)
  
  # Clean chromosome names
  peak_data$chr <- clean_chromosome_names(peak_data$chr)
  
  # Sort by genomic coordinates
  peak_data <- peak_data[order(peak_data$chr, peak_data$start), ]
  
  return(peak_data)
}

#' Standardize peak column names based on format
#'
#' @param peak_data Data frame with peak data
#' @param peak_format Format of the peak file
#' @return Data frame with standardized column names
standardize_peak_columns <- function(peak_data, peak_format) {
  
  # Get current column names
  current_cols <- colnames(peak_data)
  
  # Define expected column patterns for different formats
  if (peak_format == "csv" || peak_format == "macs") {
    # Common patterns in CSV/MACS files
    chr_patterns <- c("chr", "chrom", "chromosome", "seqnames")
    start_patterns <- c("start", "chromstart", "txstart", "pos")
    end_patterns <- c("end", "chromend", "txend", "stop")
    score_patterns <- c("score", "pileup", "fold_enrichment", "signalvalue")
    
    # Find and standardize chromosome column
    chr_idx <- find_column_by_pattern(current_cols, chr_patterns)
    if (length(chr_idx) > 0) {
      colnames(peak_data)[chr_idx[1]] <- "chr"
    }
    
    # Find and standardize start column
    start_idx <- find_column_by_pattern(current_cols, start_patterns)
    if (length(start_idx) > 0) {
      colnames(peak_data)[start_idx[1]] <- "start"
    }
    
    # Find and standardize end column
    end_idx <- find_column_by_pattern(current_cols, end_patterns)
    if (length(end_idx) > 0) {
      colnames(peak_data)[end_idx[1]] <- "end"
    }
    
    # Find and standardize score column
    score_idx <- find_column_by_pattern(current_cols, score_patterns)
    if (length(score_idx) > 0) {
      colnames(peak_data)[score_idx[1]] <- "score"
    }
  }
  
  return(peak_data)
}

#' Find column by pattern matching
#'
#' @param column_names Vector of column names
#' @param patterns Vector of patterns to match
#' @return Index of matching column
find_column_by_pattern <- function(column_names, patterns) {
  matches <- c()
  for (pattern in patterns) {
    idx <- grep(pattern, column_names, ignore.case = TRUE)
    if (length(idx) > 0) {
      matches <- c(matches, idx)
    }
  }
  return(unique(matches))
}

#' Clean chromosome names to standard format
#'
#' @param chr_names Vector of chromosome names
#' @return Vector of cleaned chromosome names
clean_chromosome_names <- function(chr_names) {
  # Remove "chr" prefix if present, then add it back consistently
  clean_names <- gsub("^chr", "", chr_names, ignore.case = TRUE)
  
  # Handle special chromosomes
  clean_names <- gsub("^MT$|^M$", "chrM", clean_names, ignore.case = TRUE)
  clean_names <- gsub("^X$", "chrX", clean_names, ignore.case = TRUE)
  clean_names <- gsub("^Y$", "chrY", clean_names, ignore.case = TRUE)
  
  # Add chr prefix to numeric chromosomes
  numeric_chr <- grepl("^[0-9]+$", clean_names)
  clean_names[numeric_chr] <- paste0("chr", clean_names[numeric_chr])
  
  # Handle already formatted chromosome names
  already_formatted <- grepl("^chr", clean_names, ignore.case = TRUE)
  clean_names[!already_formatted & !grepl("^chr", clean_names)] <- 
    paste0("chr", clean_names[!already_formatted & !grepl("^chr", clean_names)])
  
  return(clean_names)
}

#' Validate BAM files
#'
#' @param bam_files Vector of BAM file paths
#' @return Logical vector indicating which files are valid
validate_bam_files <- function(bam_files) {
  message("Validating BAM files...")
  
  valid_files <- logical(length(bam_files))
  
  for (i in seq_along(bam_files)) {
    bam_file <- bam_files[i]
    
    # Check file existence
    if (!file.exists(bam_file)) {
      warning(sprintf("BAM file not found: %s", bam_file))
      valid_files[i] <- FALSE
      next
    }
    
    # Check if BAM index exists
    index_file <- paste0(bam_file, ".bai")
    alt_index_file <- gsub("\\.bam$", ".bai", bam_file)
    
    if (!file.exists(index_file) && !file.exists(alt_index_file)) {
      warning(sprintf("BAM index not found for: %s", bam_file))
      # Try to create index
      tryCatch({
        Rsamtools::indexBam(bam_file)
        message(sprintf("Created BAM index for: %s", bam_file))
        valid_files[i] <- TRUE
      }, error = function(e) {
        warning(sprintf("Failed to create BAM index for %s: %s", bam_file, e$message))
        valid_files[i] <- FALSE
      })
    } else {
      valid_files[i] <- TRUE
    }
    
    # Basic BAM file validation
    if (valid_files[i]) {
      tryCatch({
        # Try to read header to validate BAM file
        header <- Rsamtools::scanBamHeader(bam_file)
        if (length(header) == 0) {
          warning(sprintf("Invalid BAM file (no header): %s", bam_file))
          valid_files[i] <- FALSE
        }
      }, error = function(e) {
        warning(sprintf("Error reading BAM file %s: %s", bam_file, e$message))
        valid_files[i] <- FALSE
      })
    }
  }
  
  message(sprintf("Validated %d/%d BAM files successfully", 
                  sum(valid_files), length(bam_files)))
  
  return(valid_files)
}

#' Filter peaks based on various criteria
#'
#' @param peak_data Data frame with peak data
#' @param min_width Minimum peak width
#' @param max_width Maximum peak width
#' @param min_score Minimum peak score
#' @return Filtered peak data
filter_peaks <- function(peak_data, min_width = 50, max_width = 10000, min_score = NULL) {
  initial_count <- nrow(peak_data)
  
  # Filter by width
  if (!is.null(min_width)) {
    peak_data <- peak_data[peak_data$width >= min_width, ]
  }
  
  if (!is.null(max_width)) {
    peak_data <- peak_data[peak_data$width <= max_width, ]
  }
  
  # Filter by score
  if (!is.null(min_score) && "score" %in% colnames(peak_data)) {
    peak_data <- peak_data[peak_data$score >= min_score, ]
  }
  
  # Remove invalid coordinates
  peak_data <- peak_data[peak_data$start >= 0 & peak_data$end > peak_data$start, ]
  
  # Remove peaks on unplaced/random contigs if requested
  standard_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
  peak_data <- peak_data[peak_data$chr %in% standard_chromosomes, ]
  
  final_count <- nrow(peak_data)
  
  message(sprintf("Peak filtering: %d -> %d peaks (removed %d)", 
                  initial_count, final_count, initial_count - final_count))
  
  return(peak_data)
}

#' Create genomic ranges from peak data
#'
#' @param peak_data Data frame with peak data
#' @return GRanges object
create_genomic_ranges <- function(peak_data) {
  
  # Create GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = peak_data$chr,
    ranges = IRanges::IRanges(start = peak_data$start, end = peak_data$end),
    score = peak_data$score
  )
  
  # Add metadata columns if present
  metadata_cols <- setdiff(colnames(peak_data), c("chr", "start", "end", "width", "score"))
  if (length(metadata_cols) > 0) {
    for (col in metadata_cols) {
      GenomicRanges::mcols(gr)[[col]] <- peak_data[[col]]
    }
  }
  
  return(gr)
}

#' Save batch correction information for step 2
#'
#' @param dds_lrt DESeq2 object
#' @param args Command line arguments
save_batch_correction_info <- function(dds_lrt, args) {
  if (args$batchcorrection == "limmaremovebatcheffect") {
    message("Saving batch correction information for step 2...")
    
    # Extract information needed for limma batch correction
    batch_info <- list(
      batch_factor = colData(dds_lrt)$Replicate,  # Assuming Replicate is the batch factor
      design_matrix = model.matrix(as.formula(args$design), data = as.data.frame(colData(dds_lrt)))
    )
    
    # Save to RDS file
    batch_file <- paste0(args$output, "_batch_correction_info.rds")
    saveRDS(batch_info, batch_file)
    
    message(sprintf("Batch correction information saved to: %s", batch_file))
  }
} 