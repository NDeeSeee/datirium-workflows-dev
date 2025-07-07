#!/usr/bin/env Rscript

# --- Common CLI Argument Parsing Helpers ---
# Safe, reusable functions for parsing command line arguments
# These can be used by any cli_args.R file to reduce code duplication

#' Parse multi-value arguments (like --input file1 file2 file3)
#' 
#' @param all_args Vector of all command line arguments
#' @param flag_name Flag name (e.g., "input", "name") or single char (e.g., "u", "t")
#' @return Vector of values for this flag
parse_multi_value_args <- function(all_args, flag_name) {
  # Handle both long (--flag) and short (-f) formats
  if (nchar(flag_name) == 1) {
    flag <- paste0("-", flag_name)
  } else {
    flag <- paste0("--", flag_name)
  }
  
  flag_idx <- which(all_args == flag)
  
  if (length(flag_idx) == 0) {
    return(character(0))
  }
  
  start_idx <- flag_idx[1] + 1
  if (start_idx > length(all_args)) {
    return(character(0))
  }
  
  # Find end of values (next flag or end of args)
  end_idx <- start_idx
  while (end_idx <= length(all_args) && !startsWith(all_args[end_idx], "--") && !startsWith(all_args[end_idx], "-")) {
    end_idx <- end_idx + 1
  }
  
  if (start_idx < end_idx) {
    return(all_args[start_idx:(end_idx - 1)])
  } else {
    return(character(0))
  }
}

#' Parse single-value arguments (like --meta file.csv)
#' 
#' @param all_args Vector of all command line arguments
#' @param flag_name Flag name (e.g., "meta", "design") or single char (e.g., "o", "p")
#' @param default_value Default value if flag not found
#' @return Single value for this flag
parse_single_value_arg <- function(all_args, flag_name, default_value = NULL) {
  # Handle both long (--flag) and short (-f) formats
  if (nchar(flag_name) == 1) {
    flag <- paste0("-", flag_name)
  } else {
    flag <- paste0("--", flag_name)
  }
  
  flag_idx <- which(all_args == flag)
  
  if (length(flag_idx) == 0) {
    return(default_value)
  }
  
  if (flag_idx[1] < length(all_args) && !startsWith(all_args[flag_idx[1] + 1], "--")) {
    return(all_args[flag_idx[1] + 1])
  }
  
  return(default_value)
}

#' Parse boolean flags (handles both --flag and --flag TRUE/FALSE)
#' 
#' @param all_args Vector of all command line arguments
#' @param flag_names Vector of boolean flag names
#' @return Named list of boolean values
parse_boolean_flags <- function(all_args, flag_names) {
  result <- list()
  
  for (flag_name in flag_names) {
    flag <- paste0("--", flag_name)
    flag_idx <- which(all_args == flag)
    
    if (length(flag_idx) == 0) {
      result[[flag_name]] <- FALSE
    } else {
      # Check if there's a value after the flag
      if (flag_idx[1] < length(all_args) && !startsWith(all_args[flag_idx[1] + 1], "--")) {
        # Has a value (TRUE/FALSE)
        val <- all_args[flag_idx[1] + 1]
        result[[flag_name]] <- toupper(val) == "TRUE"
      } else {
        # Just the flag present (store_true action)
        result[[flag_name]] <- TRUE
      }
    }
  }
  
  return(result)
}

#' Parse numeric arguments with type conversion
#' 
#' @param all_args Vector of all command line arguments
#' @param numeric_args Named list: arg_name -> "integer" or "double"
#' @param defaults Named list of default values
#' @return Named list of numeric values
parse_numeric_args <- function(all_args, numeric_args, defaults = list()) {
  result <- list()
  
  for (arg_name in names(numeric_args)) {
    flag <- paste0("--", arg_name)
    flag_idx <- which(all_args == flag)
    
    if (length(flag_idx) == 0 || flag_idx[1] >= length(all_args)) {
      result[[arg_name]] <- defaults[[arg_name]]
    } else {
      val <- all_args[flag_idx[1] + 1]
      
      if (numeric_args[[arg_name]] == "integer") {
        result[[arg_name]] <- as.integer(val)
      } else {
        result[[arg_name]] <- as.numeric(val)
      }
    }
  }
  
  return(result)
}

#' Validate file existence for a list of files
#' 
#' @param files Vector of file paths
#' @param file_type Description of file type for error messages
#' @return Vector of error messages (empty if all valid)
validate_file_existence <- function(files, file_type = "file") {
  errors <- character(0)
  
  for (file in files) {
    if (!file.exists(file)) {
      errors <- c(errors, paste(file_type, "does not exist:", file))
    }
  }
  
  return(errors)
}

#' Safe sourcing with fallback error handling
#' 
#' @param filepath Path to file to source
#' @return TRUE if successful, FALSE otherwise
safe_source <- function(filepath) {
  tryCatch({
    if (file.exists(filepath)) {
      source(filepath)
      return(TRUE)
    } else {
      message(paste("Warning: Helper file not found:", filepath))
      return(FALSE)
    }
  }, error = function(e) {
    message(paste("Warning: Failed to source helper file:", filepath, "-", e$message))
    return(FALSE)
  })
}

#' Source CLI helpers with multiple fallback paths
#' 
#' Tries multiple common paths for the cli_helpers.R file
#' @return TRUE if helpers were loaded, FALSE otherwise
source_cli_helpers <- function() {
  # Try multiple paths in order of preference
  paths_to_try <- c(
    # Relative path from current directory
    file.path(dirname(getwd()), "common", "cli_helpers.R"),
    "../common/cli_helpers.R",
    "common/cli_helpers.R",
    # Docker paths
    "/usr/local/bin/functions/common/cli_helpers.R",
    "/usr/local/bin/common/cli_helpers.R",
    # Absolute paths (if running from functions subdirectory)
    file.path(dirname(dirname(getwd())), "common", "cli_helpers.R")
  )
  
  for (path in paths_to_try) {
    if (safe_source(path)) {
      message(paste("Successfully loaded CLI helpers from:", path))
      return(TRUE)
    }
  }
  
  message("CLI helpers not available, using manual parsing fallbacks")
  return(FALSE)
}

# Export all functions for easy access
cli_helpers <- new.env()
cli_helpers$parse_multi_value_args <- parse_multi_value_args
cli_helpers$parse_single_value_arg <- parse_single_value_arg
cli_helpers$parse_boolean_flags <- parse_boolean_flags
cli_helpers$parse_numeric_args <- parse_numeric_args
cli_helpers$validate_file_existence <- validate_file_existence
cli_helpers$safe_source <- safe_source

# Make available globally
assign("cli_helpers", cli_helpers, envir = .GlobalEnv)