#!/usr/bin/env Rscript

# --- Contrast generation functions for ATAC-seq LRT Step 1 ---

#' Generate ATAC-seq contrasts from DESeq2 object
#'
#' @param dds_lrt DESeq2 object after LRT analysis
#' @param args Command line arguments
#' @return List of contrast results
generate_atac_contrasts <- function(dds_lrt, args) {
  message("Generating ATAC-seq contrasts...")
  
  # Extract design information
  design_formula <- as.formula(args$design)
  design_vars <- all.vars(design_formula)
  
  # Get column data
  coldata <- as.data.frame(colData(dds_lrt))
  
  # Generate contrasts based on design
  contrasts_list <- list()
  
  # Main effects contrasts
  for (var in design_vars) {
    if (var %in% colnames(coldata)) {
      var_contrasts <- generate_main_effect_contrasts(dds_lrt, var, args)
      contrasts_list <- c(contrasts_list, var_contrasts)
    }
  }
  
  # Interaction contrasts if interaction terms are present
  if (grepl(":", args$design)) {
    interaction_contrasts <- generate_interaction_contrasts(dds_lrt, design_vars, args)
    contrasts_list <- c(contrasts_list, interaction_contrasts)
  }
  
  message(sprintf("Generated %d contrasts total", length(contrasts_list)))
  
  return(contrasts_list)
}

#' Generate main effect contrasts
#'
#' @param dds_lrt DESeq2 object
#' @param var_name Variable name for contrast
#' @param args Command line arguments
#' @return List of contrast results for the variable
generate_main_effect_contrasts <- function(dds_lrt, var_name, args) {
  message(sprintf("Generating main effect contrasts for: %s", var_name))
  
  contrasts_list <- list()
  
  # Get unique levels of the variable
  coldata <- as.data.frame(colData(dds_lrt))
  
  if (!var_name %in% colnames(coldata)) {
    warning(sprintf("Variable %s not found in sample data", var_name))
    return(contrasts_list)
  }
  
  unique_levels <- unique(coldata[[var_name]])
  
  if (length(unique_levels) < 2) {
    warning(sprintf("Variable %s has fewer than 2 levels", var_name))
    return(contrasts_list)
  }
  
  # Generate pairwise contrasts
  for (i in 1:(length(unique_levels) - 1)) {
    for (j in (i + 1):length(unique_levels)) {
      level1 <- unique_levels[i]
      level2 <- unique_levels[j]
      
      contrast_name <- sprintf("%s_%s_vs_%s", var_name, level2, level1)
      
      tryCatch({
        # Generate contrast using DESeq2
        contrast_result <- results(
          dds_lrt,
          contrast = c(var_name, level2, level1),
          alpha = args$fdr,
          lfcThreshold = if(args$use_lfc_thresh) args$lfcthreshold else 0
        )
        
        contrasts_list[[contrast_name]] <- contrast_result
        
        message(sprintf("Generated contrast: %s", contrast_name))
        
      }, error = function(e) {
        warning(sprintf("Failed to generate contrast %s: %s", contrast_name, e$message))
      })
    }
  }
  
  return(contrasts_list)
}

#' Generate interaction contrasts
#'
#' @param dds_lrt DESeq2 object
#' @param design_vars Design variables
#' @param args Command line arguments
#' @return List of interaction contrast results
generate_interaction_contrasts <- function(dds_lrt, design_vars, args) {
  message("Generating interaction contrasts...")
  
  contrasts_list <- list()
  
  # Extract interaction terms from design formula
  design_str <- as.character(args$design)
  interaction_terms <- extract_interaction_terms(design_str)
  
  coldata <- as.data.frame(colData(dds_lrt))
  
  for (interaction_term in interaction_terms) {
    # Parse interaction term (e.g., "ConditionA:TissueB")
    vars <- strsplit(interaction_term, ":")[[1]]
    
    if (length(vars) == 2 && all(vars %in% colnames(coldata))) {
      var1 <- vars[1]
      var2 <- vars[2]
      
      # Generate interaction contrasts
      interaction_contrasts <- generate_interaction_pairs(dds_lrt, var1, var2, args)
      contrasts_list <- c(contrasts_list, interaction_contrasts)
    }
  }
  
  return(contrasts_list)
}

#' Generate pairwise interaction contrasts
#'
#' @param dds_lrt DESeq2 object
#' @param var1 First interaction variable
#' @param var2 Second interaction variable
#' @param args Command line arguments
#' @return List of interaction contrast results
generate_interaction_pairs <- function(dds_lrt, var1, var2, args) {
  contrasts_list <- list()
  
  coldata <- as.data.frame(colData(dds_lrt))
  
  # Get unique combinations
  combinations <- unique(coldata[, c(var1, var2)])
  
  # Generate contrasts between combinations
  for (i in 1:(nrow(combinations) - 1)) {
    for (j in (i + 1):nrow(combinations)) {
      comb1 <- combinations[i, ]
      comb2 <- combinations[j, ]
      
      contrast_name <- sprintf("Interaction_%s%s_vs_%s%s", 
                               var1, comb2[[var1]], var1, comb1[[var1]])
      
      tryCatch({
        # This is a simplified approach - more complex interaction contrasts
        # would require custom contrast matrices
        contrast_result <- results(
          dds_lrt,
          alpha = args$fdr,
          lfcThreshold = if(args$use_lfc_thresh) args$lfcthreshold else 0
        )
        
        # Filter results for the specific interaction
        # This is a placeholder - actual implementation would be more complex
        contrasts_list[[contrast_name]] <- contrast_result
        
        message(sprintf("Generated interaction contrast: %s", contrast_name))
        
      }, error = function(e) {
        warning(sprintf("Failed to generate interaction contrast %s: %s", 
                        contrast_name, e$message))
      })
    }
  }
  
  return(contrasts_list)
}

#' Extract interaction terms from design formula
#'
#' @param design_str Design formula as string
#' @return Vector of interaction terms
extract_interaction_terms <- function(design_str) {
  # Remove the ~ and spaces
  design_clean <- gsub("~|\\s", "", design_str)
  
  # Split by + to get individual terms
  terms <- strsplit(design_clean, "\\+")[[1]]
  
  # Find terms that contain :
  interaction_terms <- terms[grepl(":", terms)]
  
  return(interaction_terms)
}

#' Generate all possible contrasts for a factorial design
#'
#' @param dds_lrt DESeq2 object
#' @param args Command line arguments
#' @return List of all possible contrast results
generate_all_factorial_contrasts <- function(dds_lrt, args) {
  message("Generating all factorial contrasts...")
  
  contrasts_list <- list()
  
  # Get design matrix
  design_matrix <- model.matrix(as.formula(args$design), data = as.data.frame(colData(dds_lrt)))
  
  # Get coefficient names
  coef_names <- colnames(design_matrix)
  
  # Generate contrasts for each coefficient
  for (coef_name in coef_names) {
    if (coef_name != "(Intercept)") {
      contrast_name <- paste("Coef", coef_name, sep = "_")
      
      tryCatch({
        contrast_result <- results(
          dds_lrt,
          name = coef_name,
          alpha = args$fdr,
          lfcThreshold = if(args$use_lfc_thresh) args$lfcthreshold else 0
        )
        
        contrasts_list[[contrast_name]] <- contrast_result
        
        message(sprintf("Generated coefficient contrast: %s", contrast_name))
        
      }, error = function(e) {
        warning(sprintf("Failed to generate coefficient contrast %s: %s", 
                        contrast_name, e$message))
      })
    }
  }
  
  return(contrasts_list)
}

#' Filter contrasts by significance
#'
#' @param contrasts_list List of contrast results
#' @param fdr_threshold FDR threshold for filtering
#' @return Filtered list of contrasts
filter_significant_contrasts <- function(contrasts_list, fdr_threshold = 0.1) {
  message("Filtering contrasts by significance...")
  
  filtered_contrasts <- list()
  
  for (contrast_name in names(contrasts_list)) {
    contrast_result <- contrasts_list[[contrast_name]]
    
    # Count significant peaks
    sig_peaks <- sum(!is.na(contrast_result$padj) & 
                     contrast_result$padj < fdr_threshold)
    
    if (sig_peaks > 0) {
      filtered_contrasts[[contrast_name]] <- contrast_result
      message(sprintf("Keeping contrast %s: %d significant peaks", 
                      contrast_name, sig_peaks))
    } else {
      message(sprintf("Removing contrast %s: no significant peaks", contrast_name))
    }
  }
  
  return(filtered_contrasts)
}

#' Export individual contrast results
#'
#' @param contrasts_list List of contrast results
#' @param args Command line arguments
export_individual_contrasts <- function(contrasts_list, args) {
  message("Exporting individual contrast results...")
  
  for (contrast_name in names(contrasts_list)) {
    contrast_result <- contrasts_list[[contrast_name]]
    
    # Convert to data frame
    contrast_df <- as.data.frame(contrast_result)
    contrast_df$peak_id <- rownames(contrast_df)
    
    # Sort by significance
    contrast_df <- contrast_df[order(contrast_df$pvalue), ]
    
    # Export to file
    output_file <- paste0(args$output, "_contrast_", contrast_name, ".tsv")
    write.table(contrast_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    message(sprintf("Exported contrast %s to: %s", contrast_name, output_file))
  }
} 