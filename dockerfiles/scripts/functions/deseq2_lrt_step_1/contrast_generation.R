#!/usr/bin/env Rscript

# Load required packages
library(DESeq2)
library(dplyr)
library(purrr)

# Memory efficient function to generate main effect contrasts
generate_main_effect_contrasts <- function(dds, factors, factor_levels, args) {
  start_time <- proc.time()
  contrasts <- list()
  
  cat("Generating main effect contrasts (memory efficient version)...\n")
  
  # Get the model matrix and coefficient names
  model_matrix <- model.matrix(design(dds), data = colData(dds))
  coef_names <- resultsNames(dds)
  
  cat("Available coefficients:", paste(coef_names, collapse=", "), "\n")
  
  # For each factor, generate contrasts by comparing levels
  for (factor in factors) {
    cat(paste("Processing factor:", factor, "\n"))
    levels <- factor_levels[[factor]]
    
    if (length(levels) < 2) {
      cat(paste("Skipping factor", factor, "- less than 2 levels\n"))
      next
    }
    
    # For binary factors, create simple contrasts
    if (length(levels) == 2) {
      # Find the coefficient name for this factor (exclude interaction terms)
      factor_coef <- grep(paste0("^", factor), coef_names, value = TRUE)
      # Exclude interaction terms (those with dots)
      factor_coef <- factor_coef[!grepl("\\.", factor_coef)]
      
      if (length(factor_coef) > 0) {
        for (coef in factor_coef) {
          # Extract results for this coefficient
          contrast_res <- results(dds, name = coef, alpha = args$fdr)
          significant_genes <- sum(contrast_res$padj < args$fdr & 
                                 abs(contrast_res$log2FoldChange) > args$lfcthreshold, na.rm = TRUE)
          
          contrasts <- append(contrasts, list(list(
            effect_type = "main",
            specificity_group = factor,
            numerator = levels[2],
            denominator = levels[1], 
            contrast = coef,
            contrast_res = contrast_res,
            significant_genes = significant_genes
          )))
          
          # Force garbage collection
          gc(verbose = FALSE)
        }
      }
    } else {
      # Multi-level factor - create pairwise contrasts
      ref_level <- levels[1]
      for (i in 2:length(levels)) {
        level <- levels[i]
        
        # Try to find coefficient for this comparison (exclude interaction terms)
        coef_pattern <- paste0(factor, level)
        factor_coef <- grep(coef_pattern, coef_names, value = TRUE)
        # Exclude interaction terms (those with dots)
        factor_coef <- factor_coef[!grepl("\\.", factor_coef)]
        
        if (length(factor_coef) > 0) {
          for (coef in factor_coef) {
            contrast_res <- results(dds, name = coef, alpha = args$fdr)
            significant_genes <- sum(contrast_res$padj < args$fdr & 
                                   abs(contrast_res$log2FoldChange) > args$lfcthreshold, na.rm = TRUE)
            
            contrasts <- append(contrasts, list(list(
              effect_type = "main",
              specificity_group = factor,
              numerator = level,
              denominator = ref_level,
              contrast = paste(factor, level, "vs", ref_level, sep = "_"),
              contrast_res = contrast_res,
              significant_genes = significant_genes
            )))
            
            # Force garbage collection
            gc(verbose = FALSE)
          }
        }
      }
    }
  }
  
  end_time <- proc.time() - start_time
  cat("Main effect contrasts generation completed in", round(end_time["elapsed"], 2), "seconds\n")
  
  return(contrasts)
}

# Memory efficient function to generate interaction effect contrasts
generate_interaction_effect_contrasts <- function(dds, args) {
  start_time <- proc.time()
  contrasts <- list()
  
  cat("Generating interaction effect contrasts (memory efficient version)...\n")
  
  # Get interaction terms from coefficient names
  coef_names <- resultsNames(dds)
  # Look for interaction terms - DESeq2 uses dots (.) for interactions, not colons (:)
  interaction_terms <- grep("\\.", coef_names, value = TRUE)
  # Filter out Intercept which also contains no dots
  interaction_terms <- interaction_terms[interaction_terms != "Intercept"]
  
  cat("Found interaction terms:", paste(interaction_terms, collapse=", "), "\n")
  
  if (length(interaction_terms) == 0) {
    cat("No interaction terms found\n")
    return(contrasts)
  }
  
  # Process each interaction term
  for (interaction in interaction_terms) {
    cat(paste("Processing interaction:", interaction, "\n"))
    
    # Extract results for this interaction
    contrast_res <- results(dds, name = interaction, alpha = args$fdr)
    significant_genes <- sum(contrast_res$padj < args$fdr & 
                           abs(contrast_res$log2FoldChange) > args$lfcthreshold, na.rm = TRUE)
    
    # Parse interaction term to get factors and levels
    # For DESeq2 interaction terms like "treatmentKO.condRest"
    parts <- strsplit(interaction, "\\.")[[1]]
    if (length(parts) >= 2) {
      # Extract factor and level information from DESeq2 coefficient names
      # Format is typically: factor1Level1.factor2Level2
      factor1_part <- parts[1]
      factor2_part <- parts[2]
      
      # For DESeq2 coefficients, extract factor and level names
      # treatmentKO -> factor1="treatment", level1="KO"
      # condRest -> factor2="cond", level2="Rest"
      
      # Extract factor1 and level1
      if (grepl("treatment", factor1_part, ignore.case = TRUE)) {
        factor1 <- "treatment"
        level1 <- gsub("treatment", "", factor1_part, ignore.case = TRUE)
      } else if (grepl("cond", factor1_part, ignore.case = TRUE)) {
        factor1 <- "cond"
        level1 <- gsub("cond", "", factor1_part, ignore.case = TRUE)
      } else {
        factor1 <- "factor1"
        level1 <- factor1_part
      }
      
      # Extract factor2 and level2
      if (grepl("treatment", factor2_part, ignore.case = TRUE)) {
        factor2 <- "treatment"
        level2 <- gsub("treatment", "", factor2_part, ignore.case = TRUE)
      } else if (grepl("cond", factor2_part, ignore.case = TRUE)) {
        factor2 <- "cond"
        level2 <- gsub("cond", "", factor2_part, ignore.case = TRUE)
      } else {
        factor2 <- "factor2"
        level2 <- factor2_part
      }
      
      specificity_group <- paste(factor2, level2, sep = "_")
      numerator <- paste0(factor1, level1)
      denominator <- "baseline"
      
      contrasts <- append(contrasts, list(list(
        effect_type = "interaction",
        specificity_group = specificity_group,
        numerator = numerator,
        denominator = denominator,
        contrast = interaction,
        contrast_res = contrast_res,
        significant_genes = significant_genes
      )))
    }
    
    # Force garbage collection
    gc(verbose = FALSE)
  }
  
  end_time <- proc.time() - start_time
  cat("Interaction effect contrasts generation completed in", round(end_time["elapsed"], 2), "seconds\n")
  
  return(contrasts)
}

# Main function to generate the dataframe with all possible contrasts (memory efficient)
generate_contrasts <- function(dds, args, expression_data_df) {
  start_time <- proc.time()
  
  cat("Starting memory-efficient contrasts generation...\n")
  
  # Check available memory
  if (requireNamespace("pryr", quietly = TRUE)) {
    cat("Available memory:", round(pryr::mem_used() / 1024^3, 2), "GB\n")
  }
  
  # Extract the design formula and model matrix
  design_formula <- design(dds)
  cat("Design formula:", deparse(design_formula), "\n")
  
  # Get the levels of each factor in the design
  factors <- all.vars(design_formula)
  factors <- factors[!grepl("batch", factors)]
  cat("Analyzing factors:", paste(factors, collapse=", "), "\n")
  
  factor_levels <- lapply(factors, function(f) levels(colData(dds)[[f]]))
  names(factor_levels) <- factors
  cat("Factor levels:\n")
  for (f in names(factor_levels)) {
    cat(paste("  ", f, ":", paste(factor_levels[[f]], collapse=", "), "\n"))
  }
  
  # Set lfcthreshold from args (make sure it's available globally for helper functions)
  lfcthreshold <- args$lfcthreshold
  assign("lfcthreshold", lfcthreshold, envir = .GlobalEnv)
  
  # Generate contrasts more efficiently
  main_contrasts <- generate_main_effect_contrasts(dds, factors, factor_levels, args)
  cat("Generated", length(main_contrasts), "main effect contrasts\n")
  
  interaction_contrasts <- generate_interaction_effect_contrasts(dds, args)
  cat("Generated", length(interaction_contrasts), "interaction effect contrasts\n")
  
  all_contrasts <- c(main_contrasts, interaction_contrasts)
  cat("Total contrasts generated:", length(all_contrasts), "\n")
  
  # Create minimal results object for storage (memory efficient)
  results_for_storage <- list(
    dds = dds,
    expression_data_df = expression_data_df,
    n_contrasts = length(all_contrasts),
    design_formula = design_formula
  )
  
  # Save contrasts list to RDS (but don't store the full results in memory)
  contrasts_file <- paste0(args$output, "_contrasts.rds")
  cat("Saving contrasts to:", contrasts_file, "\n")
  saveRDS(results_for_storage, file = contrasts_file)
  
  # Create the summary dataframe (memory efficient)
  cat("Creating contrasts summary table...\n")
  contrast_df <- data.frame(
    effect = character(length(all_contrasts)),
    specificity_group = character(length(all_contrasts)),
    contrast = character(length(all_contrasts)),
    numerator = character(length(all_contrasts)),
    denominator = character(length(all_contrasts)),
    significant_genes = integer(length(all_contrasts)),
    stringsAsFactors = FALSE
  )
  
  # Fill the dataframe efficiently
  for (i in seq_along(all_contrasts)) {
    contrast <- all_contrasts[[i]]
    contrast_df[i, ] <- list(
      effect = contrast$effect_type,
      specificity_group = contrast$specificity_group,
      contrast = contrast$contrast,
      numerator = contrast$numerator,
      denominator = contrast$denominator,
      significant_genes = contrast$significant_genes
    )
    
    # Clear the contrast results from memory
    all_contrasts[[i]]$contrast_res <- NULL
  }
  
  # Remove duplicate contrasts and sort
  contrast_df <- contrast_df %>%
    distinct(contrast, .keep_all = TRUE) %>%
    arrange(desc(significant_genes)) %>%
    mutate(contrast_number = row_number()) %>%
    select(contrast_number, everything())
  
  # Final cleanup
  rm(all_contrasts, main_contrasts, interaction_contrasts)
  gc(verbose = FALSE)
  
  end_time <- proc.time() - start_time
  cat("Total contrasts generation completed in", round(end_time["elapsed"], 2), "seconds\n")
  
  if (requireNamespace("pryr", quietly = TRUE)) {
    cat("Memory used after generation:", round(pryr::mem_used() / 1024^3, 2), "GB\n")
  }
  
  return(contrast_df)
} 