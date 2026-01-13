
# Let's first check what columns actually exist in data
cat("Actual column names in your data:\n")
print(colnames(dm_data))

# Based on data structure, create a corrected version
prepare_area_data_corrected <- function(data, area_type) {
  
  cat("Preparing data for", area_type, "area...\n")
  
  # Use the exact column names from dataset
  id_col <- "ID.GMP"
  time_col <- "Time"
  indicator_col <- "Indicator"
  
  if (area_type == "segregated") {
    obs_col <- "OBS.cases.in.segregated.area"
    target_col <- "Target.group.in.segregated.area"
  } else if (area_type == "complementary") {
    obs_col <- "OBS.cases.in.complementary.area"
    target_col <- "Target.group.in.complementary.area"
  }
  
  cat("Looking for columns:\n")
  cat("  ID:", id_col, "\n")
  cat("  Time:", time_col, "\n")
  cat("  Indicator:", indicator_col, "\n")
  cat("  OBS:", obs_col, "\n")
  cat("  Target:", target_col, "\n")
  
  # Check if columns exist
  all_cols <- colnames(data)
  cat("Available columns:", paste(all_cols, collapse = ", "), "\n")
  
  required_cols <- c(id_col, time_col, indicator_col, obs_col, target_col)
  missing_cols <- required_cols[!required_cols %in% all_cols]
  
  if (length(missing_cols) > 0) {
    cat("Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Select the columns
  area_data <- data[, required_cols, drop = FALSE]
  colnames(area_data) <- c("ID.GMP", "Time", "Indicator", "OBS", "Target")
  
  # Convert to numeric
  area_data$OBS <- as.numeric(area_data$OBS)
  area_data$Target <- as.numeric(area_data$Target)
  
  # Remove NA rows
  area_data <- area_data[complete.cases(area_data[, c("OBS", "Target")]), ]
  
  # Calculate expected values
  area_data <- area_data %>%
    group_by(Indicator) %>%
    mutate(EXP = mean(OBS, na.rm = TRUE)) %>%
    ungroup()
  
  cat("Successfully prepared", nrow(area_data), "rows\n")
  return(area_data)
}

# Test with the corrected function
cat("Testing with corrected function...\n")
segregated_test <- prepare_area_data_corrected(dm_data, "segregated")
complementary_test <- prepare_area_data_corrected(dm_data, "complementary")

if (!is.null(segregated_test)) {
  cat("Segregated test successful:", nrow(segregated_test), "rows\n")
}
if (!is.null(complementary_test)) {
  cat("Complementary test successful:", nrow(complementary_test), "rows\n")
}

# Load required libraries
library(tidyverse)
library(ggplot2)
library(factoextra)
library(flextable)
library(officer)
library(robustbase)
library(MASS)
library(openxlsx)
library(ggpubr)
library(corrplot)
library(reshape2)

dm_data <- DM_Related_to_work_2

# =============================================================================
# DEBUGGING: CHECK DATA STRUCTURE FIRST
# =============================================================================

cat("=== DEBUGGING DATA STRUCTURE ===\n")
cat("Column names in dm_data:\n")
print(colnames(dm_data))

cat("\nData dimensions:", dim(dm_data), "\n")
cat("Data class:", class(dm_data), "\n")

# Check if it's a tibble and convert to data frame if needed
if (inherits(dm_data, "tbl_df") || inherits(dm_data, "tbl")) {
  dm_data <- as.data.frame(dm_data)
  cat("Converted to data frame\n")
}

# =============================================================================
# FIXED DATA PREPARATION FUNCTION - USING EXACT COLUMN NAMES FROM YOUR DATA
# =============================================================================

prepare_area_data_fixed <- function(data, area_type) {
  
  cat("Preparing data for", area_type, "area...\n")
  
  # Use exact column names from your dataset
  id_col <- "ID GMP"
  time_col <- "Time"
  indicator_col <- "Indicator"
  
  if (area_type == "segregated") {
    obs_col <- "OBS cases in segregated area"
    target_col <- "Target group in segregated area"
  } else if (area_type == "complementary") {
    obs_col <- "OBS cases in complementary area"
    target_col <- "Target group in complementary area"
  }
  
  cat("Using columns:\n")
  cat("  ID:", id_col, "\n")
  cat("  Time:", time_col, "\n")
  cat("  Indicator:", indicator_col, "\n")
  cat("  OBS:", obs_col, "\n")
  cat("  Target:", target_col, "\n")
  
  # Check if we have the necessary columns
  required_cols <- c(id_col, time_col, indicator_col, obs_col, target_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns for", area_type, "area:", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Select and rename columns
  area_data <- data[, required_cols]
  colnames(area_data) <- c("ID.GMP", "Time", "Indicator", "OBS", "Target")
  
  # Convert to numeric if needed
  area_data$OBS <- as.numeric(area_data$OBS)
  area_data$Target <- as.numeric(area_data$Target)
  
  # Remove rows where OBS or Target is NA
  area_data <- area_data[!is.na(area_data$OBS) & !is.na(area_data$Target), ]
  
  # Calculate expected as mean of OBS for each indicator
  area_data <- area_data %>%
    group_by(Indicator) %>%
    mutate(EXP = mean(OBS, na.rm = TRUE)) %>%
    ungroup()
  
  cat("Prepared", nrow(area_data), "rows for", area_type, "area\n")
  cat("Unique indicators:", length(unique(area_data$Indicator)), "\n")
  
  return(area_data)
}

# =============================================================================
# ENHANCED EMPIRICAL BAYES SHRINKAGE FUNCTIONS
# =============================================================================

empirical_bayes_shrinkage <- function(observed, expected, method = "beta-binomial") {
  if (all(observed == 0) | all(expected == 0)) {
    return(list(
      sr_raw = observed / ifelse(expected == 0, 1, expected),
      sr_shrunk = observed / ifelse(expected == 0, 1, expected),
      shrinkage_factor = 0,
      reliability = 0
    ))
  }
  
  # Calculate raw rates
  raw_rates <- observed / ifelse(expected == 0, 1, expected)
  
  # Estimate prior parameters using method of moments
  n <- length(observed)
  if (n < 2) {
    return(list(
      sr_raw = raw_rates,
      sr_shrunk = raw_rates,
      shrinkage_factor = 0,
      reliability = 0
    ))
  }
  
  # Remove extreme values for stable estimation
  finite_idx <- is.finite(raw_rates) & !is.na(raw_rates)
  if (sum(finite_idx) < 2) {
    return(list(
      sr_raw = raw_rates,
      sr_shrunk = raw_rates,
      shrinkage_factor = 0,
      reliability = 0
    ))
  }
  
  rates_finite <- raw_rates[finite_idx]
  
  # Winsorize extreme values for robust estimation
  q_low <- quantile(rates_finite, 0.05, na.rm = TRUE)
  q_high <- quantile(rates_finite, 0.95, na.rm = TRUE)
  rates_winsor <- pmin(pmax(rates_finite, q_low), q_high)
  
  # Method of moments for Beta distribution
  m <- mean(rates_winsor, na.rm = TRUE)
  v <- var(rates_winsor, na.rm = TRUE)
  
  if (v == 0 || is.na(v)) {
    alpha_prior <- 1
    beta_prior <- 1
  } else {
    alpha_prior <- m * (m * (1 - m) / v - 1)
    beta_prior <- (1 - m) * (m * (1 - m) / v - 1)
    
    # Ensure parameters are positive
    alpha_prior <- max(alpha_prior, 0.1)
    beta_prior <- max(beta_prior, 0.1)
  }
  
  # Apply shrinkage
  sr_shrunk <- raw_rates
  reliability <- numeric(length(observed))
  
  for (i in 1:length(observed)) {
    if (finite_idx[i] && expected[i] > 0) {
      alpha_post <- alpha_prior + observed[i]
      beta_post <- beta_prior + expected[i] - observed[i]
      sr_shrunk[i] <- alpha_post / (alpha_post + beta_post)
      
      # Reliability based on effective sample size
      reliability[i] <- expected[i] / (expected[i] + alpha_prior + beta_prior)
    }
  }
  
  shrinkage_factor <- 1 - (1 / (1 + alpha_prior + beta_prior))
  
  return(list(
    sr_raw = raw_rates,
    sr_shrunk = sr_shrunk,
    shrinkage_factor = shrinkage_factor,
    reliability = reliability,
    prior_alpha = alpha_prior,
    prior_beta = beta_prior
  ))
}

# =============================================================================
# COMPREHENSIVE PCA AND COMPOSITE INDICATOR FUNCTIONS
# =============================================================================

calculate_pca_components <- function(data_matrix) {
  # Remove constant columns
  data_matrix <- data_matrix[, apply(data_matrix, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  
  if (ncol(data_matrix) < 2) {
    stop("Insufficient columns for PCA after removing constant columns")
  }
  
  # Handle missing values
  data_matrix[is.na(data_matrix)] <- median(as.matrix(data_matrix), na.rm = TRUE)
  
  # Perform PCA
  pca_result <- prcomp(data_matrix, scale. = TRUE, center = TRUE)
  
  # Extract components
  variance_explained <- summary(pca_result)$importance[2, ] * 100
  loadings <- pca_result$rotation
  scores <- pca_result$x
  
  return(list(
    pca_result = pca_result,
    variance_explained = variance_explained,
    loadings = loadings,
    scores = scores
  ))
}

calculate_composite_indicator <- function(data_matrix, weights) {
  # Standardize the data
  data_standardized <- scale(data_matrix)
  
  # Calculate composite scores
  composite_scores <- data_standardized %*% weights
  
  return(as.numeric(composite_scores))
}

# =============================================================================
# FIXED RESHAPE FUNCTION AND DATA PREPARATION
# =============================================================================

reshape_for_analysis <- function(data, value_column, id_column = "ID.GMP") {
  if (is.null(data)) return(NULL)
  
  cat("Reshaping", value_column, "with", nrow(data), "rows\n")
  
  # Select only the necessary columns using base R
  data_selected <- data[, c(id_column, "Indicator", value_column)]
  
  # Remove any rows with missing ID or Indicator
  data_selected <- data_selected[complete.cases(data_selected[, c(id_column, "Indicator")]), ]
  
  # Check for duplicates
  duplicates <- data_selected %>% 
    group_by(!!sym(id_column), Indicator) %>% 
    filter(n() > 1)
  
  if (nrow(duplicates) > 0) {
    cat("Warning: Found", nrow(duplicates), "duplicate rows. Taking mean...\n")
    data_selected <- data_selected %>%
      group_by(!!sym(id_column), Indicator) %>%
      summarize(!!sym(value_column) := mean(!!sym(value_column), na.rm = TRUE),
                .groups = 'drop')
  }
  
  # Pivot to wide format
  data_wide <- data_selected %>%
    pivot_wider(
      names_from = "Indicator", 
      values_from = value_column,
      values_fill = NA
    )
  
  # Set row names
  data_wide <- as.data.frame(data_wide)
  rownames(data_wide) <- data_wide[[id_column]]
  data_wide[[id_column]] <- NULL
  
  cat("Created wide data with", nrow(data_wide), "rows and", ncol(data_wide), "columns\n")
  
  return(data_wide)
}

# =============================================================================
# INDICATOR NAME MAPPING - SIMPLIFIED NAMES
# =============================================================================

# Define shorter indicator names
indicator_mapping <- c(
  "FV30a - Prescription - Alimentary tract and metabolism (ATC A)" = "Prescription",
  "FV31a - Drug substitution - Alimentary tract and metabolism (ATC A)" = "Drug_Substitution",
  "FV32a - Drug substitution rate - Alimentary tract and metabolism (ATC A)" = "Substitution_Rate",
  "FV06 - Proportion of diabetic and/or hypertensive patients who underwent blood lipid testing" = "Lipid_Test",
  "FV08 - Proportion of diabetics who underwent hemoglobin A1c testing" = "HbA1c",
  "<br>FV09 - Proportion of diabetics who underwent ophthalmological examination" = "Eye_Exam",
  "FV13 - Proportion of people aged 40-54 who are switching from medication for diabetes" = "Med_Switch_40_54",
  "FV14 - Proportion of people aged 55-69 who are switching from medication for diabetes" = "Med_Switch_55_69",
  "FV15 - Proportion of diabetic patients who underwent serum creatinine level determination" = "Creatinine",
  "FV16a - Proportion of people under 65 years of age receiving influenza vaccination among patients with hypertension, diabetes, ischemic heart disease, or COPD" = "Flu_Vaccine",
  "PK15 - Frequency of limb amputation due to diabetes" = "Amputation",
  "PK16 - Frequency of surgery for diabetic retinopathy" = "Retinopathy_Surgery"
)

# =============================================================================
# FIXED COMPREHENSIVE ANALYSIS FUNCTION - USING ONLY 12 INDICATORS
# =============================================================================

prepare_comprehensive_analysis_fixed <- function(dm_data) {
  
  # Define the 12 indicators we want to analyze - USING EXACT NAMES FROM DATA
  target_indicators <- names(indicator_mapping)
  
  # Prepare datasets using fixed approach
  cat("Preparing area datasets...\n")
  segregated_data <- prepare_area_data_fixed(dm_data, "segregated")
  complementary_data <- prepare_area_data_fixed(dm_data, "complementary")
  
  if (is.null(segregated_data) && is.null(complementary_data)) {
    stop("No data available for analysis from either area")
  }
  
  # Filter data to only include our 12 target indicators
  if (!is.null(segregated_data)) {
    segregated_data <- segregated_data %>%
      filter(Indicator %in% target_indicators)
    cat("Segregated data after filtering:", nrow(segregated_data), "rows\n")
  }
  
  if (!is.null(complementary_data)) {
    complementary_data <- complementary_data %>%
      filter(Indicator %in% target_indicators)
    cat("Complementary data after filtering:", nrow(complementary_data), "rows\n")
  }
  
  cat("Filtered to", length(target_indicators), "target indicators\n")
  
  # Apply Bayesian shrinkage to each indicator in each area
  apply_shrinkage_by_indicator <- function(data, area_name) {
    if (is.null(data)) return(NULL)
    
    results <- list()
    
    for (ind in target_indicators) {
      ind_data <- data[data$Indicator == ind, ]
      
      if (nrow(ind_data) > 0) {
        # Remove any duplicates by taking mean
        ind_data <- ind_data %>%
          group_by(ID.GMP, Indicator) %>%
          summarize(
            OBS = mean(OBS, na.rm = TRUE),
            EXP = mean(EXP, na.rm = TRUE),
            .groups = 'drop'
          )
        
        bayes_result <- empirical_bayes_shrinkage(
          observed = ind_data$OBS,
          expected = ind_data$EXP
        )
        
        results[[ind]] <- data.frame(
          ID.GMP = ind_data$ID.GMP,
          Indicator = ind,
          SR_Raw = bayes_result$sr_raw,
          SR_Shrunk = bayes_result$sr_shrunk,
          Reliability = bayes_result$reliability,
          Shrinkage_Factor = bayes_result$shrinkage_factor,
          Area = area_name
        )
      }
    }
    
    if (length(results) == 0) return(NULL)
    return(bind_rows(results))
  }
  
  cat("Applying Bayesian shrinkage...\n")
  segregated_shrunk <- apply_shrinkage_by_indicator(segregated_data, "Segregated")
  complementary_shrunk <- apply_shrinkage_by_indicator(complementary_data, "Complementary")
  
  cat("Segregated shrunk data:", ifelse(is.null(segregated_shrunk), "NULL", paste(nrow(segregated_shrunk), "rows")), "\n")
  cat("Complementary shrunk data:", ifelse(is.null(complementary_shrunk), "NULL", paste(nrow(complementary_shrunk), "rows")), "\n")
  
  # Create combined dataset - FIXED NULL HANDLING
  if (!is.null(segregated_shrunk) && !is.null(complementary_shrunk)) {
    combined_shrunk <- segregated_shrunk %>%
      inner_join(complementary_shrunk, by = c("ID.GMP", "Indicator"), 
                 suffix = c("_seg", "_comp")) %>%
      mutate(
        SR_Shrunk_combined = (SR_Shrunk_seg + SR_Shrunk_comp) / 2,
        Reliability_combined = (Reliability_seg + Reliability_comp) / 2,
        Area_combined = "Combined"
      )
  } else if (!is.null(segregated_shrunk)) {
    combined_shrunk <- segregated_shrunk %>%
      mutate(
        SR_Shrunk_combined = SR_Shrunk,
        Reliability_combined = Reliability,
        Area_combined = "Combined"
      )
  } else if (!is.null(complementary_shrunk)) {
    combined_shrunk <- complementary_shrunk %>%
      mutate(
        SR_Shrunk_combined = SR_Shrunk,
        Reliability_combined = Reliability,
        Area_combined = "Combined"
      )
  } else {
    stop("No data available after Bayesian shrinkage")
  }
  
  # Apply shorter names to indicators
  apply_short_names <- function(data) {
    if (is.null(data)) return(NULL)
    data$Indicator_Short <- indicator_mapping[data$Indicator]
    return(data)
  }
  
  segregated_shrunk <- apply_short_names(segregated_shrunk)
  complementary_shrunk <- apply_short_names(complementary_shrunk)
  combined_shrunk <- apply_short_names(combined_shrunk)
  
  # Reshape for PCA using the fixed function
  cat("Reshaping data for PCA...\n")
  segregated_wide <- reshape_for_analysis(segregated_shrunk, "SR_Shrunk")
  complementary_wide <- reshape_for_analysis(complementary_shrunk, "SR_Shrunk")
  combined_wide <- reshape_for_analysis(combined_shrunk, "SR_Shrunk_combined")
  
  # Apply shorter names to column names in wide data
  apply_short_names_to_wide <- function(wide_data) {
    if (is.null(wide_data)) return(NULL)
    new_colnames <- sapply(colnames(wide_data), function(x) {
      if (x %in% names(indicator_mapping)) indicator_mapping[x] else x
    })
    colnames(wide_data) <- new_colnames
    return(wide_data)
  }
  
  segregated_wide <- apply_short_names_to_wide(segregated_wide)
  complementary_wide <- apply_short_names_to_wide(complementary_wide)
  combined_wide <- apply_short_names_to_wide(combined_wide)
  
  # Perform PCA for each area
  cat("Performing PCA...\n")
  pca_results <- list()
  
  if (!is.null(segregated_wide) && ncol(segregated_wide) >= 2 && nrow(segregated_wide) >= 2) {
    tryCatch({
      pca_results$segregated <- calculate_pca_components(segregated_wide)
      cat("PCA successful for segregated area\n")
    }, error = function(e) {
      cat("PCA failed for segregated area:", e$message, "\n")
    })
  }
  
  if (!is.null(complementary_wide) && ncol(complementary_wide) >= 2 && nrow(complementary_wide) >= 2) {
    tryCatch({
      pca_results$complementary <- calculate_pca_components(complementary_wide)
      cat("PCA successful for complementary area\n")
    }, error = function(e) {
      cat("PCA failed for complementary area:", e$message, "\n")
    })
  }
  
  if (!is.null(combined_wide) && ncol(combined_wide) >= 2 && nrow(combined_wide) >= 2) {
    tryCatch({
      pca_results$combined <- calculate_pca_components(combined_wide)
      cat("PCA successful for combined area\n")
    }, error = function(e) {
      cat("PCA failed for combined area:", e$message, "\n")
    })
  }
  
  # Calculate composite indicators
  cat("Calculating composite indicators...\n")
  composite_scores <- list()
  
  if (!is.null(pca_results$segregated)) {
    composite_scores$segregated <- calculate_composite_indicator(
      segregated_wide, 
      pca_results$segregated$loadings[, 1]
    )
  }
  
  if (!is.null(pca_results$complementary)) {
    composite_scores$complementary <- calculate_composite_indicator(
      complementary_wide, 
      pca_results$complementary$loadings[, 1]
    )
  }
  
  if (!is.null(pca_results$combined)) {
    composite_scores$combined <- calculate_composite_indicator(
      combined_wide, 
      pca_results$combined$loadings[, 1]
    )
  }
  
  # Create comprehensive results table
  cat("Creating comprehensive results table...\n")
  
  # First create the base results
  if (!is.null(segregated_shrunk) && !is.null(complementary_shrunk)) {
    comprehensive_results <- segregated_shrunk %>%
      inner_join(complementary_shrunk, by = c("ID.GMP", "Indicator"), 
                 suffix = c("_seg", "_comp"))
  } else if (!is.null(segregated_shrunk)) {
    comprehensive_results <- segregated_shrunk
  } else {
    comprehensive_results <- complementary_shrunk
  }
  
  # Add combined data
  if (!is.null(combined_shrunk)) {
    combined_subset <- combined_shrunk[, c("ID.GMP", "Indicator", "SR_Shrunk_combined", "Reliability_combined", "Indicator_Short")]
    comprehensive_results <- comprehensive_results %>%
      left_join(combined_subset, by = c("ID.GMP", "Indicator"))
  }
  
  # Create composite scores data frame
  composite_scores_df <- data.frame(ID.GMP = character())
  
  if (!is.null(composite_scores$segregated) && !is.null(segregated_wide)) {
    composite_scores_df <- data.frame(
      ID.GMP = rownames(segregated_wide),
      Composite_Segregated = composite_scores$segregated
    )
  }
  
  if (!is.null(composite_scores$complementary) && !is.null(complementary_wide)) {
    temp_df <- data.frame(
      ID.GMP = rownames(complementary_wide),
      Composite_Complementary = composite_scores$complementary
    )
    if (nrow(composite_scores_df) == 0) {
      composite_scores_df <- temp_df
    } else {
      composite_scores_df <- composite_scores_df %>%
        full_join(temp_df, by = "ID.GMP")
    }
  }
  
  if (!is.null(composite_scores$combined) && !is.null(combined_wide)) {
    temp_df <- data.frame(
      ID.GMP = rownames(combined_wide),
      Composite_Combined = composite_scores$combined
    )
    if (nrow(composite_scores_df) == 0) {
      composite_scores_df <- temp_df
    } else {
      composite_scores_df <- composite_scores_df %>%
        full_join(temp_df, by = "ID.GMP")
    }
  }
  
  # Create indicator-level analysis
  indicator_level_analysis <- comprehensive_results %>%
    group_by(Indicator_Short) %>%
    summarize(
      Segregated_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_Mean = mean(SR_Shrunk_comp, na.rm = TRUE),
      Combined_Mean = mean(SR_Shrunk_combined, na.rm = TRUE),
      Segregated_SD = sd(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_SD = sd(SR_Shrunk_comp, na.rm = TRUE),
      Combined_SD = sd(SR_Shrunk_combined, na.rm = TRUE),
      N_GMPs = n_distinct(ID.GMP)
    )
  
  # Variance explained
  variance_explained <- data.frame(
    Area = character(),
    PC1_Variance = numeric()
  )
  
  if (!is.null(pca_results$segregated)) {
    variance_explained <- rbind(variance_explained, 
                                data.frame(Area = "Segregated", 
                                           PC1_Variance = pca_results$segregated$variance_explained[1]))
  }
  if (!is.null(pca_results$complementary)) {
    variance_explained <- rbind(variance_explained, 
                                data.frame(Area = "Complementary", 
                                           PC1_Variance = pca_results$complementary$variance_explained[1]))
  }
  if (!is.null(pca_results$combined)) {
    variance_explained <- rbind(variance_explained, 
                                data.frame(Area = "Combined", 
                                           PC1_Variance = pca_results$combined$variance_explained[1]))
  }
  
  return(list(
    comprehensive_results = comprehensive_results,
    composite_scores = composite_scores_df,
    indicator_analysis = indicator_level_analysis,
    pca_results = pca_results,
    variance_explained = variance_explained,
    used_indicators = target_indicators,
    segregated_wide = segregated_wide,
    complementary_wide = complementary_wide,
    combined_wide = combined_wide
  ))
}

# =============================================================================
# RUN THE FIXED ANALYSIS
# =============================================================================

cat("Starting FIXED comprehensive analysis with 12 indicators...\n")
analysis_results <- prepare_comprehensive_analysis_fixed(dm_data)

cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("Used indicators:", length(analysis_results$used_indicators), "\n")
cat("Comprehensive results dimensions:", dim(analysis_results$comprehensive_results), "\n")
cat("Composite scores dimensions:", dim(analysis_results$composite_scores), "\n")

# Display key results
cat("\n=== KEY RESULTS ===\n")
cat("Variance Explained by PC1:\n")
print(analysis_results$variance_explained)

cat("\nIndicator-Level Scores:\n")
print(analysis_results$indicator_analysis)

# [REST OF THE VISUALIZATION AND EXPORT CODE REMAINS THE SAME AS BEFORE]

# =============================================================================
# COMPREHENSIVE VISUALIZATION - FOCUSED ON 12 INDICATORS
# =============================================================================

# Create output directory
output_dir <- "DM_12_Indicators_Analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# =============================================================================
# 1. INDICATOR-LEVEL MEAN COMPARISONS
# =============================================================================

cat("\n=== Creating Indicator-Level Visualizations ===\n")

# Prepare indicator data for plotting
indicator_plot_data <- analysis_results$indicator_analysis %>%
  pivot_longer(cols = c(Segregated_Mean, Complementary_Mean, Combined_Mean), 
               names_to = "Area", 
               values_to = "Score") %>%
  mutate(Area = gsub("_Mean", "", Area))

# Plot 1: Individual Indicator Means Comparison
p1 <- ggplot(indicator_plot_data, aes(x = reorder(Indicator_Short, Score), y = Score, fill = Area)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Mean Values for All 12 Diabetes Care Indicators",
       subtitle = "Comparison across Segregated, Complementary and Combined Areas",
       x = "Indicator", y = "Mean Value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "01_Indicator_Means_Comparison.png"), p1, width = 12, height = 8)

# Plot 2: Indicator Standard Deviations
p2_data <- analysis_results$indicator_analysis %>%
  pivot_longer(cols = c(Segregated_SD, Complementary_SD, Combined_SD), 
               names_to = "Area", 
               values_to = "SD") %>%
  mutate(Area = gsub("_SD", "", Area))

p2 <- ggplot(p2_data, aes(x = reorder(Indicator_Short, SD), y = SD, fill = Area)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Standard Deviation for All 12 Diabetes Care Indicators",
       x = "Indicator", y = "Standard Deviation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "02_Indicator_SD_Comparison.png"), p2, width = 12, height = 8)

# =============================================================================
# 2. HEATMAP OF INDICATOR PERFORMANCE
# =============================================================================

# Plot 3: Heatmap of indicator scores across areas
heatmap_data <- analysis_results$indicator_analysis %>%
  dplyr::select(Indicator_Short, Segregated_Mean, Complementary_Mean, Combined_Mean) %>%
  column_to_rownames("Indicator_Short") %>%
  as.matrix()

p3 <- ggplot(melt(heatmap_data), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 3)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(heatmap_data, na.rm = TRUE),
                       name = "Score") +
  labs(title = "Heatmap of Diabetes Care Indicator Performance",
       x = "Area", y = "Indicator") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "03_Indicator_Heatmap.png"), p3, width = 10, height = 8)

# =============================================================================
# 3. INDICATOR DISTRIBUTIONS BY AREA
# =============================================================================

# Plot 4: Box plots showing distribution for each indicator across areas
if ("SR_Shrunk_seg" %in% names(analysis_results$comprehensive_results) &&
    "SR_Shrunk_comp" %in% names(analysis_results$comprehensive_results)) {
  
  p4_data <- analysis_results$comprehensive_results %>%
    dplyr::select(Indicator_Short, SR_Shrunk_seg, SR_Shrunk_comp, SR_Shrunk_combined) %>%
    pivot_longer(cols = c(SR_Shrunk_seg, SR_Shrunk_comp, SR_Shrunk_combined),
                 names_to = "Area",
                 values_to = "Value") %>%
    mutate(Area = recode(Area, 
                         "SR_Shrunk_seg" = "Segregated",
                         "SR_Shrunk_comp" = "Complementary",
                         "SR_Shrunk_combined" = "Combined"))
  
  p4 <- ggplot(p4_data, aes(x = Indicator_Short, y = Value, fill = Area)) +
    geom_boxplot(alpha = 0.7) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Distribution of Values for All 12 Diabetes Care Indicators",
         x = "Indicator", y = "Value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "04_Indicator_Distributions.png"), p4, width = 12, height = 8)
}

# =============================================================================
# 4. PCA LOADINGS FOR THE 12 INDICATORS
# =============================================================================

# Plot 5: PCA Loadings comparison
if (!is.null(analysis_results$pca_results$segregated) && 
    !is.null(analysis_results$pca_results$complementary)) {
  
  loadings_seg <- data.frame(
    Indicator = rownames(analysis_results$pca_results$segregated$loadings),
    Loading = analysis_results$pca_results$segregated$loadings[, 1],
    Area = "Segregated"
  )
  
  loadings_comp <- data.frame(
    Indicator = rownames(analysis_results$pca_results$complementary$loadings),
    Loading = analysis_results$pca_results$complementary$loadings[, 1],
    Area = "Complementary"
  )
  
  loadings_combined <- rbind(loadings_seg, loadings_comp)
  
  p5 <- ggplot(loadings_combined, aes(x = reorder(Indicator, Loading), 
                                      y = Loading, fill = Area)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "PCA Loadings for Diabetes Care Indicators",
         subtitle = "First Principal Component Loadings by Area",
         x = "Indicator", y = "Loading on PC1") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir, "05_PCA_Loadings.png"), p5, width = 12, height = 8)
}

# =============================================================================
# 5. COMPOSITE SCORE COMPARISONS
# =============================================================================

# Plot 6: Composite scores comparison
if (nrow(analysis_results$composite_scores) > 0) {
  composite_plot_data <- analysis_results$composite_scores %>%
    pivot_longer(cols = starts_with("Composite_"), 
                 names_to = "Area", 
                 values_to = "Score") %>%
    mutate(Area = gsub("Composite_", "", Area))
  
  p6 <- ggplot(composite_plot_data, aes(x = Area, y = Score, fill = Area)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(alpha = 0.5, width = 0.2, size = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Composite Diabetes Care Scores by Area",
         x = "Area", y = "Composite Score") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "06_Composite_Scores.png"), p6, width = 10, height = 6)
}

# =============================================================================
# 6. SCATTER PLOT: SEGREGATED VS COMPLEMENTARY
# =============================================================================

# Plot 7: Scatter plot comparing segregated vs complementary for each indicator
if ("SR_Shrunk_seg" %in% names(analysis_results$comprehensive_results) &&
    "SR_Shrunk_comp" %in% names(analysis_results$comprehensive_results)) {
  
  # Calculate mean values for each indicator
  scatter_data <- analysis_results$comprehensive_results %>%
    group_by(Indicator_Short) %>%
    summarize(
      Segregated_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_Mean = mean(SR_Shrunk_comp, na.rm = TRUE)
    )
  
  p7 <- ggplot(scatter_data, aes(x = Segregated_Mean, y = Complementary_Mean, label = Indicator_Short)) +
    geom_point(size = 3, color = "steelblue", alpha = 0.7) +
    geom_smooth(method = "lm", color = "red", se = TRUE, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "solid") +
    geom_text(hjust = 0, vjust = 0, nudge_x = 0.01, size = 3) +
    labs(title = "Segregated vs Complementary Area Comparison",
         subtitle = "Each point represents one diabetes care indicator",
         x = "Segregated Area Mean", y = "Complementary Area Mean") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  ggsave(file.path(output_dir, "07_Segregated_vs_Complementary.png"), p7, width = 10, height = 8)
}

# =============================================================================
# 7. INDICATOR RANKING PLOTS
# =============================================================================

# Plot 8: Ranking of indicators by performance
ranking_data <- analysis_results$indicator_analysis %>%
  dplyr::select(Indicator_Short, Segregated_Mean, Complementary_Mean, Combined_Mean) %>%
  pivot_longer(cols = c(Segregated_Mean, Complementary_Mean, Combined_Mean),
               names_to = "Area",
               values_to = "Mean") %>%
  mutate(Area = gsub("_Mean", "", Area)) %>%
  group_by(Area) %>%
  mutate(Rank = rank(-Mean))

p8 <- ggplot(ranking_data, aes(x = Area, y = Rank, label = Indicator_Short, color = Indicator_Short)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.1, size = 3) +
  scale_y_reverse(breaks = 1:12) +
  labs(title = "Indicator Performance Ranking by Area",
       x = "Area", y = "Rank (1 = Best)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "08_Indicator_Rankings.png"), p8, width = 10, height = 8)

# =============================================================================
# 8. VARIANCE EXPLAINED PLOT
# =============================================================================

# Plot 9: Variance explained by PCA
p9 <- ggplot(analysis_results$variance_explained, 
             aes(x = Area, y = PC1_Variance, fill = Area)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(round(PC1_Variance, 1), "%")), 
            vjust = -0.5, size = 5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Variance Explained by First Principal Component",
       subtitle = "For Diabetes Care Indicators Across Areas",
       x = "Area", y = "Variance Explained (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "09_PCA_Variance_Explained.png"), p9, width = 10, height = 6)

# =============================================================================
# 9. RELIABILITY COMPARISON
# =============================================================================

# Plot 10: Reliability comparison
reliability_cols <- grep("Reliability", names(analysis_results$comprehensive_results), 
                         value = TRUE, ignore.case = TRUE)

if (length(reliability_cols) >= 2) {
  reliability_data <- analysis_results$comprehensive_results %>%
    dplyr::select(Indicator_Short, all_of(reliability_cols[1:2])) %>%
    pivot_longer(cols = all_of(reliability_cols[1:2]),
                 names_to = "Area",
                 values_to = "Reliability") %>%
    mutate(Area = gsub("Reliability_", "", Area))
  
  p10 <- ggplot(reliability_data, aes(x = Area, y = Reliability, fill = Area)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Reliability of Diabetes Care Indicators by Area",
         x = "Area", y = "Reliability") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "10_Reliability_Comparison.png"), p10, width = 10, height = 6)
}

# =============================================================================
# COMPREHENSIVE EXCEL EXPORT
# =============================================================================

cat("\n=== Creating Comprehensive Excel Export ===\n")

# Create workbook
wb <- createWorkbook()

# Sheet 1: Summary Statistics
addWorksheet(wb, "Summary")
summary_data <- data.frame(
  Metric = c("Total GMPs", "Total Indicators", "Segregated Area Available", 
             "Complementary Area Available", "Combined Analysis Available"),
  Value = c(
    length(unique(analysis_results$comprehensive_results$ID.GMP)),
    length(unique(analysis_results$comprehensive_results$Indicator_Short)),
    !is.null(analysis_results$segregated_wide),
    !is.null(analysis_results$complementary_wide),
    !is.null(analysis_results$combined_wide)
  )
)
writeData(wb, "Summary", summary_data)

# Sheet 2: Indicator Analysis
addWorksheet(wb, "Indicator_Analysis")
writeData(wb, "Indicator_Analysis", analysis_results$indicator_analysis)

# Sheet 3: Composite Scores
addWorksheet(wb, "Composite_Scores")
writeData(wb, "Composite_Scores", analysis_results$composite_scores)

# Sheet 4: Variance Explained
addWorksheet(wb, "Variance_Explained")
writeData(wb, "Variance_Explained", analysis_results$variance_explained)

# Sheet 5: Detailed Results by Indicator and GMP
addWorksheet(wb, "Detailed_Results")
detailed_results <- analysis_results$comprehensive_results
writeData(wb, "Detailed_Results", detailed_results)

# Sheet 6: PCA Loadings
if (!is.null(analysis_results$pca_results$segregated) && 
    !is.null(analysis_results$pca_results$complementary)) {
  
  addWorksheet(wb, "PCA_Loadings")
  
  pca_loadings <- data.frame(
    Indicator = rownames(analysis_results$pca_results$segregated$loadings),
    Segregated_PC1 = analysis_results$pca_results$segregated$loadings[, 1],
    Complementary_PC1 = analysis_results$pca_results$complementary$loadings[, 1]
  )
  
  if (!is.null(analysis_results$pca_results$combined)) {
    pca_loadings$Combined_PC1 <- analysis_results$pca_results$combined$loadings[, 1]
  }
  
  writeData(wb, "PCA_Loadings", pca_loadings)
}

# Save workbook
excel_file <- file.path(output_dir, "DM_12_Indicators_Analysis_Results.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("Excel file saved to:", excel_file, "\n")

# =============================================================================
# DISPLAY ALL PLOTS
# =============================================================================

cat("\n=== Displaying All Plots ===\n")

# Display plots
print(p1)
cat("Plot 1: Indicator Means Comparison displayed\n")

print(p2)
cat("Plot 2: Indicator SD Comparison displayed\n")

print(p3)
cat("Plot 3: Indicator Heatmap displayed\n")

if (exists("p4")) {
  print(p4)
  cat("Plot 4: Indicator Distributions displayed\n")
}

if (exists("p5")) {
  print(p5)
  cat("Plot 5: PCA Loadings displayed\n")
}

if (exists("p6")) {
  print(p6)
  cat("Plot 6: Composite Scores displayed\n")
}

if (exists("p7")) {
  print(p7)
  cat("Plot 7: Segregated vs Complementary displayed\n")
}

print(p8)
cat("Plot 8: Indicator Rankings displayed\n")

print(p9)
cat("Plot 9: PCA Variance Explained displayed\n")

if (exists("p10")) {
  print(p10)
  cat("Plot 10: Reliability Comparison displayed\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")