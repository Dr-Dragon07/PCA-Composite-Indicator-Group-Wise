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

dm_data<-DM_Related_to_work_2
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
# FLEXIBLE DATA PREPARATION FUNCTION
# =============================================================================

prepare_area_data_flexible <- function(data, area_type) {
  
  cat("Preparing data for", area_type, "area...\n")
  
  # Find columns using flexible matching
  find_column <- function(pattern) {
    cols <- colnames(data)
    matches <- cols[grepl(pattern, cols, ignore.case = TRUE)]
    if (length(matches) > 0) {
      return(matches[1])
    } else {
      return(NULL)
    }
  }
  
  # Find essential columns
  id_col <- find_column("ID|id|Id|GMP|gmp")
  time_col <- find_column("Time|time|TIME")
  group_col <- find_column("Group|group|GROUP")
  indicator_col <- find_column("Indicator|indicator|INDICATOR")
  
  cat("Found columns:\n")
  cat("  ID:", id_col, "\n")
  cat("  Time:", time_col, "\n")
  cat("  Group:", group_col, "\n")
  cat("  Indicator:", indicator_col, "\n")
  
  if (area_type == "segregated") {
    obs_col <- find_column("OBS.*segregated|segregated.*OBS|OBS.*seg|seg.*OBS")
    target_col <- find_column("Target.*segregated|segregated.*Target|Target.*seg|seg.*Target")
  } else if (area_type == "complementary") {
    obs_col <- find_column("OBS.*complementary|complementary.*OBS|OBS.*comp|comp.*OBS")
    target_col <- find_column("Target.*complementary|complementary.*Target|Target.*comp|comp.*Target")
  }
  
  cat("  OBS:", obs_col, "\n")
  cat("  Target:", target_col, "\n")
  
  # Check if we found the necessary columns
  required_cols <- c(id_col, time_col, group_col, indicator_col, obs_col, target_col)
  if (any(sapply(required_cols, is.null))) {
    warning(paste("Could not find all necessary columns for", area_type, "area"))
    cat("Missing columns:\n")
    missing <- which(sapply(required_cols, is.null))
    print(names(required_cols)[missing])
    return(NULL)
  }
  
  # Select and rename columns
  area_data <- data[, c(id_col, time_col, group_col, indicator_col, obs_col, target_col)]
  colnames(area_data) <- c("ID.GMP", "Time", "Group.of.indicators", "Indicator", "OBS", "Target")
  
  # Convert to numeric if needed
  area_data$OBS <- as.numeric(area_data$OBS)
  area_data$Target <- as.numeric(area_data$Target)
  
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
# FIXED COMPREHENSIVE ANALYSIS FUNCTION - WITH SELECT CONFLICT RESOLVED
# =============================================================================

prepare_comprehensive_analysis_fixed <- function(dm_data) {
  
  # Define indicators and groups
  indicators <- c(
    "FV30a.Prescription.Alimentary.tract.and.metabolism.ATC.A",
    "FV31a.Drug.substitution.Alimentary.tract.and.metabolism.ATC.A", 
    "FV32a.Drug.substitution.rate.Alimentary.tract.and.metabolism.ATC.A",
    "FV06.Proportion.of.diabetic.and.or.hypertensive.patients.who.underwent.blood.lipid.testing",
    "FV08.Proportion.of.diabetics.who.underwent.hemoglobin.A1c.testing",
    "FV09.Proportion.of.diabetics.who.underwent.ophthalmological.examination",
    "FV13.Proportion.of.people.aged.40.54.who.are.switching.from.medication.for.diabetes",
    "FV14.Proportion.of.people.aged.55.69.who.are.switching.from.medication.for.diabetes",
    "FV15.Proportion.of.diabetic.patients.who.underwent.serum.creatinine.level.determination",
    "FV16a.Proportion.of.people.under.65.years.of.age.receiving.influenza.vaccination.among.patients.with.hypertension.diabetes.ischemic.heart.disease.or.COPD",
    "PK15.Frequency.of.limb.amputation.due.to.diabetes",
    "PK16.Frequency.of.surgery.for.diabetic.retinopathy"
  )
  
  indicator_groups <- list(
    Pharmaceutical = c(
      "FV30a.Prescription.Alimentary.tract.and.metabolism.ATC.A",
      "FV31a.Drug.substitution.Alimentary.tract.and.metabolism.ATC.A", 
      "FV32a.Drug.substitution.rate.Alimentary.tract.and.metabolism.ATC.A"
    ),
    GP_NEFMI = c(
      "FV06.Proportion.of.diabetic.and.or.hypertensive.patients.who.underwent.blood.lipid.testing",
      "FV08.Proportion.of.diabetics.who.underwent.hemoglobin.A1c.testing",
      "FV09.Proportion.of.diabetics.who.underwent.ophthalmological.examination"
    ),
    GP_Non_NEFMI = c(
      "FV13.Proportion.of.people.aged.40.54.who.are.switching.from.medication.for.diabetes",
      "FV14.Proportion.of.people.aged.55.69.who.are.switching.from.medication.for.diabetes", 
      "FV15.Proportion.of.diabetic.patients.who.underwent.serum.creatinine.level.determination",
      "FV16a.Proportion.of.people.under.65.years.of.age.receiving.influenza.vaccination.among.patients.with.hypertension.diabetes.ischemic.heart.disease.or.COPD"
    ),
    Metabolic_Diseases = c("PK15.Frequency.of.limb.amputation.due.to.diabetes"),
    Professional_Care = c("PK16.Frequency.of.surgery.for.diabetic.retinopathy")
  )
  
  # Prepare datasets using flexible approach
  cat("Preparing area datasets...\n")
  segregated_data <- prepare_area_data_flexible(dm_data, "segregated")
  complementary_data <- prepare_area_data_flexible(dm_data, "complementary")
  
  if (is.null(segregated_data) && is.null(complementary_data)) {
    stop("No data available for analysis from either area")
  }
  
  # Get actual indicator names from the data
  actual_indicators <- unique(c(
    if(!is.null(segregated_data)) unique(segregated_data$Indicator) else NULL,
    if(!is.null(complementary_data)) unique(complementary_data$Indicator) else NULL
  ))
  
  cat("Actual indicators in data:\n")
  print(actual_indicators)
  
  # Map the actual group names to our expected structure
  group_mapping <- list(
    "Pharmaceutical Indicators" = "Pharmaceutical",
    "General practitioner indicators (included in the NEFMI decree)" = "GP_NEFMI", 
    "General practitioner indicators (not included in the NEFMI decree)" = "GP_Non_NEFMI",
    "Metabolic Diseases Indicators" = "Metabolic_Diseases",
    "Professional care utilization indicators" = "Professional_Care"
  )
  
  # Apply Bayesian shrinkage to each indicator in each area
  apply_shrinkage_by_indicator <- function(data, area_name) {
    if (is.null(data)) return(NULL)
    
    results <- list()
    
    for (ind in actual_indicators) {
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
  
  # Create combined dataset
  if (!is.null(segregated_shrunk) && !is.null(complementary_shrunk)) {
    combined_shrunk <- segregated_shrunk %>%
      inner_join(complementary_shrunk, by = c("ID.GMP", "Indicator"), 
                 suffix = c("_seg", "_comp"), relationship = "many-to-many") %>%
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
  } else {
    combined_shrunk <- complementary_shrunk %>%
      mutate(
        SR_Shrunk_combined = SR_Shrunk,
        Reliability_combined = Reliability,
        Area_combined = "Combined"
      )
  }
  
  # Reshape for PCA using the fixed function
  cat("Reshaping data for PCA...\n")
  segregated_wide <- reshape_for_analysis(segregated_shrunk, "SR_Shrunk")
  complementary_wide <- reshape_for_analysis(complementary_shrunk, "SR_Shrunk")
  combined_wide <- reshape_for_analysis(combined_shrunk, "SR_Shrunk_combined")
  
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
  
  # Create comprehensive results table - FIXED SELECT CONFLICT
  cat("Creating comprehensive results table...\n")
  
  # First create the base results
  if (!is.null(segregated_shrunk) && !is.null(complementary_shrunk)) {
    comprehensive_results <- segregated_shrunk %>%
      inner_join(complementary_shrunk, by = c("ID.GMP", "Indicator"), 
                 suffix = c("_seg", "_comp"), relationship = "many-to-many")
  } else if (!is.null(segregated_shrunk)) {
    comprehensive_results <- segregated_shrunk
  } else {
    comprehensive_results <- complementary_shrunk
  }
  
  # Add combined data using base R to avoid select conflict
  if (!is.null(combined_shrunk)) {
    # Use base R to select the columns we need
    combined_subset <- combined_shrunk[, c("ID.GMP", "Indicator", "SR_Shrunk_combined", "Reliability_combined")]
    comprehensive_results <- comprehensive_results %>%
      left_join(combined_subset, by = c("ID.GMP", "Indicator"))
  }
  
  # Add group mapping
  comprehensive_results <- comprehensive_results %>%
    mutate(
      Group = case_when(
        Indicator == "Pharmaceutical Indicators" ~ "Pharmaceutical",
        Indicator == "General practitioner indicators (included in the NEFMI decree)" ~ "GP_NEFMI",
        Indicator == "General practitioner indicators (not included in the NEFMI decree)" ~ "GP_Non_NEFMI", 
        Indicator == "Metabolic Diseases Indicators" ~ "Metabolic_Diseases",
        Indicator == "Professional care utilization indicators" ~ "Professional_Care",
        TRUE ~ "Other"
      )
    )
  
  # Add PCA weights if available - using base R to avoid select conflict
  if (!is.null(pca_results$segregated)) {
    pca_weights_seg <- data.frame(
      Indicator = colnames(segregated_wide),
      PCA_Weight_Segregated = pca_results$segregated$loadings[, 1]
    )
    comprehensive_results <- comprehensive_results %>%
      left_join(pca_weights_seg, by = "Indicator")
  }
  
  if (!is.null(pca_results$complementary)) {
    pca_weights_comp <- data.frame(
      Indicator = colnames(complementary_wide),
      PCA_Weight_Complementary = pca_results$complementary$loadings[, 1]
    )
    comprehensive_results <- comprehensive_results %>%
      left_join(pca_weights_comp, by = "Indicator")
  }
  
  if (!is.null(pca_results$combined)) {
    pca_weights_comb <- data.frame(
      Indicator = colnames(combined_wide),
      PCA_Weight_Combined = pca_results$combined$loadings[, 1]
    )
    comprehensive_results <- comprehensive_results %>%
      left_join(pca_weights_comb, by = "Indicator")
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
  
  # Create group-level analysis
  group_level_analysis <- comprehensive_results %>%
    group_by(Group, ID.GMP) %>%
    summarize(
      Segregated_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_Mean = mean(SR_Shrunk_comp, na.rm = TRUE),
      Combined_Mean = mean(SR_Shrunk_combined, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    group_by(Group) %>%
    summarize(
      Segregated_Group_Score = mean(Segregated_Mean, na.rm = TRUE),
      Complementary_Group_Score = mean(Complementary_Mean, na.rm = TRUE),
      Combined_Group_Score = mean(Combined_Mean, na.rm = TRUE)
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
    group_analysis = group_level_analysis,
    pca_results = pca_results,
    variance_explained = variance_explained,
    used_indicators = actual_indicators,
    segregated_wide = segregated_wide,
    complementary_wide = complementary_wide,
    combined_wide = combined_wide
  ))
}

# =============================================================================
# RUN THE FIXED ANALYSIS
# =============================================================================

cat("Starting FIXED comprehensive analysis...\n")
analysis_results <- prepare_comprehensive_analysis_fixed(dm_data)

cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("Used indicators:", length(analysis_results$used_indicators), "\n")
cat("Comprehensive results dimensions:", dim(analysis_results$comprehensive_results), "\n")
cat("Composite scores dimensions:", dim(analysis_results$composite_scores), "\n")

# Display key results
cat("\n=== KEY RESULTS ===\n")
cat("Variance Explained by PC1:\n")
print(analysis_results$variance_explained)

cat("\nGroup-Level Scores:\n")
print(analysis_results$group_analysis)

# Show first few rows of comprehensive results
cat("\nFirst few rows of comprehensive results:\n")
print(head(analysis_results$comprehensive_results))

# Show first few rows of composite scores
cat("\nFirst few rows of composite scores:\n")
print(head(analysis_results$composite_scores))

# Show PCA weights if available
if (!is.null(analysis_results$pca_results$segregated)) {
  cat("\nPCA Weights for Segregated Area:\n")
  print(analysis_results$pca_results$segregated$loadings[, 1])
}

if (!is.null(analysis_results$pca_results$complementary)) {
  cat("\nPCA Weights for Complementary Area:\n")
  print(analysis_results$pca_results$complementary$loadings[, 1])
}

if (!is.null(analysis_results$pca_results$combined)) {
  cat("\nPCA Weights for Combined Area:\n")
  print(analysis_results$pca_results$combined$loadings[, 1])
}

# =============================================================================
# =============================================================================
# COMPREHENSIVE VISUALIZATION AND EXCEL EXPORT
# =============================================================================

# Create output directory
output_dir <- "DM_Analysis_Output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# =============================================================================
# 1. COMPOSITE INDICATOR COMPARISONS
# =============================================================================

cat("\n=== Creating Composite Indicator Visualizations ===\n")

# Prepare composite scores for plotting
composite_plot_data <- analysis_results$composite_scores %>%
  pivot_longer(cols = starts_with("Composite_"), 
               names_to = "Area", 
               values_to = "Score") %>%
  mutate(Area = gsub("Composite_", "", Area))

# Box plot comparison
p1 <- ggplot(composite_plot_data, aes(x = Area, y = Score, fill = Area)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(alpha = 0.3, width = 0.2) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Composite Indicator Comparison Across Areas",
       x = "Area", y = "Composite Score") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Violin plot comparison
p2 <- ggplot(composite_plot_data, aes(x = Area, y = Score, fill = Area)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of Composite Indicators by Area",
       x = "Area", y = "Composite Score") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Density plot comparison
p3 <- ggplot(composite_plot_data, aes(x = Score, fill = Area)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Density Distribution of Composite Scores",
       x = "Composite Score", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Scatter plot: Segregated vs Complementary
if ("Composite_Segregated" %in% names(analysis_results$composite_scores) &&
    "Composite_Complementary" %in% names(analysis_results$composite_scores)) {
  p4 <- ggplot(analysis_results$composite_scores, 
               aes(x = Composite_Segregated, y = Composite_Complementary)) +
    geom_point(alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    labs(title = "Segregated vs Complementary Composite Scores",
         x = "Segregated Area Score", y = "Complementary Area Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Calculate correlation
  cor_val <- cor(analysis_results$composite_scores$Composite_Segregated,
                 analysis_results$composite_scores$Composite_Complementary,
                 use = "complete.obs")
  cat("Correlation between Segregated and Complementary:", round(cor_val, 3), "\n")
}

# Save composite plots
ggsave(file.path(output_dir, "01_Composite_Boxplot.png"), p1, width = 10, height = 6)
ggsave(file.path(output_dir, "02_Composite_Violin.png"), p2, width = 10, height = 6)
ggsave(file.path(output_dir, "03_Composite_Density.png"), p3, width = 10, height = 6)
if (exists("p4")) {
  ggsave(file.path(output_dir, "04_Composite_Scatter.png"), p4, width = 10, height = 6)
}

# =============================================================================
# 2. PCA VISUALIZATION
# =============================================================================

cat("\n=== Creating PCA Visualizations ===\n")

# Variance explained bar plot
p5 <- ggplot(analysis_results$variance_explained, 
             aes(x = Area, y = PC1_Variance, fill = Area)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = paste0(round(PC1_Variance, 1), "%")), 
            vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Variance Explained by First Principal Component",
       x = "Area", y = "Variance Explained (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "05_PCA_Variance_Explained.png"), p5, width = 10, height = 6)

# PCA loadings comparison
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
  
  p6 <- ggplot(loadings_combined, aes(x = reorder(Indicator, Loading), 
                                      y = Loading, fill = Area)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "PCA Loadings Comparison: Segregated vs Complementary",
         x = "Indicator", y = "Loading on PC1") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 8))
  
  ggsave(file.path(output_dir, "06_PCA_Loadings_Comparison.png"), p6, 
         width = 12, height = 8)
}

# =============================================================================
# 3. GROUP-LEVEL ANALYSIS PLOTS
# =============================================================================

cat("\n=== Creating Group-Level Visualizations ===\n")

# Prepare group data for plotting
group_plot_data <- analysis_results$group_analysis %>%
  pivot_longer(cols = ends_with("_Group_Score"), 
               names_to = "Area", 
               values_to = "Score") %>%
  mutate(Area = gsub("_Group_Score", "", Area))

# Group comparison bar plot
p7 <- ggplot(group_plot_data, aes(x = Group, y = Score, fill = Area)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Group-Level Scores by Area",
       x = "Indicator Group", y = "Mean Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "07_Group_Comparison.png"), p7, width = 12, height = 6)

# Heatmap of group scores
group_matrix <- analysis_results$group_analysis %>%
  column_to_rownames("Group") %>%
  as.matrix()

p8 <- ggplot(melt(group_matrix), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(group_matrix, na.rm = TRUE)) +
  labs(title = "Heatmap of Group-Level Scores",
       x = "Area", y = "Group", fill = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "08_Group_Heatmap.png"), p8, width = 10, height = 6)

# =============================================================================
# 4. INDICATOR-LEVEL COMPARISONS
# =============================================================================

cat("\n=== Creating Indicator-Level Visualizations ===\n")

# Prepare indicator comparison data
if ("SR_Shrunk_seg" %in% names(analysis_results$comprehensive_results) &&
    "SR_Shrunk_comp" %in% names(analysis_results$comprehensive_results)) {
  
  indicator_comparison <- analysis_results$comprehensive_results %>%
    group_by(Indicator) %>%
    summarize(
      Segregated_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_Mean = mean(SR_Shrunk_comp, na.rm = TRUE),
      Segregated_SD = sd(SR_Shrunk_seg, na.rm = TRUE),
      Complementary_SD = sd(SR_Shrunk_comp, na.rm = TRUE)
    ) %>%
    pivot_longer(cols = c(Segregated_Mean, Complementary_Mean),
                 names_to = "Area",
                 values_to = "Mean") %>%
    mutate(Area = gsub("_Mean", "", Area))
  
  p9 <- ggplot(indicator_comparison, aes(x = reorder(Indicator, Mean), 
                                         y = Mean, fill = Area)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Mean Indicator Values: Segregated vs Complementary",
         x = "Indicator", y = "Mean Value") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "09_Indicator_Comparison.png"), p9, 
         width = 12, height = 10)
}

# =============================================================================
# 5. RELIABILITY AND SHRINKAGE ANALYSIS
# =============================================================================

cat("\n=== Creating Reliability and Shrinkage Visualizations ===\n")

# Check which reliability columns exist
reliability_cols <- grep("Reliability", names(analysis_results$comprehensive_results), 
                         value = TRUE, ignore.case = TRUE)
cat("Available reliability columns:", paste(reliability_cols, collapse = ", "), "\n")

if (length(reliability_cols) >= 2) {
  
  # Use the actual column names found
  reliability_data <- analysis_results$comprehensive_results %>%
    dplyr::select(Indicator, all_of(reliability_cols[1:2])) %>%
    pivot_longer(cols = all_of(reliability_cols[1:2]),
                 names_to = "Area",
                 values_to = "Reliability") %>%
    mutate(Area = gsub("Reliability_", "", Area))
  
  p10 <- ggplot(reliability_data, aes(x = Area, y = Reliability, fill = Area)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Reliability Distribution by Area",
         x = "Area", y = "Reliability") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(output_dir, "10_Reliability_Comparison.png"), p10, 
         width = 10, height = 6)
  
  cat("Reliability plot created successfully\n")
} else {
  cat("Warning: Could not find sufficient reliability columns for comparison\n")
}

# =============================================================================
# 6. COMPREHENSIVE EXCEL EXPORT
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
    length(analysis_results$used_indicators),
    !is.null(analysis_results$segregated_wide),
    !is.null(analysis_results$complementary_wide),
    !is.null(analysis_results$combined_wide)
  )
)
writeData(wb, "Summary", summary_data)

# Sheet 2: Composite Scores
addWorksheet(wb, "Composite_Scores")
writeData(wb, "Composite_Scores", analysis_results$composite_scores)

# Sheet 3: Group Analysis
addWorksheet(wb, "Group_Analysis")
writeData(wb, "Group_Analysis", analysis_results$group_analysis)

# Sheet 4: Variance Explained
addWorksheet(wb, "Variance_Explained")
writeData(wb, "Variance_Explained", analysis_results$variance_explained)

# Sheet 5: Detailed Results by Indicator and GMP
addWorksheet(wb, "Detailed_Results")

# Create comprehensive table with all metrics
detailed_results <- analysis_results$comprehensive_results

# Calculate additional statistics per indicator
indicator_stats <- detailed_results %>%
  group_by(Indicator, Group) %>%
  summarize(
    # Segregated area statistics
    Seg_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
    Seg_SD = sd(SR_Shrunk_seg, na.rm = TRUE),
    Seg_Median = median(SR_Shrunk_seg, na.rm = TRUE),
    Seg_Mean_Reliability = mean(Reliability_seg, na.rm = TRUE),
    Seg_Mean_Shrinkage = mean(Shrinkage_Factor_seg, na.rm = TRUE),
    
    # Complementary area statistics
    Comp_Mean = mean(SR_Shrunk_comp, na.rm = TRUE),
    Comp_SD = sd(SR_Shrunk_comp, na.rm = TRUE),
    Comp_Median = median(SR_Shrunk_comp, na.rm = TRUE),
    Comp_Mean_Reliability = mean(Reliability_comp, na.rm = TRUE),
    Comp_Mean_Shrinkage = mean(Shrinkage_Factor_comp, na.rm = TRUE),
    
    # T-test comparison
    P_Value = tryCatch({
      t.test(SR_Shrunk_seg, SR_Shrunk_comp)$p.value
    }, error = function(e) NA),
    
    N = n(),
    .groups = 'drop'
  )

writeData(wb, "Detailed_Results", indicator_stats)

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

# Sheet 7: Raw Comprehensive Results (all GMPs, all indicators)
addWorksheet(wb, "All_GMP_Results")
writeData(wb, "All_GMP_Results", detailed_results)

# Sheet 8: Group-by-GMP Analysis
addWorksheet(wb, "Group_by_GMP")
group_by_gmp <- detailed_results %>%
  group_by(ID.GMP, Group) %>%
  summarize(
    Segregated_Mean = mean(SR_Shrunk_seg, na.rm = TRUE),
    Complementary_Mean = mean(SR_Shrunk_comp, na.rm = TRUE),
    Combined_Mean = mean(SR_Shrunk_combined, na.rm = TRUE),
    Segregated_Reliability = mean(Reliability_seg, na.rm = TRUE),
    Complementary_Reliability = mean(Reliability_comp, na.rm = TRUE),
    .groups = 'drop'
  )
writeData(wb, "Group_by_GMP", group_by_gmp)

# Sheet 9: Statistical Tests Summary
addWorksheet(wb, "Statistical_Tests")
statistical_tests <- data.frame(
  Test = character(),
  Comparison = character(),
  Statistic = numeric(),
  P_Value = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

# Overall comparison test
if ("Composite_Segregated" %in% names(analysis_results$composite_scores) &&
    "Composite_Complementary" %in% names(analysis_results$composite_scores)) {
  
  t_result <- t.test(analysis_results$composite_scores$Composite_Segregated,
                     analysis_results$composite_scores$Composite_Complementary,
                     paired = TRUE)
  
  statistical_tests <- rbind(statistical_tests, data.frame(
    Test = "Paired t-test",
    Comparison = "Segregated vs Complementary Composite",
    Statistic = t_result$statistic,
    P_Value = t_result$p.value,
    Interpretation = ifelse(t_result$p.value < 0.05, 
                            "Significant difference", 
                            "No significant difference")
  ))
}

writeData(wb, "Statistical_Tests", statistical_tests)

# Save workbook
excel_filename <- file.path(output_dir, "DM_Comprehensive_Analysis.xlsx")
saveWorkbook(wb, excel_filename, overwrite = TRUE)

cat("\n=== EXPORT COMPLETED ===\n")
cat("All visualizations saved to:", output_dir, "\n")
cat("Excel file saved as:", excel_filename, "\n")

# Print summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total GMPs analyzed:", length(unique(analysis_results$comprehensive_results$ID.GMP)), "\n")
cat("Total indicators:", length(analysis_results$used_indicators), "\n")
cat("Segregated area variance explained:", 
    analysis_results$variance_explained$PC1_Variance[analysis_results$variance_explained$Area == "Segregated"], "%\n")
cat("Complementary area variance explained:", 
    analysis_results$variance_explained$PC1_Variance[analysis_results$variance_explained$Area == "Complementary"], "%\n")
cat("\nFiles created:\n")
cat("  - 10 visualization plots (PNG)\n")
cat("  - 1 comprehensive Excel workbook with 9 sheets\n")

# =============================================================================
# DISPLAY ALL PLOTS IN R STUDIO
# =============================================================================

cat("\n=== Displaying All Plots ===\n")

# Display plots one by one
if (exists("p1")) {
  print(p1)
  cat("Plot 1: Composite Boxplot displayed\n")
}

if (exists("p2")) {
  print(p2)
  cat("Plot 2: Composite Violin displayed\n")
}

if (exists("p3")) {
  print(p3)
  cat("Plot 3: Composite Density displayed\n")
}

if (exists("p4")) {
  print(p4)
  cat("Plot 4: Composite Scatter displayed\n")
}

if (exists("p5")) {
  print(p5)
  cat("Plot 5: PCA Variance Explained displayed\n")
}

if (exists("p6")) {
  print(p6)
  cat("Plot 6: PCA Loadings Comparison displayed\n")
}

if (exists("p7")) {
  print(p7)
  cat("Plot 7: Group Comparison displayed\n")
}

if (exists("p8")) {
  print(p8)
  cat("Plot 8: Group Heatmap displayed\n")
}

if (exists("p9")) {
  print(p9)
  cat("Plot 9: Indicator Comparison displayed\n")
}

if (exists("p10")) {
  print(p10)
  cat("Plot 10: Reliability Comparison displayed\n")
}

cat("\n=== All Available Plots Displayed ===\n")


cat("\n=== Creating Individual Indicator Analysis ===\n")