# Load required libraries without startup messages
suppressPackageStartupMessages({
  library(survival)
  library(readxl)
  library(dplyr)
  library(riskRegression)
  library(ggplot2)
  library(caret)  # For cross-validation folds
  library(pec)    # For cindex and calPlot functions
})

#--------------------------------------------------------------------------
# Step 1: Load the Data
#--------------------------------------------------------------------------
file_path <- "//Volumes/PhD Data/Final Datasets/Merged File_Final_Corrected_DVHcc_Prediction.xlsx"
df <- read_excel(file_path, sheet = "Factors")

#--------------------------------------------------------------------------
# Step 2: Prepare the data
#--------------------------------------------------------------------------
# Convert categorical variables to factors
categorical_vars <- c('Insurance_Type',
                      'Node',
                      'Periodontal_Grading', 
                      'Disease_Site_Merged_2')
df[categorical_vars] <- lapply(df[categorical_vars], as.factor)

# Define continuous variables
continuous_vars <- c('Age', 'Smoking_Pack_per_Year', 'Income_1000',
                     'Number_Teeth_after_Extraction', 'RT_Dose',
                      'D20')

# Survival-related variables
df$time  <- as.numeric(as.character(df$ClinRad_Time_Indicator_M...8))
df$delta <- as.numeric(as.character(df$ClinRad_M_Competing))  # 0: censored, 1: event, 2: competing risk

# Create the event indicator (ignoring competing events by treating them as censored)
df$delta <- ifelse(df$delta == 1, 1, 0)

# Cap time at 114
df <- df %>% mutate(time = ifelse(time > 114, 114, time))

#--------------------------------------------------------------------------
# Step 3: Set up repeated 5-fold cross-validation (5×5 = 25 splits)
#--------------------------------------------------------------------------
set.seed(123)  # For reproducibility
k_folds   <- 5
n_repeats <- 5

multi_folds <- createMultiFolds(factor(df$delta), k = k_folds, times = n_repeats)
all_indices <- seq_len(nrow(df))

#--------------------------------------------------------------------------
# Step 4: Initialize lists to store results
#--------------------------------------------------------------------------
auc_df_list       <- list()
brier_df_list     <- list()
ibs_values_list   <- list()
cindex_df_list    <- list()
cif_list          <- list()

# NEW: We also store each Cox model and each test index vector for retrieval
cox_models        <- list()
test_indices_list <- list()

# Determine the maximum observed time for global time horizon
all_test_indices <- lapply(multi_folds, function(train_i) setdiff(all_indices, train_i))
fold_max_times   <- sapply(all_test_indices, function(idx) max(df$time[idx]))
global_max_time  <- max(fold_max_times)
cat("Using global maximum time horizon:", global_max_time, "\n")

# Specify time horizons of interest
time_horizons <- c(36, 60, 113)

#--------------------------------------------------------------------------
# Helper function for computing mean + 95% CI
#--------------------------------------------------------------------------
compute_mean_ci <- function(x, ci_level = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) {
    return(c(Mean = mean(x), Lower = NA, Upper = NA))
  }
  m  <- mean(x)
  sd <- sd(x)
  se <- sd / sqrt(n)
  alpha <- 1 - ci_level
  z     <- qnorm(1 - alpha/2)  # ~1.96 for 95% CI
  lower <- m - z * se
  upper <- m + z * se
  c(Mean = m, Lower = lower, Upper = upper)
}

#--------------------------------------------------------------------------
# Step 5: Repeated CV Loop (25 unique splits)
#--------------------------------------------------------------------------
fold_counter <- 1
for (fold_name in names(multi_folds)) {
  cat("\nProcessing Split:", fold_name, "(", fold_counter, "of", length(multi_folds), ")\n")
  
  # Get training indices
  train_indices <- multi_folds[[fold_name]]
  # Invert to get test indices
  test_indices  <- setdiff(all_indices, train_indices)
  
  train_data <- df[train_indices, ]
  test_data  <- df[test_indices, ]
  
  # Build the formula for the Cox model
  predictor_vars <- paste(c(categorical_vars, continuous_vars), collapse = " + ")
  cox_formula <- as.formula(paste0("Surv(time, delta) ~ ", predictor_vars))
  
  # Fit the Cox model on training data
  cox_model <- coxph(cox_formula, data = train_data, x = TRUE, y = TRUE)
  
  # Store the model & test indices for later retrieval
  cox_models[[fold_counter]]        <- cox_model   # <-- Fixed: use `cox_model`, not `fg_model`
  test_indices_list[[fold_counter]] <- test_indices
  
  #----------------------------------------------------------------
  # Compute the predicted risk at each time for the test set
  #----------------------------------------------------------------
  time_points <- seq(0, 114, by = 1)
  # predictRisk() from {riskRegression} calculates survival probabilities for Cox 
  # so 1 - survival_probability = cumulative incidence (assuming only 1 type of event).
  test_surv <- predictRisk(cox_model, newdata = test_data, times = time_points)
  
  # Convert survival to cumulative incidence (CIF) = 1 - S(t)
  test_cif <- test_surv
  
  # Average CIF across test patients for each time point
  cif_df <- data.frame(
    Time = time_points,
    CIF  = colMeans(test_cif),
    Fold = fold_counter
  )
  cif_list[[fold_counter]] <- cif_df
  
  #----------------------------------------------------------------
  # Evaluate AUC, Brier scores, IBS, C-index
  #----------------------------------------------------------------
  for (time_horizon in time_horizons) {
    cat("  Evaluating up to time horizon:", time_horizon, "\n")
    
    time_seq <- seq(1, time_horizon, by = 1)
    
    scores <- Score(
      object       = list("Cox" = cox_model),  
      formula      = Surv(time, delta) ~ 1,
      data         = test_data,
      summary      = "ibs",   # to compute IBS
      times        = time_seq,
      metrics      = c("AUC", "Brier"),
      cens.model   = "cox",
      split.method = "none",
      verbose      = FALSE
    )
    
    # AUC
    auc_values <- scores$AUC$score %>%
      filter(model == "Cox") %>%
      select(times, AUC)
    auc_values$Fold         <- fold_counter
    auc_values$Time_Horizon <- time_horizon
    
    # Brier
    brier_values <- scores$Brier$score %>%
      filter(model == "Cox") %>%
      select(times, Brier)
    brier_values$Fold         <- fold_counter
    brier_values$Time_Horizon <- time_horizon
    
    # IBS for the last time in time_seq
    ibs_value <- scores$Brier$score %>%
      filter(model == "Cox", times == max(time_seq)) %>%
      pull(IBS)
    if (length(ibs_value) == 0) ibs_value <- NA
    
    ibs_value_df <- data.frame(
      Fold         = fold_counter,
      Time_Horizon = time_horizon,
      IBS          = ibs_value
    )
    
    # Store
    auc_df_list[[length(auc_df_list) + 1]]         <- auc_values
    brier_df_list[[length(brier_df_list) + 1]]     <- brier_values
    ibs_values_list[[length(ibs_values_list) + 1]] <- ibs_value_df
    
    # C-index
    cat("  Computing C-index at time horizon:", time_horizon, "\n")
    cindex_result <- pec::cindex(
      object      = list("Cox" = cox_model),
      formula     = Surv(time, delta) ~ 1,
      data        = test_data,
      eval.times  = time_horizon,
      cens.model  = "cox",
      splitMethod = "none",
      keep.matrix = FALSE,
      verbose     = FALSE
    )
    cindex_value <- cindex_result$AppCindex$Cox
    cindex_df <- data.frame(
      Fold         = fold_counter,
      Time_Horizon = time_horizon,
      Cindex       = cindex_value
    )
    cindex_df_list[[length(cindex_df_list) + 1]] <- cindex_df
  }
  
  fold_counter <- fold_counter + 1
}

#--------------------------------------------------------------------------
# Step 6: Combine results (all 25 splits) into data frames
#--------------------------------------------------------------------------
auc_df    <- do.call(rbind, auc_df_list)
brier_df  <- do.call(rbind, brier_df_list)
ibs_df    <- do.call(rbind, ibs_values_list)
cindex_df <- do.call(rbind, cindex_df_list)
cif_all   <- do.call(rbind, cif_list)

#--------------------------------------------------------------------------
# Step 7: Identify the best model based on C-index at max time (114)
#--------------------------------------------------------------------------
# Filter to time_horizon = 114
cindex_113 <- cindex_df %>%
  filter(Time_Horizon == 113)

# Find the row with the highest Cindex
best_fold <- cindex_113 %>%
  arrange(desc(Cindex)) %>%
  slice(1) %>%
  pull(Fold)

cat("\nBest fold (model) is Fold #:", best_fold, "with C-index at 113 =",
    cindex_113 %>% filter(Fold == best_fold) %>% pull(Cindex), "\n")

#--------------------------------------------------------------------------
# Step 8: Compute summary statistics (Mean + 95% CI) for AUC, Brier, IBS, C-index
#--------------------------------------------------------------------------
library(tidyr)

## AUC summary
auc_summary <- auc_df %>%
  group_by(times, Time_Horizon) %>%
  summarise(
    AUC_Mean    = compute_mean_ci(AUC)[1],
    AUC_LowerCI = compute_mean_ci(AUC)[2],
    AUC_UpperCI = compute_mean_ci(AUC)[3],
    .groups     = "drop"
  )

## Brier summary
brier_summary <- brier_df %>%
  group_by(times, Time_Horizon) %>%
  summarise(
    Brier_Mean    = compute_mean_ci(Brier)[1],
    Brier_LowerCI = compute_mean_ci(Brier)[2],
    Brier_UpperCI = compute_mean_ci(Brier)[3],
    .groups       = "drop"
  )

## IBS summary
ibs_summary <- ibs_df %>%
  group_by(Time_Horizon) %>%
  summarise(
    IBS_Mean    = compute_mean_ci(IBS)[1],
    IBS_LowerCI = compute_mean_ci(IBS)[2],
    IBS_UpperCI = compute_mean_ci(IBS)[3],
    .groups     = "drop"
  )

## C-index summary
cindex_summary <- cindex_df %>%
  group_by(Time_Horizon) %>%
  summarise(
    Cindex_Mean    = compute_mean_ci(Cindex)[1],
    Cindex_LowerCI = compute_mean_ci(Cindex)[2],
    Cindex_UpperCI = compute_mean_ci(Cindex)[3],
    .groups        = "drop"
  )

# Categorize time ranges in AUC and Brier data frames (for plotting)
auc_summary <- auc_summary %>%
  mutate(Time_Range = case_when(
    times > 0 & times <= 60 ~ "0-60",
    times > 60             ~ "60-Max"
  ))

brier_summary <- brier_summary %>%
  mutate(Time_Range = case_when(
    times > 0 & times <= 60 ~ "0-60",
    times > 60             ~ "60-Max"
  ))

#--------------------------------------------------------------------------
# Step 9: Plotting with 95% CI ribbons
#--------------------------------------------------------------------------

# A) Time-dependent AUC
ggplot(auc_summary, aes(x = times, y = AUC_Mean, color = Time_Range)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_ribbon(aes(ymin = AUC_LowerCI, ymax = AUC_UpperCI, fill = Time_Range),
              alpha = 0.2, color = NA) +
  ggtitle("Repeated 5×5 CV Time-dependent AUC with Time Ranges (95% CI)") +
  xlab("Time (months)") +
  ylab("AUC") +
  scale_x_continuous(
    breaks = seq(0, max(auc_summary$times, na.rm = TRUE), by = 10),
    limits = c(0, max(auc_summary$times, na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  scale_fill_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  theme_minimal()

# B) Time-dependent Brier score
ggplot(brier_summary, aes(x = times, y = Brier_Mean, color = Time_Range)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_ribbon(aes(ymin = Brier_LowerCI, ymax = Brier_UpperCI, fill = Time_Range),
              alpha = 0.2, color = NA) +
  ggtitle("Repeated 5×5 CV Brier Scores with Time Ranges (95% CI)") +
  xlab("Time (months)") +
  ylab("Brier Score") +
  scale_x_continuous(
    breaks = seq(0, max(brier_summary$times, na.rm = TRUE), by = 10),
    limits = c(0, max(brier_summary$times, na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  scale_fill_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  theme_minimal()

#--------------------------------------------------------------------------
# Print IBS summary (95% CI)
#--------------------------------------------------------------------------
cat("\nRepeated (5×5) Cross-Validated Integrated Brier Scores (IBS) for each time horizon:\n")
print(ibs_summary)

#--------------------------------------------------------------------------
# Print C-index summary (95% CI)
#--------------------------------------------------------------------------
cat("\nRepeated (5×5) Cross-Validated C-index for each time horizon:\n")
print(cindex_summary)

#--------------------------------------------------------------------------
# Combine & summarize CIF
#--------------------------------------------------------------------------
cif_summary <- cif_all %>%
  group_by(Time) %>%
  summarise(
    CIF_Mean    = compute_mean_ci(CIF)[1],
    CIF_LowerCI = compute_mean_ci(CIF)[2],
    CIF_UpperCI = compute_mean_ci(CIF)[3],
    .groups     = "drop"
  )

# Extract CIF values for months 60 and 114 (if present)
cif_at_60  <- cif_summary %>% filter(Time == 60)  %>% pull(CIF_Mean)
cif_at_114 <- cif_summary %>% filter(Time == 114) %>% pull(CIF_Mean)
cat("\nAverage CIF at month 60:", cif_at_60, "\n")
cat("Average CIF at month 114:", cif_at_114, "\n")

#--------------------------------------------------------------------------
# Plot the Average CIF with 95% CI
#--------------------------------------------------------------------------
ggplot(cif_summary, aes(x = Time, y = CIF_Mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = CIF_LowerCI, ymax = CIF_UpperCI),
              fill = "blue", alpha = 0.2) +
  geom_point(size = 1, color = "blue") +
  ggtitle("Repeated (5×5) CV Average CIF Across Folds (95% CI)") +
  xlab("Time (months)") +
  ylab("Average CIF") +
  scale_x_continuous(
    breaks = seq(0, 114, by = 10),
    limits = c(0, 114)
  ) +
  theme_minimal()


##############################################################################
# ADDITIONAL CODE: FOLD-SPECIFIC CALIBRATION AT TIME=60 AND TIME=114 USING calPlot()
##############################################################################
# We'll assume the rest of your script above is already run,
# so you have "fg_models", "test_indices_list", etc.

library(dplyr)
library(ggplot2)

# Define the calibration horizons
time_horizons <- c(60, 114)

# We'll store fold-level calibration in a list
# It's a list of lists: outer list by time_horizon, inner list by fold
cal_list <- vector("list", length(time_horizons))
names(cal_list) <- paste0("Time_", time_horizons)

# Function to compute mean and 95% CI
compute_mean_ci <- function(x, conf.level = 0.95) {
  mean_val <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(x))
  alpha <- 1 - conf.level
  lower <- mean_val - qnorm(1 - alpha/2) * se
  upper <- mean_val + qnorm(1 - alpha/2) * se
  return(c(mean_val, lower, upper))
}

# Iterate over each time horizon
for (t in seq_along(time_horizons)) {
  time_horizon <- time_horizons[t]
  
  # Initialize list for this time_horizon
  cal_list[[t]] <- vector("list", length(cox_models))
  
  # For each fold
  for (f in seq_along(cox_models)) {
    
    # Get the test indices and data for this fold
    test_idx  <- test_indices_list[[f]]
    test_data <- df[test_idx, ]
    
    # Use calPlot() to get calibration data at `time_horizon`
    cal_result <- calPlot(
      object    = cox_models[[f]],
      time      = time_horizon,
      formula   = Hist(time, delta) ~ 1,
      data      = test_data,
      cause     = 1,
      method    = "quantile",
      q         = 10,
      pseudo    = TRUE,  # typical for survival/competing-risks calibration
      returnData= TRUE
    )
    
    # Extract the data for our single model from calPlot's output
    # Usually located in cal_result$plotFrames$Model.1
    cal_data <- as.data.frame(cal_result$plotFrames$Model.1)
    
    # Rename/organize for clarity
    fold_df <- data.frame(
      predicted_risk = 1 - cal_data$Pred,
      observed_risk  = 1 - cal_data$Obs,
      bin            = seq_len(nrow(cal_data)),  # each row is one bin
      Time_Horizon   = time_horizon,
      Fold           = f
    )
    
    # Store in the list for this time_horizon
    cal_list[[t]][[f]] <- fold_df
  }
}

# Combine all calibration data into one data frame
cal_all <- bind_rows(
  lapply(seq_along(cal_list), function(t) {
    bind_rows(cal_list[[t]]) %>% mutate(Time_Horizon = time_horizons[t])
  })
)

# Now, cal_all has columns:
# predicted_risk, observed_risk, bin, Time_Horizon, Fold

# Compute summary statistics: mean and 95% CI per bin and time_horizon
cal_summary <- cal_all %>%
  group_by(Time_Horizon, bin) %>%
  summarise(
    PredRisk_Mean    = compute_mean_ci(predicted_risk)[1],
    PredRisk_LowerCI = compute_mean_ci(predicted_risk)[2],
    PredRisk_UpperCI = compute_mean_ci(predicted_risk)[3],
    
    ObsRisk_Mean     = compute_mean_ci(observed_risk)[1],
    ObsRisk_LowerCI  = compute_mean_ci(observed_risk)[2],
    ObsRisk_UpperCI  = compute_mean_ci(observed_risk)[3],
    .groups = "drop"
  )

# Convert Time_Horizon to a factor for plotting
cal_summary$Time_Horizon <- factor(cal_summary$Time_Horizon, 
                                   levels = time_horizons,
                                   labels = paste0(time_horizons, " Months"))

#--------------------------------------------------------------------------
# 1. Combined Plot: Both Time Horizons on One Figure with Different Colors
#--------------------------------------------------------------------------
combined_plot <- ggplot(cal_summary, aes(x = PredRisk_Mean, y = ObsRisk_Mean, color = Time_Horizon)) +
  # Points + lines for each time horizon
  geom_point(size = 2) +
  geom_line() +
  
  # Vertical error bars for observed risk
  geom_errorbar(
    aes(ymin = ObsRisk_LowerCI, ymax = ObsRisk_UpperCI),
    width = 0.01, alpha = 0.6
  ) +
  
  # (Optional) horizontal error bars for predicted risk
  # geom_errorbarh(
  #   aes(xmin = PredRisk_LowerCI, xmax = PredRisk_UpperCI),
  #   height = 0.01, alpha = 0.6
  # ) +
  
  # Perfect calibration lines for each time horizon
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  
  coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4)) +  # Adjust as needed
  xlab("Mean Predicted Probability (decile)") +
  ylab("Mean Observed Probability (decile)") +
  ggtitle("Calibration at Multiple Time Horizons (Fold-by-Fold Deciles)") +
  theme_minimal() +
  scale_color_manual(values = c("60 Months" = "blue", "114 Months" = "red")) +
  theme(legend.title = element_blank())

# Display the combined plot
print(combined_plot)

#--------------------------------------------------------------------------
# 2. Separate Plots: One Figure per Time Horizon
#--------------------------------------------------------------------------

# Define colors for separate plots
time_colors <- c("60 Months" = "blue", "114 Months" = "red")

# Loop over each time horizon to create separate plots
for (t in levels(cal_summary$Time_Horizon)) {
  
  # Filter data for the current time horizon
  cal_summary_subset <- cal_summary %>% filter(Time_Horizon == t)
  
  # Create the plot
  separate_plot <- ggplot(cal_summary_subset, aes(x = PredRisk_Mean, y = ObsRisk_Mean)) +
    # Points + line
    geom_point(size = 2, color = time_colors[t]) +
    geom_line(color = time_colors[t]) +
    
    # Vertical error bars for observed risk
    geom_errorbar(
      aes(ymin = ObsRisk_LowerCI, ymax = ObsRisk_UpperCI),
      width = 0.01, color = time_colors[t], alpha = 0.6
    ) +
    
    # (Optional) horizontal error bars for predicted risk
    # geom_errorbarh(
    #   aes(xmin = PredRisk_LowerCI, xmax = PredRisk_UpperCI),
    #   height = 0.01, color = time_colors[t], alpha = 0.6
    # ) +
    
    # Perfect calibration line (45°)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    
    coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4)) +  # Adjust as needed
    xlab("Mean Predicted Probability (decile)") +
    ylab("Mean Observed Probability (decile)") +
    ggtitle(paste0("Calibration at ", t)) +
    theme_minimal()
  
  # Display the separate plot
  print(separate_plot)
}

#--------------------------------------------------------------------------
# Optional: Saving the Plots
#--------------------------------------------------------------------------
# Uncomment and modify the following lines if you wish to save the plots

# Save Combined Plot
# ggsave("Calibration_Combined_Plot.png", plot = combined_plot, width = 8, height = 6, dpi = 300)

# Save Separate Plots
# for (t in levels(cal_summary$Time_Horizon)) {
#   # Filter data for the current time horizon
#   cal_summary_subset <- cal_summary %>% filter(Time_Horizon == t)
  
#   # Create the plot
#   separate_plot <- ggplot(cal_summary_subset, aes(x = PredRisk_Mean, y = ObsRisk_Mean)) +
#     geom_point(size = 2, color = time_colors[t]) +
#     geom_line(color = time_colors[t]) +
#     geom_errorbar(
#       aes(ymin = ObsRisk_LowerCI, ymax = ObsRisk_UpperCI),
#       width = 0.01, color = time_colors[t], alpha = 0.6
#     ) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     coord_cartesian(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
#     xlab("Mean Predicted Probability (decile)") +
#     ylab("Mean Observed Probability (decile)") +
#     ggtitle(paste0("Calibration at ", t)) +
#     theme_minimal()
  
#   # Define the filename
#   filename <- paste0("Calibration_", gsub(" ", "_", t), ".png")
  
#   # Save the plot
#   ggsave(filename, plot = separate_plot, width = 8, height = 6, dpi = 300)
# }
