# Machine Learning Framework using Cox Proportioal Hazards 

## Overview

This repository contains an R script designed to perform survival model development and evaluation using Cox Proportioal Hazards (CPH). The script facilitates data preprocessing, model training using the CPH approach, and extensive performance evaluation through various metrics and visualizations. It employs repeated 5-fold cross-validation to ensure robust and reliable assessment of the predictive models.

## Features

- **Data Loading and Preparation:**
  - Reads data from an Excel file.
  - Converts specified categorical variables to factors.
  - Scales continuous variables as needed.
  - Handles missing values and caps time-related data to a specified maximum.

- **Cross-Validation Setup:**
  - Implements repeated 5-fold cross-validation (totaling 25 unique splits) to validate model performance consistently.

- **Model Training:**
  - Fits Fine-Gray models to account for competing risks in survival analysis.
  
- **Performance Evaluation:**
  - Calculates key performance metrics including time-dependent Area Under the Curve (AUC), Brier Scores, Integrated Brier Score (IBS), and Concordance Index (C-index).

- **Visualization:**
  - Generates time-dependent AUC and Brier Score plots with 95% confidence intervals.
  - Creates Calibration plots at specified time horizons to assess the agreement between predicted and observed risks.
  - Plots the Cumulative Incidence Function (CIF) across different folds.

- **Calibration Analysis:**
  - Performs fold-specific calibration at multiple time horizons using decile-based calibration plots.

- **Output Management:**
  - Organizes results and plots systematically, with options to save visualizations for further review.

## Installations

### Prerequisites

Ensure that you have **R** installed on your system. You can download R from the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/).

### Required R Packages

The script relies on the following R packages. You can install them using the commands provided below.

- **survival**: Analysis of survival data.
- **readxl**: Reading Excel files.
- **dplyr**: Data manipulation.
- **riskRegression**: Regression models for risk analysis.
- **ggplot2**: Data visualization.
- **caret**: Classification and Regression Training.
- **pec**: Prediction Error Curves.
- **tidyr**: Tidy Messy Data.

#### Installation Instructions

Open your R console or RStudio and run the following command to install all required packages:

```r
install.packages(c(
  "survival",
  "readxl",
  "dplyr",
  "riskRegression",
  "ggplot2",
  "caret",
  "pec",
  "tidyr"
))
```
