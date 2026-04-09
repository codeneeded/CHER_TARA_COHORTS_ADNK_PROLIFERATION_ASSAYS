################################################################################
# ADNK Assay QC and Normalization Script
# 
# Purpose: Quality control and normalize ADNK flow cytometry data
# Input: Excel file with batch/antigen sheets containing FlowJo exported data
# Output: QC'd, normalized dataset ready for analysis
#
# Cohorts: CHER (longitudinal, treatment interruption)
#          TARA (cross-sectional, sup vs unsup)
#
# Author: Generated for CHER/TARA ADNK analysis
# Date: April 2026
################################################################################

# =============================================================================
# SETUP
# =============================================================================

required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "writexl", 
                       "stringr", "patchwork", "scales")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(stringr)
library(patchwork)
library(scales)

theme_set(theme_minimal(base_size = 12) +
            theme(
              panel.grid.minor = element_blank(),
              strip.background = element_rect(fill = "#f0f0f0", color = NA),
              legend.position = "bottom"
            ))

# =============================================================================
# CONFIGURATION - EDIT THESE PATHS
# =============================================================================

BASE_DIR <- "/home/akshay-iyer/Documents/CHER_TARA_COHORTS_ADNK_PROLIFERATION_ASSAYS"
INPUT_FILE <- file.path(BASE_DIR, "ADNK_Data_CHER_and_TARA_with_metadata_040826.xlsx")

# Output directory structure
OUTPUT_DIR <- file.path(BASE_DIR, "ADNK_QC_Output")

DIRS <- list(
  root = OUTPUT_DIR,
  qc_reports = file.path(OUTPUT_DIR, "01_QC_Reports"),
  qc_plots = file.path(OUTPUT_DIR, "02_QC_Plots"),
  batch_correction = file.path(OUTPUT_DIR, "03_Batch_Correction"),
  normalized_data = file.path(OUTPUT_DIR, "04_Normalized_Data"),
  metadata = file.path(OUTPUT_DIR, "05_Cleaned_Metadata"),
  summary = file.path(OUTPUT_DIR, "06_Summary")
)

for (dir_path in DIRS) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created: ", dir_path)
  }
}

# =============================================================================
# CONTROL DEFINITIONS
# =============================================================================

CONTROLS <- list(
  PBS = c("PBS control", "PBS control_1", "PBS control_2", "PBS control_3", "PBS_control"),
  HIV_NEG = c("HIV neg plasma"),
  N6_BNAB = c("N6 bnab"),
  STAINING = c("NK stained", "Unstained")
)

ALL_CONTROLS <- unlist(CONTROLS, use.names = FALSE)

# =============================================================================
# FUNCTIONAL MARKER DEFINITIONS
# =============================================================================

FUNCTIONAL_MARKERS <- c(
  "NK 56bright/CD107a | Freq. of Parent",
  "NK 56bright/FasL | Freq. of Parent",
  "NK 56bright/IFNg | Freq. of Parent",
  "NK 56bright/MIP1b+ | Freq. of Parent",
  "NK CD56dim/CD107a | Freq. of Parent",
  "NK CD56dim/FasL | Freq. of Parent",
  "NK CD56dim/IFNgamma | Freq. of Parent",
  "NK CD56dim/MIP1b+ | Freq. of Parent",
  "non-NK/CD107a | Freq. of Parent",
  "non-NK/FasL | Freq. of Parent",
  "non-NK/IFNg | Freq. of Parent",
  "non-NK/MIP1b+ | Freq. of Parent"
)

MARKER_SHORT_NAMES <- c(
  "NK 56bright/CD107a | Freq. of Parent" = "CD56br_CD107a",
  "NK 56bright/FasL | Freq. of Parent" = "CD56br_FasL",
  "NK 56bright/IFNg | Freq. of Parent" = "CD56br_IFNg",
  "NK 56bright/MIP1b+ | Freq. of Parent" = "CD56br_MIP1b",
  "NK CD56dim/CD107a | Freq. of Parent" = "CD56dim_CD107a",
  "NK CD56dim/FasL | Freq. of Parent" = "CD56dim_FasL",
  "NK CD56dim/IFNgamma | Freq. of Parent" = "CD56dim_IFNg",
  "NK CD56dim/MIP1b+ | Freq. of Parent" = "CD56dim_MIP1b",
  "non-NK/CD107a | Freq. of Parent" = "nonNK_CD107a",
  "non-NK/FasL | Freq. of Parent" = "nonNK_FasL",
  "non-NK/IFNg | Freq. of Parent" = "nonNK_IFNg",
  "non-NK/MIP1b+ | Freq. of Parent" = "nonNK_MIP1b"
)

PRIMARY_MARKERS <- c(
  "NK CD56dim/CD107a | Freq. of Parent",
  "NK CD56dim/IFNgamma | Freq. of Parent"
)

# =============================================================================
# QC THRESHOLDS
# =============================================================================

QC_THRESHOLDS <- list(
  MIN_N6_CD107A = 0.05,
  MAX_PBS_BACKGROUND = 0.15,
  MIN_EVENT_COUNT = 5000,
  MIN_LYMPHOCYTE_FREQ = 0.50,
  MAX_BATCH_FOLD_DIFF = 3.0
)

# =============================================================================
# FUNCTION: Load and combine all sheets
# =============================================================================

load_adnk_data <- function(filepath) {
  
  sheets <- excel_sheets(filepath)
  message("Found sheets (plates): ", paste(sheets, collapse = ", "))
  
  all_data <- lapply(sheets, function(sheet) {
    df <- read_excel(filepath, sheet = sheet)
    names(df) <- gsub("Well Position", "Well", names(df))
    
    # Keep original timepoint column for processing
    if ("Timepoint relative to phase" %in% names(df)) {
      df$Timepoint_Raw <- df$`Timepoint relative to phase`
    }
    
    # Convert all numeric-ish columns to character first to avoid bind_rows type conflicts
    # This handles cases where some sheets have "N/A" or "-" in numeric columns
    numeric_cols <- names(df)[sapply(df, function(x) is.numeric(x) | is.character(x))]
    for (col in numeric_cols) {
      if (grepl("Freq|Mean|Count|Median", col)) {
        df[[col]] <- as.character(df[[col]])
      }
    }
    
    df$Sheet <- sheet
    df$Plate <- sheet
    df$Batch <- ifelse(grepl("batch 1", sheet, ignore.case = TRUE), "Batch 1", "Batch 2")
    df$Antigen <- ifelse(grepl("p24", sheet, ignore.case = TRUE), "p24", "gp120")
    return(df)
  })
  
  combined <- bind_rows(all_data)
  
  # Now convert measurement columns back to numeric (non-numeric values become NA)
  measurement_cols <- names(combined)[grepl("Freq|Mean|Count|Median", names(combined))]
  for (col in measurement_cols) {
    combined[[col]] <- suppressWarnings(as.numeric(combined[[col]]))
  }
  
  combined$PID <- as.character(combined$PID)
  combined$Cohort <- ifelse(grepl("^CP", combined$PID), "TARA", "CHER")
  combined$IsControl <- combined$PID %in% ALL_CONTROLS
  combined$ControlType <- case_when(
    combined$PID %in% CONTROLS$PBS ~ "PBS",
    combined$PID %in% CONTROLS$HIV_NEG ~ "HIV_neg",
    combined$PID %in% CONTROLS$N6_BNAB ~ "N6_bnAb",
    combined$PID %in% CONTROLS$STAINING ~ "Staining",
    TRUE ~ NA_character_
  )
  
  # -------------------------------------------------------------------------
  # TIMEPOINT HANDLING
  # -------------------------------------------------------------------------
  # TARA: "Timepoint relative to phase" is actually Age (months) - already numeric
  # CHER: Extract numeric value from formats like "0 B", "12 WK", "14 VST", "0WK"
  # -------------------------------------------------------------------------
  
  combined <- combined %>%
    mutate(
      # For TARA: This column represents Age in months
      `Age (months)` = ifelse(
        Cohort == "TARA" & !IsControl,
        as.numeric(gsub("[^0-9.]", "", as.character(Timepoint_Raw))),
        NA_real_
      ),
      
      # For CHER: Clean timepoint to just the number
      # Handles: "0 B", "0 WK", "0WK", "12 WK", "14 VST", "32 WK", etc.
      Timepoint_Numeric = as.numeric(gsub("[^0-9.]", "", as.character(Timepoint_Raw))),
      
      # Keep original for reference
      Timepoint_Original = Timepoint_Raw
    )
  
  # Log timepoint cleaning
  message("\nTimepoint cleaning summary:")
  
  cher_tp <- combined %>% 
    filter(Cohort == "CHER", !IsControl) %>%
    select(Timepoint_Original, Timepoint_Numeric) %>%
    distinct() %>%
    arrange(Timepoint_Numeric)
  
  message("  CHER timepoints (original -> numeric):")
  for (i in 1:min(nrow(cher_tp), 10)) {
    message(sprintf("    '%s' -> %s", cher_tp$Timepoint_Original[i], cher_tp$Timepoint_Numeric[i]))
  }
  if (nrow(cher_tp) > 10) message("    ... and ", nrow(cher_tp) - 10, " more")
  
  tara_ages <- combined %>%
    filter(Cohort == "TARA", !IsControl) %>%
    pull(`Age (months)`) %>%
    unique() %>%
    sort()
  
  message("  TARA Age (months) values: ", paste(tara_ages, collapse = ", "))
  
  message("\nLoaded ", nrow(combined), " rows")
  message("  - Samples: ", sum(!combined$IsControl))
  message("  - Controls: ", sum(combined$IsControl))
  message("  - CHER: ", sum(!combined$IsControl & combined$Cohort == "CHER"), " samples")
  message("  - TARA: ", sum(!combined$IsControl & combined$Cohort == "TARA"), " samples")
  
  return(combined)
}

# =============================================================================
# FUNCTION: Calculate control statistics
# =============================================================================

calculate_control_stats <- function(data) {
  
  # Each sheet = one plate, each plate has its own controls
  # Sheets are: batch 1 p24, batch 1 gp120, batch 2 p24, batch 2 gp120
  # So grouping by (Batch, Antigen) or by Sheet gives the same result
  # We group by Sheet explicitly to be clear that each plate uses its own controls
  
  controls <- data %>%
    filter(IsControl) %>%
    filter(ControlType %in% c("PBS", "HIV_neg", "N6_bnAb"))
  
  markers_present <- FUNCTIONAL_MARKERS[FUNCTIONAL_MARKERS %in% names(data)]
  
  # Group by Sheet (plate) to ensure each plate uses its own controls
  control_stats <- controls %>%
    group_by(Sheet, Batch, Antigen, ControlType) %>%
    summarise(
      across(all_of(markers_present),
             list(mean = ~mean(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE),
                  n = ~sum(!is.na(.x))),
             .names = "{.col}___{.fn}"),
      .groups = "drop"
    )
  
  message("Control stats calculated per plate (sheet):")
  message("  ", paste(unique(control_stats$Sheet), collapse = ", "))
  
  return(control_stats)
}

# =============================================================================
# FUNCTION: QC Check 1 - Control hierarchy
# =============================================================================

qc_control_hierarchy <- function(data, control_stats, output_dir) {
  
  message("\n", strrep("=", 60))
  message("QC CHECK 1: Control Hierarchy Validation")
  message(strrep("=", 60))
  
  marker <- "NK CD56dim/CD107a | Freq. of Parent"
  marker_col <- paste0(marker, "___mean")
  
  hierarchy_check <- control_stats %>%
    select(Batch, Antigen, ControlType, all_of(marker_col)) %>%
    pivot_wider(names_from = ControlType, values_from = all_of(marker_col)) %>%
    mutate(
      PBS_lt_N6 = PBS < N6_bnAb,
      PBS_lt_HIVneg = PBS <= HIV_neg,
      Hierarchy_OK = PBS_lt_N6 & PBS_lt_HIVneg
    )
  
  results_text <- character()
  for (i in 1:nrow(hierarchy_check)) {
    row <- hierarchy_check[i, ]
    status <- ifelse(row$Hierarchy_OK, "PASS", "FAIL")
    line <- sprintf("%s %s: PBS=%.4f, HIV_neg=%.4f, N6=%.4f  [%s]",
                    row$Batch, row$Antigen, row$PBS, row$HIV_neg, row$N6_bnAb, status)
    message(line)
    results_text <- c(results_text, line)
  }
  
  writeLines(results_text, file.path(output_dir, "QC1_Control_Hierarchy.txt"))
  write.csv(hierarchy_check, file.path(output_dir, "QC1_Control_Hierarchy.csv"), row.names = FALSE)
  
  return(hierarchy_check)
}

# =============================================================================
# FUNCTION: QC Check 2 - N6 performance
# =============================================================================

qc_positive_control <- function(data, control_stats, output_dir) {
  
  message("\n", strrep("=", 60))
  message("QC CHECK 2: N6 Positive Control Performance")
  message(strrep("=", 60))
  
  marker_col <- "NK CD56dim/CD107a | Freq. of Parent___mean"
  
  n6_check <- control_stats %>%
    filter(ControlType == "N6_bnAb") %>%
    select(Batch, Antigen, all_of(marker_col)) %>%
    rename(N6_CD107a = all_of(marker_col)) %>%
    mutate(
      N6_OK = N6_CD107a >= QC_THRESHOLDS$MIN_N6_CD107A,
      Status = ifelse(N6_OK, "PASS", "LOW_RESPONSE")
    )
  
  for (i in 1:nrow(n6_check)) {
    row <- n6_check[i, ]
    message(sprintf("%s %s: N6 CD107a = %.2f%%  [%s]",
                    row$Batch, row$Antigen, row$N6_CD107a * 100, row$Status))
  }
  
  write.csv(n6_check, file.path(output_dir, "QC2_N6_Performance.csv"), row.names = FALSE)
  
  return(n6_check)
}

# =============================================================================
# FUNCTION: QC Check 3 - Batch effects
# =============================================================================

qc_batch_effects <- function(data, control_stats, output_dir) {
  
  message("\n", strrep("=", 60))
  message("QC CHECK 3: Batch Effects (PBS Background)")
  message(strrep("=", 60))
  
  batch_results <- list()
  
  for (marker in PRIMARY_MARKERS) {
    marker_col <- paste0(marker, "___mean")
    marker_short <- MARKER_SHORT_NAMES[marker]
    
    if (!marker_col %in% names(control_stats)) next
    
    batch_compare <- control_stats %>%
      filter(ControlType == "PBS") %>%
      select(Batch, Antigen, all_of(marker_col)) %>%
      pivot_wider(names_from = Batch, values_from = all_of(marker_col))
    
    batch_compare <- batch_compare %>%
      mutate(
        Fold_Diff = pmax(`Batch 1`, `Batch 2`) / pmin(`Batch 1`, `Batch 2`),
        Batch_OK = Fold_Diff < QC_THRESHOLDS$MAX_BATCH_FOLD_DIFF,
        Marker = marker_short
      )
    
    message(sprintf("\n%s:", marker_short))
    for (i in 1:nrow(batch_compare)) {
      row <- batch_compare[i, ]
      status <- ifelse(row$Batch_OK, "OK", "BATCH_EFFECT")
      message(sprintf("  %s: Batch1=%.4f, Batch2=%.4f, Fold=%.2f [%s]",
                      row$Antigen, row$`Batch 1`, row$`Batch 2`, row$Fold_Diff, status))
    }
    
    batch_results[[marker_short]] <- batch_compare
  }
  
  batch_df <- bind_rows(batch_results)
  write.csv(batch_df, file.path(output_dir, "QC3_Batch_Effects.csv"), row.names = FALSE)
  
  return(batch_df)
}

# =============================================================================
# FUNCTION: QC Check 4 - Sample-level
# =============================================================================

qc_sample_level <- function(data, output_dir) {
  
  message("\n", strrep("=", 60))
  message("QC CHECK 4: Sample-Level Quality")
  message(strrep("=", 60))
  
  samples <- data %>% filter(!IsControl)
  
  sample_qc <- samples %>%
    mutate(
      QC_LowEvents = Count < QC_THRESHOLDS$MIN_EVENT_COUNT,
      QC_LowLymphocytes = `Lymphocytes | Freq. of Parent` < QC_THRESHOLDS$MIN_LYMPHOCYTE_FREQ,
      QC_Pass = !QC_LowEvents & !QC_LowLymphocytes
    )
  
  qc_summary <- sample_qc %>%
    group_by(Cohort, Batch, Antigen) %>%
    summarise(
      Total = n(),
      Passed = sum(QC_Pass),
      Failed_Events = sum(QC_LowEvents),
      Failed_Lymph = sum(QC_LowLymphocytes),
      Pass_Rate = mean(QC_Pass) * 100,
      .groups = "drop"
    )
  
  message(sprintf("\nOverall: %d/%d samples pass QC (%.1f%%)",
                  sum(sample_qc$QC_Pass), nrow(sample_qc), mean(sample_qc$QC_Pass) * 100))
  
  write.csv(qc_summary, file.path(output_dir, "QC4_Sample_Summary.csv"), row.names = FALSE)
  
  failed <- sample_qc %>% filter(!QC_Pass)
  if (nrow(failed) > 0) {
    write.csv(failed %>% select(PID, Timepoint_Original, Timepoint_Numeric, Cohort, Batch, Antigen, Count, 
                                `Lymphocytes | Freq. of Parent`, QC_LowEvents, QC_LowLymphocytes),
              file.path(output_dir, "QC4_Failed_Samples.csv"), row.names = FALSE)
  }
  
  return(sample_qc)
}

# =============================================================================
# FUNCTION: Background subtraction
# =============================================================================

subtract_background <- function(data, control_stats) {
  
  message("\n", strrep("=", 60))
  message("NORMALIZATION: Background Subtraction")
  message(strrep("=", 60))
  message("Using PBS control from each plate (sheet) for that plate's samples")
  
  markers_present <- FUNCTIONAL_MARKERS[FUNCTIONAL_MARKERS %in% names(data)]
  
  # Get PBS means per plate (Sheet)
  pbs_means <- control_stats %>%
    filter(ControlType == "PBS") %>%
    select(Sheet, Batch, Antigen, ends_with("___mean")) %>%
    rename_with(~gsub("___mean", "_PBS", .x), ends_with("___mean"))
  
  samples <- data %>% filter(!IsControl)
  
  # Join by Sheet to ensure each plate uses its own PBS control
  samples <- samples %>% left_join(pbs_means, by = c("Sheet", "Batch", "Antigen"))
  
  for (marker in markers_present) {
    pbs_col <- paste0(marker, "_PBS")
    bgsub_col <- paste0(marker, "_BGsub")
    
    if (pbs_col %in% names(samples)) {
      samples[[bgsub_col]] <- samples[[marker]] - samples[[pbs_col]]
      samples[[bgsub_col]] <- pmax(samples[[bgsub_col]], 0)
    }
  }
  
  message("Created background-subtracted columns (*_BGsub)")
  return(samples)
}

# =============================================================================
# FUNCTION: N6 normalization
# =============================================================================

normalize_to_n6 <- function(data, control_stats) {
  
  message("\n", strrep("=", 60))
  message("NORMALIZATION: % of N6 Maximum Response")
  message(strrep("=", 60))
  message("Using N6 and PBS controls from each plate (sheet) for that plate's samples")
  
  markers_present <- FUNCTIONAL_MARKERS[FUNCTIONAL_MARKERS %in% names(data)]
  
  # Get N6 means per plate (Sheet)
  n6_means <- control_stats %>%
    filter(ControlType == "N6_bnAb") %>%
    select(Sheet, Batch, Antigen, ends_with("___mean")) %>%
    rename_with(~gsub("___mean", "_N6", .x), ends_with("___mean"))
  
  # Data already has _PBS columns from background subtraction, so just join N6
  samples <- data %>%
    left_join(n6_means, by = c("Sheet", "Batch", "Antigen"))
  
  # Create normalized columns
  for (marker in markers_present) {
    n6_col <- paste0(marker, "_N6")
    pbs_col <- paste0(marker, "_PBS")  # Already exists from BG subtraction
    norm_col <- paste0(marker, "_NormN6")
    
    if (n6_col %in% names(samples) & pbs_col %in% names(samples)) {
      denominator <- samples[[n6_col]] - samples[[pbs_col]]
      denominator <- ifelse(denominator < 0.001, NA, denominator)
      samples[[norm_col]] <- ((samples[[marker]] - samples[[pbs_col]]) / denominator) * 100
      samples[[norm_col]] <- pmin(pmax(samples[[norm_col]], 0), 500)
    } else {
      message("  WARNING: Missing columns for ", marker)
      message("    N6 col exists: ", n6_col %in% names(samples))
      message("    PBS col exists: ", pbs_col %in% names(samples))
    }
  }
  
  message("Created N6-normalized columns (*_NormN6)")
  return(samples)
}

# =============================================================================
# FUNCTION: QC Plots
# =============================================================================

generate_qc_plots <- function(data, control_stats, output_dir) {
  
  message("\n", strrep("=", 60))
  message("GENERATING QC PLOTS")
  message(strrep("=", 60))
  
  control_data <- data %>%
    filter(IsControl, ControlType %in% c("PBS", "HIV_neg", "N6_bnAb")) %>%
    mutate(ControlType = factor(ControlType, levels = c("PBS", "HIV_neg", "N6_bnAb")))
  
  samples <- data %>% filter(!IsControl)
  marker <- "NK CD56dim/CD107a | Freq. of Parent"
  
  # Plot 1: Control performance
  p1 <- ggplot(control_data, aes(x = ControlType, y = .data[[marker]] * 100, fill = ControlType)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
    facet_grid(Batch ~ Antigen) +
    scale_fill_manual(values = c("PBS" = "#888888", "HIV_neg" = "#5DCAA5", "N6_bnAb" = "#7F77DD")) +
    labs(title = "Control Performance by Batch and Antigen",
         subtitle = "CD56dim CD107a - Expected: PBS < HIV neg < N6",
         y = "CD107a+ (% of parent)", x = "") +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "01_Control_Performance.png"), p1, width = 10, height = 8, dpi = 150)
  
  # Plot 2: Batch distribution RAW
  p2 <- ggplot(samples, aes(x = Batch, y = .data[[marker]] * 100, fill = Batch)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
    facet_grid(Cohort ~ Antigen) +
    scale_fill_manual(values = c("Batch 1" = "#3B8BD4", "Batch 2" = "#D85A30")) +
    labs(title = "RAW Sample Distribution by Batch (Before Normalization)",
         subtitle = "CD56dim CD107a - Note batch differences",
         y = "CD107a+ (% of parent)", x = "") +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "02_Batch_Distribution_RAW.png"), p2, width = 10, height = 8, dpi = 150)
  
  # Plot 3: Event counts
  p3 <- ggplot(samples, aes(x = interaction(Cohort, Batch), y = Count, fill = Cohort)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = QC_THRESHOLDS$MIN_EVENT_COUNT, linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = c("CHER" = "#7F77DD", "TARA" = "#D4537E")) +
    labs(title = "Event Counts by Cohort and Batch",
         subtitle = sprintf("Red line = QC threshold (%d events)", QC_THRESHOLDS$MIN_EVENT_COUNT),
         y = "Event Count", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "03_Event_Counts.png"), p3, width = 8, height = 6, dpi = 150)
  
  # Plot 4: Lymphocyte frequency
  p4 <- ggplot(samples, aes(x = interaction(Cohort, Batch), 
                            y = `Lymphocytes | Freq. of Parent` * 100, fill = Cohort)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = QC_THRESHOLDS$MIN_LYMPHOCYTE_FREQ * 100, 
               linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = c("CHER" = "#7F77DD", "TARA" = "#D4537E")) +
    labs(title = "Lymphocyte Frequency by Cohort and Batch",
         subtitle = sprintf("Red line = QC threshold (%.0f%%)", QC_THRESHOLDS$MIN_LYMPHOCYTE_FREQ * 100),
         y = "Lymphocytes (% of parent)", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "04_Lymphocyte_Freq.png"), p4, width = 8, height = 6, dpi = 150)
  
  message("Saved QC plots to: ", output_dir)
}

# =============================================================================
# FUNCTION: Batch correction confirmation plots
# =============================================================================

generate_batch_correction_plots <- function(raw_data, normalized_data, output_dir) {
  
  message("\n", strrep("=", 60))
  message("GENERATING BATCH CORRECTION CONFIRMATION PLOTS")
  message(strrep("=", 60))
  
  marker_raw <- "NK CD56dim/CD107a | Freq. of Parent"
  
  # Find the actual column names (they may have been modified)
  bgsub_cols <- names(normalized_data)[grepl("CD56dim.*CD107a.*BGsub", names(normalized_data))]
  norm_cols <- names(normalized_data)[grepl("CD56dim.*CD107a.*NormN6", names(normalized_data))]
  
  if (length(bgsub_cols) == 0 | length(norm_cols) == 0) {
    message("WARNING: Could not find normalized columns. Available columns:")
    message(paste(names(normalized_data)[grepl("CD107a", names(normalized_data))], collapse = "\n"))
    return(NULL)
  }
  
  marker_bgsub <- bgsub_cols[1]
  marker_norm <- norm_cols[1]
  
  message("Using columns:")
  message("  Raw: ", marker_raw)
  message("  BG-sub: ", marker_bgsub)
  message("  Normalized: ", marker_norm)
  
  samples_raw <- raw_data %>% filter(!IsControl)
  
  plot_data <- bind_rows(
    samples_raw %>%
      select(PID, Sheet, Batch, Antigen, Cohort, all_of(marker_raw)) %>%
      mutate(Method = "1_Raw (%)", Value = .data[[marker_raw]] * 100),
    normalized_data %>%
      select(PID, Sheet, Batch, Antigen, Cohort, all_of(marker_bgsub)) %>%
      mutate(Method = "2_BG_Subtracted (%)", Value = .data[[marker_bgsub]] * 100),
    normalized_data %>%
      select(PID, Sheet, Batch, Antigen, Cohort, all_of(marker_norm)) %>%
      mutate(Method = "3_N6_Normalized (% of max)", Value = .data[[marker_norm]])
  )
  
  # Main comparison plot
  p_compare <- ggplot(plot_data, aes(x = Batch, y = Value, fill = Batch)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
    facet_grid(Method ~ Antigen, scales = "free_y") +
    scale_fill_manual(values = c("Batch 1" = "#3B8BD4", "Batch 2" = "#D85A30")) +
    labs(title = "BATCH EFFECT CORRECTION: Before vs After Normalization",
         subtitle = "CD56dim CD107a - Batches should overlap after normalization",
         y = "Value (see row label)", x = "") +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "01_Batch_Correction_Comparison.png"), p_compare, 
         width = 10, height = 12, dpi = 150)
  
  # Statistical test
  batch_stats <- plot_data %>%
    group_by(Method, Antigen) %>%
    summarise(
      Batch1_mean = mean(Value[Batch == "Batch 1"], na.rm = TRUE),
      Batch2_mean = mean(Value[Batch == "Batch 2"], na.rm = TRUE),
      Batch1_median = median(Value[Batch == "Batch 1"], na.rm = TRUE),
      Batch2_median = median(Value[Batch == "Batch 2"], na.rm = TRUE),
      Fold_Diff_Mean = max(Batch1_mean, Batch2_mean) / min(Batch1_mean, Batch2_mean),
      Wilcox_p = tryCatch(
        wilcox.test(Value[Batch == "Batch 1"], Value[Batch == "Batch 2"])$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Batch_Effect_Significant = ifelse(!is.na(Wilcox_p) & Wilcox_p < 0.05, "YES", "NO"),
      Fold_Diff_Mean = round(Fold_Diff_Mean, 2),
      Wilcox_p = round(Wilcox_p, 4)
    )
  
  write.csv(batch_stats, file.path(output_dir, "02_Batch_Effect_Statistics.csv"), row.names = FALSE)
  
  # Density plot
  p_density <- ggplot(plot_data %>% filter(Method == "3_N6_Normalized (% of max)"), 
                      aes(x = Value, fill = Batch)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ Antigen, scales = "free") +
    scale_fill_manual(values = c("Batch 1" = "#3B8BD4", "Batch 2" = "#D85A30")) +
    labs(title = "Distribution After N6 Normalization",
         subtitle = "Batches should largely overlap - confirms batch correction",
         x = "% of N6 maximum response", y = "Density") +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "03_Normalized_Density_Overlay.png"), p_density, 
         width = 10, height = 5, dpi = 150)
  
  message("Saved batch correction plots to: ", output_dir)
  
  return(batch_stats)
}

# =============================================================================
# FUNCTION: Export cleaned metadata
# =============================================================================

export_cleaned_metadata <- function(data, output_dir) {
  
  message("\n", strrep("=", 60))
  message("EXPORTING CLEANED METADATA")
  message(strrep("=", 60))
  
  samples <- data %>% filter(!IsControl)
  
  # -------------------------------------------------------------------------
  # CHER: Age (weeks), Phase, Timepoint_Numeric (cleaned from Timepoint_Original)
  # -------------------------------------------------------------------------
  cher_meta <- samples %>%
    filter(Cohort == "CHER") %>%
    select(PID, Timepoint_Original, Timepoint_Numeric, `Age (weeks)`, Phase, Gender, Cohort) %>%
    distinct() %>%
    arrange(PID, `Age (weeks)`)
  
  write.csv(cher_meta, file.path(output_dir, "CHER_Metadata.csv"), row.names = FALSE)
  message("Saved: CHER_Metadata.csv (", nrow(cher_meta), " unique PID-timepoint combinations)")
  
  # -------------------------------------------------------------------------
  # TARA: Age (months) extracted from Timepoint column, Gender needed
  # -------------------------------------------------------------------------
  tara_meta <- samples %>%
    filter(Cohort == "TARA") %>%
    select(PID, `Age (months)`, Gender, Cohort) %>%
    distinct() %>%
    mutate(
      Age_Missing = is.na(`Age (months)`),
      Gender_Missing = is.na(Gender)
    ) %>%
    arrange(PID)
  
  write.csv(tara_meta, file.path(output_dir, "TARA_Metadata.csv"), row.names = FALSE)
  message("Saved: TARA_Metadata.csv (", nrow(tara_meta), " unique PID-age combinations)")
  
  # -------------------------------------------------------------------------
  # CHER timepoint cleaning log - show what was converted
  # -------------------------------------------------------------------------
  cher_tp_map <- samples %>%
    filter(Cohort == "CHER") %>%
    select(Timepoint_Original, Timepoint_Numeric) %>%
    distinct() %>%
    arrange(Timepoint_Numeric)
  
  write.csv(cher_tp_map, file.path(output_dir, "CHER_Timepoint_Mapping.csv"), row.names = FALSE)
  message("Saved: CHER_Timepoint_Mapping.csv (", nrow(cher_tp_map), " unique timepoint formats)")
  
  # -------------------------------------------------------------------------
  # Missing data summary
  # -------------------------------------------------------------------------
  missing_df <- data.frame(
    Cohort = c("CHER", "TARA"),
    N_Subjects = c(n_distinct(cher_meta$PID), n_distinct(tara_meta$PID)),
    N_Sample_Rows = c(nrow(samples %>% filter(Cohort == "CHER")), 
                      nrow(samples %>% filter(Cohort == "TARA"))),
    Age_Column = c("Age (weeks)", "Age (months)"),
    Missing_Age = c(sum(is.na(cher_meta$`Age (weeks)`)), sum(tara_meta$Age_Missing)),
    Missing_Gender = c(sum(is.na(cher_meta$Gender)), sum(tara_meta$Gender_Missing)),
    Missing_Phase = c(sum(is.na(cher_meta$Phase)), NA)
  )
  
  write.csv(missing_df, file.path(output_dir, "Missing_Data_Summary.csv"), row.names = FALSE)
  message("Saved: Missing_Data_Summary.csv")
  
  return(list(cher = cher_meta, tara = tara_meta, tp_map = cher_tp_map, summary = missing_df))
}

# =============================================================================
# FUNCTION: Summary report (HTML) - Self-contained with embedded images
# =============================================================================

# Helper function to encode image as base64
encode_image_base64 <- function(image_path) {
  if (file.exists(image_path)) {
    img_data <- base64enc::base64encode(image_path)
    return(paste0("data:image/png;base64,", img_data))
  } else {
    return("")
  }
}

generate_summary_report <- function(data, hierarchy_check, n6_check, batch_check, 
                                    sample_qc, batch_correction_stats, control_stats,
                                    output_dir, plots_dir, batch_corr_dir) {
  
  # Check for base64enc package
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    install.packages("base64enc")
  }
  library(base64enc)
  
  samples <- data %>% filter(!IsControl)
  
  # Calculate summary statistics
  n_cher_subjects <- n_distinct(samples$PID[samples$Cohort == "CHER"])
  n_cher_samples <- sum(samples$Cohort == "CHER")
  n_tara_subjects <- n_distinct(samples$PID[samples$Cohort == "TARA"])
  n_tara_samples <- sum(samples$Cohort == "TARA")
  
  # Encode all images as base64
  message("Embedding images in HTML report...")
  img_control_perf <- encode_image_base64(file.path(plots_dir, "01_Control_Performance.png"))
  img_batch_raw <- encode_image_base64(file.path(plots_dir, "02_Batch_Distribution_RAW.png"))
  img_event_counts <- encode_image_base64(file.path(plots_dir, "03_Event_Counts.png"))
  img_lymph_freq <- encode_image_base64(file.path(plots_dir, "04_Lymphocyte_Freq.png"))
  img_batch_corr <- encode_image_base64(file.path(batch_corr_dir, "01_Batch_Correction_Comparison.png"))
  img_density <- encode_image_base64(file.path(batch_corr_dir, "03_Normalized_Density_Overlay.png"))
  
  # Build HTML report
  html <- c('
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ADNK Assay QC Report</title>
  <style>
    * { box-sizing: border-box; }
    body { 
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      line-height: 1.6; 
      max-width: 1200px; 
      margin: 0 auto; 
      padding: 20px;
      background: #f8f9fa;
      color: #333;
    }
    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
    h2 { color: #34495e; border-bottom: 2px solid #bdc3c7; padding-bottom: 8px; margin-top: 40px; }
    h3 { color: #7f8c8d; margin-top: 25px; }
    .card {
      background: white;
      border-radius: 8px;
      padding: 20px;
      margin: 20px 0;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .summary-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 15px;
      margin: 20px 0;
    }
    .stat-box {
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      padding: 20px;
      border-radius: 8px;
      text-align: center;
    }
    .stat-box.green { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }
    .stat-box.orange { background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); }
    .stat-box.blue { background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); }
    .stat-number { font-size: 2.5em; font-weight: bold; }
    .stat-label { font-size: 0.9em; opacity: 0.9; }
    table { 
      width: 100%; 
      border-collapse: collapse; 
      margin: 15px 0;
      background: white;
    }
    th, td { 
      padding: 12px 15px; 
      text-align: left; 
      border-bottom: 1px solid #ddd; 
    }
    th { 
      background: #34495e; 
      color: white; 
      font-weight: 500;
    }
    tr:hover { background: #f5f5f5; }
    .pass { color: #27ae60; font-weight: bold; }
    .fail { color: #e74c3c; font-weight: bold; }
    .warn { color: #f39c12; font-weight: bold; }
    .info-box {
      background: #e8f4f8;
      border-left: 4px solid #3498db;
      padding: 15px;
      margin: 15px 0;
      border-radius: 0 8px 8px 0;
    }
    .warning-box {
      background: #fef9e7;
      border-left: 4px solid #f1c40f;
      padding: 15px;
      margin: 15px 0;
      border-radius: 0 8px 8px 0;
    }
    .success-box {
      background: #eafaf1;
      border-left: 4px solid #27ae60;
      padding: 15px;
      margin: 15px 0;
      border-radius: 0 8px 8px 0;
    }
    .formula {
      background: #f4f4f4;
      padding: 10px 15px;
      border-radius: 5px;
      font-family: "Courier New", monospace;
      margin: 10px 0;
      overflow-x: auto;
    }
    img { 
      max-width: 100%; 
      height: auto; 
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      margin: 10px 0;
    }
    .img-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
      gap: 20px;
    }
    .toc {
      background: white;
      padding: 20px;
      border-radius: 8px;
      margin: 20px 0;
    }
    .toc ul { list-style: none; padding-left: 0; }
    .toc li { padding: 5px 0; }
    .toc a { color: #3498db; text-decoration: none; }
    .toc a:hover { text-decoration: underline; }
    footer {
      margin-top: 40px;
      padding-top: 20px;
      border-top: 1px solid #ddd;
      color: #7f8c8d;
      font-size: 0.9em;
    }
  </style>
</head>
<body>

<h1>ADNK Assay QC Report</h1>
<p><strong>Generated:</strong> ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
<p><em>This is a self-contained report - all images are embedded and no external files are needed.</em></p>

<div class="toc card">
  <h3>Contents</h3>
  <ul>
    <li><a href="#overview">1. Data Overview</a></li>
    <li><a href="#controls">2. Control Performance</a></li>
    <li><a href="#qc-checks">3. QC Checks</a></li>
    <li><a href="#normalization">4. Normalization Strategy</a></li>
    <li><a href="#batch-correction">5. Batch Correction Results</a></li>
    <li><a href="#outputs">6. Output Files</a></li>
  </ul>
</div>

<!-- ===== SECTION 1: DATA OVERVIEW ===== -->
<h2 id="overview">1. Data Overview</h2>

<div class="card">
  <h3>Study Design</h3>
  <p>This dataset contains ADNK (Antibody-Dependent NK cell) functional assay data from two pediatric HIV cohorts:</p>
  
  <table>
    <tr>
      <th>Cohort</th>
      <th>Design</th>
      <th>Subjects</th>
      <th>Samples</th>
      <th>Key Question</th>
    </tr>
    <tr>
      <td><strong>CHER</strong></td>
      <td>Longitudinal, treatment interruption</td>
      <td>', n_cher_subjects, '</td>
      <td>', n_cher_samples, '</td>
      <td>How does ADCC change across treatment phases?</td>
    </tr>
    <tr>
      <td><strong>TARA</strong></td>
      <td>Cross-sectional</td>
      <td>', n_tara_subjects, '</td>
      <td>', n_tara_samples, '</td>
      <td>Does viral suppression status affect ADCC?</td>
    </tr>
  </table>
  
  <div class="info-box">
    <p><strong>Note on data structure:</strong></p>
    <ul>
      <li><strong>Batch 1:</strong> CHER samples only (42 per antigen)</li>
      <li><strong>Batch 2:</strong> CHER (20) + TARA (23) samples per antigen</li>
      <li><strong>CHER age:</strong> Recorded in weeks in the source data</li>
      <li><strong>TARA age:</strong> Derived from "Timepoint relative to phase" column (values are months: 36-71)</li>
    </ul>
  </div>
</div>

<div class="summary-grid">
  <div class="stat-box">
    <div class="stat-number">', nrow(samples), '</div>
    <div class="stat-label">Total Samples</div>
  </div>
  <div class="stat-box green">
    <div class="stat-number">', sum(data$IsControl), '</div>
    <div class="stat-label">Controls</div>
  </div>
  <div class="stat-box orange">
    <div class="stat-number">4</div>
    <div class="stat-label">Plates (Sheets)</div>
  </div>
  <div class="stat-box blue">
    <div class="stat-number">2</div>
    <div class="stat-label">Antigens (p24, gp120)</div>
  </div>
</div>

<div class="card">
  <h3>Assay Structure</h3>
  <p>The ADNK assay measures NK cell functional responses to HIV antigens (p24 and gp120) in the presence of participant plasma antibodies.</p>
  
  <h4>Gating Hierarchy (FlowJo)</h4>
  <p>Total Events → Lymphocytes → Single Cells → Live CD3- → NK subsets:</p>
  <ul>
    <li><strong>CD56bright</strong> (~1.8% of CD3-) - Cytokine producers</li>
    <li><strong>CD56dim</strong> (~84% of CD3-) - Main ADCC effectors</li>
    <li><strong>non-NK</strong> (~11%) - CD3-CD56- cells</li>
  </ul>
  
  <h4>Functional Markers Measured</h4>
  <table>
    <tr>
      <th>Marker</th>
      <th>Function</th>
      <th>Interpretation</th>
    </tr>
    <tr><td>CD107a</td><td>Degranulation</td><td>Cytotoxic granule release (primary ADCC readout)</td></tr>
    <tr><td>IFNγ</td><td>Cytokine</td><td>Antiviral inflammatory response</td></tr>
    <tr><td>MIP1β</td><td>Chemokine</td><td>Immune cell recruitment</td></tr>
    <tr><td>FasL</td><td>Death ligand</td><td>Alternative killing mechanism</td></tr>
  </table>
</div>

<!-- ===== SECTION 2: CONTROLS ===== -->
<h2 id="controls">2. Control Performance</h2>

<div class="card">
  <h3>Control Definitions</h3>
  <div class="info-box">
    <p><strong>Why controls matter:</strong> Controls establish the dynamic range of the assay and allow us to normalize across batches. Each plate has its own set of controls.</p>
    <p><strong>Note:</strong> Batch 1 had triplicate PBS controls (averaged for analysis); Batch 2 had single PBS controls.</p>
  </div>
  
  <table>
    <tr>
      <th>Control</th>
      <th>Purpose</th>
      <th>Expected Response</th>
    </tr>
    <tr>
      <td><strong>PBS</strong></td>
      <td>No antibody baseline</td>
      <td>LOWEST - Background/spontaneous NK activation only</td>
    </tr>
    <tr>
      <td><strong>HIV neg plasma</strong></td>
      <td>Negative control</td>
      <td>LOW - No HIV-specific antibodies present</td>
    </tr>
    <tr>
      <td><strong>N6 bnAb</strong></td>
      <td>Positive control</td>
      <td>HIGHEST - Known broadly neutralizing antibody, strong ADCC</td>
    </tr>
  </table>
  
  <p><strong>Expected hierarchy:</strong> PBS &lt; HIV neg plasma ≤ N6 bnAb</p>
</div>

<div class="card">
  <h3>Control Hierarchy Check (CD56dim CD107a)</h3>
  <table>
    <tr>
      <th>Plate</th>
      <th>PBS</th>
      <th>HIV neg</th>
      <th>N6 bnAb</th>
      <th>Hierarchy OK?</th>
    </tr>')
  
  # Add hierarchy check rows
  for (i in 1:nrow(hierarchy_check)) {
    row <- hierarchy_check[i, ]
    status_class <- ifelse(row$Hierarchy_OK, "pass", "fail")
    status_text <- ifelse(row$Hierarchy_OK, "✓ PASS", "✗ FAIL")
    html <- c(html, sprintf('
    <tr>
      <td>%s %s</td>
      <td>%.2f%%</td>
      <td>%.2f%%</td>
      <td>%.2f%%</td>
      <td class="%s">%s</td>
    </tr>', row$Batch, row$Antigen, row$PBS*100, row$HIV_neg*100, row$N6_bnAb*100, status_class, status_text))
  }
  
  html <- c(html, '
  </table>')
  
  # Add embedded image
  if (nchar(img_control_perf) > 0) {
    html <- c(html, paste0('<img src="', img_control_perf, '" alt="Control Performance Plot">'))
  }
  
  html <- c(html, '
</div>

<!-- ===== SECTION 3: QC CHECKS ===== -->
<h2 id="qc-checks">3. QC Checks</h2>

<div class="card">
  <h3>QC Check 1: Control Hierarchy</h3>
  <div class="info-box">
    <p><strong>What we check:</strong> PBS should have the lowest response (background), and N6 bnAb should have the highest (maximum ADCC).</p>
    <p><strong>Why it matters:</strong> If this hierarchy is violated, the assay may have failed (e.g., contamination, pipetting error, reagent issues).</p>
  </div>
  <p><strong>Result:</strong> <span class="pass">', sum(hierarchy_check$Hierarchy_OK), '/', nrow(hierarchy_check), ' plates pass</span></p>
</div>

<div class="card">
  <h3>QC Check 2: N6 Positive Control Performance</h3>
  <div class="info-box">
    <p><strong>What we check:</strong> N6 bnAb should trigger strong degranulation (CD107a &gt; 5%).</p>
    <p><strong>Why it matters:</strong> A weak positive control suggests the assay dynamic range is compressed, making it harder to detect real differences.</p>
  </div>
  
  <table>
    <tr>
      <th>Plate</th>
      <th>N6 CD107a (%)</th>
      <th>Threshold</th>
      <th>Status</th>
    </tr>')
  
  for (i in 1:nrow(n6_check)) {
    row <- n6_check[i, ]
    status_class <- ifelse(row$N6_OK, "pass", "warn")
    status_text <- ifelse(row$N6_OK, "✓ PASS", "⚠ LOW")
    html <- c(html, sprintf('
    <tr>
      <td>%s %s</td>
      <td>%.2f%%</td>
      <td>&gt; 5%%</td>
      <td class="%s">%s</td>
    </tr>', row$Batch, row$Antigen, row$N6_CD107a*100, status_class, status_text))
  }
  
  html <- c(html, '
  </table>')
  
  # Add warning if Batch 2 is low
  if (any(!n6_check$N6_OK)) {
    html <- c(html, '
  <div class="warning-box">
    <strong>⚠ Note:</strong> Batch 2 shows lower N6 response than Batch 1. This could indicate:
    <ul>
      <li>N6 stock degradation or concentration issue</li>
      <li>Target cell coating variability</li>
      <li>NK cell viability differences between batches</li>
    </ul>
    <p>The normalization strategy (below) accounts for this by expressing results as % of each plate\'s N6 maximum.</p>
  </div>')
  }
  
  html <- c(html, '
</div>

<div class="card">
  <h3>QC Check 3: Batch Effects</h3>
  <div class="info-box">
    <p><strong>What we check:</strong> Background (PBS) levels should be similar across batches. Large differences indicate technical variability.</p>
    <p><strong>Threshold:</strong> &lt; 3-fold difference is acceptable.</p>
  </div>
  
  <table>
    <tr>
      <th>Marker</th>
      <th>Antigen</th>
      <th>Batch 1 PBS</th>
      <th>Batch 2 PBS</th>
      <th>Fold Difference</th>
      <th>Status</th>
    </tr>')
  
  for (i in 1:nrow(batch_check)) {
    row <- batch_check[i, ]
    status_class <- ifelse(row$Batch_OK, "pass", "warn")
    status_text <- ifelse(row$Batch_OK, "✓ OK", "⚠ BATCH EFFECT")
    html <- c(html, sprintf('
    <tr>
      <td>%s</td>
      <td>%s</td>
      <td>%.2f%%</td>
      <td>%.2f%%</td>
      <td>%.1fx</td>
      <td class="%s">%s</td>
    </tr>', row$Marker, row$Antigen, row$`Batch 1`*100, row$`Batch 2`*100, row$Fold_Diff, status_class, status_text))
  }
  
  html <- c(html, '
  </table>')
  
  # Add embedded image
  if (nchar(img_batch_raw) > 0) {
    html <- c(html, paste0('<img src="', img_batch_raw, '" alt="Batch Distribution RAW">'))
  }
  
  html <- c(html, '
</div>

<div class="card">
  <h3>QC Check 4: Sample-Level Quality</h3>
  <div class="info-box">
    <p><strong>What we check:</strong></p>
    <ul>
      <li>Minimum event count: ', format(QC_THRESHOLDS$MIN_EVENT_COUNT, big.mark=","), ' events (statistical power)</li>
      <li>Minimum lymphocyte frequency: ', QC_THRESHOLDS$MIN_LYMPHOCYTE_FREQ*100, '% (sample quality)</li>
    </ul>
  </div>
  
  <div class="success-box">
    <strong>Result:</strong> ', sum(sample_qc$QC_Pass), ' / ', nrow(sample_qc), ' samples pass QC (', round(mean(sample_qc$QC_Pass)*100, 1), '%)
  </div>
  
  <div class="img-grid">')
  
  # Add embedded images
  if (nchar(img_event_counts) > 0) {
    html <- c(html, paste0('<img src="', img_event_counts, '" alt="Event Counts">'))
  }
  if (nchar(img_lymph_freq) > 0) {
    html <- c(html, paste0('<img src="', img_lymph_freq, '" alt="Lymphocyte Frequency">'))
  }
  
  html <- c(html, '
  </div>
</div>

<!-- ===== SECTION 4: NORMALIZATION ===== -->
<h2 id="normalization">4. Normalization Strategy</h2>

<div class="card">
  <h3>Why Normalize?</h3>
  <div class="info-box">
    <p>Raw percentages cannot be directly compared across plates/batches because:</p>
    <ul>
      <li>Background activation varies (PBS levels differ)</li>
      <li>Maximum response varies (N6 levels differ)</li>
      <li>Each plate may have different NK cell responsiveness</li>
    </ul>
  </div>
  
  <h3>Step 1: Background Subtraction</h3>
  <p>Remove spontaneous NK activation (PBS background) from each sample:</p>
  <div class="formula">
    Sample_BGsub = Sample_raw - PBS_plate
  </div>
  <p>This isolates the antibody-dependent response.</p>
  
  <h3>Step 2: N6 Normalization</h3>
  <p>Express each sample as a percentage of the maximum possible response (N6):</p>
  <div class="formula">
    Sample_normalized = (Sample_raw - PBS_plate) / (N6_plate - PBS_plate) × 100
  </div>
  <p>This gives "% of maximum ADCC response" which is comparable across all plates.</p>
  
  <div class="success-box">
    <strong>Key point:</strong> Each plate uses its OWN controls for normalization. Batch 1 samples are normalized to Batch 1 controls, and Batch 2 samples to Batch 2 controls.
  </div>
</div>

<!-- ===== SECTION 5: BATCH CORRECTION ===== -->
<h2 id="batch-correction">5. Batch Correction Results</h2>

<div class="card">
  <h3>Did Normalization Work?</h3>
  <p>We compare batch distributions before and after normalization. Successful normalization should:</p>
  <ul>
    <li>Reduce the fold-difference between batches</li>
    <li>Make batch distributions overlap</li>
    <li>Eliminate statistically significant batch effects</li>
  </ul>')
  
  # Add embedded image
  if (nchar(img_batch_corr) > 0) {
    html <- c(html, paste0('<img src="', img_batch_corr, '" alt="Batch Correction Comparison">'))
  }
  
  # Add batch correction statistics if available
  if (!is.null(batch_correction_stats) && nrow(batch_correction_stats) > 0) {
    html <- c(html, '
  <h3>Statistical Confirmation</h3>
  <table>
    <tr>
      <th>Stage</th>
      <th>Antigen</th>
      <th>Batch 1 Mean</th>
      <th>Batch 2 Mean</th>
      <th>Fold Diff</th>
      <th>p-value</th>
      <th>Significant?</th>
    </tr>')
    
    for (i in 1:nrow(batch_correction_stats)) {
      row <- batch_correction_stats[i, ]
      sig_class <- ifelse(row$Batch_Effect_Significant == "YES", "warn", "pass")
      html <- c(html, sprintf('
    <tr>
      <td>%s</td>
      <td>%s</td>
      <td>%.2f</td>
      <td>%.2f</td>
      <td>%.2fx</td>
      <td>%.4f</td>
      <td class="%s">%s</td>
    </tr>', row$Method, row$Antigen, row$Batch1_mean, row$Batch2_mean, 
                              row$Fold_Diff_Mean, row$Wilcox_p, sig_class, row$Batch_Effect_Significant))
    }
    
    html <- c(html, '
  </table>')
    
    # Add embedded density image
    if (nchar(img_density) > 0) {
      html <- c(html, paste0('<img src="', img_density, '" alt="Normalized Density Overlay">'))
    }
  }
  
  html <- c(html, '
</div>

<!-- ===== SECTION 6: OUTPUTS ===== -->
<h2 id="outputs">6. Output Files</h2>

<div class="card">
  <table>
    <tr>
      <th>Folder</th>
      <th>Contents</th>
      <th>Description</th>
    </tr>
    <tr>
      <td><strong>01_QC_Reports/</strong></td>
      <td>CSV files</td>
      <td>Detailed QC check results for each plate</td>
    </tr>
    <tr>
      <td><strong>02_QC_Plots/</strong></td>
      <td>PNG images</td>
      <td>Control performance, batch distributions, event counts</td>
    </tr>
    <tr>
      <td><strong>03_Batch_Correction/</strong></td>
      <td>PNG images + CSV</td>
      <td>Before/after normalization comparisons</td>
    </tr>
    <tr>
      <td><strong>04_Normalized_Data/</strong></td>
      <td>XLSX files</td>
      <td>Analysis-ready datasets with normalized values</td>
    </tr>
    <tr>
      <td><strong>05_Cleaned_Metadata/</strong></td>
      <td>CSV files</td>
      <td>Subject metadata with issues flagged</td>
    </tr>
    <tr>
      <td><strong>06_Summary/</strong></td>
      <td>This report</td>
      <td>QC summary and methodology documentation</td>
    </tr>
  </table>
  
  <h3>Key Output Files for Analysis</h3>
  <ul>
    <li><strong>ADNK_Analysis_Ready.xlsx</strong> - All QC-passed samples with raw + normalized columns</li>
    <li><strong>CHER_Analysis_Ready.xlsx</strong> - CHER cohort only</li>
    <li><strong>TARA_Analysis_Ready.xlsx</strong> - TARA cohort only</li>
  </ul>
  
  <h3>Column Naming Convention</h3>
  <table>
    <tr>
      <th>Suffix</th>
      <th>Meaning</th>
      <th>Example</th>
    </tr>
    <tr>
      <td>(none)</td>
      <td>Raw value from FlowJo</td>
      <td>NK CD56dim/CD107a | Freq. of Parent</td>
    </tr>
    <tr>
      <td>_PBS</td>
      <td>PBS control value for this plate</td>
      <td>...CD107a | Freq. of Parent_PBS</td>
    </tr>
    <tr>
      <td>_BGsub</td>
      <td>Background-subtracted</td>
      <td>...CD107a | Freq. of Parent_BGsub</td>
    </tr>
    <tr>
      <td>_N6</td>
      <td>N6 control value for this plate</td>
      <td>...CD107a | Freq. of Parent_N6</td>
    </tr>
    <tr>
      <td>_NormN6</td>
      <td>Normalized to N6 (% of max)</td>
      <td>...CD107a | Freq. of Parent_NormN6</td>
    </tr>
  </table>
</div>

<footer>
  <p>Report generated by ADNK QC Pipeline | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
  <p><em>This is a self-contained HTML file with all images embedded - no external files required.</em></p>
</footer>

</body>
</html>')
  
  # Write HTML file
  html_file <- file.path(output_dir, "ADNK_QC_Report.html")
  writeLines(html, html_file)
  message("Saved self-contained HTML report to: ", html_file)
  
  # Also save a simple CSV summary
  qc_table <- data.frame(
    Check = c("Control Hierarchy", "N6 Positive Control", "Batch Effects (raw)", "Sample QC"),
    Passed = c(sum(hierarchy_check$Hierarchy_OK), sum(n6_check$N6_OK), 
               sum(batch_check$Batch_OK), sum(sample_qc$QC_Pass)),
    Total = c(nrow(hierarchy_check), nrow(n6_check), nrow(batch_check), nrow(sample_qc)),
    Pass_Rate_Pct = round(c(mean(hierarchy_check$Hierarchy_OK), mean(n6_check$N6_OK),
                            mean(batch_check$Batch_OK), mean(sample_qc$QC_Pass)) * 100, 1)
  )
  write.csv(qc_table, file.path(output_dir, "QC_Summary_Table.csv"), row.names = FALSE)
  
  return(html_file)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

run_adnk_qc <- function(input_file = INPUT_FILE, dirs = DIRS) {
  
  message("\n", strrep("#", 70))
  message("  ADNK ASSAY QC AND NORMALIZATION PIPELINE")
  message(strrep("#", 70))
  message("\nInput:  ", input_file)
  message("Output: ", dirs$root)
  
  # 1. Load
  message("\n[1/9] Loading data...")
  data <- load_adnk_data(input_file)
  
  # 2. Control stats
  message("\n[2/9] Calculating control statistics...")
  control_stats <- calculate_control_stats(data)
  write.csv(control_stats, file.path(dirs$qc_reports, "Control_Statistics_Full.csv"), row.names = FALSE)
  
  # 3. QC checks
  message("\n[3/9] Running QC checks...")
  hierarchy_check <- qc_control_hierarchy(data, control_stats, dirs$qc_reports)
  n6_check <- qc_positive_control(data, control_stats, dirs$qc_reports)
  batch_check <- qc_batch_effects(data, control_stats, dirs$qc_reports)
  sample_qc <- qc_sample_level(data, dirs$qc_reports)
  
  # 4. QC plots
  message("\n[4/9] Generating QC plots...")
  generate_qc_plots(data, control_stats, dirs$qc_plots)
  
  # 5. Background subtraction
  message("\n[5/9] Background subtraction...")
  samples_bgsub <- subtract_background(data, control_stats)
  
  # 6. N6 normalization
  message("\n[6/9] N6 normalization...")
  samples_normalized <- normalize_to_n6(samples_bgsub, control_stats)
  
  # Add QC flags - join by Sheet as well to avoid many-to-many
  samples_normalized <- samples_normalized %>%
    left_join(
      sample_qc %>% select(PID, Timepoint_Original, Sheet, Batch, Antigen, 
                           QC_LowEvents, QC_LowLymphocytes, QC_Pass),
      by = c("PID", "Timepoint_Original", "Sheet", "Batch", "Antigen")
    )
  
  # 7. Batch correction plots
  message("\n[7/9] Generating batch correction confirmation...")
  batch_correction_stats <- generate_batch_correction_plots(data, samples_normalized, dirs$batch_correction)
  
  # 8. Metadata
  message("\n[8/9] Exporting cleaned metadata...")
  metadata <- export_cleaned_metadata(data, dirs$metadata)
  
  # 9. Save data
  message("\n[9/9] Saving normalized datasets...")
  
  write_xlsx(samples_normalized, file.path(dirs$normalized_data, "ADNK_All_Normalized.xlsx"))
  
  analysis_ready <- samples_normalized %>% filter(QC_Pass)
  write_xlsx(analysis_ready, file.path(dirs$normalized_data, "ADNK_Analysis_Ready.xlsx"))
  
  write_xlsx(analysis_ready %>% filter(Cohort == "CHER"), 
             file.path(dirs$normalized_data, "CHER_Analysis_Ready.xlsx"))
  write_xlsx(analysis_ready %>% filter(Cohort == "TARA"), 
             file.path(dirs$normalized_data, "TARA_Analysis_Ready.xlsx"))
  
  message(sprintf("Saved: ADNK_Analysis_Ready.xlsx (%d samples)", nrow(analysis_ready)))
  
  # Summary HTML report
  report <- generate_summary_report(data, hierarchy_check, n6_check, batch_check,
                                    sample_qc, batch_correction_stats, control_stats,
                                    dirs$summary, dirs$qc_plots, dirs$batch_correction)
  
  message("\n", strrep("#", 70))
  message("  PIPELINE COMPLETE")
  message(strrep("#", 70))
  message("\nOutput: ", dirs$root)
  
  return(list(
    raw_data = data,
    control_stats = control_stats,
    normalized_data = samples_normalized,
    analysis_ready = analysis_ready,
    qc = list(hierarchy = hierarchy_check, n6 = n6_check, batch = batch_check, samples = sample_qc),
    batch_correction_stats = batch_correction_stats,
    metadata = metadata
  ))
}

# =============================================================================
# RUN
# =============================================================================

message("\n=== ADNK QC Pipeline ===")
message("Working directory: ", getwd())

if (file.exists(INPUT_FILE)) {
  results <- run_adnk_qc()
} else {
  message("\nERROR: Input file not found!")
  message("Expected: ", INPUT_FILE)
  message("\nPlease check the path in the CONFIGURATION section.")
}

