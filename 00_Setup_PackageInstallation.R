# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is searched for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  (leave "" if unsure)

local({
  cand = c(if (nzchar(CODE_FOLDER)) file.path(path.expand(CODE_FOLDER), "config_paths.R"),
           "config_paths.R", "../config_paths.R", "Code/config_paths.R")
  hit  = cand[file.exists(cand)]
  if (!length(hit))
    stop("config_paths.R not found.\n",
         "Fix either way:\n",
         "  1. setwd(\"/path/to/AlphaBurstRhythm/Code\")   then re-run, or\n",
         "  2. set CODE_FOLDER at the top of this script to that same path.\n",
         "Currently looking from: ", getwd(), call. = FALSE)
  source(hit[1], local = FALSE)
})

# =============================================================================
# -----------------------------------------------------------------------------
# Script: 00_Setup_PackageInstallation.R
# Purpose: Installs and verifies all R packages required across the project.
#          Defines ALL shared global constants (thresholds, colours, themes)
#          that every other script sources by calling:
#              source("00_Setup_PackageInstallation.R")
#          You don't need to run the script per se. It is called every time in the other code. 
# =============================================================================
# Outputs: None
# =============================================================================


# =============================================================================
# SECTION 1: CRAN PACKAGES
# =============================================================================
# Complete list of all packages used across the project, organised by purpose.

required_packages = c(

  # --- Core Tidyverse & Data Manipulation ---
  "tidyverse",    # Meta-package: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats
  "dplyr",        # Data manipulation (explicit for namespace clarity)
  "tidyr",        # Data reshaping
  "readr",        # Read CSV
  "readxl",       # Reading Excel 
  "haven",        # Reading SPSS / Stata / SAS files
  "reshape2",     # Legacy data reshaping
  "stringr",      # String manipulation
  "purrr",        # 
  "tibble",       # 

  # --- Plotting ---
  "ggplot2",      # Core plotting system
  "ggpubr",       # Publication-ready ggplot2 themes and helpers
  "patchwork",    # Combining multiple ggplot2 plots into one figure
  "ggnewscale",   # Color scales
  "colorspace",   # Colour manipulation in the plots
  "scales",       # Scale functions for ggplot2 

  # --- Mixed Models & Statistics ---
  "lme4",         # MLM models
  "lmerTest",     # p-values for lme4 via Satterthwaite approximation in lme4
  "nlme",         # Linear/nonlinear mixed-effects models used in GAMMs
  "mgcv",         # GAMMs
  "gratia",       # Tidy tools for GAM/GAMM
  "AICcmodavg",   # AICc-based model comparison and selection
  "MuMIn",        # Multi-model inference; R² for mixed models
  "performance",  # Model diagnostics: R², VIF, normality, heteroscedasticity
  "broom.mixed",  # Tidy methods for mixed-effects models

  # --- Partial Correlations, Bootstrapping & SEM ---
  "ppcor",        # Partial and semi-partial correlations
  "boot",         # Bootstrap resampling

  # --- Classification & ROC Analysis ---
  "caret",        # Cross-validation and data partitioning for classification
  "pROC",         # ROC curve analysis and AUC computation

  # --- Descriptive Statistics & Tables ---
  "psych",        # Descriptive statistics, correlation matrices
  "gtsummary",    # Summary Tables
  "flextable",    # Formatted tables for Word / HTML output
  "officer",      # Creating and modifying Word (.docx) documents

  # --- Additional Statistics ---
  "VGAM",         # Vector GAM models (multinomial logistic regression)

  # --- Utilities ---
  "glue",         # String interpolation
  "rempsyc",      # Utility functions for psychology research
  "rtf"           # RTF file output
)


# =============================================================================
# SECTION 2: INSTALLATION
# =============================================================================

# Identify which packages are already installed while keeping those still needing installation 
# for the next step.
already_installed   = installed.packages()[, "Package"]
packages_to_install = setdiff(required_packages, already_installed)

if (length(packages_to_install) == 0) {
  cat("All required packages are already installed.\n\n")
} else {
  cat(sprintf("Installing %d missing package(s):\n", length(packages_to_install)))
  cat(paste(" -", packages_to_install, collapse = "\n"), "\n\n")
  install.packages(packages_to_install, dependencies = TRUE)
}

# =============================================================================
# SECTION 3: LOADING & VERIFICATION
# =============================================================================

# Initialize data frame to track installation/load status for each package
load_results = data.frame(
  Package   = required_packages,
  Installed = NA,
  Version   = NA,
  Loaded    = NA,
  stringsAsFactors = FALSE
)

# For each required package, check if installed and attempt to load
for (i in seq_along(required_packages)) {
  pkg = required_packages[i]

  if (pkg %in% installed.packages()[, "Package"]) {
    load_results$Installed[i] = TRUE
    load_results$Version[i]   = as.character(packageVersion(pkg))

    # Attempt to load the package; suppress non-critical messages/warnings
    load_ok = suppressMessages(suppressWarnings(
      require(pkg, character.only = TRUE, quietly = TRUE)
    ))
    load_results$Loaded[i] = load_ok

  } else {
    load_results$Installed[i] = FALSE
    load_results$Version[i]   = NA
    load_results$Loaded[i]    = FALSE
  }
}

# --- Print verification table with status for each package ---
for (i in seq_len(nrow(load_results))) {
  r             = load_results[i, ]
  installed_str = ifelse(r$Installed, "YES", "NO  *** MISSING ***")
  loaded_str    = ifelse(r$Loaded,    "YES", "FAILED")
  version_str   = ifelse(is.na(r$Version), "---", r$Version)
  cat(sprintf("%-30s %-12s %-12s %-8s\n",
              r$Package, installed_str, version_str, loaded_str))
}

# Count and report on failed package loads
n_failed = sum(!load_results$Loaded)
cat(strrep("=", 65), "\n")
if (n_failed == 0) {
  cat("SUCCESS: All", nrow(load_results), "packages loaded successfully.\n")
} else {
  cat(sprintf("WARNING: %d package(s) failed to load:\n", n_failed))
  cat(paste(" -", load_results$Package[!load_results$Loaded], collapse = "\n"), "\n")
  cat("\nTry re-running install.packages() for the listed packages above.\n")
}
cat(strrep("=", 65), "\n\n")

# =============================================================================
# SECTION 4: SESSION INFO
# =============================================================================

cat("Session Information (for reproducibility):\n")
cat(strrep("-", 65), "\n")
cat(sprintf("R version: %s\n", R.version$version.string))
cat(sprintf("Platform:  %s\n", R.version$platform))
cat(sprintf("Date:      %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(strrep("-", 65), "\n\n")


# =============================================================================
# SECTION 5: SHARED PROJECT CONSTANTS
# =============================================================================
# These global parameters are used identically across ALL SCRIPTS.
# Sourcing this file into other scripts ensures uniform parameter values.
# Do NOT redefine these in individual scripts. INSTEAD MODIFY THEM HERE IF YOU WANT TO CHANGE THE PARAMETERS. 

# ---------------------------------------------------------------------------
# 5a. EEG Quality Thresholds
# ---------------------------------------------------------------------------
# Specparam model quality thresholds applied at the electrode level before
# computing any summary statistics or running models.
#

EPOCHS_THRESHOLD = 5      # Minimum number of clean epochs required (15s of data)
R2_THRESH        = 0.900  # Minimum Specparam model R-squared (goodness of fit)
MAE_THRESH       = 0.10   # Maximum acceptable Specparam model MAE (Mean Absolute Error)
CH_THRESHOLD     = 12     # Minimum number of good electrodes for subject inclusion

# ---------------------------------------------------------------------------
# 5b. Study Visits (in months). Note: 15 to 18 are called 18 in the manuscript.
# ---------------------------------------------------------------------------
# All visit ages (in months) present in the study.
STUDY_VISITS = c(1, 6, 12, 15, 18, 30, 36, 40, 48)

# ---------------------------------------------------------------------------
# 5c. Participants to Exclude
# ---------------------------------------------------------------------------
# These four subjects are excluded from all analyses due to data quality or because they were pilots
EXCLUDED_SUBJECTS = c("SUB-X2THXQ", "SUB-YVEB6H", "SUB-D9QSB5")

# ---------------------------------------------------------------------------
# 5d. Reproducibility Seed
# ---------------------------------------------------------------------------
RANDOM_SEED = 42  # Used for all bootstrap / random-sampling operations

# ---------------------------------------------------------------------------
# 5e. Colour Palettes
# ---------------------------------------------------------------------------
# COLORS_MAIN      : Two-colour palette (blue and burgundy) for primary comparisons (e.g., burst vs. non-burst containing epochs)
# COLORS_GRADIENT  : 8-colour gradient for visit-group panels
# COLORS_NATURE    : Three-colour palette for derivative significance maps
#                     (Increase = purple, Decrease = orange, No Change = grey)
# COLORS_NATURE_BURST : Same concept but for burst × no-burst interaction models
#                     (Diverging (+) / Diverging (-) / Parallel) - It is equal to the above
# COLORS_COHORT    : Three-colour palette for cohort-level plots 

COLORS_MAIN    = c("#7A1B61", "#A9E0EE")  


COLORS_GRADIENT = rev(c(
  "#3D0A2E",  
  "#7A1B61",  
  "#B23386",  
  "#D364A0", 
  "#B07BD4",   
  "#6E8FE0",  
  "#2DA8C8",  
  "#A9E0EE"   
))


COLORS_NATURE = c(
  "Decrease"   = "#A9E0EE",   
  "Increase"   = "#7A1B61",   
  "No Change"  = "gray90"   
)


COLORS_NATURE_BURST = c(
  "Diverging (-)" = "#A9E0EE",  
  "Diverging (+)" = "#7A1B61",  
  "Parallel"      = "gray90"    
)

COLORS_COHORT = c("brown", "#7687AB", "lightblue")  # Cohorts 1, 2, 3 (Respectively)

# ---------------------------------------------------------------------------
# 5f. Standard ggplot2 Themes
# ---------------------------------------------------------------------------
# THEME_BASE        : Publication theme for all x–y trajectory / bar plots
# THEME_TEXT        : Axis-text sizes that supplement THEME_BASE
# THEME_COMPACT     : Compact theme for MLM / classification heatmap plots
# THEME_TEXT_COMPACT : Axis-text sizes (with 45° x-axis labels) for THEME_COMPACT

THEME_BASE = ggpubr::theme_pubr() +
  theme(
    legend.position       = "bottom",
    legend.direction      = "horizontal",
    legend.box            = "vertical",
    legend.margin         = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin     = margin(t = -5, r = 0, b = 0, l = 0),
    legend.spacing.y      = unit(0, "pt"),
    legend.title          = element_text(face = "bold", size = 12, vjust = 0.5),
    legend.text           = element_text(size = 14),
    legend.key            = element_rect(fill = "transparent", color = NA),
    panel.border          = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks            = element_line(color = "black"),
    axis.ticks.length     = unit(3, "pt"),
    panel.grid            = element_blank(),
    strip.background      = element_blank(),
    strip.text            = element_text(face = "bold", size = 12),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )


THEME_BASE_NOTICKS = ggpubr::theme_pubr() +
  theme(
    legend.position       = "bottom",
    legend.direction      = "horizontal",
    legend.box            = "vertical",
    legend.margin         = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin     = margin(t = -5, r = 0, b = 0, l = 0),
    legend.spacing.y      = unit(0, "pt"),
    legend.title          = element_text(face = "bold", size = 12, vjust = 0.5),
    legend.text           = element_text(size = 14),
    legend.key            = element_rect(fill = "transparent", color = NA),
    panel.border          = element_rect(color = "black", fill = NA, linewidth = 1.2),
    panel.grid            = element_blank(),
    strip.background      = element_blank(),
    strip.text            = element_text(face = "bold", size = 12),
    plot.margin           = margin(t = 5, r = 10, b = 5, l = 5)
  )


THEME_TEXT = theme(
  axis.title.x = element_text(size = 13, face = "bold"),
  axis.title.y = element_text(size = 13, face = "bold"),
  axis.text.x  = element_text(size = 12),
  axis.text.y  = element_text(size = 12)
)

THEME_COMPACT = ggpubr::theme_pubr() +
  theme(
    legend.title          = element_text(angle = 0),
    legend.title.align    = 0.5,
    legend.frame          = element_rect(color = "black", linewidth = 0.5),
    legend.text.position  = "right",
    legend.ticks          = element_line(linewidth = 0.5, color = "black"),
    panel.border          = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks            = element_line(color = "black"),
    axis.ticks.length     = unit(5, "pt"),
    panel.grid            = element_blank(),
    strip.background      = element_blank(),
    strip.text            = element_text(face = "bold", size = 14)
  )

THEME_TEXT_COMPACT = theme(
  axis.title.x     = element_text(size = 14, face = "bold"),
  axis.title.y     = element_text(size = 14, face = "bold"),
  axis.text.x      = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
  axis.text.y      = element_text(size = 13),
  legend.spacing.x = unit(-5, "pt"),
  legend.text      = element_text(size = 11.5)
)


# 95% CI function to use with gtsummary and create the confidence intervals. 
p2.5  <- function(x) quantile(x, probs = 0.025, na.rm = TRUE)
p97.5 <- function(x) quantile(x, probs = 0.975, na.rm = TRUE)

