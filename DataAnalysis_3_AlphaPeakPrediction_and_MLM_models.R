# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term to columns
#   These on-disk CSV column identifiers are intentionally NOT renamed to preserve
#   pipeline integrity and reproducibility. Display labels map here to manuscript text.
#
#   DATA columns to DISPLAY LABELS (used in plot y-axis/legend labels):
#     slope                   → "Slope" (Aperiodic slope)
#     offset                  → "Offset" (Aperiodic offset)
#     peak_freq               → "Peak Freq." (Alpha Peak Frequency, Hz)
#     peak_ampl               → "Peak Amp." (Alpha Peak Amplitude, μV²)
#     peak_prop                → "Prop. of Peaks" (Alpha Peak Proportion, 0–1)
#     osc_ampl               → "Band Power" (Oscillatory Alpha Band Power, μV²)
#     volt_amp / *_corrected  → "Volt. Amp." (Peak-to-peak voltage, absolute/corrected)
#     band_amp / *_corrected  → "Band Amp." (Burst band amplitude, absolute/corrected)
#     prop_bursty_epochs     → "Prop. of Epochs\nw/ Burst" (Proportion of epochs with alpha burst)
#     prop_bursty_cycles_burst → "Prop. of Cycles\nw/ Burst" (Proportion of cycles with alpha burst)
#     avg_burst_duration      → "Burst Duration" (Consecutive cycles per burst)
#     alpha_LAcH              → "Lifespan" (Alpha lifespan in cycles; LAcH cumsum ≥90%)
#     is_burst / burst_type   → "Cycle Type" (categorical: "Burst" vs. "Non-Burst")
# -----------------------------------------------------------------------------
# Script: DataAnalysis_3_AlphaPeakPrediction_and_MLM_models.R
# Purpose: Iterative mixed-effects logistic regression predicting alpha peak
#          presence at 1 month from burst metrics (corrected band power).
#          Section 1: GLMER classification (proc_iter; 1000 iterations)
#          Section 2: LMM extraction (heatmap of beta estimates by visit)
# =============================================================================
# Inputs:
#   - Data/Aperiodic_Oscillatory_ByCycle_Long.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long.csv
#   - Data/BurstProperties_ByCycle_Long.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Scripts/electrodes.csv
# Outputs (saved to path2save and path2figs):
#   - Complementary_Fig_5_Alpha_Peak_Classification_Beta_Coeffs_Table.html
#   - Complementary_Fig_5_Alpha_Peak_Classification_Table.html
#.  - Complementary_Fig_5_Contributions_Burst_Rhythm_to_ParametrizedPSD_**visit**.html
#   - Fig_5_Contributions_Burst_Rhythm_to_ParametrizedPSD.jpeg
# Dependencies: 00_Setup_PackageInstallation.R, 01_Utils_ProcFunctions.R
#   NOTE: proc_graph_ggplot2() and median_se() are defined in
#   01_Utils_ProcFunctions.R and sourced from there. They are no longer
#   defined locally in this script.
# =============================================================================

# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  If you have open the code from the project, you don't need to modify this line. Otherwise, select where the code folder that contains the config_paths.R is

local({
  cand = c(if (nzchar(CODE_FOLDER)) file.path(path.expand(CODE_FOLDER), "config_paths.R"),
           "config_paths.R", "../config_paths.R", "Code/config_paths.R")
  hit  = cand[file.exists(cand)]
  if (!length(hit)) # If the code cannot find the config_path.R it will stop the execution avoiding crashing.
    stop("config_paths.R not found.\n",
         "Fix either way:\n",
         "  1. setwd(\"/path/to/AlphaBurstRhythm/Code\")   then re-run, or\n", # Solution proposed 1: just add the directory of the code
         "  2. set CODE_FOLDER at the top of this script to that same path.\n", #Solution proposed 2: set the CODE_FOLDER variable to the directory of the code
         "Currently looking from: ", getwd(), call. = FALSE)
  source(hit[1], local = FALSE)
})

source(file.path(path2code, "00_Setup_PackageInstallation.R"))

source(file.path(path2code, "SupplementaryTables_Helper.R"))
source(file.path(path2code, "01_Utils_ProcFunctions.R"))


# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================

r2_thresh        = R2_THRESH         # Minimum Specparam R-squared (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)

# pred_model = TRUE:  run GLMER classification loop.
# pred_model = FALSE: skip the classification loop but it will crash.
pred_model   = TRUE
n_bootstraps = 1000  # Iterations for bootstrapped CIs and proc_iter classification

set.seed(RANDOM_SEED)  # Reproducibility seed (42)


# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================
# IMPORTANT: Update only `path2root`. Sub-folders mirror the manuscript (README).
# path2data : ready input data
# path2save : MainText/Tables/AlphaPrediction/ -> peak-prediction + MLM tables
# path2figs : MainText/Figures                 -> Fig 5

path2data = path2sets   # EDIT   # = Data/; this script prepends "" to each file name
# path2root comes from config_paths.R

path2save = file.path(path2root, "MainText", "Tables", "AlphaPrediction")  # Alpha-prediction results folder
path2figs = file.path(path2root, "MainText", "Figures")                    # Fig 5 folder

for (p in c(path2save, path2figs)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
set_flextable_defaults(digits = 3)

# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- EEG cleaning descriptives ---
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    # Recode visit ages: 15 mo. to 18 mo.; 40 mo. to 42 mo. (standardize age labels across cohorts)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic data ---
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS, sujid %in% eeg_desc$sujid) |>
  mutate(
    # Recode visit ages to match eeg_desc harmonization
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  # Normalize GestationalAge_weeks
  mutate(across(where(is.numeric), ~ scale(.x)[, 1])) 

# --- Electrode map ---
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  # Expand electrode region abbreviations to full ROI names for display/analysis
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE           ~ "Central"
  )) |>
  filter(chinclu == 1) # Keep the 60 electrodes of interest

# Combined descriptives frame (left anchor for merging)
descriptives = left_join(
    desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
    eeg_desc
  ) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs)) #Keep only participants with EEG data (prop_epochs is from the EEG descriptive information)

# --- Aperiodic / oscillatory data ---
# 'mae' = Specparam mean absolute error (MAE). goodch (R2 > .900, MAE < .10)
aper_voi = c("mae", "sujid", "session_age", "ch", "region", "chinclu", "epochs",
              "r2value", "alpha_peak", "inclusion_final_dummy",
              "alpha_ampl", "alpha_osc", "alpha_freq", "goodch")

aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv")) |>
  dplyr::select(all_of(aper_voi)) |>
  filter(r2value > r2_thresh, mae < mae_thresh, chinclu == 1,
         epochs >= epochs_threshold, goodch >= CH_THRESHOLD) |>
  mutate(
    # Recode visit ages to match downstream harmonization
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Lagged coherence data ---
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv")) |>
  # Merge with electrodes to filter by channel inclusion flag (cross-validate electrode consistency)
  merge(electrodes |> dplyr::select(label, chinclu), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  mutate(
    # Recode visit ages to match other datasets
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Burst data: corrected amplitude ratios (Burst / NoBurst) ---
# Corrected band amplitude removes non-specific changes from burst amplitude estimate.
burst_data = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv")) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1) |>
  mutate(
    # Recode visit ages
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  group_by(sujid, session_age, ch, region, is_burst) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  # Wide format: separate columns for Burst vs. NoBurst values, enables ratio calculation
  pivot_wider(
    names_from  = is_burst,
    values_from = c(prop_bursty_epochs, avg_burst_duration, volt_amp, band_amp,
                    frequency, sd_frequency)
  ) |>
  rowwise() |>
  mutate(
    # Compute corrected burst metrics as Burst / NoBurst ratios
    band_amp_corrected = band_amp_Burst / band_amp_NoBurst,
    volt_amp_corrected = volt_amp_Burst / volt_amp_NoBurst
  ) |>
  ungroup() |>
  dplyr::select(-any_of(c('avg_burst_duration_NoBurst', 'prop_bursty_epochs_NoBurst'))) |>
  # Rename Burst-specific columns to remove "_Burst" suffix for clean reference in models
  rename(
    avg_burst_duration  = avg_burst_duration_Burst,
    prop_bursty_epochs = prop_bursty_epochs_Burst
  )

# Merge all datasets
data_merged = left_join(descriptives, aperiodic_data) |>
  left_join(hlcoh_data) |>
  left_join(burst_data) |>
  mutate(avg_burst_duration = if_else(is.na(avg_burst_duration), 0, avg_burst_duration)) |>
  filter(!is.na(alpha_peak))

# Variables of interest for this script (band power version)
voi = c("alpha_LAcH", "prop_bursty_epochs", "avg_burst_duration", "band_amp_corrected")

# Age-1 data (primary prediction target: alpha peak presence at 1 month)
data_plot_age1 = data_merged |> filter(session_age == 1)


# =============================================================================
# SECTION 4: EXPLORATORY PLOT (Burst Metrics by Alpha Peak Status at 1 Month)
# =============================================================================

plot50 = data_plot_age1 |>
  pivot_longer(
    cols    = c(alpha_LAcH, band_amp_corrected, avg_burst_duration, prop_bursty_epochs),
    names_to  = "variable",
    values_to = "value"
  ) |>
  # Z-score within electrode to control for topographical variance and enable cross-variable comparison
  mutate(value = scale(value)[, 1], .by = c(variable, ch)) |>
  # Recode variable names to publication-ready labels
  mutate(variable = factor(variable,
    levels = c("alpha_LAcH", "band_amp_corrected", "avg_burst_duration", "prop_bursty_epochs"),
    labels = c("Lifespan", "Corrected\nBand Amp.", "Burst\nDuration", "Prop. of Epochs\nw/ Burst")
  )) |>
  ggplot(aes(x     = variable,
             y     = value,
             color = factor(alpha_peak, labels = c("No Peak", "Peak")),
             fill  = factor(alpha_peak, labels = c("No Peak", "Peak")))) +
  geom_point(alpha = 0.33,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  stat_summary(geom = "bar",      fun = "mean",       position = position_dodge(0.8),
               width = 0.7, alpha = 0.66, linewidth = 1, color = "black") +
  stat_summary(geom = "point",    fun = "median",     position = position_dodge(0.8),
               size = 2, color = "black") +
  stat_summary(geom = "errorbar", fun.data = median_se, position = position_dodge(0.8),
               width = 0.2, linewidth = 0.8, color = "black") +
  THEME_COMPACT + THEME_TEXT_COMPACT +
  scale_color_manual(values = COLORS_MAIN) +
  scale_fill_manual(values  = COLORS_MAIN) +
  labs(y = "Z-Score", color = "Alpha Peak", fill = "Alpha Peak") +
  theme(legend.position = c(0.85, 0.85),
        axis.title.x    = element_blank())


# =============================================================================
# SECTION 5: GLMER CLASSIFICATION (pred_model = TRUE branch)
# =============================================================================
# proc_iter() (from 01_Utils_ProcFunctions.R) runs n_iter train/test splits,
# fits a GLMER at each, and returns AUC, accuracy, sensitivity, and coefficients.

if (pred_model) {

  # --- Prepare age-1 model data (z-scored) ---
  data_model = data_merged |>
    filter(session_age == 1) |>
    dplyr::select(sujid, session_age, alpha_peak, all_of(voi),
                  r2value, prop_epochs, age_months, GestationalAge_weeks) |>
    mutate(across(all_of(c(voi, "r2value", "prop_epochs", "age_months")),
                  ~ scale(.x)[, 1])) |>
    drop_na(alpha_peak) # We need to standarize the rest of covariates. We control by age at 1 month because how fast EEG develops in this period. 

  model_alpha_peak = proc_iter(
    dataset   = data_model,
    n_iter    = n_bootstraps,
    dep_var   = "alpha_peak",
    fixed_eff = voi,
    rand_eff  = "(1|sujid)",
    covariates = "r2value + prop_epochs + age_months + GestationalAge_weeks"
  ) 

  # --- ROC curve plot ---
  plot5a1 = proc_graph_ggplot2(model_alpha_peak$roc) +
    ggpubr::theme_pubr() + THEME_COMPACT + THEME_TEXT_COMPACT +
    labs(title = "") +
    scale_fill_manual(values = COLORS_MAIN[2]) # Call the function to generate a proc curve

  # --- Coefficient bar chart ---
  # Recode predictor names to publication-ready labels for display
  model_estimates = as.data.frame(model_alpha_peak$results) |>
    mutate(term = factor(term,
      levels = c("alpha_LAcH", "band_amp_corrected", "avg_burst_duration", "prop_bursty_epochs"),
      labels = c("Lifespan", "Corrected\nBand Amp.", "Burst\nDuration", "Prop. of Epochs\nw/ Burst")
    ))

  plot5a2 = ggplot(model_estimates, aes(x = term, y = beta)) +
    geom_point(color = "grey50", position = position_jitter(0.15), alpha = 0.33) +
    stat_summary(geom = "bar", fun = "mean", color = "black",
                 position = position_dodge(0.8), width = 0.6, alpha = 0.66,
                 linewidth = 1, fill = COLORS_MAIN[1]) +
    stat_summary(geom = "errorbar", fun.data = median_se, color = "black",
                 position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
    stat_summary(geom = "point", fun = "median", color = "black",
                 position = position_dodge(0.8), size = 2) +
    labs(x = "", y = "Coefficient (\u03B2)") +
    THEME_COMPACT + THEME_TEXT_COMPACT +
    theme(legend.position = "none", axis.title.x = element_blank()) # Generate the code with the estimates

  plot5_full = ggarrange(plot5a1, plot5a2, plot50, ncol = 3,
                          labels = c("a", "b", "c"), widths = c(1, 0.75, 1.05))

  # --- Export classification accuracy table ---
  doc = read_docx()
  model_acc_df = as.data.frame(model_alpha_peak[c("accuracy", "sensitivity",
                                                    "precision", "f1", "brier", "mcc")])
  colnames(model_acc_df) = c("Accuracy", "Sensitivity", "Precision", "F1", "Brier", "MCC")
  
  
  model_acc_df |> 
    tbl_summary(statistic = all_continuous() ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]", digits = all_continuous() ~ 3) |>
    as_flex_table() |> autofit() |> theme_booktabs() |> fontsize(size = 10) |>
    set_header_labels(
      Characteristics = "Classificator Performance Metrics")|>
     add_footer_lines(paste0(
        "Proc_iter classification results for alpha peak presence at 1 month. ",
        "Metrics averaged across ", n_bootstraps, " iterations with random train/test splits. ",
        "95% CIs from iteration distribution. See Methods for model formula and covariates."
      )) |>
      theme_booktabs() |> fix_border_issues() |> autofit() |>
      set_table_properties(layout = "fixed") |>
      align(align = "center", part = "all") |>
      bold(part = "header")|>
      padding(padding = 1, part = "all") |>
      flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
      align(i = 1, part = "header", align = "center") |> 
    save_as_html(path = file.path(path2save, "Complementary_Fig_5_Alpha_Peak_Classification_Table.html"))
  
  # Recode and rename model estimates for table export
  model_alpha_peak$results |>
    mutate(term = factor(term,
      levels = c("alpha_LAcH", "band_amp_corrected", "avg_burst_duration", "prop_bursty_epochs"),
      labels = c("Lifespan", "Corrected band amplitude", "Burst duration",
                 "Prop. epochs with burst"))) |>
    group_by(Predictor = term) |>
    summarise(
      `Beta (mean)`   = mean(beta,  na.rm = TRUE),
      `Beta (s.d.)`   = sd(beta,    na.rm = TRUE),
      `95% CI lower`  = quantile(beta, 0.025, na.rm = TRUE),
      `95% CI upper`  = quantile(beta, 0.975, na.rm = TRUE),
      `s.e. (mean)`   = mean(SE,    na.rm = TRUE),
      `z (mean)`      = mean(z_val, na.rm = TRUE),
      `P (mean)`      = mean(p_val, na.rm = TRUE),
      `Iterations`    = dplyr::n(),
      .groups = "drop") |>
    nn_supp_table(id = "SR5", note = paste(
      "Logistic mixed model classifying the presence of an alpha peak per electrode",
      "at the 1-month visit: Alpha peak (0/1) ~ Corrected band amplitude + Prop.",
      "epochs with burst + Burst duration + Alpha lifespan + covariates + (1|child).",
      "Values summarise 1,000 bootstrap iterations (80/20 train/test split with equal",
      "condition sampling). A predictor was considered reliable when its bootstrap",
      "95% CI excluded zero. This is a Wald z test, so no degrees of freedom apply."))

  # Recode and rename model estimates for table export
  alpha_estimates = model_alpha_peak$results |>
    mutate(term = factor(term,
      levels = c("alpha_LAcH", "band_amp_corrected", "avg_burst_duration", "prop_bursty_epochs"),
      labels = c("Lifespan", "Corrected Band Amp.", "Burst Duration", "Prop. Epochs\nw/ Burst")
    ))|>
    dplyr::select(beta, SE, z_val, p_val, term) |>
    rename(Beta = beta, SE = SE, Z = z_val, P = p_val) |>
    tbl_summary(by = term, statistic = all_continuous() ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]", digits = all_continuous() ~ 3) |>
    as_flex_table() |> autofit() |> theme_booktabs() |> fontsize(size = 10) |>
    set_header_labels(
      Characteristics = "Predictor Coefficients")|>
     add_footer_lines(paste0(
        "Proc_iter classification results for alpha peak presence at 1 month. ",
        "Coefficients averaged across ", n_bootstraps, " iterations with random train/test splits. ",
        "95% CIs from iteration distribution. See Methods for model formula and covariates."
      )) |>
      theme_booktabs() |> fix_border_issues() |> autofit() |>
      set_table_properties(layout = "fixed") |>
      align(align = "center", part = "all") |>
      bold(part = "header")|>
      padding(padding = 1, part = "all") |>
      flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
      align(i = 1, part = "header", align = "center") |>
    save_as_html(path = file.path(path2save, "Complementary_Fig_5_Alpha_Peak_Classification_Beta_Coeffs_Table.html"))

}

  # =============================================================================
  # SECTION 6: LMM EXTRACTION ACROSS ALL VISITS (Heatmap of Beta Estimates)
  # =============================================================================
  # Fits a separate LMM per visit × dependent variable combination; bootstraps CIs.
  # Produces a heatmap summarising burst predictors of alpha parameters over age.

  voi_lmm  = c("alpha_ampl", "alpha_peak", "alpha_osc")
  predictors_lmm = c("band_amp_corrected", "avg_burst_duration",
                      "prop_bursty_epochs", "alpha_LAcH")

  data_merged_region = data_merged |>
    group_by(sujid, session_age, region) |>
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
    group_by(session_age) |>
    mutate(across(where(is.numeric), ~ scale(.x)[, 1])) |>
    ungroup()

  results_list_lmm = list()

  for (age in unique(data_merged_region$session_age)) {
    age_data = filter(data_merged_region, session_age == age)

    for (f in voi_lmm) {
      fixed_formula = as.formula(paste(
        f, "~ band_amp_corrected + avg_burst_duration + prop_bursty_epochs +
             alpha_LAcH + region + r2value + prop_epochs + age_months +
             GestationalAge_weeks + (1|sujid)"
      ))

      m = tryCatch(lmerTest::lmer(fixed_formula, data = age_data),
                   error = function(e) {
                     warning(paste("Model failed for", f, "at age", age)); NULL
                   })

      if (!is.null(m)) {
        model_r2 = performance::r2(m)
        boot_ci  = confint(m, parm = "beta_", method = "boot",
                           nsim = n_bootstraps, oldNames = FALSE)
        boot_ci_df = as.data.frame(boot_ci) |>
          tibble::rownames_to_column("term") |>
          rename(Boot_CI_Lower = `2.5 %`, Boot_CI_Upper = `97.5 %`)

        # ---------------------------------------------------------------
        # Satterthwaite denominator degrees of freedom (lmerTest).
        # Pulled straight from coef(summary(m)). For an `lmerModLmerTest` object the
        # summary coefficient matrix is guaranteed to carry the columns
        # Estimate | Std. Error | df | t value | Pr(>|t|).
        # ---------------------------------------------------------------
        coef_tab = as.data.frame(coef(summary(m)))
        satt_df  = tibble::tibble(
          term = rownames(coef_tab),
          df   = if ("df" %in% names(coef_tab)) coef_tab[["df"]] else NA_real_
        )
        if (all(is.na(satt_df$df))) {
          warning(paste("Satterthwaite df unavailable for", f, "at age", age,
                        "- check that lmerTest::lmer (not lme4::lmer) fitted the model."))
        }

        broom.mixed::tidy(m) |>
          filter(effect == "fixed", term %in% predictors_lmm) |>
          # Drop any `df` broom.mixed may already supply, so the join below
          # cannot silently produce df.x / df.y.
          dplyr::select(-dplyr::any_of("df")) |>
          left_join(boot_ci_df, by = "term") |>
          left_join(satt_df,    by = "term") |>
          mutate(session_age    = age,
                 Dependent_Var  = f,
                 R2_Marginal    = model_r2$R2_marginal,
                 R2_Conditional = model_r2$R2_conditional) -> tidy_m

        results_list_lmm[[paste(age, f, sep = "_")]] = tidy_m
      }
    }
  }

  results_df = bind_rows(results_list_lmm) |>
    group_by(session_age) |>
    mutate(
      # Apply FDR correction (Benjamini-Hochberg) within each visit; create significance label
      p_fdr     = p.adjust(p.value, method = "fdr"),
      sig_label = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**",
                             p_fdr < 0.05 ~ "*",   TRUE          ~ "")
    ) |>
    ungroup() |>
    ungroup() |>
    mutate(
      # Recode predictor names to publication-ready labels
      term_clean    = case_when(
        term == "band_amp_corrected"  ~ "Corrected Band Amp.",
        term == "avg_burst_duration" ~ "Burst Duration",
        term == "prop_bursty_epochs" ~ "Prop. of Epoch\nw/ Burst",
        term == "alpha_LAcH"         ~ "Lifespan",
        TRUE                         ~ gsub("_", " ", term)
      ),
      # Recode dependent variables to publication-ready labels and ensure consistent factor levels
      Dependent_Var = case_when(
        Dependent_Var == "alpha_peak" ~ "Prop. of Peak",
        Dependent_Var == "alpha_ampl" ~ "Peak Amp.",
        Dependent_Var == "alpha_osc"  ~ "Band Power",
        TRUE                          ~ Dependent_Var
      ),
      Dependent_Var = factor(Dependent_Var, levels = c("Prop. of Peak", "Peak Amp.", "Band Power"))
    )

  # Heatmap: beta estimates per predictor × visit, faceted by dependent variable
  heatmap_fig5 = ggplot(results_df, aes(x = factor(session_age), y = term_clean,
                                      fill = estimate)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(estimate, 2), sig_label)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = COLORS_MAIN[2], mid = "white", high = COLORS_MAIN[1],
                         midpoint = 0, name = "\u03B2 Estimate") +
    facet_wrap(~ Dependent_Var, ncol = 3) +
    THEME_COMPACT + THEME_TEXT_COMPACT +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          strip.text  = element_text(size = 12, face = "bold"),
          legend.position = "right") +
    labs(x = "Age (months)", y = "Predictor")

  # Export LMM table
  # Export LMM results per dependent variable with recoded predictor names
  supp_rows = list()   # accumulates the numbered supplementary table

  # NOTE: iterate over UNIQUE dependent variables. Looping over the raw column
  # re-wrote each of the three files once per row of results_df.
  for (d in unique(results_df$Dependent_Var)) {

    tbl_data = results_df |>
      filter(Dependent_Var == d) |>
      # Recode predictor names for table display
      mutate(Predictor = factor(term, levels = c("band_amp_corrected", "avg_burst_duration",
                                                 "prop_bursty_epochs", "alpha_LAcH"),
                                labels = c("Corrected Band Amp.", "Burst Duration",
                                           "Prop. of Epochs\nw/ Burst", "Lifespan"))) |>
      # Create significance symbol column based on FDR-corrected p-value
      mutate(sig = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**",
                             p_fdr < 0.05 ~ "*",   TRUE          ~ "")) |>
      # SE, t and Satterthwaite df are reported so that each contrast can be
      # cited as t(df) = ..., P = ... 
      dplyr::select(session_age, Predictor, Beta = estimate, SE = std.error,
                    t = statistic, df,
                    Boot_CI_Lower, Boot_CI_Upper, p.value, p_fdr, sig,
                    R2_Marginal, R2_Conditional) |>
      arrange(session_age) |>
      mutate(session_age = paste(session_age, 'mo.'))

    # Numbered supplementary table: one file for all three dependent variables,
    # accumulated across loop iterations then written on the last pass.
    supp_rows[[as.character(d)]] = tbl_data |> mutate(`Dependent variable` = as.character(d))

    flextable(tbl_data) |>
      set_header_labels(
        session_age    = "Visit",
        Predictor      = "Predictor",
        Beta           = "\u03B2",
        SE             = "SE",
        t              = "t",
        df             = "df",
        Boot_CI_Lower  = "95% CI Lower",
        Boot_CI_Upper  = "95% CI Upper",
        p.value        = "p-value",
        p_fdr          = "p (FDR)",
        sig            = "",
        R2_Marginal    = "R\u00b2 Marginal",
        R2_Conditional = "R\u00b2 Conditional"
      ) |>
      colformat_double(j = c("Beta", "SE", "Boot_CI_Lower", "Boot_CI_Upper",
                             "R2_Marginal", "R2_Conditional"), digits = 3) |>
      colformat_double(j = "t",  digits = 3) |>
      # Satterthwaite df are fractional, so keep one decimal.
      colformat_double(j = "df", digits = 1) |>
      colformat_double(j = c("p.value", "p_fdr"), digits = 4) |>
      bold(j = "p_fdr", i = ~ p_fdr < 0.05, part = "body") |>
      theme_booktabs() |> autofit() |>
      add_footer_lines(paste0(
        "Model:", d,  " ~  Corrected Band Amplitude + Burst Duration + Prop. of Epoch w/ Burst + Alpha Lifespan + 
        Region + Parametrization Model Fit  + Prop. of Epoch + Age at Visit (months) + Gestational Age (weeks) + (1|sujid). ",
        "Degrees of freedom from Satterthwaite's approximation (lmerTest). ",
        "FDR correction (Benjamini-Hochberg) applied within each visit, across predictors and dependent variables. ",
        "95% CIs from parametric bootstrapping (", n_bootstraps, " iterations). ",
        "*p < .05; **p < .01; ***p < .001."
      )) |>
      theme_booktabs() |> fix_border_issues() |> autofit() |>
      set_table_properties(layout = "fixed") |>
      merge_v(j = "session_age") |>
      align(align = "left", part = "all") |>
      # Centre every column except the first; robust to added/removed columns.
      align(j = 2:ncol(tbl_data), align = "center", part = "all") |>
      bold(part = "header")|>
      padding(padding = 1, part = "all") |>
      flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
      align(i = 1, part = "header", align = "center") |> # Centers the new grouped headers
      save_as_html(path = file.path(path2save, paste0("Complementary_Fig_5_Contributions_Burst_Rhythm_to_ParametrizedPSD_", d, ".html")))
  }

  bind_rows(supp_rows) |>
    dplyr::relocate(`Dependent variable`, .before = 1) |>
    nn_supp_table(id = "SR6", cols = c(
      `Dependent variable` = "Dependent variable",
      session_age          = "Visit (months)",
      Predictor            = "Predictor",
      Beta                 = "Beta",
      SE                   = "s.e.",
      t                    = "t",
      df                   = "df",
      Boot_CI_Lower        = "95% CI lower",
      Boot_CI_Upper        = "95% CI upper",
      p.value              = "P",
      p_fdr                = "P (FDR)",
      R2_Marginal          = "Marginal R2",
      R2_Conditional       = "Conditional R2"),
      note = paste0(
        "Visit-stratified linear mixed models: Alpha metric ~ Corrected Band Amp. + Burst duration + ",
        "Prop. of epochs with burst + Alpha lifespan + ROI + covariates + (1|child), ",
        "fitted separately per visit. All predictors z-scored within visit, so betas are ",
        "standardised and comparable across predictors. Degrees of freedom from ",
        "Satterthwaite's approximation. 95% CIs from parametric bootstrapping ",
        "(1,000 iterations). P values FDR-corrected (Benjamini-Hochberg) within visit ",
        "across predictors and dependent variables."))

  # Combine classification and LMM panels
  plot_full = ggarrange(plot5_full, heatmap_fig5, labels = c("", "d"),
                         nrow = 2, heights = c(1, 0.75))
  ggsave(plot = plot_full,
         filename = file.path(path2figs, "Fig_5_Contributions_Burst_Rhythm_to_ParametrizedPSD.jpeg"),
         width = 30, height = 20, dpi = 300, units = "cm")


