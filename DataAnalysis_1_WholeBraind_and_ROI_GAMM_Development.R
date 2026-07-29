# =============================================================================
# CODE Notes
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
#     osc_ampl                → "Band Power" (Oscillatory Alpha Band Power, μV²)
#     volt_amp / *_corrected  → "Volt. Amp." (Peak-to-peak voltage, absolute/corrected)
#     band_amp / *_corrected  → "Band Amp." (Burst band amplitude, absolute/corrected)
#     prop_bursty_epochs      → "Prop. of Epochs\nw/ Burst" (Proportion of epochs with alpha burst)
#     prop_bursty_cycles_burst → "Prop. of Cycles\nw/ Burst" (Proportion of cycles with alpha burst)
#     avg_burst_duration      → "Burst Duration" (Consecutive cycles per burst)
#     alpha_LAcH              → "Lifespan" (Alpha lifespan in cycles; LAcH cumsum ≥90%)
#     is_burst / burst_type   → "Cycle Type" (categorical: "Burst" vs. "Non-Burst")
# -----------------------------------------------------------------------------
# Script: DataAnalysis_1_WholeBraind_and_ROI_GAMM_Development.R
# Purpose: Fits GAMM models for developmental trajectories of all EEG metrics
#          (whole brain and regional). Computes simultaneous derivative-based
#          "moments of significant change", and exports summary tables and plots.
# =============================================================================
# Inputs:
#   - Data/Aperiodic_Oscillatory_ByCycle_Long.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long.csv
#   - Data/BurstProperties_ByCycle_Long.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Data/electrodes.csv
# Outputs (all saved to path2tabs either main or supplement):
#   - Complementary_Fig2_3_Table_WholeBrain_GAMM_Developtment_**variabletype**.html
#   - Complementary_Fig2_3_Table_WholeBrain_GAMM_Developtment_BurstInteraction.html
#   - Supplementary_Results_FigR1_Table_ROIs_Trajectories_**variabletype**.html
#   - FigSR1_GAMM_ROIs_Development_all_metrics.jpeg
#   - WholeBrain_MomentsOfChange_**metric**.csv # These csv files contain the derivatives of the GAMM smooths and identify periods of significant change.
#   - WholeBrain_Interaction_Moments_**metric**.csv These csv files contain the derivatives of the difference smooths and identify periods of significant divergence between Burst and NoBurst trajectories.
# Dependencies: 00_Setup_PackageInstallation.R
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



# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================
# NOTE: Themes and colours are loaded from 00_Setup_PackageInstallation.R.
# Do NOT redefine THEME_BASE, THEME_TEXT, COLORS_* here.

r2_thresh        = R2_THRESH         # Minimum FOOOF R-squared (0.900)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)
mae_thresh       = MAE_THRESH        # Max Specparam MAE (mean absolute error)
set.seed(RANDOM_SEED)
# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================
# IMPORTANT: Update only `path2root` to match your local directory.
# Sub-folders derive from it and mirror the manuscript organisation (see README).
# path2data : merged/analysis-ready input data (read-only here)
# path2tabs : Tables/MainText/Development/   -> whole-brain GAMM dev. table
# path2figs : Figures/SupplementaryResults/  -> Fig SR1 (ROI GAMM trajectories)

# path2data comes from config_paths.R
# path2root comes from config_paths.R
# path2code comes from config_paths.R
source(file.path(path2code, "SupplementaryTables_Helper.R"))

path2tabs    = file.path(path2root, "MainText",  "Tables", "Development")  # GAMM development results
path2suppres = file.path(path2root, "SupplementaryInformation", "Results", "Tables")   # ROI-stratified GAMM tables
path2figs    = file.path(path2root, "SupplementaryInformation", "Results", "Figures")     # Fig SR1 lives in Supp. Results

if (!dir.exists(path2tabs)) dir.create(path2tabs, recursive = TRUE)
if (!dir.exists(path2suppres)) dir.create(path2suppres, recursive = TRUE)
if (!dir.exists(path2figs)) dir.create(path2figs, recursive = TRUE)

# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- EEG cleaning descriptives (defines the analysis sample) ---
# prop_epochs is z-scored so it has a comparable scale when used as a GAMM covariate.
# inclusion_lmm = 1 marks participants that pass the minimum number of visits to be included in the analysis. 
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic / longitudinal age data ---
# age_months = exact chronological age at the visit (used as the GAMM smooth).
# SES covariates are z-scored
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS, sujid %in% eeg_desc$sujid) |>
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  # ITN_mean is mean-imputed only if present. It is not a covariate in any
  # model here; it is swept in by select(contains("mean")) and is absent
  # from the public data release (see prepare_public_data.R). any_of()
  # makes this a no-op when the column is not there.
  mutate(across(any_of("ITN_mean"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x))) |>
  mutate(across(where(is.numeric), ~ scale(.x)[, 1]))

# --- Electrode map ---
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  # Harmonize electrode region abbreviations to full region names for consistent labeling across pipeline
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1) # Chinclu == 1 retains only the 60 electrodes of interest

# Combined descriptives frame (left anchor for all data merges)
descriptives = left_join(
    desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
    eeg_desc
  ) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs))

# --- Aperiodic / oscillatory data ---
# Filters applied: R² > threshold, MAE < threshold, included channels,
aper_voi = c("sujid", "session_age", "ch", "region", "chinclu", "epochs", "r2value",
              "mae", "goodch", "offset", "slope", "alpha_freq", "alpha_ampl",
              "alpha_osc", "alpha_peak", "inclusion_final_dummy") #Variables of interst

aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv")) |>
  dplyr::select(all_of(aper_voi)) |>
  group_by(session_age, sujid, region) |>
  filter(r2value > r2_thresh, mae < mae_thresh, chinclu == 1,
         epochs >= epochs_threshold, goodch > CH_THRESHOLD) |>
  summarise(
    # Average all aperiodic/oscillatory metrics within each participant-visit-region cell
    slope     = mean(slope,                       na.rm = TRUE),
    offset    = mean(offset,                      na.rm = TRUE),
    peak_ampl = mean(alpha_ampl,                  na.rm = TRUE),
    peak_freq = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),
    osc_ampl  = mean(alpha_osc,                   na.rm = TRUE),
    peak_prop  = mean(alpha_peak,                  na.rm = TRUE),
    epochs    = mean(epochs,                      na.rm = TRUE),
    r2value   = mean(r2value,                     na.rm = TRUE),
    goodch    = mean(goodch,                      na.rm = TRUE),
    .groups   = "drop"
  ) |> # Harmonize session age coding: correct labeling inconsistencies (15 to 18, 40 to 42 months)
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Lagged coherence data ---
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv")) |>
  merge(electrodes |> dplyr::select(label, chinclu, region), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  group_by(session_age, sujid, region) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop") |>
  # Harmonize session age coding to align with scanner protocol ages
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Burst properties (whole-brain level for Part A and B) ---
# Both Burst and NoBurst rows are retained; is_burst identifies the cycle type.
burst_data = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv")) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1) |>
  group_by(session_age, sujid, is_burst, region) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  # Rename 'region' to 'region' for consistent column naming; is_burst layer already differentiates Burst vs NoBurst
  rename(region = region) |>
  # Harmonize session age coding
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Burst properties (regional level for Part C) ---
# Computes corrected amplitude ratios (Burst / NoBurst) for regional models.
burst_voi = c("volt_amp", "band_amp", "frequency",
              "prop_bursty_cycles_burst", "prop_bursty_epochs", "avg_burst_duration")

burst_data_regional = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv")) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1) |>
  dplyr::select(sujid, session_age, ch, region, is_burst, all_of(burst_voi)) |>
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  group_by(sujid, session_age, is_burst, region) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  # Pivot is_burst wide to separate Burst and NoBurst measurements side-by-side for ratio computation
  pivot_wider(names_from = is_burst,
              values_from = c(volt_amp, band_amp, frequency,
                               prop_bursty_cycles_burst, prop_bursty_epochs, avg_burst_duration)) |>
  # Drop NoBurst proportions (undefined ratios; keep only Burst-exclusive metrics)
  dplyr::select(-c(prop_bursty_cycles_burst_NoBurst, prop_bursty_epochs_NoBurst,
                   avg_burst_duration_NoBurst)) |>
  # Recode Burst-conditional columns to root names (prop_bursty_cycles_burst_Burst -> prop_bursty_cycles_burst)
  # for use in regional GAMM models as "burst properties conditional on burst presence"
  rename(prop_bursty_cycles_burst = prop_bursty_cycles_burst_Burst,
         prop_bursty_epochs     = prop_bursty_epochs_Burst,
         avg_burst_duration      = avg_burst_duration_Burst) |>
  # Compute amplitude correction factors (Burst / NoBurst ratios) to normalize burst energy by cycle type
  mutate(
    volt_amp_corrected = volt_amp_Burst / volt_amp_NoBurst,
    band_amp_corrected = band_amp_Burst / band_amp_NoBurst
  )

# --- Merge EEG datasets with descriptives ---
# Excluded subjects and those not meeting the inclusion criteria are removed.
aperiodic_data      = left_join(descriptives, aperiodic_data)      |> filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
hlcoh_data          = left_join(descriptives, hlcoh_data)          |> filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
burst_data          = left_join(descriptives, burst_data)          |> filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
burst_data_regional = left_join(descriptives, burst_data_regional) |> filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)


# =============================================================================
# PART A: WHOLE-BRAIN MAIN AGE EFFECTS
# =============================================================================
# Fits GAMM: pow ~ s(age_months, k=4) + region + prop_epochs + Cohort +
#                  GestationalAge_weeks [+ model_fit if available] - Except for the offset (K = 3)
# Compares against a reduced model without the age smooth to compute partial R².
# Derivatives computed via gratia::derivatives() flag periods of significant change.

metrics_config = list(
  list(metric = "slope",                   type = "Aperiodic",   data = "aperiodic_data",  label = 'Slope', tablesub = 'powerspectrum'),
  list(metric = "offset",                  type = "Aperiodic",   data = "aperiodic_data",  label = 'Offset', tablesub = 'powerspectrum'),
  list(metric = "peak_freq",               type = "Oscillatory", data = "aperiodic_data",  label = 'Peak Freq.', tablesub = 'powerspectrum'),
  list(metric = "osc_ampl",                type = "Oscillatory", data = "aperiodic_data",  label = 'Band Power', tablesub = 'powerspectrum'),
  list(metric = "peak_prop",                type = "Oscillatory", data = "aperiodic_data",  label = 'Prop. of Peaks',      tablesub = 'powerspectrum'),
  list(metric = "peak_ampl",               type = "Oscillatory", data = "aperiodic_data",  label = 'Peak Amp.',  tablesub = 'powerspectrum'),
  list(metric = "alpha_LAcH",              type = "Oscillatory", data = "hlcoh_data",      label = 'Lifespan', tablesub = 'alphalifespan'),
  list(metric = "prop_bursty_cycles_burst", type = "Oscillatory", data = "burst_data",      label = 'Prop. of Cycles\nw/ Burst', tablesub = 'burst'),
  list(metric = "prop_bursty_epochs",     type = "Oscillatory", data = "burst_data",      label = 'Prop. of Epochs\nw/ Burst', tablesub = 'burst'),
  list(metric = "avg_burst_duration",      type = "Oscillatory", data = "burst_data",      label = 'Burst Duration', tablesub = 'burst')
)

results_list_main = list()
diag_list           = list()   # gam.check / k.check diagnostics -> Table SR12  # Stores one stats row per metric
main_plots        = list()  # Stores one ggplot per metric

print("--- PART A: WHOLE-BRAIN MAIN AGE EFFECTS ---")

for (item in metrics_config) {
  f      = item$metric
  p      = item$type
  d_name = item$data
  flabel = item$label
  stable = item$tablesub

  print(paste("Running model for:", p, "-", f))

  # Select and prepare the correct dataset for this metric
  if (d_name == "aperiodic_data") {
    model_data = aperiodic_data |> rename(pow = !!sym(f)) |>
      mutate(model_fit = scale(r2value)[, 1]) |> drop_na(pow)
  } else if (d_name == "hlcoh_data") {
    model_data = hlcoh_data    |> rename(pow = !!sym(f)) |> drop_na(pow)
  } else if (d_name == "burst_data") {
    # For burst presence and duration metrics: keep only Burst rows; replace NA with 0
    model_data = burst_data    |> rename(pow = !!sym(f)) |>
      filter(is_burst == "Burst") |>
      mutate(pow = if_else(is.na(pow), 0, pow)) |> drop_na(pow)
  }

  # Drop rows with missing numeric covariates
  numeric_covars = c("age_months", "prop_epochs", "GestationalAge_weeks")
  if ("model_fit" %in% names(model_data)) numeric_covars = c(numeric_covars, "model_fit")
  model_data = model_data |> drop_na(any_of(numeric_covars))

  n_unique_ids = length(unique(model_data$sujid))
  total_obs    = model_data |>
    group_by(sujid, session_age) |> summarise(.groups = "drop") |> nrow()

  if (n_unique_ids > 1 && nrow(model_data) > 2 * n_unique_ids) {

    if (f == "offset") { # Offset with K > 3 produced negative r2, which indicates a poor model.
      fixed_effects = "pow ~ s(age_months, k=3) + region + prop_epochs + Cohort + GestationalAge_weeks"
    } else { 
      fixed_effects = "pow ~ s(age_months, k=4) + region + prop_epochs + Cohort + GestationalAge_weeks"
      }
    
    # Build GAMM formula; add model_fit covariate only when available (aperiodic data)
    if ("model_fit" %in% names(model_data)) fixed_effects = paste(fixed_effects, "+ model_fit")

    fixed_effects_red = "pow ~ region + prop_epochs + Cohort + GestationalAge_weeks"
    if ("model_fit" %in% names(model_data)) fixed_effects_red = paste(fixed_effects_red, "+ model_fit")

    # Full model (with age smooth)
    m = gamm(
      as.formula(fixed_effects),
      random  = list(sujid =~ 1),
      correlation = nlme::corCAR1(form = ~ age_months | sujid / region),
      data    = model_data, method = "REML",
      control = nlme::lmeControl(opt = "optim", niterEM = 100, maxIter = 10000, msMaxIter = 10000)
    )

    # Reduced model (age smooth removed) — used to compute partial R² for age
    m_red = gamm(
      as.formula(fixed_effects_red),
      random  = list(sujid =~ 1),
      correlation = nlme::corCAR1(form = ~ age_months | sujid / region), 
      data    = model_data, method = "REML",
      control = nlme::lmeControl(opt = "optim", niterEM = 100, maxIter = 10000, msMaxIter = 10000)
    )

    # Partial R²: (RSS_reduced − RSS_full) / RSS_reduced
    partial_r2  = (sum(residuals(m_red$lme)^2) - sum(residuals(m$lme)^2)) /
                   sum(residuals(m_red$lme)^2)

    gam_results = summary(m$gam)

    # Compute simultaneous derivatives to identify periods of significant change
    moments_of_change = derivatives(m$gam, interval = "simultaneous") |>
      mutate(inc_or_dec = case_when(
        `.upper_ci` >= 0 & `.lower_ci` >= 0 ~ "Increase",
        `.upper_ci` <= 0 & `.lower_ci` <= 0 ~ "Decrease",
        TRUE                             ~ "No Change"
      ))
    write_csv(moments_of_change,
              paste0(path2data, "/WholeBrain_MomentsOfChange_", p, "_", f, ".csv"))

    # Store results row
    results_list_main[[f]] = data.frame(
      stable          = stable,
      type_of_pow     = p,
      metric          = flabel,
      n               = n_unique_ids,
      total_obs       = total_obs,
      age_edf         = gam_results$s.table["s(age_months)", "edf"],
      age_refdf       = gam_results$s.table["s(age_months)", "Ref.df"],
      age_F           = gam_results$s.table["s(age_months)", "F"],
      age_pval        = gam_results$s.table["s(age_months)", "p-value"],
      epochs_estimate = gam_results$p.table["prop_epochs", "Estimate"],
      epochs_SE       = gam_results$p.table["prop_epochs", "Std. Error"],
      epochs_t        = gam_results$p.table["prop_epochs", "t value"],
      epochs_df       = nn_gam_df(m$gam),
      epochs_pval     = gam_results$p.table["prop_epochs", "Pr(>|t|)"],
      partial_r2      = partial_r2
    )
  } else {
    print(paste("SKIPPING:", f, "(Insufficient data)"))
  }
}


# =============================================================================
# PART B: BURST × AGE INTERACTION EFFECTS
# =============================================================================
# Fits a difference-smooth GAMM: s(age_months) + s(age_months, by=is_burst_ord)
# The difference smooth tests whether Burst and NoBurst trajectories diverge over age.

inter_voi          = c("volt_amp", "band_amp")
results_list_inter = list()

print("--- PART B: BURST × AGE INTERACTION EFFECTS ---")

for (f in inter_voi) {
  
  if (f == "volt_amp") {flabel = "Burst Peak-to-Peak Amp."} else {flabel = "Burst Band Amp."}  # Table 1 absolute burst-energy metrics
  stable = 'burst'

  print(paste("Running interaction for:", f))

  model_data = burst_data |> rename(pow = !!sym(f))

  # Create ordered factor for is_burst to enable difference smooth testing burst trajectory divergence
  model_data = model_data |> mutate(is_burst = factor(is_burst))
  model_data$is_burst_ord = as.ordered(model_data$is_burst)
  contrasts(model_data$is_burst_ord) = "contr.treatment"

  model_data = model_data |> drop_na(any_of(c("age_months", "prop_epochs", "GestationalAge_weeks")))

  
  n_unique_ids = length(unique(model_data$sujid))
  total_obs    = model_data |>
    group_by(sujid, session_age) |> summarise(.groups = "drop") |> nrow()
  
  
  # Full model with difference smooth
  m = gamm(
    pow ~ is_burst_ord + s(age_months, k=4) + s(age_months, k=4, by = is_burst_ord) +
          region + prop_epochs + Cohort + GestationalAge_weeks,
    random  = list(sujid =~ 1),
    correlation = nlme::corCAR1(form = ~ age_months  | sujid / region / is_burst_ord), # We need to specify all the clusters, so we go with burst nested within region nested within participant
    data    = model_data, method = "REML",
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100)
  )

  # Reduced model (no difference smooth) for partial R²
  m_red = gamm(
    pow ~ is_burst_ord + region + prop_epochs + Cohort + GestationalAge_weeks,
    random  = list(sujid =~ 1),
    correlation = nlme::corCAR1(form = ~ age_months | sujid / region / is_burst_ord),
    data    = model_data, method = "REML",
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100)
  )

  partial_r2  = (sum(residuals(m_red$lme)^2) - sum(residuals(m$lme)^2)) /
                 sum(residuals(m_red$lme)^2)
  gam_results = summary(m$gam)

  # Extract name of the by-smooth term (difference smooth s(age_months, by=is_burst_ord)) for derivative testing
  term_diff_name = grep("is_burst_ord", rownames(gam_results$s.table), value = TRUE)

  # Derivatives of the difference smooth
  moments_diff = derivatives(m$gam, term = term_diff_name, interval = "simultaneous") |>
    mutate(inc_or_dec = case_when(
      .upper_ci >= 0 & .lower_ci >= 0 ~ "Diverging (+)",
      .upper_ci <= 0 & .lower_ci <= 0 ~ "Diverging (-)",
      TRUE                             ~ "Parallel"
    ))
  write_csv(moments_diff, paste0(path2data, "/WholeBrain_Interaction_Moments_", f, ".csv"))

  # Plot the difference smooth (departure from 0 = growing divergence)
  plot_smooth_data = smooth_estimates(m$gam, select = term_diff_name) |> add_confint()

  # Store interaction results
  results_list_inter[[f]] = data.frame(
    metric          = flabel,
    stable          = stable,
    n               = n_unique_ids,
    total_obs       = total_obs,
    diff_edf        = gam_results$s.table[term_diff_name, "edf"],
    diff_refdf      = gam_results$s.table[term_diff_name, "Ref.df"],
    diff_F          = gam_results$s.table[term_diff_name, "F"],
    diff_pval       = gam_results$s.table[term_diff_name, "p-value"],
    burst_estimate  = gam_results$p.table["is_burst_ordNoBurst", "Estimate"],
    burst_SE        = gam_results$p.table["is_burst_ordNoBurst", "Std. Error"],
    burst_t         = gam_results$p.table["is_burst_ordNoBurst", "t value"],
    burst_df        = nn_gam_df(m$gam),
    burst_pval      = gam_results$p.table["is_burst_ordNoBurst", "Pr(>|t|)"],
    partial_r2      = partial_r2
  )

  # Also add a main age-effect row for burst metrics in the main results table
  results_list_main[[f]] = data.frame(
    type_of_pow     = "Oscillatory",
    metric          = flabel,
    stable          = stable, 
    n = n_unique_ids,
    total_obs       = total_obs,
    age_edf         = gam_results$s.table["s(age_months)", "edf"],
    age_refdf       = gam_results$s.table["s(age_months)", "Ref.df"],
    age_F           = gam_results$s.table["s(age_months)", "F"],
    age_pval        = gam_results$s.table["s(age_months)", "p-value"],
    epochs_estimate = gam_results$p.table["prop_epochs", "Estimate"],
    epochs_SE       = gam_results$p.table["prop_epochs", "Std. Error"],
    epochs_t        = gam_results$p.table["prop_epochs", "t value"],
    epochs_df       = nn_gam_df(m$gam),
    epochs_pval     = gam_results$p.table["prop_epochs", "Pr(>|t|)"],
    partial_r2      = partial_r2
  )
}


# =============================================================================
# PART C: TABLE EXPORT — WHOLE-BRAIN
# =============================================================================

# --- Table 1: Whole-Brain Main Age Effects ---
if (length(results_list_main) > 0) {
  
  results_list_main = bind_rows(results_list_main)
  # --- Numbered supplementary table SR1 (one file covering all model families) ---
  results_list_main |>
    group_by(stable) |>
    mutate(age_pval_fdr    = p.adjust(age_pval,    method = "fdr"),
           epochs_pval_fdr = p.adjust(epochs_pval, method = "fdr")) |>
    ungroup() |>
    arrange(stable, type_of_pow, metric) |>
    nn_supp_table(id = "SR1", cols = c(
      stable          = "Model family",
      type_of_pow     = "Component",
      metric          = "Metric",
      n               = "n (children)",
      total_obs       = "n (observations)",
      age_edf         = "EDF",
      age_refdf       = "Ref. df",
      age_F           = "F",
      age_pval        = "P (age)",
      age_pval_fdr    = "P (age, FDR)",
      epochs_estimate = "Beta (prop. epochs)",
      epochs_SE       = "s.e.",
      epochs_t        = "t",
      epochs_df       = "df",
      epochs_pval     = "P (epochs)",
      epochs_pval_fdr = "P (epochs, FDR)",
      partial_r2      = "Partial R2"),
      note = paste(
        "Generalized additive mixed models: Metric ~ s(Age) + ROI + Cohort +",
        "proportion of retained epochs + gestational age (+ model fit for",
        "parameterized metrics), with a continuous-time AR1 correlation structure",
        "nested by region within child. The age smooth is tested with an F statistic",
        "on the reference degrees of freedom; the proportion-of-epochs covariate is",
        "a t test on the model residual degrees of freedom. P values are",
        "FDR-corrected (Benjamini-Hochberg) across metrics within each model family.",
        "Partial R2 = (RSS_reduced - RSS_full)/RSS_reduced.",
        "EDF, effective degrees of freedom."))

  for (t in unique(results_list_main$stable)) {
    results_list_main|>
      filter(stable == t) |>
      dplyr::select(-stable)|>
    mutate(
      age_pval_fdr    = p.adjust(age_pval,    method = "fdr"),
      epochs_pval_fdr = p.adjust(epochs_pval, method = "fdr")
    ) |>
      dplyr::select(type_of_pow, metric, n, total_obs, age_edf, age_F, age_pval, age_pval_fdr,
                    epochs_estimate, epochs_SE, epochs_t, epochs_pval, epochs_pval_fdr, partial_r2) |>
    flextable() |>
    set_header_labels(
      type_of_pow     = "Type",        metric          = "Metric",
      n               = "N",           total_obs       = "N Obs.",
      age_edf         = "EDF",         age_F           = "F",
      age_pval        = "p (age)",     age_pval_fdr    = "p (age, FDR)",
      epochs_estimate = "\u03B2 (epochs)", epochs_SE   = "SE",
      epochs_t        = "t",           epochs_pval     = "p (epochs)",
      epochs_pval_fdr = "p (epochs, FDR)", partial_r2  = "Partial R\u00b2"
    ) |>
    set_caption("GAMM Whole-Brain Age Effects — Main Developmental Trajectories") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("age_pval", "age_pval_fdr", "epochs_pval", "epochs_pval_fdr"), digits = 4) |>
    bold(j = "age_pval_fdr", i = ~ age_pval_fdr < 0.05, part = "body") |>
    add_footer_lines("EDF = Estimated degrees of freedom of the age smooth (>1 = non-linear). Partial R\u00b2 = (RSS_red \u2212 RSS_full) / RSS_red. FDR = Benjamini-Hochberg correction applied across all metrics.") |>
    theme_booktabs() |> autofit() |>
      merge_v(j = "type_of_pow") |>
      set_table_properties(layout = "fixed") |>
      valign(j = "metric", valign = "top") |>
      align(align = "left", part = "all") |>
      align(j = 3:12, align = "center", part = "all") |>
      bold(part = "header")|>
      padding(padding = 1, part = "all") |>
      flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    save_as_html(path = file.path(path2tabs, paste0("Complementary_Fig2_3_Table_WholeBrain_GAMM_Developtment_", t, ".html")))
  }
}

# --- Table 2: Burst Interaction Effects ---
if (length(results_list_inter) > 0) {
  # --- Numbered supplementary table SR2 (full frame: includes Ref. df and df) ---
  bind_rows(results_list_inter) |>
    mutate(diff_pval_fdr = p.adjust(diff_pval, method = "fdr")) |>
    nn_supp_table(id = "SR2", cols = c(
      metric         = "Metric",
      diff_edf       = "EDF (difference smooth)",
      diff_refdf     = "Ref. df",
      diff_F         = "F",
      diff_pval      = "P",
      diff_pval_fdr  = "P (FDR)",
      burst_estimate = "Beta (burst vs non-burst)",
      burst_SE       = "s.e.",
      burst_t        = "t",
      burst_df       = "df",
      burst_pval     = "P (burst)",
      partial_r2     = "Partial R2"),
      note = paste(
        "Difference-smooth GAMMs: Amplitude ~ Cycle type + s(Age, by = Cycle type)",
        "+ ROI + Cohort + covariates. The difference smooth tests whether burst and",
        "non-burst trajectories diverge with age and is evaluated with an F statistic",
        "on the reference degrees of freedom; the cycle-type contrast is a t test on",
        "the model residual degrees of freedom. P values FDR-corrected",
        "(Benjamini-Hochberg) across metrics. EDF, effective degrees of freedom."))

  bind_rows(results_list_inter) |>
    mutate(diff_pval_fdr = p.adjust(diff_pval, method = "fdr")) |>
    dplyr::select(metric, diff_edf, diff_F, diff_pval, diff_pval_fdr,
                  burst_estimate, burst_SE, burst_t, burst_pval, partial_r2) |>
    flextable() |>
    set_header_labels(
      metric          = "Metric",
      diff_edf        = "EDF (Difference)",    diff_F         = "F",
      diff_pval       = "p (difference)",      diff_pval_fdr  = "p (FDR)",
      burst_estimate  = "\u03B2 (Burst vs. NoBurst)", burst_SE  = "SE",
      burst_t         = "t",                   burst_pval     = "p (burst)",
      partial_r2      = "Partial R\u00b2"
    ) |>
    set_caption("GAMM Burst × Age Interaction — Trajectory Divergence (Burst vs. Non-Burst)") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("diff_pval", "diff_pval_fdr", "burst_pval"), digits = 4) |>
    bold(j = "diff_pval_fdr", i = ~ diff_pval_fdr < 0.05, part = "body") |>
    add_footer_lines("Difference smooth = s(age, by=is_burst_ord): tests whether Burst and NoBurst trajectories diverge over development. FDR = Benjamini-Hochberg correction applied across both metrics.") |>
    theme_booktabs() |> autofit() |>
    merge_v(j = "metric") |>
    set_table_properties(layout = "fixed") |>
    valign(j = "metric", valign = "top") |>
    align(align = "left", part = "all") |>
    align(j = 2:9, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    save_as_html(path = file.path(path2tabs, "Complementary_Fig2_3_Table_WholeBrain_GAMM_Developtment_BurstInteraction.html"))
}

print("--- WHOLE-BRAIN MODELS COMPLETE ---")


# =============================================================================
# PART D: REGIONAL AGE TRAJECTORIES
# =============================================================================
# Fits per-region age smooths: s(age_months, k=4, by=region)
# Extracts regional intercept differences and regional smooth significance.

metrics_config_reg = list(
  list(metric = "slope",                   type = "Aperiodic",   data = "aperiodic_data",    plot_label = "Slope", stab = "powerspectrum"),
  list(metric = "offset",                  type = "Aperiodic",   data = "aperiodic_data",    plot_label = "Offset", stab = "powerspectrum"),
  list(metric = "peak_freq",               type = "Oscillatory", data = "aperiodic_data",    plot_label = "Peak Freq.", stab = "powerspectrum"),
  list(metric = "osc_ampl",                type = "Oscillatory", data = "aperiodic_data",    plot_label = "Band Power", stab = "powerspectrum"),
  list(metric = "peak_prop",                type = "Oscillatory", data = "aperiodic_data",    plot_label = "Prop. of Peaks", stab = "powerspectrum"),
  list(metric = "peak_ampl",               type = "Oscillatory", data = "aperiodic_data",    plot_label = "Peak Amp.", stab = "powerspectrum"),
  list(metric = "alpha_LAcH",              type = "Oscillatory", data = "hlcoh_data",        plot_label = "Lifespan", stab = "alphalifespan"),
  list(metric = "prop_bursty_cycles_burst", type = "Oscillatory", data = "burst_data_regional", plot_label = "Prop. of Cycles\nw/ Burst", stab = "burst"),
  list(metric = "prop_bursty_epochs",     type = "Oscillatory", data = "burst_data_regional", plot_label = "Prop. of Epochs\nw/ Burst", stab = "burst"),
  list(metric = "volt_amp_corrected",      type = "Oscillatory", data = "burst_data_regional", plot_label = "Corrected Volt. Amp.", stab = "burst"),
  list(metric = "band_amp_corrected",      type = "Oscillatory", data = "burst_data_regional", plot_label = "Corrected Band Amp.", stab = "burst"),
  list(metric = "avg_burst_duration",      type = "Oscillatory", data = "burst_data_regional", plot_label = "Burst Duration", stab = "burst")
)

results_list_region = list()
results_list_reg_smooth = list()
region_plots        = list()

print("--- PART D: REGIONAL TRAJECTORY MODELS ---")

for (item in metrics_config_reg) {
  f       = item$metric
  p       = item$type
  d_name  = item$data
  plabel  = item$plot_label
  stable  = item$stab

  # Select and prepare dataset
  if (d_name == "aperiodic_data") {
    model_data = aperiodic_data      |> rename(pow = !!sym(f)) |>
      mutate(model_fit = scale(r2value)[, 1]) |> drop_na(pow)
  } else if (d_name == "hlcoh_data") {
    model_data = hlcoh_data          |> rename(pow = !!sym(f)) |> drop_na(pow)
  } else if (d_name == "burst_data_regional") {
    model_data = burst_data_regional |> rename(pow = !!sym(f))
    if (f %in% c("prop_bursty_cycles_burst", "prop_bursty_epochs", "avg_burst_duration")) {
      model_data = model_data |> mutate(pow = if_else(is.na(pow), 0, pow)) |> drop_na(pow)
    } else {
      model_data = model_data |> drop_na(pow)
    }
  }

  model_data$region = as.factor(model_data$region)
  numeric_covars    = c("age_months", "prop_epochs", "GestationalAge_weeks")
  if ("model_fit" %in% names(model_data)) numeric_covars = c(numeric_covars, "model_fit")
  model_data = model_data |> drop_na(any_of(numeric_covars))

  n_unique_ids = length(unique(model_data$sujid))

  if (nrow(model_data) > 100) {

    #We are interested in fitting 5 different smooth curves, so we do not order the factor
    if ( f == "offset") { # Offset with K > 3 produced negative r2, which indicates a poor model.
      fixed_effects = "pow ~ region + s(age_months, k=3, by = region) + prop_epochs + Cohort + GestationalAge_weeks"
      fixed_effects_red = "pow ~ region + s(age_months, k=3) + prop_epochs + Cohort + GestationalAge_weeks"
    } else {
      # We want to examine if individually ploting the smooths for each region is better than a single smooth across all regions, so we use the by argument to fit separate smooths for each region.
      # This way r2 < 0 or small r2 would represent no apparent benefit allowing us to disentangle in which variable regional trajectories may be of relevance. 
      fixed_effects = "pow ~ region + s(age_months, k=4, by = region) + prop_epochs + Cohort + GestationalAge_weeks"
      fixed_effects_red = "pow ~ region + s(age_months, k=4) + prop_epochs + Cohort + GestationalAge_weeks"
    }

    if ("model_fit" %in% names(model_data)) {
      fixed_effects     = paste(fixed_effects,     "+ model_fit")
      fixed_effects_red = paste(fixed_effects_red, "+ model_fit")
    }

    m = gamm(
      as.formula(fixed_effects),
      random  = list(sujid =~ 1),
      correlation = nlme::corCAR1(form = ~ age_months | sujid / region),
      data    = model_data, method = "REML",
      control = nlme::lmeControl(opt = "optim", niterEM = 100, maxIter = 10000, msMaxIter = 10000)
    )

    m_red = gamm(
      as.formula(fixed_effects_red),
      random  = list(sujid =~ 1),
      correlation = nlme::corCAR1(form = ~ age_months | sujid / region),
      data    = model_data, method = "REML",
      control = nlme::lmeControl(opt = "optim", niterEM = 100, maxIter = 10000, msMaxIter = 10000)
    )

    partial_r2  = (sum(residuals(m_red$lme)^2) - sum(residuals(m$lme)^2)) /
                   sum(residuals(m_red$lme)^2)
    gam_results = summary(m$gam)


    # Extract regional intercept differences (parametric table: region coefficients)
    reg_rows <- as.data.frame(gam_results$p.table) %>%
      filter(grepl("region", rownames(.))) %>%
      mutate(stable = stable, metric = plabel, region_level = rownames(.),
             df = nn_gam_df(m$gam), effect_size = partial_r2)
    results_list_region[[f]] <- reg_rows

    # Extract regional smooth significance (s.table: s(age_months):region rows)
    reg_age_rows <- as.data.frame(gam_results$s.table) %>%
      filter(grepl("s\\(age_months\\)", rownames(.))) %>%
      mutate(stable = stable, metric = plabel, region_smooth = rownames(.))
    results_list_reg_smooth[[f]] = reg_age_rows

    # Generate predictions on a fine age grid for plotting
    # Extract Predictions and Calculate Derivatives for Significance
    pred_grid <- expand.grid(
      age_months = seq(min(model_data$age_months, na.rm = TRUE), max(model_data$age_months, na.rm = TRUE), length.out = 100),
      region = unique(model_data$region),
      prop_epochs = mean(model_data$prop_epochs, na.rm = TRUE),
      Cohort = model_data$Cohort[[1]], # Holds cohort constant
      GestationalAge_weeks = mean(model_data$GestationalAge_weeks, na.rm = TRUE)
    )
    if("model_fit" %in% names(model_data)){
      pred_grid$model_fit <- mean(model_data$model_fit, na.rm = TRUE)
    }
    
    # Predict values ignoring the random intercept (population level)
    preds <- predict(m$gam, newdata = pred_grid, se.fit = TRUE)
    pred_grid$.estimate <- preds$fit
    pred_grid$.lower_ci <- preds$fit - (1.96 * preds$se.fit)
    pred_grid$.upper_ci <- preds$fit + (1.96 * preds$se.fit)
    
    # --- GRATIA DERIVATIVE CALCULATION START ---
    # Compute pointwise derivatives (slopes) of the age smooth to identify developmental transitions
    # CI of derivative not crossing 0 indicates significant change (sig_change = TRUE)
    derivs <- gratia::derivatives(m$gam, term = "s(age_months)", data = pred_grid, partial_match = T)
    
    # Identify where the 95% CI of the derivative does not cross 0
    # and keep only the necessary columns to merge back
    sig_derivs <- derivs%>%
      mutate(sig_change = !(.lower_ci <= 0 & .upper_ci >= 0))%>%
      dplyr::select(age_months, region, sig_change)
    
    # Merge the significance flag back into the main prediction grid
    pred_grid <- pred_grid |>
      left_join(sig_derivs, by = c("age_months", "region"))
    
    # --- GRATIA DERIVATIVE CALCULATION END ---
    
    # Position the significance lines dynamically at the bottom of the plot
    y_min <- min(pred_grid$.lower_ci, na.rm = TRUE)
    y_range <- max(pred_grid$.upper_ci, na.rm = TRUE) - y_min
    
    pred_grid <- pred_grid |>
      group_by(region) |>
      mutate(
        # Stagger the lines by region so they don't overlap
        region_num = as.numeric(as.factor(region)),
        sig_y_pos = y_min - (0.05 * y_range * region_num),
        # If not significant, set to NA to break the geom_line
        sig_y_plot = ifelse(sig_change, sig_y_pos, NA)
      ) |> ungroup()
    
    # Plotting Regional Trajectories
    reg_plot <- ggplot(pred_grid, aes(x = age_months, group = region)) +
      geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = region), alpha = 0.15, color = NA) +
      geom_line(aes(y = .estimate, color = region), linewidth = 1.2) +
      # Add the significance line segments at the bottom
      geom_line(aes(y = sig_y_plot, color = region), linewidth = 2.5, na.rm = TRUE) +
      labs(x = "Age (months)", 
           y = plabel,
           color = "Region", fill = "Region") + 
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") + 
      THEME_BASE + THEME_TEXT
    
    
    region_plots[[f]] <- reg_plot 
    

  }
}

# --- Table 3: Regional Intercept Differences ---
if (length(results_list_region) > 0) {
  results_list_region = bind_rows(results_list_region)
  # --- Numbered supplementary table SR7 store in an R document ---
  results_list_region |>
    mutate(p_fdr = p.adjust(`Pr(>|t|)`, method = "fdr"), .by = c(stable, metric)) |>
    mutate(region_level = case_when(
      grepl('Central',   region_level) ~ "Central (reference)",
      grepl('Parietal',  region_level) ~ "Parietal",
      grepl('Occipital', region_level) ~ "Occipital",
      grepl('Temporal',  region_level) ~ "Temporal",
      grepl('Frontal',   region_level) ~ "Frontal",
      TRUE ~ region_level)) |>
    arrange(stable, metric, region_level) |>
    nn_supp_table(id = "SR7", cols = c(
      stable        = "Model family",
      metric        = "Metric",
      region_level  = "Region (vs central)",
      Estimate      = "Beta",
      `Std. Error`  = "s.e.",
      `t value`     = "t",
      df            = "df",
      `Pr(>|t|)`    = "P",
      p_fdr         = "P (FDR)",
      effect_size   = "Partial R2"),
      note = paste(
        "Parametric region contrasts from the ROI GAMMs (Metric ~ s(Age, by = ROI)",
        "+ ROI + Cohort + covariates). Each contrast is a t test against the central",
        "ROI on the model residual degrees of freedom. P values FDR-corrected",
        "(Benjamini-Hochberg) within metric across regions."))

  for (t in unique(results_list_region$stable)) {
    results_list_region |>
      filter(stable == t) |>
      mutate(p_fdr = p.adjust(`Pr(>|t|)`, method = "fdr"), .by = metric) |>
      dplyr::select(metric, Estimate, region_level, `Std. Error`, `t value`, `Pr(>|t|)`, p_fdr, effect_size) |>
      # Standardize region coefficient names: extract region labels and mark reference group
      mutate(region_level = case_when(grepl('Central', region_level) ~ "Central (Reference)",
                                      grepl('Parietal', region_level) ~ "Parietal",
                                      grepl('Occipital', region_level) ~ "Occipital",
                                      grepl('Temporal', region_level) ~ "Temporal",
                                      grepl('Frontal', region_level) ~ "Frontal",
                                      TRUE ~ region_level)) |>
    flextable() |>
    set_header_labels(
      metric       = "Metric",  region_level = "Region (vs. Central)",
      Estimate     = "\u03B2",  `Std. Error` = "SE",
      `t value`    = "t",       `Pr(>|t|)`   = "p",
      effect_size   = "Partial R\u00b2", p_fdr = "p (FDR)"
    ) |>
    set_caption("GAMM Regional Intercept Differences (Central Region as Reference)") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("Pr(>|t|)", "p_fdr"), digits = 4) |>
    bold(j = "p_fdr", i = ~ p_fdr < 0.05, part = "body") |>
    add_footer_lines("FDR correction (Benjamini-Hochberg) applied within each metric across regions.") |>
    theme_booktabs() |> autofit() |>
      merge_v(j = "metric") |>
      set_table_properties(layout = "fixed") |>
      valign(j = "metric", valign = "top") |>
      align(align = "left", part = "all") |>
      align(j = 3:8, align = "center", part = "all") |>
      bold(part = "header")|>
      padding(padding = 1, part = "all") |>
      flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    save_as_html(path = file.path(path2suppres, paste0("Supplementary_Results_FigR1_Table_ROIs_Differences_", t, ".html")))
  }
}

# --- Table 4: Regional Age Smooth Summary ---
if (length(results_list_reg_smooth) > 0) {
  results = bind_rows(results_list_reg_smooth) 

  # --- Numbered supplementary table SR8 ---
  results |>
    mutate(p_fdr = p.adjust(`p-value`, method = "fdr"), .by = c(stable, metric)) |>
    mutate(region_smooth = case_when(
      grepl('Central',   region_smooth) ~ "Central",
      grepl('Parietal',  region_smooth) ~ "Parietal",
      grepl('Occipital', region_smooth) ~ "Occipital",
      grepl('Temporal',  region_smooth) ~ "Temporal",
      grepl('Frontal',   region_smooth) ~ "Frontal",
      TRUE ~ region_smooth)) |>
    arrange(stable, metric, region_smooth) |>
    nn_supp_table(id = "SR8", cols = c(
      stable        = "Model family",
      metric        = "Metric",
      region_smooth = "Region",
      edf           = "EDF",
      Ref.df        = "Ref. df",
      `F`           = "F",
      `p-value`     = "P",
      p_fdr         = "P (FDR)"),
      note = paste(
        "Region-specific age smooths from the ROI GAMMs (Metric ~ s(Age, by = ROI)",
        "+ ROI + Cohort + covariates). Each smooth is tested with an F statistic on",
        "the reference degrees of freedom. P values FDR-corrected",
        "(Benjamini-Hochberg) within metric across regions.",
        "EDF, effective degrees of freedom (>1 indicates a non-linear trajectory)."))

  
  for (t in unique(results_list_region$stable)) {
  results|>
      filter(stable == t)|>
      mutate(p_fdr = p.adjust(`p-value`, method = "fdr"), .by = metric) |>
      dplyr::select(metric, region_smooth, edf, F, `p-value`, p_fdr) |>
      # Standardize region smooth term names for display (extract region labels from s(age_months):region terms)
      mutate(region_smooth = case_when(grepl('Central', region_smooth) ~ "Central",
                                       grepl('Parietal', region_smooth) ~ "Parietal",
                                       grepl('Occipital', region_smooth) ~ "Occipital",
                                       grepl('Temporal', region_smooth) ~ "Temporal",
                                       grepl('Frontal', region_smooth) ~ "Frontal",
                                       TRUE ~ region_smooth)) |>
    flextable() |>
    set_header_labels(
      metric        = "Metric",    region_smooth = "Region Smooth",
      edf           = "EDF",
      F             = "F",         `p-value`     = "p",
      p_fdr         = "p (FDR)"
    ) |>
    set_caption("GAMM Regional Age Smooth Summary") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("p-value", "p_fdr"), digits = 4) |>
    bold(j = "p_fdr", i = ~ p_fdr < 0.05, part = "body") |>
    add_footer_lines("EDF > 1 indicates a non-linear developmental trajectory. FDR correction applied within each metric.") |>
    theme_booktabs() |> autofit() |>
    merge_v(j = "metric") |>
    set_table_properties(layout = "fixed") |>
    valign(j = "metric", valign = "top") |>
    align(align = "left", part = "all") |>
    align(j = 3:6, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    save_as_html(path = file.path(path2suppres, paste0("Supplementary_Results_FigR1_Table_ROIs_Trajectories_", t, ".html")))
  }
}

# --- Regional Trajectory Figure ---
if (length(region_plots) > 0) {
  regional_fig = ggpubr::ggarrange(
    plotlist = region_plots,
    ncol = 3, nrow = ceiling(length(region_plots) / 3),
    common.legend = TRUE,
    legend = "bottom",
    labels = letters[seq_along(region_plots)]
  )
  ggsave(filename = file.path(path2figs, "FigSR1_GAMM_ROIs_Development_all_metrics.jpeg"),
         plot = regional_fig, width = 10, height = 13, units = "in", dpi = 300)
}
