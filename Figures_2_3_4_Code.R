# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - code identifier -> Table 1 terms
#
#   DATA columns -> DISPLAY LABELS (used in plot y-axis/legend labels):
#     offset                   -> "Offset"           (Aperiodic offset)
#     slope                    -> "Slope"            (Aperiodic exponent)
#     peak_freq                -> "Peak Freq."       (Peak Frequency, Hz)
#     peak_ampl                -> "Peak Amp."        (Peak Amplitude)
#     peak_prop                -> "Prop. of Peaks"   (Prop. of electrodes w/ oscillatory peak)
#     osc_ampl                 -> "Band Power"       (Aperiodic-corrected alpha band power)
#     volt_amp / *_corrected   -> "Volt. Amp."       (Voltage Amplitude, absolute/corrected)
#     band_amp / *_corrected   -> "Band Amp."        (Band Amplitude, absolute/corrected)
#     prop_bursty_epochs       -> "Prop. of Epochs w/ Burst"
#     prop_bursty_cycles_burst -> "Prop. of Cycles w/ Burst"
#     avg_burst_duration       -> "Burst Duration"   (Consecutive cycles per burst)
#     alpha_LAcH               -> "Lifespan"         (Alpha lifespan in cycles; LAcH cumsum >=90%)
#     is_burst / burst_type    -> "Cycle Type"       ("Burst" vs. "NoBurst")
#     prop_epochs              -> "Clean Epochs"     (EEG-quality covariate, z-scored)
#     mae (raw specparam col)  -> reported as MAE (mean absolute error)
# -----------------------------------------------------------------------------
# Script: Figures_2_3_4_Code.R
# Purpose: Generates part of main-text figures and some of the supplements. This code needs to be run after analysis.
# =============================================================================
# Inputs:
#   - Data/Aperiodic_Oscillatory_ByCycle_Long.csv
#   - Data/Aperiodic_Oscillatory_ByCycle_Long_Burst.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long_Burst.csv
#   - Data/BurstProperties_ByCycle_Long.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Data/WholeBrain_MomentsOfChange_*.csv  (from script 11)
#   - Data/WholeBrain_Interaction_Moments_*.csv (from script 11)
#   - Data/MLM_WholeBrain_Summary_Emmeans_PerSession_BurstDiffs.csv
#   - Data/lcohhilb_tileplot_region*.csv
#   - Data/lcohhilb_cumplot*.csv
#   - Data/psds_aperosc_long_region*.csv
#   - Data/electrodes.csv
# Outputs saved to path2figs:
#   - Fig_2_ParametrizedPSD_Development_WholeBrain_FullGAMMModeled.jpeg   (Parametrized Power-Spectrum GAMM curves and aperiodic/oscillatory power-spectrum)
#   - Fig_3_AlphaBurst_and_Lifespan_Development_WholeBrain_FullGAMMModeled.jpeg (Burst and Lifespan GAMM curves and cumulative plot/heatmap of lagged coherence)
#   - Fig_4_BurstVsNonBurst_AperOsc_WholeBrain.jpeg (burst vs. non-burst segments barplots and power-spectrum/alpha cumulative plot/heatmap)               
#   - Extended_Data_Fig_3_BurstDuration_Development_WholeBrain_FullGAMMModeled.jpeg    (GAMM curve of burst duration development)
#   - Extended_Data_Fig_4_BurstImpact_AlphaLifespan_by_Visit.jpeg. (burst vs. non-burst segments alpha lifespan barplot)
# Dependencies: 00_Setup_PackageInstallation.R, 01_Utils_ProcFunctions.R
# =============================================================================

# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
# somewhere else, set CODE_FOLDER on the next line to this script's folder.
CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  If you have opened the code from the project, you don't need to modify this line. Otherwise, select where the code folder that contains the config_paths.R is

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
source(file.path(path2code, "01_Utils_ProcFunctions.R"))

# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================

r2_thresh        = R2_THRESH         # Minimum FOOOF R-squared (0.900)
mae_thresh       = MAE_THRESH        # Maximum FOOOF MAE (0.10)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)

# Colour scheme for derivative/age-effect significance bars beneath GAMM curves.
# moments_colors_main: main age-effect bars (Increase/Decrease/No Change mapped to COLORS_GRADIENT)
# moments_colors_burst: burst-interaction bars (Diverging modes/Parallel mapped separately)

moments_colors_main  = c("Increase" = COLORS_GRADIENT[2], "Decrease" = COLORS_GRADIENT[7], "No Change" = "gray90")
moments_colors_burst = c("Diverging (-)" = COLORS_GRADIENT[2], "Diverging (+)" = COLORS_GRADIENT[7], "Parallel" = "gray90")


# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================

path2data = path2sets   # EDIT — root data folder   # = Data/; this script prepends "" to each file name
# path2root comes from config_paths.R

# This script emits BOTH main-text (Fig 2, 3, 4) and supplementary (Extended Data Figs. 3 and 4)
path2figs     = file.path(path2root, "MainText", "Figures")       # Fig 2, 3, 4
path2suppfig  = file.path(path2root, "ExtendedData", "Figures")  # Extended Data Fig. 3, Extended Data Fig. 4

for (p in c(path2figs, path2suppfig)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- List of moments-of-change CSV files (written by script 11) ---
change_files = list.files(
  path       = file.path(path2data, ""),
  pattern    = "WholeBrain.*\\.csv$",
  full.names = TRUE
)

change_files = change_files[!grepl('djEpoch', change_files)]  # This figure is generated in a different script

# --- EEG cleaning descriptives ---
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv"),
                     show_col_types = FALSE) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic data ---
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv"),
                          show_col_types = FALSE) |>
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
electrodes = read_csv(file.path(path2data, "electrodes.csv"),
                       show_col_types = FALSE) |>
  dplyr::select(label, region, hemis, chinclu) |>
  # Recode electrode region codes to more legible labels
  # region: Fr = Frontal, P = Parietal, T = Temporal, O = Occipital, C = Central
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1)

# Combined descriptives frame (left anchor for merges)
descriptives = left_join(
    desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
    eeg_desc, by = c("sujid", "session_age")
  ) |>
  left_join(desc_and_ages_wide, by = "sujid") |>
  filter(!is.na(prop_epochs))

# --- Aperiodic/oscillatory data (whole-brain averaged, not burst-conditioned) ---
aper_voi = c("sujid", "session_age", "mae", "ch", "region", "chinclu", "epochs",
              "r2value", "goodch", "offset", "slope",
              "alpha_freq", "alpha_ampl", "alpha_osc", "alpha_peak", "inclusion_final_dummy") # Variables of interest of the parametrized power-spectrum

aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv"),
                           show_col_types = FALSE) |>
  dplyr::select(all_of(aper_voi)) |>
  group_by(session_age, sujid) |>
  # Quality filters: FOOOF R², MAE, electrode inclusion, epoch count, channel quality
  filter(r2value >= r2_thresh, mae <= mae_thresh, chinclu == 1,
         epochs >= epochs_threshold, goodch > CH_THRESHOLD,
         !(sujid == "SUB-XCM45B" & session_age == 48)) |>  # Exclude a participant who crashed during processing 
  # Collapse creating the single whole-brain metric
  summarise(
    slope     = mean(slope,                       na.rm = TRUE),
    offset    = mean(offset,                      na.rm = TRUE),
    peak_ampl = mean(alpha_ampl,                  na.rm = TRUE),  # Peak amplitude was set to zero during processing, so it is safe to average
    peak_freq = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),  # Only when alpha_peak==1
    osc_ampl  = mean(alpha_osc,                   na.rm = TRUE),
    peak_prop  = mean(alpha_peak,                  na.rm = TRUE),  # Proportion with peak
    epochs    = mean(epochs,                      na.rm = TRUE),
    r2value   = mean(r2value,                     na.rm = TRUE),
    .groups   = "drop"
  )

# --- Lagged coherence (alpha lifespan) data ---
# alpha_LAcH: estimated cycle at which cumulative LAcH ≥ 90%; aggregated to subject-age
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv"),
                       show_col_types = FALSE) |>
  merge(electrodes |> dplyr::select(label, chinclu), by.x = "ch", by.y = "label") |>
  group_by(session_age, sujid) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop")

# --- Burst properties (amplitude, duration) stratified by cycle type (is_burst) ---
# Input has one row per cycle per electrode; aggregated to subject-age-burst_type level
burst_data = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv"),
                       show_col_types = FALSE) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  group_by(session_age, sujid, is_burst) |>  # Keep is_burst for comparison
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Merge with descriptives and apply final inclusion filter
aperiodic_data = left_join(aperiodic_data, descriptives, by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
hlcoh_data     = left_join(hlcoh_data,     descriptives, by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
burst_data     = left_join(burst_data,     descriptives, by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)

# =============================================================================
# SECTION 4: PSD DATA
# =============================================================================
# Region-averaged PSDs aggregated across subjects to one row per session_age × freq.
# IQR ribbons (Q1–Q3) computed to estimate the variability of the signal. 

psd_data_files = list.files(path2data, pattern = "psds_aperosc_long_region.*\\.csv$",
                              full.names = TRUE)
psd_data_files = psd_data_files[!grepl("burst", psd_data_files)]

psd_data = lapply(psd_data_files, read_csv, show_col_types = FALSE) |> bind_rows() |>
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  filter(inclusion_final_dummy == TRUE, sujid %in% descriptives$sujid) |>
  group_by(session_age, freq, sujid) |>
  summarise(aperiodic = mean(aperiodic, na.rm = TRUE),
            oscillatory = mean(oscillatory, na.rm = TRUE), .groups = "drop") |>
  group_by(session_age, freq) |>
  summarise(
    aper_mean  = mean(aperiodic,   na.rm = TRUE),
    aper_first = quantile(aperiodic,   0.25, na.rm = TRUE),
    aper_third = quantile(aperiodic,   0.75, na.rm = TRUE),
    osc_mean   = mean(oscillatory, na.rm = TRUE),
    osc_first  = quantile(oscillatory, 0.25, na.rm = TRUE),
    osc_third  = quantile(oscillatory, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )


# =============================================================================
# SECTION 4B: FULL-GAMM TRAJECTORY
# =============================================================================
# For comparability with the analysis, in which we included covariates, GAMM curves are drawn from the full model (DataAnalysis_1).

# Prediction grid resolution
GRID_N        = 200
CI_MULT       = 1.96
OVERLAY_SHIFT = TRUE

# --- Aperiodic/oscillatory data REGION-STRATIFIED (for GAMM with region random effects) ---
aper_voi_reg = c("sujid", "session_age", "mae", "ch", "region", "chinclu", "epochs",
                 "r2value", "goodch", "offset", "slope",
                 "alpha_freq", "alpha_ampl", "alpha_osc", "alpha_peak", "inclusion_final_dummy")

aperiodic_reg = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv"),
                          show_col_types = FALSE) |>
  dplyr::select(all_of(aper_voi_reg)) |>
  group_by(session_age, sujid, region) |>
  filter(r2value > r2_thresh, mae < mae_thresh, chinclu == 1,
         epochs >= epochs_threshold, goodch > CH_THRESHOLD) |>
  # Collapse electrodes across regions but keep region for model stratification
  summarise(
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
  )

# --- Lagged coherence REGION-STRATIFIED (for GAMM) ---
hlcoh_reg = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv"),
                     show_col_types = FALSE) |>
  merge(electrodes |> dplyr::select(label, chinclu, region), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  group_by(session_age, sujid, region) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop")

# --- Burst properties REGION and CYCLE-TYPE STRATIFIED (for GAMM interaction model) ---
# Retains both is_burst and region for model: allows region×burst_type variation
burst_reg = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long.csv"),
                     show_col_types = FALSE) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1) |>
  group_by(session_age, sujid, is_burst, region) |>  # Keep both is_burst and region
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Merge with descriptives and apply the same inclusion filter as the models.
aperiodic_reg = left_join(descriptives, aperiodic_reg, by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
hlcoh_reg     = left_join(descriptives, hlcoh_reg,     by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)
burst_reg     = left_join(descriptives, burst_reg,     by = c("session_age", "sujid")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)

# GAMM trajectory function: fits full model per metric type and predicts developmental curve
# f: variable name (column in model_data); metric_type: "Aperiodic"/"Oscillatory"/"Lagged Coherence"/"Burst"

gamm_full_trajectory = function(f, metric_type, align_to = NULL) {

  # Select appropriate region-stratified dataset and rename target column to 'pow' 
  if (metric_type %in% c("Aperiodic", "Oscillatory")) {
    # Aperiodic/oscillatory: use aperiodic_reg, add scaled FOOOF R² as covariate
    model_data = aperiodic_reg |> rename(pow = !!sym(f)) |>
      mutate(model_fit = scale(r2value)[, 1]) |> tidyr::drop_na(pow)
  } else if (metric_type == "Lagged Coherence") {
    # Lagged coherence: use hlcoh_reg (alpha_LAcH)
    model_data = hlcoh_reg |> rename(pow = !!sym(f)) |> tidyr::drop_na(pow)
  } else if (metric_type == "Burst") {
    # Burst metrics: use burst_reg, filter to only Burst cycles (exclude Non-Burst)
    model_data = burst_reg |> rename(pow = !!sym(f)) |>
      filter(is_burst == "Burst") |>
      mutate(pow = if_else(is.na(pow), 0, pow)) |> tidyr::drop_na(pow)
  } else return(NULL)

  # Prepare model data: scale numeric covariates, convert region/Cohort to factors
  numeric_covars = c("age_months", "prop_epochs", "GestationalAge_weeks")
  if ("model_fit" %in% names(model_data)) numeric_covars = c(numeric_covars, "model_fit")
  model_data = model_data |> tidyr::drop_na(any_of(numeric_covars)) |>
    mutate(region = factor(region), Cohort = factor(Cohort))

  # Minimum sample size check for stable model fit
  n_unique_ids = length(unique(model_data$sujid))
  if (!(n_unique_ids > 1 && nrow(model_data) > 2 * n_unique_ids)) return(NULL)

  # ---- Model formula (matches DataAnalysis_1 scripts for reproducibility) ----
  # Fixed: smooth age term + region + covariates; Random: subject intercept; Correlation: CAR(1)
  if (f != "offset") {
  fixed_effects = "pow ~ s(age_months, k = 4) + region + prop_epochs + Cohort + GestationalAge_weeks"
  } else { 
    # Offset model had worse fit when the spline was included, signalling a worse fit; use k=3.
    fixed_effects = "pow ~ s(age_months, k = 3) + region + prop_epochs + Cohort + GestationalAge_weeks"
    }
  if ("model_fit" %in% names(model_data)) fixed_effects = paste(fixed_effects, "+ model_fit")

  # Fit GAMM with error handling
  m = tryCatch(
    gamm(as.formula(fixed_effects),
         random  = list(sujid =~ 1),
         correlation = nlme::corCAR1(form = ~ age_months | sujid / region),
         data    = model_data, method = "REML",
         control = nlme::lmeControl(opt = "optim", niterEM = 100, maxIter = 10000, msMaxIter = 10000)),
    error = function(e) { message("  ! GAMM fit failed (", f, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)

  # ---- Prediction grid: reference region, mean-centered covariates (z-score = 0) ----
  # This isolates the smooth age term; regional/cohort intercepts are factored out
  agerng = range(model_data$age_months, na.rm = TRUE)
  grid = tibble(
    age_months           = seq(agerng[1], agerng[2], length.out = GRID_N),
    region               = factor(levels(model_data$region)[1], levels = levels(model_data$region)),
    prop_epochs           = 0,  # Mean value (scaled)
    Cohort               = factor(levels(model_data$Cohort)[1], levels = levels(model_data$Cohort)),
    GestationalAge_weeks = 0   # Mean value (scaled)
  )
  if ("model_fit" %in% names(model_data)) grid$model_fit = 0

  # Predict with 95% CI
  pr = predict(m$gam, newdata = grid, se.fit = TRUE)
  grid$fit = as.numeric(pr$fit); grid$se = as.numeric(pr$se.fit)

  # Display-only vertical shift to overlay curve on raw data (preserves shape and CI width)
  if (OVERLAY_SHIFT) {
    ref = if (!is.null(align_to)) mean(align_to, na.rm = TRUE) else mean(model_data$pow, na.rm = TRUE)
    grid$fit = grid$fit + (ref - mean(grid$fit, na.rm = TRUE))
  }
  # Confidence intervals
  grid$lo = grid$fit - CI_MULT * grid$se
  grid$hi = grid$fit + CI_MULT * grid$se
  grid
}

# ---- GAMM with Burst×Age interaction (DataAnalysis_1 Part B) ----
# Fits separate age smooths for Burst vs. Non-Burst, returns both trajectories
gamm_full_interaction = function(f, align_to = NULL) {

  # Prepare data: rename metric to 'pow', convert is_burst to ordered factor (treatment contrast)
  model_data = burst_reg |> rename(pow = !!sym(f)) |>
    mutate(is_burst = factor(is_burst))
  model_data$is_burst_ord = as.ordered(model_data$is_burst)
  contrasts(model_data$is_burst_ord) = "contr.treatment"

  # Drop NA in covariates and convert to factors
  model_data = model_data |>
    tidyr::drop_na(any_of(c("age_months", "prop_epochs", "GestationalAge_weeks"))) |>
    mutate(region = factor(region), Cohort = factor(Cohort))

  # Fit interaction model: is_burst_ord main effect + smooth age + smooth age by is_burst_ord
  m = tryCatch(
    gamm(pow ~ is_burst_ord + s(age_months, k=4) + s(age_months, k=4, by = is_burst_ord) +
               region + prop_epochs + Cohort + GestationalAge_weeks,
         random  = list(sujid =~ 1),
         correlation = nlme::corCAR1(form = ~ age_months | sujid / region / is_burst_ord),
         data    = model_data, method = "REML",
         control = nlme::lmeControl(maxIter = 100, msMaxIter = 100)),
    error = function(e) { message("  ! GAMM interaction fit failed (", f, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)

  # ---- Prediction: separate trajectories for Burst and Non-Burst across age ----
  agerng = range(model_data$age_months, na.rm = TRUE)
  ages   = seq(agerng[1], agerng[2], length.out = GRID_N)
  levs   = levels(model_data$is_burst_ord)

  # Map over is_burst_ord levels, create grid row per level-age combination
  grid = purrr::map_dfr(levs, function(lv) {
    tibble(
      age_months           = ages,
      is_burst_ord         = ordered(lv, levels = levs),
      is_burst             = factor(lv, levels = levels(model_data$is_burst)),  # For colour mapping in plot
      region               = factor(levels(model_data$region)[1], levels = levels(model_data$region)),
      prop_epochs           = 0,  # Mean (scaled)
      Cohort               = factor(levels(model_data$Cohort)[1], levels = levels(model_data$Cohort)),
      GestationalAge_weeks = 0   # Mean (scaled)
    )
  })
  # Predict both curves
  pr = predict(m$gam, newdata = grid, se.fit = TRUE)
  grid$fit = as.numeric(pr$fit); grid$se = as.numeric(pr$se.fit)

  # Display-only shift (COMMON for both curves to preserve gap)
  if (OVERLAY_SHIFT) {
    ref = if (!is.null(align_to)) mean(align_to, na.rm = TRUE) else mean(model_data$pow, na.rm = TRUE)
    grid$fit = grid$fit + (ref - mean(grid$fit, na.rm = TRUE))
  }
  grid$lo = grid$fit - CI_MULT * grid$se
  grid$hi = grid$fit + CI_MULT * grid$se
  grid
}


# =============================================================================
# SECTION 5: FIGURE A — DEVELOPMENTAL TRAJECTORY PLOTS WITH SIGNIFICANCE BARS
# =============================================================================
# Generates one panel per metric; bars below each curve indicate periods of
# significant positive/negative change from the simultaneous derivatives (from change_files CSVs).

# Internal column names (from CSVs; never renamed)
variables      = c("offset", "slope", "peak_freq", "peak_ampl", "osc_ampl", "peak_prop",
                   "prop_bursty_epochs", "prop_bursty_cycles_burst", "alpha_LAcH", "avg_burst_duration")
# Data source type for each variable (selects aperiodic_data, burst_data, or hlcoh_data)
type_of_metric = c("Aperiodic", "Aperiodic", "Oscillatory", "Oscillatory", "Oscillatory", "Oscillatory",
                   "Burst", "Burst", "Lagged Coherence", "Burst")

# DISPLAY LABELS: manuscript-friendly names mapped to y-axis in plots
# Must match variable[i] to label[i] correspondence. 
labels         = c("Offset", "Slope", "Peak Freq.", "Peak Amp.", "Band Power", "Prop. of Peaks",
                   "Prop. of Epochs\nw/ Burst", "Prop. of Cycles\nw/ Burst", "Lifespan", "Burst Duration")

plot_list = list()

for (i in seq_along(variables)) {
  var_name    = variables[i]
  metric_type = type_of_metric[i]
  y_label     = labels[i]

  # Select the correct source dataset
  if (metric_type %in% c("Aperiodic", "Oscillatory")) {
    plot_data = aperiodic_data
    if (var_name == "peak_freq") plot_data = plot_data |> filter(peak_freq > 0)
  } else if (metric_type == "Burst") {
    plot_data = burst_data
  } else if (metric_type == "Lagged Coherence") {
    plot_data = hlcoh_data
  } else next

  if (!var_name %in% names(plot_data)) next

  # Match the relevant moments-of-change CSV (exclude burst-conditioned files for non-burst metrics)
  file_matches = change_files[grepl(var_name, change_files)]

  if (length(file_matches) == 0) next

  change_subfile_data = read_csv(file_matches[1], show_col_types = FALSE)

  # Position significance bar just below the lowest data point
  y_bar_pos = min(plot_data[[var_name]], na.rm = TRUE) -
              0.15 * sd(plot_data[[var_name]], na.rm = TRUE)

  bar_data = change_subfile_data |>
    filter(!is.na(age_months), !is.na(inc_or_dec)) |>
    arrange(age_months) |>
    mutate(status_change = inc_or_dec != lag(inc_or_dec, default = first(inc_or_dec)),
           segment_id    = cumsum(status_change),
           y_location    = y_bar_pos)

  # Helper function to plot the GAMM trajectory above
  gamm_fit = gamm_full_trajectory(var_name, metric_type,
                                  align_to = plot_data[[var_name]])

  p = ggplot(plot_data, aes(x = age_months, y = .data[[var_name]])) +
    geom_point(alpha = 0.2, color = "grey60", size = 1,
               position = position_jitter(width = 0.2, height = 0)) +
    geom_line(data = bar_data, inherit.aes = FALSE,
              aes(x = age_months, y = y_location, color = inc_or_dec, group = segment_id),
              linewidth = 5, alpha = 0.5, lineend = "butt") +
    scale_color_manual(values = COLORS_NATURE) +
    THEME_BASE + THEME_TEXT +
    labs(x = "Age (months)", y = y_label, color = "Age-related change")

  if (!is.null(gamm_fit)) {
    p = p +
      geom_ribbon(data = gamm_fit, inherit.aes = FALSE,
                  aes(x = age_months, ymin = lo, ymax = hi),
                  fill = colorspace::lighten(COLORS_MAIN[1], 0.3), alpha = 0.4) +
      geom_line(data = gamm_fit, inherit.aes = FALSE,
                aes(x = age_months, y = fit),
                color = COLORS_MAIN[1], linewidth = 1.2)
  }

  plot_list[[var_name]] = p
}

ggsave(plot_list$avg_burst_duration,
       filename = file.path(path2suppfig, "Extended_Data_Fig_3_BurstDuration_Development_WholeBrain_FullGAMMModeled.jpeg"),
       width = 5, height = 4.5, dpi = 300)


# =============================================================================
# SECTION 6: BURST × AGE INTERACTION PANELS (Volt/Band Amp)
# =============================================================================

voi_burst    = c("volt_amp", "band_amp")
labels_burst = c("Volt. Amp.", "Band Amp.")
plot_list_isburst = list()

for (i in seq_along(voi_burst)) {
  var          = voi_burst[i]
  file_matches = change_files[grepl(var, change_files)]
  if (length(file_matches) == 0) next

  change_inter = read_csv(file_matches[1], show_col_types = FALSE)
  y_bar_pos    = min(burst_data[[var]], na.rm = TRUE) -
                 0.15 * sd(burst_data[[var]], na.rm = TRUE)

  bar_data = change_inter |>
    filter(!is.na(age_months), !is.na(inc_or_dec)) |>
    arrange(age_months) |>
    mutate(status_change = inc_or_dec != lag(inc_or_dec, default = first(inc_or_dec)),
           segment_id    = cumsum(status_change),
           y_location    = y_bar_pos)

  gamm_fit = gamm_full_interaction(var, align_to = burst_data[[var]])

  p = ggplot(burst_data, aes(x = age_months, y = .data[[var]])) +
    geom_point(aes(color = is_burst), alpha = 0.2, size = 1,
               position = position_jitter(width = 0.2, height = 0))

  if (!is.null(gamm_fit)) {
    p = p +
      geom_ribbon(data = gamm_fit, inherit.aes = FALSE,
                  aes(x = age_months, ymin = lo, ymax = hi, fill = is_burst),
                  alpha = 0.33) +
      geom_line(data = gamm_fit, inherit.aes = FALSE,
                aes(x = age_months, y = fit, color = is_burst), linewidth = 1.2)
  }

  p = p +
    scale_color_manual(values = COLORS_MAIN, name = "Cycle Type",
                       guide = guide_legend(override.aes = list(alpha = 1, size = 3,
                                                                 linewidth = 2, fill = NA))) +
    scale_fill_manual(values = COLORS_MAIN, name = "Cycle Type", guide = "none") +
    ggnewscale::new_scale_color() +
    geom_line(data = bar_data, inherit.aes = FALSE,
              aes(x = age_months, y = y_location, color = inc_or_dec, group = segment_id),
              linewidth = 5, alpha = 0.5, lineend = "butt") +
    scale_color_manual(values = COLORS_NATURE_BURST, name = "Age-related change") +
    THEME_BASE + THEME_TEXT +
    labs(x = "Age (months)", y = labels_burst[i])

  if (i == 2) p = p + theme(legend.position = "none")
  plot_list_isburst[[var]] = p
}


# =============================================================================
# SECTION 7: PSD FIGURE PANELS
# =============================================================================

plot_aper_psd = ggplot(psd_data, aes(x = freq, y = aper_mean, color = factor(session_age))) +
  geom_ribbon(aes(ymin = aper_first, ymax = aper_third, fill = factor(session_age)),
              alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = COLORS_GRADIENT) +
  scale_fill_manual(values  = COLORS_GRADIENT) +
  labs(x = "Frequency (Hz)", y = "log10(mV\u00b2/Hz)", color = "Age (months)", fill = "Age (months)") +
  THEME_BASE + THEME_TEXT

plot_osc_psd = ggplot(psd_data, aes(x = freq, y = osc_mean, color = factor(session_age))) +
  geom_ribbon(aes(ymin = osc_first, ymax = osc_third, fill = factor(session_age)),
              alpha = 0.20, color = NA) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = COLORS_GRADIENT) +
  scale_fill_manual(values  = COLORS_GRADIENT) +
  labs(x = "Frequency (Hz)", y = "log10(mV\u00b2/Hz)", color = "Age (months)", fill = "Age (months)") +
  THEME_BASE + THEME_TEXT

plot_psd_combined = ggpubr::ggarrange(plot_aper_psd, plot_osc_psd,
                                       ncol = 1, nrow = 2, labels = c('a', 'b'),
                                       common.legend = TRUE, legend = "bottom")

# Aperiodic/oscillatory trajectory panel array (panels C–H)
plot_list_aperosc = plot_list[c("offset", "slope", "peak_freq", "peak_ampl", "osc_ampl", "peak_prop")]
plot_trajectories = ggpubr::ggarrange(plotlist = plot_list_aperosc, ncol = 2, nrow = 3,
                                       common.legend = TRUE, legend = "bottom",
                                       labels = c('c', 'd', 'e', 'f', 'g', 'h'))

ggsave(
  filename = file.path(path2figs, "Fig2_ParametrizedPSD_Development_WholeBrain_FullGAMMModeled.jpeg"),
  plot = ggpubr::ggarrange(plot_psd_combined, plot_trajectories, ncol = 2, widths = c(1, 1.5)),
  width = 30, height = 18, dpi = 300, units = "cm"
)


# =============================================================================
# SECTION 8: BURST & LCoH COMBINED FIGURE (PANELS A–G)
# =============================================================================

# --- Panel E: LCoH tile heatmap (cycle × frequency, faceted by age) ---
lcoh_tileplot_files = list.files(path2data, pattern = "lcohhilb_tileplot_region.*\\.csv$",
                                  recursive = TRUE, full.names = TRUE)
lcoh_tileplot_files = lcoh_tileplot_files[!grepl("burst", lcoh_tileplot_files)]

lcoh_tiledata = lapply(lcoh_tileplot_files, read.csv, header = TRUE) |> bind_rows() |>
  filter(inclusion == 1) |> # Precomputed during the dataset generation for the plot
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  group_by(session_age, sujid, freq, cycle) |>
  summarise(LAcH = mean(LAcH, na.rm = TRUE), .groups = "drop")

lcoh_tileplot = ggplot(lcoh_tiledata, aes(x = cycle, y = freq, fill = LAcH)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       guide = guide_colorbar(barwidth = unit(12, "lines"),
                                              barheight = unit(0.5, "lines"))) +
  facet_grid(. ~ session_age) +
  labs(x = "Cycle", y = "Frequency (Hz)", fill = "LAcH") +
  THEME_BASE + THEME_TEXT +
  theme(legend.title = element_text(size = 9, face = "bold", vjust = 0.8),
        legend.text  = element_text(size = 8))

# --- Panel G: Cumulative lifespan curve by age ---
lcoh_lifeplot_files = list.files(path2data, pattern = "lcohhilb_cumplot.*\\.csv$",
                                  recursive = TRUE, full.names = TRUE)
lcoh_lifeplot_files = lcoh_lifeplot_files[grepl("region", lcoh_lifeplot_files) &
                                           !grepl("burst", lcoh_lifeplot_files)]

lcoh_lifedata = lapply(lcoh_lifeplot_files, read.csv, header = TRUE) |> bind_rows() |>
  filter(inclusion == 1) |> # Precomputed during the dataset generation for the plot
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  group_by(sujid, session_age, cycle) |>
  summarise(LAcH = mean(cumLAcH, na.rm = TRUE), .groups = "drop") |>
  group_by(session_age, cycle) |>
  summarise(m_lcoh    = mean(LAcH, na.rm = TRUE),
            lcoh_first = quantile(LAcH, 0.25, na.rm = TRUE),
            lcoh_third = quantile(LAcH, 0.75, na.rm = TRUE), .groups = "drop")

lcoh_lifeplot = ggplot(lcoh_lifedata,
                        aes(x = cycle, y = m_lcoh, color = factor(session_age),
                            fill = after_scale(color))) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lcoh_first, ymax = lcoh_third), alpha = 0.20, linewidth = 0) +
  scale_color_manual(values = COLORS_GRADIENT) +
  THEME_BASE + THEME_TEXT +
  labs(x = "Cycle", y = "Cumulative LAcH",
       color = "Age (months)", fill = "Age (months)") +
  theme(legend.position  = c(0.80, 0.35), legend.direction = "vertical",
        legend.key        = element_blank(), legend.key.size = unit(0.7, "lines"),
        legend.title      = element_text(size = 9, face = "bold"),
        legend.text       = element_text(size = 8))

# --- Burst-conditioned LCoH tile and lifespan (difference: burst − no-burst) ---
lcoh_sets_tile  = list.files(path2data, pattern = "tile.*\\.csv$", recursive = TRUE, full.names = TRUE)
lcoh_set_tile_plot = lcoh_sets_tile[grepl("burst", lcoh_sets_tile) & grepl("region", lcoh_sets_tile)]

lcoh_data_tile = lapply(lcoh_set_tile_plot, read.csv, header = TRUE) |>
  bind_rows() |>
  filter(inclusion == 1, sujid %in% descriptives$sujid) |>
  group_by(session_age, burst_type, cycle, freq) |>
  summarise(LAcH = mean(LAcH, na.rm = TRUE), .groups = "drop") |>
  ungroup() |>
  mutate(lcoh_diff = LAcH[burst_type == "burst"] - LAcH[burst_type == "noburst"],
         .by = c(session_age, cycle, freq))

lcoh_tileplot_burst = ggplot(lcoh_data_tile, aes(x = cycle, y = freq, fill = lcoh_diff)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       guide = guide_colorbar(barwidth = unit(12, "lines"),
                                              barheight = unit(0.5, "lines"))) +
  facet_grid(. ~ session_age) +
  labs(x = "Cycle", y = "Frequency (Hz)", fill = "\u0394LAcH") +
  THEME_BASE + THEME_TEXT +
  theme(legend.title = element_text(size = 9, face = "bold", vjust = 0.8),
        legend.text  = element_text(size = 8))

# Burst-conditioned cumulative lifespan (No Burst − Burst)
lcoh_set_alpha_lifecycle = list.files(path2data, pattern = "lcohhilb_cumplot_allch.*\\.csv$",
                                       recursive = TRUE, full.names = TRUE)
lcoh_set_alpha_lifecycle = lcoh_set_alpha_lifecycle[grepl("burst", lcoh_set_alpha_lifecycle)]

lcoh_alphacycle_diff = lapply(lcoh_set_alpha_lifecycle, read.csv, header = TRUE) |>
  bind_rows() |>
  filter(ch %in% electrodes$label, inclusion == 1, sujid %in% descriptives$sujid) |>
  mutate(goodch = n_distinct(ch), .by = c(sujid, burst_type, session_age)) |>
  filter(goodch > CH_THRESHOLD) |>
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  group_by(sujid, session_age, burst_type, cycle, region) |>
  summarise(LAcH = mean(cumLAcH, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(values_from = LAcH, names_from = burst_type) |>
  rowwise() |>
  mutate(cumlcoh_diff = noburst - burst) |>
  group_by(session_age, cycle) |>
  summarise(m_lcoh     = mean(cumlcoh_diff, na.rm = TRUE),
            lcoh_first = quantile(cumlcoh_diff, 0.25, na.rm = TRUE),
            lcoh_third = quantile(cumlcoh_diff, 0.75, na.rm = TRUE), .groups = "drop")

lcoh_lifeplot_diff = ggplot(lcoh_alphacycle_diff,
                             aes(x = cycle, y = m_lcoh, color = factor(session_age),
                                 fill = after_scale(color))) +
  geom_ribbon(aes(ymin = lcoh_first, ymax = lcoh_third), alpha = 0.20, linewidth = 0) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = COLORS_GRADIENT) +
  THEME_BASE + THEME_TEXT +
  labs(x = "Cycle", y = "Cumulative LAcH", color = "Age (months)") +
  theme(legend.position = "none")


# Assemble the burst + LCoH combined figure
plot_burst_combined1 = ggpubr::ggarrange(plotlist = plot_list_isburst, ncol = 2,
                                          labels = c('a', 'b'), common.legend = TRUE,
                                          legend = "bottom")
plot_burst_combined2 = ggpubr::ggarrange(
  plotlist = plot_list[c("prop_bursty_cycles_burst", "prop_bursty_epochs", "alpha_LAcH")],
  ncol = 3, common.legend = TRUE, legend = "bottom", labels = c('c', 'd', 'f')
)
lcoh_plot_FG         = ggpubr::ggarrange(plot_burst_combined2, lcoh_lifeplot,
                                          ncol = 2, widths = c(3, 1), labels = c("", 'g'))
lcoh_combined        = ggpubr::ggarrange(plot_burst_combined1, lcoh_tileplot,
                                          ncol = 2, labels = c("", 'e'), common.legend = FALSE)
lcoh_and_burst_plot  = ggpubr::ggarrange(lcoh_combined, lcoh_plot_FG, ncol = 1, nrow = 2)

ggsave(
  filename = file.path(path2figs, "Fig_3_AlphaBurst_and_Lifespan_Development_WholeBrain_FullGAMMModeled.jpeg"),
  plot     = lcoh_and_burst_plot,
  width = 30, height = 18, dpi = 300, units = "cm"
)


# =============================================================================
# SECTION 9: BURST VS. NON-BURST APERIODIC / OSCILLATORY PANELS (Bar charts)
# =============================================================================
# Bar charts comparing aperiodic/oscillatory metrics between burst and
# non-burst cycle segments. Significance stars overlaid from pre-computed emmeans.

# --- Burst-conditioned aperiodic data ---
aper_voi_b = c("sujid", "session_age", "ch", "region", "chinclu", "r2value", "goodch",
                "offset", "slope", "alpha_freq", "alpha_ampl", "alpha_osc",
                "alpha_peak", "inclusion_final_dummy", "epochs", "burst")

aperiodic_data_b = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long_Burst.csv"),
                              show_col_types = FALSE) |>
  dplyr::select(all_of(aper_voi_b)) |>
  group_by(session_age, sujid, region, burst) |>
  filter(r2value > r2_thresh, chinclu == 1, epochs >= epochs_threshold,
         goodch >= CH_THRESHOLD, sujid != "SUB-XCM45B") |> # This participant crashed when parametrizing the PSD with only bursts/non-bursts 
  group_by(session_age, sujid, burst) |>
  summarise(
    slope     = mean(slope,                       na.rm = TRUE),
    offset    = mean(offset,                      na.rm = TRUE),
    peak_ampl = mean(alpha_ampl,                  na.rm = TRUE),
    peak_freq = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),
    osc_ampl  = mean(alpha_osc,                   na.rm = TRUE),
    peak_prop  = mean(alpha_peak,                  na.rm = TRUE),
    epochs    = mean(epochs,                 na.rm = TRUE),
    r2value   = mean(r2value,                     na.rm = TRUE),
    .groups   = "drop"
  ) |>
  filter(sujid %in% descriptives$sujid) |>
  mutate(session_age = if_else(session_age == 15, 18, session_age))

# --- Burst-conditioned coherence data ---
hlcoh_data_b = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long_Burst.csv"),
                         show_col_types = FALSE) |>
  merge(electrodes |> dplyr::select(label, chinclu), by.x = "ch", by.y = "label") |>
  filter(epochs >= epochs_threshold, chinclu == 1, sujid %in% descriptives$sujid)|>
  mutate(goodch = n_distinct(ch), .by = c(sujid, session_age, burst)) |>
  filter(goodch >= CH_THRESHOLD) |>
  group_by(session_age, sujid, burst) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop") |>
  mutate(session_age = if_else(session_age == 15, 18, session_age))

# Significance overlay data from pre-computed emmeans contrasts
emmeans_results = read.csv(file.path(path2data, "Extended_Results_MLM_Emmeans_Stratified_by_visit.csv"), # This file is generated in DataAnalysis_2 Code
                            header = TRUE) |>
  mutate(stars = Burst_sig) |>
  rename(p.value.fdr = Burst_sig)|>
  #we extract now the number in visit 
  mutate(visit = as.numeric(gsub(" mo.", "", visit))) |>
  dplyr::select(visit, metric, p.value.fdr, stars) |>
  rename(session_age = visit)

# --- Build bar chart panels ---
variables_b      = c("offset", "slope", "peak_freq", "peak_ampl", "osc_ampl",
                      "peak_prop", "alpha_LAcH")
type_of_metric_b = c("Aperiodic", "Aperiodic", "Oscillatory", "Oscillatory",
                      "Oscillatory", "Oscillatory", "Lagged Coherence")
labels_b         = c("Offset", "Slope", "Peak Freq.", "Peak Amp.", "Band Power",
                      "Prop. of Peaks", "Lifespan")

plot_list_b = list()
for (i in seq_along(variables_b)) {
  var_name    = variables_b[i]
  metric_type = type_of_metric_b[i]
  y_label     = labels_b[i]

  if (metric_type %in% c("Aperiodic", "Oscillatory")) {
    plot_data = aperiodic_data_b |>
      mutate(session_age = if_else(session_age == 15, 18, session_age))
    if (var_name == "peak_freq") plot_data = plot_data |> filter(peak_freq > 0)
  } else if (metric_type == "Lagged Coherence") {
    plot_data = hlcoh_data_b |>
      mutate(session_age = if_else(session_age == 15, 18, session_age))
  } else next

  if (!var_name %in% names(plot_data)) next

  p = ggplot(plot_data |> mutate(burst = if_else(burst == "burst", "Burst", "Non-Burst")),
             aes(x = factor(session_age), y = .data[[var_name]], fill = burst)) +
    stat_summary(geom = "bar", fun = "mean", color = "black",
                 position = position_dodge(0.8), alpha = 0.66, linewidth = 1) +
    geom_point(aes(color = burst),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
               alpha = 0.5) +
    stat_summary(geom = "errorbar", fun.data = median_se, color = "black",
                 position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
    stat_summary(geom = "point", fun = "median", color = "black",
                 position = position_dodge(0.8), size = 2) +
    stat_summary(geom = "line", fun = "median", aes(group = burst), color = "black",
                 position = position_dodge(0.8), linewidth = 1)

  # Add significance stars if this variable has emmeans contrasts
  if (any(emmeans_results$metric == y_label)) {
    emmean_info = emmeans_results |> filter(metric == y_label)
    max_y       = max(plot_data[[var_name]], na.rm = TRUE)
    emmean_info = emmean_info |> mutate(interval_y = max_y * 1.2)
    p = p + geom_text(data = emmean_info,
                      aes(x = factor(session_age), y = interval_y, label = stars),
                      inherit.aes = FALSE, size = 5, vjust = 1)
  }

  p = p +
    THEME_BASE + THEME_TEXT +
    scale_fill_manual(values = COLORS_MAIN) +
    scale_color_manual(values = COLORS_MAIN) +
    labs(x = "Visit", y = y_label, color = "Dataset Type", fill = "Dataset Type")

  plot_list_b[[var_name]] = p
}


psd_data_files <- list.files(path = path2data, pattern = "psds_aperosc_long_.*\\.csv$", full.names = TRUE)
psd_data_files <- psd_data_files[grepl('burst', psd_data_files)] 
psd_data_files <- psd_data_files[(grepl('region', psd_data_files) | grepl('Region', psd_data_files))] 

psd_data = lapply(psd_data_files, read_csv, show_col_types = FALSE) |> bind_rows() |> 
  mutate(session_age = if_else(session_age == 15, 18, session_age))|>
  filter(inclusion_final_dummy == T, sujid %in% descriptives$sujid) |>
  group_by(session_age, sujid, burst, freq)|>
  summarise(aperiodic = mean(aperiodic, na.rm = TRUE),
            oscillatory = mean(oscillatory, na.rm = TRUE), .groups = 'drop')|>
  pivot_wider(values_from = c('aperiodic', 'oscillatory'), names_from = burst)|>
  mutate(aperiodic_diff = aperiodic_burst - aperiodic_noburst,
         oscillatory_diff = oscillatory_burst - oscillatory_noburst) |>
  group_by(session_age, freq) |> 
  summarise(aper_mean = mean(aperiodic_diff, na.rm = TRUE),
            aper_first = quantile(aperiodic_diff, 0.25, na.rm = TRUE),
            aper_third = quantile(aperiodic_diff, 0.75, na.rm = TRUE),
            osc_mean = mean(oscillatory_diff, na.rm = TRUE),
            osc_first = quantile(oscillatory_diff, 0.25, na.rm = TRUE),
            osc_third = quantile(oscillatory_diff, 0.75, na.rm = TRUE), .groups = 'drop')

plot_aper_psd = ggplot(psd_data, aes(x = freq, y = aper_mean, color = factor(session_age))) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = aper_first, ymax = aper_third, fill = factor(session_age)), alpha = 0.20, color = NA) +
  scale_color_manual(values = COLORS_GRADIENT) +
  scale_fill_manual(values = COLORS_GRADIENT) +
  labs(x = "Frequency (Hz)", y = "log10(mV^2/Hz)", color = "Age (months)", fill = "Age (months)") +
  THEME_BASE + THEME_TEXT

plot_osc_psd = ggplot(psd_data, aes(x = freq, y = osc_mean, color = factor(session_age))) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = osc_first, ymax = osc_third, fill = factor(session_age)), alpha = 0.20, color = NA) +
  scale_color_manual(values = COLORS_GRADIENT) +
  scale_fill_manual(values = COLORS_GRADIENT) +
  labs(x = "Frequency (Hz)", y = "log10(mV^2/Hz)", color = "Age (months)", fill = "Age (months)") +
  THEME_BASE + THEME_TEXT

plot_psd_combined = ggpubr::ggarrange(plot_aper_psd, plot_osc_psd, ncol = 1, nrow = 2,
                                       labels = c('a', 'b'), common.legend = TRUE, legend = "bottom")



plot_list_aperosc = plot_list_b[c('offset', 'slope', 'peak_freq', 'peak_ampl', 'osc_ampl', 'peak_prop')]
plot_trajectories = ggpubr::ggarrange(plotlist = plot_list_aperosc, 
                                      ncol = 2, nrow = 3, 
                                      common.legend = TRUE, legend = "bottom", 
                                      labels = c('c', 'd', 'e', 'f', 'g', 'h'))


plot_trajectories_comb = ggpubr::ggarrange(plot_psd_combined, plot_trajectories, 
                                           ncol =2, 
                                           widths = c(1, 1.5)) 


plot_lcoh = ggpubr::ggarrange(lcoh_lifeplot_diff, lcoh_tileplot_burst,
                              ncol = 2, nrow = 1, labels = c('i', 'j'), 
                              widths = c(1, 1.5))

ggsave(
  filename = file.path(path2figs, "Fig_4_BurstImpact_ParametrizedPSD_LAcH.jpeg"),
  plot = ggpubr::ggarrange(plot_trajectories_comb, plot_lcoh, ncol = 1, nrow = 2,
                           heights = c(2, 1)),
  width = 30, height = 27, dpi = 300, units = "cm"
)

ggsave(filename = file.path(path2suppfig, "Extended_Data_Fig_4_BurstImpact_AlphaLifespan_by_Visit.jpeg"),
       plot = plot_list_b$alpha_LAcH, 
       width = 5, height = 4.5, dpi = 300)
