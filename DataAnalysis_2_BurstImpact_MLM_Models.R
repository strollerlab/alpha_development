# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term
#   Column identifiers below are an on-disk CSV / cross-script DATA CONTRACT and
#   are intentionally NOT renamed (would break the pipeline + OSF data). Display
#   labels and this map carry the manuscript terminology.
#     slope                   -> Aperiodic slope
#     offset                  -> Aperiodic offset
#     peak_freq               -> Alpha Peak Frequency
#     peak_ampl               -> Alpha Peak Amplitude
#     peak_prop                -> Alpha Peak Proportion
#     osc_ampl                -> Oscillatory Alpha Band Power
#     alpha_LAcH              -> Alpha lifespan (Lagged Auto coherence Hilbert - cycle in which cumulative sum == 90% of the total)
# -----------------------------------------------------------------------------
# Script: DataAnalysis_2_BurstImpact_MLM_Models.R
# Purpose: Tests whether EEG metrics differ between burst and non-burst cycle
#          epochs across development using Linear Mixed-Effects Models (LMER).
#
#          For each metric (aperiodic slope/offset, oscillatory amplitude/
#          frequency/percentage, alpha lifespan), fits:
#
#            pow ~ burst * factor(session_age) + prop_epochs + Cohort +
#                  GestationalAge_weeks +  [+ model_fit] + (1|sujid)
# =============================================================================
# Inputs:
#   - Data/Aperiodic_Oscillatory_ByCycle_Long_Burst.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long_Burst.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Scripts/electrodes.csv
# Outputs (to path2tabs):
#  - Complementary_Fig4_Table_MLM_BurstImpact_by_Visit.html
#  - Complementary_Fig4_Table_MLM_BurstImpact_by_Visit_**variable**.html
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

source(file.path(path2code, "SupplementaryTables_Helper.R"))


# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================
# Local copies of global thresholds from 00_Setup_PackageInstallation.R.

r2_thresh        = R2_THRESH         # Minimum Specparam R-squared (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)
set.seed(RANDOM_SEED)


n_bootstraps = 1000  # Iterations for bootstrapped confidence intervals

# Additional single-visit exclusion: it crashed during Specparam by burst. 
ADDITIONAL_EXCLUSIONS = c("SUB-XCM45B")


# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================
# IMPORTANT: Update these to match your local directory structure.

# IMPORTANT: Update only `path2root` to match your local directory structure.
# Sub-folders derive from it and mirror the manuscript organisation (README).
# path2data : merged/analysis-ready input data (read-only here)
# path2tabs : MLM burst-impact result tables
# path2desc : Extended Data Table 2 (by-burst descriptives)

path2data = path2sets   # EDIT   # = Data/; this script prepends "" to each file name
# path2root comes from config_paths.R

path2tabs = file.path(path2root, "MainText",  "Tables", "BurstImpact")  # burst-impact MLM tables
path2desc = file.path(path2root, "ExtendedData",  "Tables")             # Extended Data Table 2 (descriptive)

for (p in c(path2tabs, path2desc)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- Electrode map ---
# chinclu == 1 marks the 60 pre-selected analysis channels.
# Remap electrode region abbreviations to full display labels for manuscript tables
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1)

# --- EEG cleaning descriptives ---
# prop_epochs is z-scored so it is on a comparable scale to other covariates.
# inclusion_lmm = 1 marks participants meeting minimum epoch count.
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(
    clean_epochs_rest >= epochs_threshold,
    !sujid %in% EXCLUDED_SUBJECTS
  ) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic / longitudinal age data ---
# GestationalAge_weeks and z-scored SES covariates are used as confound controls.
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(
    dev_filter == 1,
    !sujid %in% EXCLUDED_SUBJECTS,
    sujid %in% eeg_desc$sujid
  ) |>
  mutate(
    # Harmonize session_age labels: raw CSV uses 15to18mo and 40to42mo
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# Wide format for covariate extraction (one row per subject, across visits)
desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  # Mean-impute ITN_mean so no participant is dropped due to a single missing SES value
  # ITN_mean is mean-imputed only if present. It is not a covariate in any
  # model here; it is swept in by select(contains("mean")) and is absent
  # from the public data release (see prepare_public_data.R). any_of()
  # makes this a no-op when the column is not there.
  mutate(across(any_of("ITN_mean"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x))) |>
  mutate(across(where(is.numeric), ~ scale(.x)[, 1]))

# Combined descriptives frame (left anchor for all data merges)
descriptives = left_join(
    desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
    eeg_desc
  ) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs))

# --- Burst-conditioned aperiodic / oscillatory data ---
# This file contains separate Specparam estimates for burst-cycle and non-burst-cycle
# epochs (column: burst = "burst" | "noburst").
# NOTE: 'mae' is the Specparam MAE column
# It is filtered here under its raw column name, and then renamed to 'mae'
#'epochs' is renamed from 'epochs'.

aper_voi = c("sujid", "session_age", "ch", "region", "chinclu", "r2value",
              "goodch", "offset", "slope", "alpha_freq", "alpha_ampl",
              "alpha_osc", "alpha_peak", "inclusion_final_dummy",
              "epochs", "burst", "mae")

aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long_Burst.csv")) |>
  dplyr::select(all_of(aper_voi)) |>
  # Rename epochc to epochs for consistency across datasets
  rename(epochs = epochs) |>
  # Quality filters: R², MAE, analysis channels, minimum epochs, min good channels
  filter(
    r2value > r2_thresh,
    mae     < mae_thresh,
    chinclu == 1,
    epochs  >= epochs_threshold,
    goodch  >= CH_THRESHOLD,
    !sujid %in% ADDITIONAL_EXCLUSIONS
  ) |>
  # Collapse to one row per subject × visit × burst condition
  group_by(session_age, sujid, burst) |>
  summarise(
    slope     = mean(slope,                       na.rm = TRUE),
    offset    = mean(offset,                      na.rm = TRUE),
    peak_ampl = mean(alpha_ampl,                  na.rm = TRUE),
    peak_freq = mean(alpha_freq[alpha_peak == 1], na.rm = TRUE),
    osc_ampl  = mean(alpha_osc,                   na.rm = TRUE),
    peak_prop  = mean(alpha_peak,                  na.rm = TRUE),
    epochs    = mean(epochs,                      na.rm = TRUE),
    r2value   = mean(r2value,                     na.rm = TRUE),
    mae       = mean(mae,                         na.rm = TRUE),  # mae to mae rename
    .groups   = "drop"
  ) |>
  mutate(
    # Harmonize session_age labels: raw CSV uses 15to18mo and 40to42mo (see PROJECT_NAMING_CONVENTIONS.md)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age),
    # Standardize burst labels to match all other scripts ("Burst" / "Non-Burst") and set factor levels
    burst       = if_else(burst == "burst", "Burst", "Non-Burst"),
    burst       = factor(burst, levels = c("Non-Burst", "Burst"))
  )

# --- Burst-conditioned lagged coherence data ---
# goodch counts channels per subject × session × burst condition.
# Subjects below CH_THRESHOLD good channels are excluded from this dataset.
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long_Burst.csv")) |>
  merge(electrodes |> dplyr::select(label, chinclu), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  mutate(goodch = n(), .by = c(session_age, sujid, burst)) |>
  filter(goodch >= CH_THRESHOLD) |>
  group_by(session_age, sujid, burst) |>
  summarise(
    alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE),
    epochs     = mean(epochs,     na.rm = TRUE),
    .groups    = "drop"
  ) |>
  mutate(
    # Harmonize session_age labels: raw CSV uses 15to18mo and 40to42mo (see PROJECT_NAMING_CONVENTIONS.md)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age),
    # Standardize burst labels to match all other scripts ("Burst" / "Non-Burst") and set factor levels
    burst       = if_else(burst == "burst", "Burst", "Non-Burst"),
    burst       = factor(burst, levels = c("Non-Burst", "Burst"))
  )

# Merge with descriptives and apply final inclusion/exclusion criteria
aperiodic_data = left_join(descriptives, aperiodic_data, by = c("sujid", "session_age")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)

hlcoh_data = left_join(descriptives, hlcoh_data, by = c("sujid", "session_age")) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)

# Generate the descriptive tables
# Define canonical visit order and map variable names to display labels (Table 1 nomenclature)
visit_levels = paste0(c(1, 6, 12, 18, 30, 36, 42, 48), " mo.")
var_levels = c("slope", "offset", "peak_freq", "osc_ampl", "peak_prop", "peak_ampl", "alpha_LAcH")
var_labels = c("slope" = "Slope", "offset" = "Offset", "peak_freq" = "Peak Freq.",
               "osc_ampl" = "Band Power", "peak_prop" = "Prop. of Peaks",
               "peak_ampl" = "Peak Amp.", "alpha_LAcH" = "Lifespan")

aperiodic_data_table = dplyr::select(aperiodic_data, sujid, session_age, burst, slope, offset, peak_freq,
              osc_ampl, peak_prop, peak_ampl)
hlcoh_data_table = dplyr::select(hlcoh_data, sujid, session_age, burst, alpha_LAcH)

data_table = left_join(aperiodic_data_table, hlcoh_data_table)|>
  pivot_longer(cols = c(slope, offset, peak_freq, osc_ampl, peak_prop, peak_ampl, alpha_LAcH),
               names_to = "Variable", values_to = "Value") |>
  # Format session_age labels and apply canonical visit and variable ordering
  mutate(session_age = paste0(session_age, " mo."),
         session_age = factor(session_age, levels = visit_levels),
         Variable    = factor(Variable,    levels = var_levels, labels = var_labels)) |>
  group_by(burst, session_age, Variable) |>
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(burst))|>
  mutate(Display_Stat = sprintf("%.3f\n(%.3f)", Mean, SD)) |>
  dplyr::select(Variable, burst, session_age, Display_Stat) |>
  pivot_wider(names_from = session_age, values_from = Display_Stat) |>
  arrange(Variable, burst) |>
  rename(Metric = Variable)|>
  flextable() |>
  merge_v(j = "Metric") |>
  valign(j = "Metric", valign = "top") |>
  set_caption("**EEG Metric Descriptives by Region and Visit (Mean (SD))**") |>
  add_footer_lines("Values are Mean (SD). Rows grouped by metric; columns are visit ages.") |>
  theme_booktabs() |> fix_border_issues() |> autofit() |>
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(j = 2:9, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2desc, "Extended_Data_Table_2_ParametrizedPSDandLifespan_Descriptives_by_burstdataset.docx"))

# =============================================================================
# SECTION 4: LMER MODELS — BURST DISTINCTION ACROSS DEVELOPMENT
# =============================================================================
# For each metric, the model tests:
#   (1) Main effect of burst: do Burst cycles differ from NoBurst cycles on average?
#   (2) Main effect of visit
#   (3) Burst × Visit interaction: does the Burst – NoBurst difference vary across visits?
#
# Bootstrapped 95% CIs are computed for all fixed-effect coefficients.
# When the interaction is significant (FDR-corrected p < .05), per-visit follow-up
# LMERs estimate the burst effect separately at each study visit.

metrics_config = list(
  list(metric = "slope",      type = "Aperiodic",   data = "aperiodic_data", label = 'Slope'),
  list(metric = "offset",     type = "Aperiodic",   data = "aperiodic_data", label = "Offset"),
  list(metric = "peak_freq",  type = "Oscillatory", data = "aperiodic_data", label = "Peak Freq."),
  list(metric = "osc_ampl",   type = "Oscillatory", data = "aperiodic_data", label = "Band Power"),
  list(metric = "peak_prop",   type = "Oscillatory", data = "aperiodic_data", label = "Prop. of Peaks"),
  list(metric = "peak_ampl",  type = "Oscillatory", data = "aperiodic_data", label = "Peak Amp."),
  list(metric = "alpha_LAcH", type = "Oscillatory", data = "hlcoh_data", label = "Lifespan")
)

results_list_main  = list()  # Main effects (age + burst) per metric
results_list_inter = list()  # Interaction effects (burst × age) per metric
emmeans_list       = list()  # Per-age follow-up results (when interaction significant)
inter_plots        = list()  # One plot per metric


for (item in metrics_config) {
  f      = item$metric
  p      = item$type
  d_name = item$data
  ftable = item$label

  # Select and prepare the appropriate dataset for this metric
  if (d_name == "aperiodic_data") {
    model_data = aperiodic_data |>
      rename(pow = !!sym(f)) |>
      # Z-score Specparam model fit (r2value) to use as a quality covariate
      mutate(model_fit = scale(r2value)[, 1]) |>
      # Z-score trial count within visit to remove age-related differences
      # in epoch length from the covariate (older children have longer clean epochs)
      mutate(epochs    = scale(epochs)[, 1], .by = session_age) |>
      drop_na(pow)

  } else if (d_name == "hlcoh_data") {
    model_data = hlcoh_data |>
      rename(pow = !!sym(f)) |>
      mutate(epochs = scale(epochs)[, 1], .by = session_age) |>
      drop_na(pow)
    # LAcH does not have a model_fit (Specparam R²) covariate
  }

  # Count unique subjects and total observations for the results tables
  n_unique_ids = length(unique(model_data$sujid))
  total_obs    = nrow(model_data)

  if (n_unique_ids <= 1 || total_obs <= 2 * n_unique_ids) {
    cat(sprintf("  SKIPPING %s: insufficient data (N=%d, obs=%d).\n",
                f, n_unique_ids, total_obs))
    next
  }

  # --- Build model formula ---
  # Core formula tests burst × age interaction; model_fit (R²) included for aperiodic
  # data only. All fixed effects are z-scored where applicable (epochs, prop_epochs, model_fit).
  fixed_effects = paste(
    "pow ~ burst * session_age + prop_epochs + Cohort +",
    "GestationalAge_weeks + (1|sujid)"
  )
  if ("model_fit" %in% names(model_data)) {
    fixed_effects = paste(fixed_effects, "+ model_fit")
  }

  m = lmer(as.formula(fixed_effects), data = model_data, REML = TRUE)

  # Extract model performance metrics (marginal and conditional R²)
  perf = performance::r2(m)

  model_summary = summary(m)
  coefs         = as.data.frame(model_summary$coefficients)

  # Identify coefficient row names robustly using grep to extract burst main effect and interaction
  burst_row = grep("^burstBurst$", rownames(coefs), value = TRUE)[1]
  inter_row = grep(":",             rownames(coefs), value = TRUE)[1]

  # Compute bootstrapped 95% CIs for all fixed-effect coefficients (drift-correction method)
  boot_ci = confint(m, parm = "beta_", method = "boot",
                    nsim = n_bootstraps, oldNames = FALSE)

  # --- Store main effects results ---
  # SE columns use the standardised suffix (not the misspelled 'Ste')
  results_list_main[[f]] = data.frame(
    metric         = ftable,
    type_of_pow    = p,
    n              = n_unique_ids,
    total_obs      = total_obs,
    # Age smooth coefficients
    Visit_estimate   = coefs[inter_row,   "Estimate"],
    Visit_SE         = coefs[inter_row,   "Std. Error"],
    Visit_CI_lower   = boot_ci[inter_row,  "2.5 %"],
    Visit_CI_upper   = boot_ci[inter_row,  "97.5 %"],
    Visit_t          = coefs[inter_row,   "t value"],
    Visit_df         = coefs[inter_row,   "df"],
    Visit_pval       = coefs[inter_row,   "Pr(>|t|)"],
    # Burst main effect coefficients (Burst vs. NoBurst reference)
    Burst_estimate = coefs[burst_row, "Estimate"],
    Burst_SE       = coefs[burst_row, "Std. Error"],
    Burst_CI_lower = boot_ci[burst_row, "2.5 %"],
    Burst_CI_upper = boot_ci[burst_row, "97.5 %"],
    Burst_t        = coefs[burst_row, "t value"],
    Burst_df       = coefs[burst_row, "df"],
    Burst_pval     = coefs[burst_row, "Pr(>|t|)"],
    # Effect size 
    Marginal_R2 = perf$R2_marginal,
    Conditional_R2 = perf$R2_conditional
  )

  inter_pval = coefs[inter_row, "Pr(>|t|)"]

  # --- Per-visit follow-up LMERs (run only when interaction is significant) ---
  # When burst × visit interaction p < 0.05, fit separate LMER at each visit to test
  # burst effect independently (without visit as a factor). Data collapsed to one row
  # per subject × burst condition before fitting.
  if (!is.na(inter_pval) && inter_pval < 0.05) {

    # Simplified formula for follow-up (no age term; already stratified by visit)
    sub_fixed_effects = "pow ~ burst + prop_epochs + Cohort + GestationalAge_weeks + (1|sujid)"
    if ("model_fit" %in% names(model_data)) {
      sub_fixed_effects = paste(sub_fixed_effects, "+ model_fit")
    }

    sub_results = list()

    for (age in sort(unique(model_data$session_age))) {

      # Collapse to one row per subject × burst condition
      # Preserves scalar covariates (Cohort, GestationalAge_weeks) with first().
      sub_model_data = model_data |>
        filter(session_age == age) |>
        group_by(sujid, burst) |>
        summarise(
          pow                  = mean(pow,                  na.rm = TRUE),
          prop_epochs           = mean(prop_epochs,           na.rm = TRUE),
          Cohort               = first(Cohort),
          GestationalAge_weeks = first(GestationalAge_weeks),
          model_fit            = if ("model_fit" %in% names(model_data))
                                   mean(model_fit, na.rm = TRUE) else NA_real_,
          .groups = "drop"
        )

      # Remove model_fit column if not needed (avoids NA-in-formula error for hlcoh)
      if (!"model_fit" %in% names(model_data)) {
        sub_model_data = dplyr::select(sub_model_data, -any_of("model_fit"))
      }

      # Descriptive stats for Burst and NoBurst at this visit
      burst_descs = sub_model_data |>
        group_by(burst) |>
        summarise(
          m_pow = mean(pow, na.rm = TRUE),
          s_pow = sd(pow,   na.rm = TRUE),
          .groups = "drop"
        )

      # Fit per-visit LMER
      sub_m = tryCatch(
        lmer(as.formula(sub_fixed_effects), data = sub_model_data, REML = TRUE),
        error = function(e) {
          cat(sprintf("  Sub-model failed at age %d: %s\n", age, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(sub_m)) next

      sub_perf    = performance::r2(sub_m)
      sub_coefs   = as.data.frame(summary(sub_m)$coefficients)
      sub_brow    = grep("^burstBurst$", rownames(sub_coefs), value = TRUE)[1]

      cat(sprintf("  Running %d bootstrap iterations for age %d mo...\n",
                  n_bootstraps, age))
      boot_ci_sub = confint(sub_m, parm = "beta_", method = "boot",
                            nsim = n_bootstraps, oldNames = FALSE)

      # Format descriptive statistics as "Mean (SD)" strings
      fmt_msd = function(df, burst_label) {
        r = df[df$burst == burst_label, ]
        if (nrow(r) == 0) return("—")
        sprintf("%.3f (%.3f)", r$m_pow, r$s_pow)
      }

      sub_results[[paste0(f, "_", age)]] = data.frame(
        metric         = ftable,
        visit          = age,
        Burst_MSD      = fmt_msd(burst_descs, "Burst"),
        NoBurst_MSD    = fmt_msd(burst_descs, "Non-Burst"),
        Burst_estimate = sub_coefs[sub_brow, "Estimate"],
        Burst_SE       = sub_coefs[sub_brow, "Std. Error"],
        Burst_CI_lower = boot_ci_sub[sub_brow, "2.5 %"],
        Burst_CI_upper = boot_ci_sub[sub_brow, "97.5 %"],
        Burst_t        = sub_coefs[sub_brow, "t value"],
        Burst_df       = sub_coefs[sub_brow, "df"],
        Burst_pval     = sub_coefs[sub_brow, "Pr(>|t|)"],
        Marginal_R2    = sub_perf$R2_marginal,
        Conditional_R2 = sub_perf$R2_conditional
      )
    }

    # Combine per-visit results and apply FDR correction within this metric (across all visits)
    emm_df = bind_rows(sub_results) |>
      mutate(Burst_pval_fdr = p.adjust(Burst_pval, method = "fdr"))

    emmeans_list[[f]] = emm_df
  }
}


# =============================================================================
# SECTION 5: TABLE GENERATION & EXPORT
# =============================================================================
# All tables use standardised formatting:
#   - Greek β symbol (Unicode \u03B2) for coefficient columns
#   - R² symbol (\u00b2)
#   - Bold rows where FDR-corrected p < .05
#   - theme_booktabs() + fix_border_issues() + autofit()
#   - Footer notes describing the model, FDR correction, and CI method

# ---------------------------------------------------------------------------
# TABLE: Main Effect of Burst and Interaction (Burst*Session_Age (i.e., Visit))
# ---------------------------------------------------------------------------
# One row per metric. Reports the main effect of age_months and the main effect
# of burst (Burst vs. NoBurst reference) with bootstrapped 95% CIs.

if (length(results_list_main) > 0) {

  # --- Numbered supplementary table SR3 ---
  bind_rows(results_list_main) |>
    mutate(Visit_pval_fdr = p.adjust(Visit_pval, method = "fdr"),
           Burst_pval_fdr = p.adjust(Burst_pval, method = "fdr")) |>
    nn_supp_table(id = "SR3", cols = c(
      metric         = "Metric",
      n              = "n (children)",
      total_obs      = "n (observations)",
      Burst_estimate = "Beta (burst vs non-burst)",
      Burst_SE       = "s.e.",
      Burst_CI_lower = "95% CI lower",
      Burst_CI_upper = "95% CI upper",
      Burst_t        = "t",
      Burst_df       = "df",
      Burst_pval     = "P",
      Burst_pval_fdr = "P (FDR)",
      Visit_estimate = "Beta (burst x visit)",
      Visit_SE       = "s.e. (interaction)",
      Visit_CI_lower = "95% CI lower (interaction)",
      Visit_CI_upper = "95% CI upper (interaction)",
      Visit_t        = "t (interaction)",
      Visit_df       = "df (interaction)",
      Visit_pval     = "P (interaction)",
      Visit_pval_fdr = "P (interaction, FDR)",
      Marginal_R2    = "Marginal R2",
      Conditional_R2 = "Conditional R2"),
      note = paste(
        "Linear mixed models: Metric ~ Epoch type (burst vs non-burst) x Visit +",
        "ROI + Cohort + covariates + (1|child). Degrees of freedom from",
        "Satterthwaite's approximation. 95% CIs from parametric bootstrapping",
        "(1,000 iterations). P values FDR-corrected (Benjamini-Hochberg) across",
        "metrics. Positive beta indicates a higher value in burst-containing epochs."))

  bind_rows(results_list_main) |>
    mutate(
      # FDR correction (Benjamini-Hochberg) applied separately to age and burst p-value vectors
      Visit_pval_fdr   = p.adjust(Visit_pval,   method = "fdr"),
      Burst_pval_fdr = p.adjust(Burst_pval, method = "fdr")
    ) |>
    # Create significance star columns (*, **, ***, or empty) for table display
    mutate(
      Age_sig   = case_when(Visit_pval_fdr   < 0.001 ~ "***",
                             Visit_pval_fdr   < 0.01  ~ "**",
                             Visit_pval_fdr   < 0.05  ~ "*",  TRUE ~ ""),
      Burst_sig = case_when(Burst_pval_fdr < 0.001 ~ "***",
                             Burst_pval_fdr < 0.01  ~ "**",
                             Burst_pval_fdr < 0.05  ~ "*",  TRUE ~ "")
    ) |>
    # Format bootstrapped 95% CI bounds as [lower, upper] strings for table display
    mutate(
      Visit_CI   = sprintf("[%.3f, %.3f]", Visit_CI_lower,   Visit_CI_upper),
      Burst_CI = sprintf("[%.3f, %.3f]", Burst_CI_lower, Burst_CI_upper)
    ) |>
    dplyr::select(
      metric, n, total_obs,
      Burst_estimate, Burst_SE, Burst_CI, Burst_df, Burst_pval, Burst_pval_fdr, Burst_sig,
      Visit_estimate, Visit_SE, Visit_CI, Visit_df, Visit_pval, Visit_pval_fdr, Age_sig, Marginal_R2, Conditional_R2
    ) |>
    # Format as flextable: one metric per row with burst main effect (β, SE, CI, p) and interaction (β, SE, CI, p) columns
    flextable() |>
    set_header_labels(
      metric         = "Metric",
      n              = "N",
      total_obs      = "N Obs.",
      Visit_estimate   = "\u03B2 (Burst*Visit Age)",
      Visit_SE         = "SE",
      Visit_CI         = "95% CI",
      Visit_df         = "df",
      Visit_pval       = "p",
      Visit_pval_fdr   = "p (FDR)",
      Age_sig        = "",
      Burst_estimate = "\u03B2 (Burst)",
      Burst_SE       = "SE",
      Burst_CI       = "95% CI",
      Burst_df       = "df",
      Burst_pval     = "p",
      Burst_pval_fdr = "p (FDR)",
      Burst_sig      = "",
      Marginal_R2    = "Marginal\n R\u00b2",
      Conditional_R2 = "Conditional\nR\u00b2"
    ) |>
    add_header_row(
      values = c("", "Burst Effect", "Interaction (Burst x Visit)", ""),
      colwidths = c(3, 7, 7, 2) # 3 general info columns, 7 burst columns, 7 age columns
    ) |>
    # END OF ADDED BLOCK
    set_caption("Table 1: LMER Main Effects — Visit and Burst on EEG Metrics (Whole Brain)") |>
    set_caption("LMER Visit*Burst on EEG Metrics (Whole Brain)") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("Visit_pval", "Visit_pval_fdr", "Burst_pval", "Burst_pval_fdr"),
                     digits = 4) |>
    # Bold the entire row when FDR-corrected significance is met
    bold(i = ~ Visit_pval_fdr   < 0.05, j = c("Visit_pval",   "Visit_pval_fdr",   "Age_sig"),   part = "body") |>
    bold(i = ~ Burst_pval_fdr < 0.05, j = c("Burst_pval", "Burst_pval_fdr", "Burst_sig"), part = "body") |>
    add_footer_lines(paste0(
      "Model: Metric ~ Burst * Visit Age + Prop. of Epochs Retained +  # Trials + Cohort + ",
      "Gestational Age (weeks) [+ Parametrization Model Fit] + (1|sujid). ",
      "Reference level for burst: NoBurst. ",
      "FDR correction (Benjamini-Hochberg) applied across all metrics separately for Age and Burst. ",
      "95% CIs from parametric bootstrapping (", n_bootstraps, " iterations). ",
      "*p < .05; **p < .01; ***p < .001."
    )) |>
    theme_booktabs() |> fix_border_issues() |> autofit() |>
    set_table_properties(layout = "fixed") |>
    align(align = "left", part = "all") |>
    align(j = 2:14, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    align(i = 1, part = "header", align = "center") |> # Centers the new grouped headers
    save_as_html(path = file.path(path2tabs, "Complementary_Fig4_Table_MLM_BurstImpact_by_Visit.html"))

}

# ---------------------------------------------------------------------------
# TABLE: Per-Visit Follow-Up (Emmeans)
# ---------------------------------------------------------------------------
# Rows per metric × visit (only for metrics with a significant interaction).
# Reports burst–NoBurst differences at each visit with bootstrapped CIs,
# descriptive means, and FDR-corrected p-values within each metric.

if (length(emmeans_list) > 0) {

  # Combine per-age follow-up results from significant interactions; format CI and significance stars
  # --- Numbered supplementary table SR4 ---
  bind_rows(emmeans_list) |>
    arrange(metric, visit) |>
    nn_supp_table(id = "SR4", cols = c(
      metric         = "Metric",
      visit          = "Visit (months)",
      Burst_MSD      = "Burst, mean (s.d.)",
      NoBurst_MSD    = "Non-burst, mean (s.d.)",
      Burst_estimate = "Beta",
      Burst_SE       = "s.e.",
      Burst_CI_lower = "95% CI lower",
      Burst_CI_upper = "95% CI upper",
      Burst_t        = "t",
      Burst_df       = "df",
      Burst_pval     = "P",
      Burst_pval_fdr = "P (FDR)",
      Marginal_R2    = "Marginal R2",
      Conditional_R2 = "Conditional R2"),
      note = paste(
        "Visit-stratified linear mixed models, fitted where the epoch type x visit",
        "interaction was significant: Metric ~ Epoch type + ROI + Cohort +",
        "covariates + (1|child). Degrees of freedom from Satterthwaite's",
        "approximation. 95% CIs from parametric bootstrapping (1,000 iterations).",
        "P values FDR-corrected (Benjamini-Hochberg) within metric across visits."))

  final_emmeans_df = bind_rows(emmeans_list)|>
    dplyr::select(
      metric, visit, Burst_MSD, NoBurst_MSD,
      Burst_estimate, Burst_SE, Burst_CI_lower, Burst_CI_upper,
      Burst_df, Burst_pval, Burst_pval_fdr, Marginal_R2, Conditional_R2
    )|> mutate(
      Burst_sig = case_when(Burst_pval_fdr < 0.001 ~ "***",
                            Burst_pval_fdr < 0.01  ~ "**",
                            Burst_pval_fdr < 0.05  ~ "*",  TRUE ~ "NS"),
      Burst_CI  = sprintf("[%.3f, %.3f]", Burst_CI_lower, Burst_CI_upper),
      visit     = paste0(visit, " mo.")
    ) |>
    dplyr::select(
      metric, visit, Burst_MSD, NoBurst_MSD,
      Burst_estimate, Burst_SE, Burst_CI, Burst_df, Burst_pval, Burst_pval_fdr, Burst_sig,
      Marginal_R2, Conditional_R2
    ) 

  # Save as CSV for use in script 16 significance overlays
  write.csv(final_emmeans_df,
            file      = paste0(path2data,
                                "/Extended_Results_MLM_Emmeans_Stratified_by_visit.csv"),
            row.names = FALSE)

  
  for (v in unique(final_emmeans_df$metric)) {
    final_emmeans_df|>
      filter(metric == v) |>
    flextable() |>
    # Merge repeated metric cells vertically for clean reading
    merge_v(j = "metric") |>
    set_header_labels(
      metric         = "Metric",
      visit          = "Visit",
      Burst_MSD      = "Burst Mean (SD)",
      NoBurst_MSD    = "NoBurst Mean (SD)",
      Burst_estimate = "\u03B2 (Burst)",
      Burst_SE       = "SE",
      Burst_CI       = "95% CI",
      Burst_df       = "df",
      Burst_pval     = "p",
      Burst_pval_fdr = "p (FDR)",
      Burst_sig      = "",
      Marginal_R2    = "R\u00b2 Marginal",
      Conditional_R2 = "R\u00b2 Conditional"
    ) |>
    set_caption("Visit stratified MLM Follow-Up — Burst vs. NoBurst Differences at Each Visit") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("Burst_pval", "Burst_pval_fdr"), digits = 4) |>
    colformat_double(j = c("Marginal_R2", "Conditional_R2"), digits = 3) |>
    # Bold the p-values and stars for significant rows
    bold(i = ~ Burst_pval_fdr < 0.05,
         j = c("Burst_pval", "Burst_pval_fdr", "Burst_sig"), part = "body") |>
    # set the first row bold 
    add_footer_lines(paste0(
      "Visit stratified follow-up models fitted only for metrics with a significant Burst \u00d7 Age interaction ",
      "(FDR-corrected p < .05). Model: Metric ~ Burst + Prop. of Epochs Retained + Cohort + ",
      "Gestational Age (weeks) [+ Parametrization Model Fit] + (1|sujid), averaged across brain. ",
      "FDR correction (Benjamini-Hochberg) applied within each metric across all visits. ",
      "95% CIs from parametric bootstrapping (", n_bootstraps, " iterations). ",
      "*p < .05; **p < .01; ***p < .001."
    )) |>
    flextable::font(fontname = "Arial", part = "all") |>
    fontsize(size = 10, part = "all") |>
    theme_booktabs() |> fix_border_issues() |> autofit() |>
      bold(i = 1, part = "header") |>
    save_as_html(path = file.path(path2tabs, paste0("Complementary_Fig4_Table_MLM_BurstImpact_by_Visit_", v, ".html")))

}
} else {
  cat("No significant interactions found — per-age follow-up tables not generated.\n")
}
