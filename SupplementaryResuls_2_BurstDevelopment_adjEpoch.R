# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term
#   Column identifiers below are an on-disk CSV / cross-script names can be found in Crosswalk Naming file
#   are intentionally NOT renamed (would break the pipeline + Github data). Display
#   labels and this map carry the manuscript terminology.
#     prop_bursty_epochs     -> Prop. epochs w/ Alpha burst
#     prop_bursty_cycles_burst -> Prop. cycles w/ alpha burst
#     is_burst / burst_type   -> cycle type: burst vs. non-burst
# -----------------------------------------------------------------------------
# Script: SupplementaryResults_2_BurstDevelopment_adjEpoch.R
# Purpose: Sensitivity GAMM analysis validating burst density trajectories using
#          epoch-duration-adjusted epochs. Mirrors the structure of script DataAnalysis_1....R
# =============================================================================
# Inputs: - Inputs (all read from path2data):
#   - Data/BurstProperties_ByCycle_Long_AdjEpoch.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Data/electrodes.csv
# Outputs (all saved to path2tabs / path2figs):
#   - Supplementary_Results_Table_WholeBrain_GAMM_Development_burst_adjEpoch.html
#   - Supplementary_Results_Table_WholeBrain_GAMM_Development_burstInteraction_adjEpoch.html
#   - FigSR2_WholeBrain_GAMM_BurstDevelopment_AdjEpoch.jpeg
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

epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)

# --- Full-GAMM trajectory display options (Fig SR2 curves) -------------------
# Prediction grid resolution, 95% Gaussian CI multiplier, and a display-only
# vertical alignment so the fitted curve overlays the raw scatter
GRID_N        = 200
CI_MULT       = 1.96
OVERLAY_SHIFT = TRUE
set.seed(RANDOM_SEED)

# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================

# path2data comes from config_paths.R
# path2root comes from config_paths.R

# Supplementary Results outputs (cycle-controlled adj-epoch GAMMs). See README.
path2tabs = file.path(path2root, "SupplementaryInformation", "Results", "Tables")  # range-comparison tables
path2figs = file.path(path2root, "SupplementaryInformation", "Results", "Figures")  # Fig SM1

if (!dir.exists(path2tabs)) dir.create(path2tabs, recursive = TRUE)
if (!dir.exists(path2figs)) dir.create(path2figs, recursive = TRUE)


# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================
# --- EEG cleaning descriptives ---
# Filter subjects by minimum clean epochs, standardize prop_epochs (z-score),
# and harmonize age labels (15 to 18, 40 to 42 months to align cohorts).
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],  # z-score normalization for model covariate
    session_age = if_else(session_age == 15, 18, session_age),  # harmonize age cohorts
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic / longitudinal age data ---
# Harmonize age cohorts (15 to 18, 40 to 42 months) in longitudinal data
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS, sujid %in% eeg_desc$sujid) |>
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# Wide format for covariates: extract all "mean" columns (e.g., PNC_mean, ITN_mean, etc.)
# and z-score them for model entry. Impute missing ITN_mean with sample mean.
desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  # ITN_mean is mean-imputed only if present. It is not a covariate in any
  # model here; it is swept in by select(contains("mean")) and is absent
  # from the public data release (see prepare_public_data.R). any_of()
  # makes this a no-op when the column is not there.
  mutate(across(any_of("ITN_mean"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x))) |>
  mutate(across(where(is.numeric), ~ scale(.x)[, 1]))  # z-score all numeric covariates

# --- Electrode map ---
# Standardize region labels from abbreviations (Fr, P, T, O) to full names
# and filter to included electrodes (chinclu == 1).
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",  # Recode region abbreviations
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1)

# Combined descriptives frame
descriptives = left_join(
  desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
  eeg_desc
) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs))

# --- Epoch-adjusted burst data ---
# Load adjusted-epoch burst properties, aggregate to region level (mean per subject×session×burst_state×region),
# standardize region labels (region to region), and harmonize age cohorts.
burst_data = read_csv(file.path(path2data, "BurstProperties_ByCycle_Long_AdjEpoch.csv")) |>
  merge(electrodes, by.x = "ch", by.y = "label") |>
  filter(chinclu == 1) |>
  group_by(session_age, sujid, is_burst, region) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  rename(region = region) |>  # Rename region -> region for consistency with descriptives merge
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# Merge with descriptives and apply final inclusion filter
burst_data = left_join(descriptives, burst_data) |>
  filter(inclusion_lmm == 1, dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS)

# =============================================================================
# PART A: WHOLE-BRAIN MAIN AGE EFFECTS (Burst presence metrics only)
# =============================================================================

metrics_config = list(
  list(metric = "prop_bursty_cycles_burst", type = "Oscillatory", data = "burst_data",    label = 'Prop. of Cycles\nw/ Burst', tablesub = 'burst'),
  list(metric = "prop_bursty_epochs",     type = "Oscillatory", data = "burst_data",      label = 'Prop. of Epochs\nw/ Burst', tablesub = 'burst'))

results_list_main = list()  # Stores one stats row per metric
main_plots        = list()  # Stores one ggplot per metric
main_models       = list()  # Stores the fitted GAMM (+ its data) per metric, reused for Fig SR2 curves

print("--- PART A: WHOLE-BRAIN MAIN AGE EFFECTS ---")

for (item in metrics_config) {
  f      = item$metric
  p      = item$type
  d_name = item$data
  flabel = item$label
  stable = item$tablesub
  
  print(paste("Running model for:", p, "-", f))
  
  # Select and prepare the correct dataset for this metric.
  # Standardize outcome variable name to 'pow' (power/amplitude proxy) for consistent model formula.
  if (d_name == "aperiodic_data") {
    model_data = aperiodic_data |> rename(pow = !!sym(f)) |>
      mutate(model_fit = scale(r2value)[, 1]) |> drop_na(pow)
  } else if (d_name == "hlcoh_data") {
    model_data = hlcoh_data    |> rename(pow = !!sym(f)) |> drop_na(pow)
  } else if (d_name == "burst_data") {
    # For burst presence metrics (prop_bursty_cycles_burst, prop_bursty_epochs):
    # retain only Burst cycles (is_burst == "Burst"), replace missing values with 0 (no burst detected).
    model_data = burst_data    |> rename(pow = !!sym(f)) |>
      filter(is_burst == "Burst") |>
      mutate(pow = if_else(is.na(pow), 0, pow)) |> drop_na(pow)
  }
  
  # Drop rows with missing numeric covariates to ensure valid GAMM fit.
  numeric_covars = c("age_months", "prop_epochs", "GestationalAge_weeks")
  if ("model_fit" %in% names(model_data)) numeric_covars = c(numeric_covars, "model_fit")  # aperiodic models only
  model_data = model_data |> drop_na(any_of(numeric_covars))
  
  n_unique_ids = length(unique(model_data$sujid))
  total_obs    = model_data |>
    group_by(sujid, session_age) |> summarise(.groups = "drop") |> nrow()
  
  if (n_unique_ids > 1 && nrow(model_data) > 2 * n_unique_ids) {
    
    # Build GAMM formula; conditionally add model_fit covariate for aperiodic models only.
    # Full model includes age smooth s(age_months, k=4); reduced model uses age as fixed effect.
    fixed_effects = "pow ~ s(age_months, k=4) + region + prop_epochs + Cohort + GestationalAge_weeks"

    # Reduced model (no age smooth) used to compute partial R² attributable to age trajectory.
    fixed_effects_red = "pow ~ region + prop_epochs + Cohort + GestationalAge_weeks"

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
    
    # Keep the fitted model + its data so Fig SR2 can draw the model's OWN
    # predicted trajectory (reused, not re-fit).
    main_models[[f]] = list(m = m, data = model_data)
    
    # Compute simultaneous derivatives to identify periods of significant change
    moments_of_change = derivatives(m$gam, interval = "simultaneous") |>
      mutate(inc_or_dec = case_when(
        .upper_ci >= 0 & .lower_ci >= 0 ~ "Increase",
        .upper_ci <= 0 & .lower_ci <= 0 ~ "Decrease",
        TRUE                            ~ "No Change"
      ))
    write_csv(moments_of_change,
              paste0(path2data, "/WholeBrain_MomentsOfChange_", p, "_", f, "_adjEpoch.csv"))
    
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
# The difference smooth tests whether Burst and No Burst trajectories diverge over age.

inter_voi          = c("volt_amp", "band_amp")
results_list_inter = list()
inter_plots        = list()
inter_models       = list()  # Stores the fitted interaction GAMM (+ its data) per metric, reused for Fig SR2 curves

print("--- PART B: BURST × AGE INTERACTION EFFECTS ---")

for (f in inter_voi) {
  
  if (f == "volt_amp") {flabel = "Volt. Amp."} else {flabel = "Band Amp."}  #  absolute burst-energy metrics
  stable = 'burst'
  
  print(paste("Running interaction for:", f))
  
  model_data = burst_data |> rename(pow = !!sym(f))

  # Create ordered factor (is_burst_ord) with treatment contrasts for difference smooth.
  # The difference smooth s(age_months, by = is_burst_ord) tests whether Burst and NoBurst
  # trajectories diverge over age (interaction effect).
  model_data = model_data |> mutate(is_burst = factor(is_burst))
  model_data$is_burst_ord = as.ordered(model_data$is_burst)
  contrasts(model_data$is_burst_ord) = "contr.treatment"  # Reference level (Burst) vs. NoBurst difference
  
  model_data = model_data |> drop_na(any_of(c("age_months", "prop_epochs", "GestationalAge_weeks")))
  
  
  n_unique_ids = length(unique(model_data$sujid))
  total_obs    = model_data |>
    group_by(sujid, session_age) |> summarise(.groups = "drop") |> nrow()
  
  
  # Full interaction model: includes both main age smooth s(age, k=4) and difference smooth
  # s(age, by=is_burst_ord). The difference smooth tests trajectory divergence between Burst vs. NoBurst.
  m = gamm(
    pow ~ is_burst_ord + s(age_months, k=4) + s(age_months, k=4, by = is_burst_ord) +
      region + prop_epochs + Cohort + GestationalAge_weeks,
    random  = list(sujid =~ 1),
    correlation = nlme::corCAR1(form = ~ age_months | sujid / region / is_burst_ord),  # CAR1: closer ages more correlated
    data    = model_data, method = "REML",
    control = nlme::lmeControl(maxIter = 100, msMaxIter = 100)
  )

  # Reduced model (no difference smooth) for computing partial R² of the interaction effect.
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
  
  # Keep the fitted interaction model + its data so Fig SR2 can draw the model's
  inter_models[[f]] = list(m = m, data = model_data)
  
  # Identify difference smooth term in model summary for derivative analysis.
  term_diff_name = grep("is_burst_ord", rownames(gam_results$s.table), value = TRUE)

  # Compute simultaneous derivatives of the difference smooth to identify periods where
  # Burst and NoBurst trajectories diverge (non-zero derivative, CI excludes 0).
  moments_diff = derivatives(m$gam, term = term_diff_name, interval = "simultaneous") |>
    mutate(inc_or_dec = case_when(
      .upper_ci >= 0 & .lower_ci >= 0 ~ "Diverging (+)",  # Burst increasing faster
      .upper_ci <= 0 & .lower_ci <= 0 ~ "Diverging (-)",  # NoBurst increasing faster
      TRUE                             ~ "Parallel"       # No significant divergence
    ))
  write_csv(moments_diff,
            paste0(path2data, "/WholeBrain_MomentsOfChange_", p, "_", f, "_adjEpoch.csv"))
  
  # Plot difference smooth: predicted gap between Burst and No Burst trajectories over age.
  # Departure from 0 indicates divergence; CI excludes 0 = significant difference.
  plot_smooth_data = smooth_estimates(m$gam, select = term_diff_name) |> add_confint()
  inter_plots[[f]] = ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_ribbon(data = plot_smooth_data,
                aes(x = age_months, ymin = .lower_ci, ymax = .upper_ci),
                fill = scales::alpha("black", 0.2)) +
    geom_line(data = plot_smooth_data, aes(x = age_months, y = .estimate),
              color = "black", linewidth = 1.2) +
    geom_rug(data = model_data, aes(x = age_months), sides = "b",
             length = grid::unit(0.02, "npc")) +
    labs(title    = paste(f, "— Difference (Burst − NoBurst)"),
         subtitle = "Departure from 0 = trajectories diverging",
         x = "Age (months)", y = "Difference") +
    THEME_BASE + THEME_TEXT
  
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
# Export main age trajectories and epoch-adjustment effects for each metric type.
if (length(results_list_main) > 0) {

  results_list_main = bind_rows(results_list_main)
  # --- Numbered supplementary table SR9 ---
  results_list_main |>
    group_by(stable) |>
    mutate(age_pval_fdr    = p.adjust(age_pval,    method = "fdr"),
           epochs_pval_fdr = p.adjust(epochs_pval, method = "fdr")) |>
    ungroup() |>
    arrange(stable, type_of_pow, metric) |>
    nn_supp_table(id = "SR9", cols = c(
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
        "Robustness check repeating the whole-brain burst GAMMs with epoch length",
        "rescaled per visit, so that each epoch can contain a constant number of",
        "alpha cycles regardless of the developmental shift in alpha peak frequency.",
        "Age smooths tested with F on the reference degrees of freedom; parametric",
        "terms are t tests on the model residual degrees of freedom. P values",
        "FDR-corrected (Benjamini-Hochberg) across metrics.",
        "EDF, effective degrees of freedom."))

  for (t in unique(results_list_main$stable)) {  # Filter by metric type (e.g., 'burst')
    results_list_main|>
      filter(stable == t) |>
      dplyr::select(-stable) |>
      # Apply FDR correction across all metrics for multiple-comparison control.
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
      save_as_html(path = file.path(path2tabs, paste0("Supplementary_Results_Table_WholeBrain_GAMM_Development_", t, "_adjEpoch.html")))
  }
}

# --- Table 2: Burst Interaction Effects ---
# Export difference smooth statistics testing Burst vs. NoBurst trajectory divergence.
if (length(results_list_inter) > 0) {
  # --- Numbered supplementary table SR10 ---
  bind_rows(results_list_inter) |>
    mutate(diff_pval_fdr = p.adjust(diff_pval, method = "fdr")) |>
    nn_supp_table(id = "SR10", cols = c(
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
        "Robustness check repeating the burst x age difference smooths with epoch",
        "length rescaled per visit. Difference smooths tested with F on the",
        "reference degrees of freedom; the cycle-type contrast is a t test on the",
        "model residual degrees of freedom. P values FDR-corrected",
        "(Benjamini-Hochberg) across metrics."))

  bind_rows(results_list_inter) |>
    mutate(diff_pval_fdr = p.adjust(diff_pval, method = "fdr")) |>  # FDR correction across metrics
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
    save_as_html(path = file.path(path2tabs, "Supplementary_Results_Table_WholeBrain_GAMM_Development_BurstInteraction_adjEpochLength.html"))
}
print("--- WHOLE-BRAIN MODELS (ADJUSTED EPOCH) COMPLETE ---")


# =============================================================================
# PART D: FIGURE SR2
# =============================================================================
# Plot adjusted-epoch burst development trajectories with fitted GAMM curves and change annotations.

change_files = list.files(path = path2data, pattern = "*adjEpoch\\.csv$",
                          full.names = TRUE)

# Define metric variables, display labels, and variable-to-label mappings for each plot panel.
variables      = c("prop_bursty_epochs", "prop_bursty_cycles_burst")
labels_main    = c("Prop. of Epochs\nw/ Burst", "Prop. of Cycles\nw/ Burst")
voi_burst      = c("volt_amp",            "band_amp")  # Amplitude metrics in interaction model
labels_burst   = c("Volt. Amp.", "Band Amp.")

# Prepare data for raw scatter plots: aggregate to one row per subject × age × burst_state,
# averaging across regions for cleaner visualization.
burst_data_plot = burst_data |>
  group_by(session_age, sujid, is_burst) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  mutate(is_adj = "Adjusted")  # Label for sensitivity analysis identification

# --- Full-GAMM trajectory predictors (reuse the models fit in Part A / Part B)
# Reference level function: use mean for numeric, first level for factor
ref_level = function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else levels(factor(x))[1]

# Generate predicted trajectories for main age-effect models.
# Prediction grid spans observed age range at reference levels of all other covariates.
# Optional OVERLAY_SHIFT aligns curve to raw data mean for visual comparison (shape unchanged).
predict_main_curve = function(fit_obj, align_to = NULL) {
  m  = fit_obj$m
  md = fit_obj$data
  agerng = range(md$age_months, na.rm = TRUE)
  grid = tibble(
    age_months           = seq(agerng[1], agerng[2], length.out = GRID_N),
    region               = ref_level(md$region),
    prop_epochs           = 0,   # z-scored -> 0 is sample mean
    Cohort               = ref_level(md$Cohort),
    GestationalAge_weeks = 0
  )
  pr = predict(m$gam, newdata = grid, se.fit = TRUE)
  grid$fit = as.numeric(pr$fit); grid$se = as.numeric(pr$se.fit)
  if (OVERLAY_SHIFT) {   # Display-only vertical shift to overlay with raw data (shape unchanged)
    ref = if (!is.null(align_to)) mean(align_to, na.rm = TRUE) else mean(md$pow, na.rm = TRUE)
    grid$fit = grid$fit + (ref - mean(grid$fit, na.rm = TRUE))
  }
  grid$lo = grid$fit - CI_MULT * grid$se  # 95% CI bounds
  grid$hi = grid$fit + CI_MULT * grid$se
  grid
}

# Generate predicted trajectories for interaction models (separate curves per burst state).
# Creates prediction grid for each level of is_burst_ord, preserving Burst/No Burst gap
predict_inter_curve = function(fit_obj, align_to = NULL) {
  m  = fit_obj$m
  md = fit_obj$data
  agerng = range(md$age_months, na.rm = TRUE)
  ages   = seq(agerng[1], agerng[2], length.out = GRID_N)
  levs   = levels(md$is_burst_ord)
  # Create grid with separate rows for each burst_state level
  grid = purrr::map_dfr(levs, function(lv) tibble(
    age_months           = ages,
    is_burst_ord         = ordered(lv, levels = levs),
    is_burst             = factor(lv, levels = levels(factor(md$is_burst))),
    region               = ref_level(md$region),
    prop_epochs           = 0,
    Cohort               = ref_level(md$Cohort),
    GestationalAge_weeks = 0
  ))
  pr = predict(m$gam, newdata = grid, se.fit = TRUE)
  grid$fit = as.numeric(pr$fit); grid$se = as.numeric(pr$se.fit)
  if (OVERLAY_SHIFT) {   # Single common shift preserves the Burst/NoBurst divergence
    ref = if (!is.null(align_to)) mean(align_to, na.rm = TRUE) else mean(md$pow, na.rm = TRUE)
    grid$fit = grid$fit + (ref - mean(grid$fit, na.rm = TRUE))
  }
  grid$lo = grid$fit - CI_MULT * grid$se  # 95% CI bounds
  grid$hi = grid$fit + CI_MULT * grid$se
  grid
}

# --- Main age-effect panels (burst density trajectories) ---
# For each burst-density metric, plot raw scatter, GAMM fit curve, and age-related-change annotations.
plot_list = list()
for (i in seq_along(variables)) {
  var_name    = variables[i]
  y_label     = labels_main[i]

  # Load derivative/change annotations from Part A output files
  change_match = change_files[grepl(var_name, change_files)]
  if (length(change_match) == 0) next

  change_data = read_csv(change_match[1], show_col_types = FALSE)
  # Position change annotation bars below the raw data range
  y_bar_pos   = min(burst_data_plot[[var_name]], na.rm = TRUE) -
    0.15 * sd(burst_data_plot[[var_name]], na.rm = TRUE)
  
  # Create segments of change annotations: identify transitions between Increase/Decrease/NoChange states
  bar_data = change_data |>
    filter(!is.na(age_months), !is.na(inc_or_dec)) |>
    arrange(age_months) |>
    mutate(status_change = inc_or_dec != lag(inc_or_dec, default = first(inc_or_dec)),
           segment_id    = cumsum(status_change),  # Unique ID per contiguous segment
           y_location    = y_bar_pos)

  # Generate predicted curve from main model if available
  gamm_fit = if (!is.null(main_models[[var_name]]))
    predict_main_curve(main_models[[var_name]],
                       align_to = burst_data_plot[[var_name]]) else NULL
  
  p_main = ggplot(burst_data_plot,
                  aes(x = age_months, y = .data[[var_name]])) +
    geom_point(alpha = 0.2, color = "grey60", size = 1,
               position = position_jitter(width = 0.2, height = 0)) +
    geom_line(data = bar_data, inherit.aes = FALSE,
              aes(x = age_months, y = y_location, color = inc_or_dec, group = segment_id),
              linewidth = 5, alpha = 0.5, lineend = "butt") +
    scale_color_manual(values = COLORS_NATURE) +
    THEME_BASE + THEME_TEXT +
    labs(x = "Age (months)", y = y_label, color = "Age-related change")
  
  # Overlay fitted GAMM curve with CI ribbon if model exists
  if (!is.null(gamm_fit)) {
    p_main = p_main +
      geom_ribbon(data = gamm_fit, inherit.aes = FALSE,
                  aes(x = age_months, ymin = lo, ymax = hi),
                  fill = colorspace::lighten(COLORS_MAIN[1], 0.3), alpha = 0.4) +
      geom_line(data = gamm_fit, inherit.aes = FALSE,
                aes(x = age_months, y = fit),
                color = COLORS_MAIN[1], linewidth = 1.2)
  }

  plot_list[[var_name]] = p_main
}

# --- Burst × age interaction panels (amplitude trajectories) ---
# For each amplitude metric, plot Burst vs. NoBurst divergence with change annotations.
plot_list_burst = list()
for (i in seq_along(voi_burst)) {
  var     = voi_burst[i]
  y_label = labels_burst[i]

  # Load divergence/change annotations from Part B output files
  change_match = change_files[grepl(var, change_files)]
  if (length(change_match) == 0) next

  change_data_inter = read_csv(change_match[1], show_col_types = FALSE)
  # Position divergence annotation bars below the raw data range
  y_bar_pos         = min(burst_data_plot[[var]], na.rm = TRUE) -
    0.15 * sd(burst_data_plot[[var]], na.rm = TRUE)

  # Create segments of divergence state: Diverging(+) vs. Diverging(-) vs. Parallel
  bar_data = change_data_inter |>
    filter(!is.na(age_months), !is.na(inc_or_dec)) |>
    arrange(age_months) |>
    mutate(status_change = inc_or_dec != lag(inc_or_dec, default = first(inc_or_dec)),
           segment_id    = cumsum(status_change),
           y_location    = y_bar_pos)

  # Generate predicted curves for each burst_state from interaction model
  gamm_fit = if (!is.null(inter_models[[var]]))
    predict_inter_curve(inter_models[[var]],
                        align_to = burst_data_plot[[var]]) else NULL
  
  # Create scatter plot colored by burst_state (Burst vs. NoBurst)
  p_burst = ggplot(burst_data_plot, aes(x = age_months, y = .data[[var]])) +
    geom_point(aes(color = is_burst), alpha = 0.2, size = 1,
               position = position_jitter(width = 0.2, height = 0))

  # Overlay fitted curves with separate lines and CIs per burst_state
  if (!is.null(gamm_fit)) {
    p_burst = p_burst +
      geom_ribbon(data = gamm_fit, inherit.aes = FALSE,
                  aes(x = age_months, ymin = lo, ymax = hi, fill = is_burst),
                  alpha = 0.33) +
      geom_line(data = gamm_fit, inherit.aes = FALSE,
                aes(x = age_months, y = fit, color = is_burst), linewidth = 1.2)
  }

  # Add divergence annotations (separate scale for second color axis)
  p_burst = p_burst +
    scale_color_manual(values = COLORS_MAIN, name = "Cycle Type",
                       guide = guide_legend(override.aes = list(alpha = 1, size = 3,
                                                                linewidth = 2, fill = NA))) +
    scale_fill_manual(values = COLORS_MAIN, name = "Cycle Type", guide = "none") +
    ggnewscale::new_scale_color() +  # New color scale for divergence annotations
    geom_line(data = bar_data, inherit.aes = FALSE,
              aes(x = age_months, y = y_location, color = inc_or_dec, group = segment_id),
              linewidth = 5, alpha = 0.5, lineend = "butt") +
    scale_color_manual(values = COLORS_NATURE_BURST, name = "Age-related change") +
    THEME_BASE + THEME_TEXT +
    labs(x = "Age (months)", y = y_label)

  if (i == 2) p_burst = p_burst + theme(legend.position = "none")  # Hide legend on second panel
  plot_list_burst[[var]] = p_burst
}

# Arrange figure panels: top row = main age effects (burst density), bottom row = interactions (amplitude)
plot_top    = ggpubr::ggarrange(plotlist = rev(plot_list),       ncol = 2, common.legend = TRUE, legend = "bottom", labels = c('a', 'b'))
plot_bottom = ggpubr::ggarrange(plotlist = plot_list_burst, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c('c', 'd'))
plot_final  = ggpubr::ggarrange(plot_top, plot_bottom, ncol = 1, nrow = 2)

# Export Figure SR2: GAMM trajectories for adjusted-epoch analysis
ggsave(filename = file.path(path2figs, "FigSR2_WholeBrain_GAMM_BurstDevelopment_AdjEpoch.jpeg"),
       plot = plot_final, width = 6.5, height = 6.5, dpi = 300)

