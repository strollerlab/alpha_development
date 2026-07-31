# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term to columns
#   DATA columns to DISPLAY LABELS (used in plot y-axis/legend labels):
#     alpha_LAcH              → "Lifespan" (Alpha lifespan in cycles; LAcH cumsum ≥90%)
# -----------------------------------------------------------------------------
# Script: SupplementaryMethods_3_Lifespan_RangeSelection.R
# Purpose: Supplementary analysis comparing two alpha frequency band definition
#          strategies to compute alpha Lifespan:
#            OB (Only-Burst): band range estimated from burst-cycle peak frequencies
#            AC (All-Cycles): band range estimated from all-cycle peak frequencies
# =============================================================================
# Inputs:
#   - Data/frequency_bands_onlyburst_adapted_*.csv   (one per age)
#   - Data/frequency_bands_allcycles_adapted_*.csv   (one per age)
#   - Data/LaggedCoh_Hilb_ByCycle_Long.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long_allcycles.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Scripts/electrodes.csv
# Outputs (all saved to path2tabs or path2figs):
#   - Supplement_FrequencyBands_OB_vs_AC.docx
#   - lagged_coherence_summary_alpharange_selection.docx
#   - lmm_results_alpha_range_selection.docx
#   - emm_posthoc_results_alpha_range_selection.docx
#   - alpha_lifespan_by_freq_selection.jpeg
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

source(file.path(path2code, "SupplementaryTables_Helper.R"))
source(file.path(path2code, "01_Utils_ProcFunctions.R"))  # Provides median_se()


# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================

epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode 
n_bootstraps     = 1000              # Bootstrapped CI iterations
set.seed(RANDOM_SEED)

# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================

path2freqs = path2sets        # EDIT — root folder for frequency band CSVs
# path2data comes from config_paths.R
# path2root comes from config_paths.R

# Supplementary Methods outputs (LAcH alpha-range selection). See README taxonomy.
path2tabs = file.path(path2root, "SupplementaryInformation",  "Methods", "Tables")  # range-comparison tables
path2figs = file.path(path2root, "SupplementaryInformation", "Methods", "Figures")  # Fig SM1

if (!dir.exists(path2tabs)) dir.create(path2tabs, recursive = TRUE)
if (!dir.exists(path2figs)) dir.create(path2figs, recursive = TRUE)

# =============================================================================
# SECTION 3: FREQUENCY BAND COMPARISON (OB vs. AC)
# =============================================================================
# OB and AC files share the same structure. session_age is inferred from the
# filename.

# Helper: extract session_age integer from filename.
extract_age_from_filename = function(filepath) {
  dplyr::case_when(
    grepl("1mo",  filepath) ~ 1,
    grepl("6mo",  filepath) ~ 6,
    grepl("12mo", filepath) ~ 12,
    grepl("15mo", filepath) ~ 18,  # Remapped
    grepl("18mo", filepath) ~ 18,
    grepl("30mo", filepath) ~ 30,
    grepl("36mo", filepath) ~ 36,
    grepl("42mo", filepath) ~ 42,
    grepl("48mo", filepath) ~ 48)
}

# Load OB (Only-Burst) frequency band definitions: ranges estimated from burst-cycle peak frequencies
freqs_files    = list.files(path2freqs, pattern = "frequency_bands", full.names = TRUE)
freqs_files_ob = freqs_files[grepl("onlyburst_adapted", freqs_files)]

data_freq = lapply(freqs_files_ob, function(f) {
  # Label freq_group as "OB" to track band definition strategy downstream
  read.csv(f) |> mutate(session_age = extract_age_from_filename(f), freq_group = "OB")
}) |> bind_rows()

# Load AC (All-Cycles) frequency band definitions: ranges estimated from all-cycle peak frequencies
freqs_files_ac = freqs_files[grepl("allcycles_adapted", freqs_files)]

data_freq_allcycles = lapply(freqs_files_ac, function(f) {
  # Label freq_group as "AC" to track band definition strategy downstream
  read.csv(f) |> mutate(session_age = extract_age_from_filename(f), freq_group = "AC")
}) |> bind_rows()

# Compute summary statistics per visit × group
data_combined = bind_rows(data_freq, data_freq_allcycles) |>
  filter(session_age != -99) |>
  mutate(
    freq_group  = factor(freq_group, levels = c("OB", "AC")),
    session_age = factor(paste(session_age, 'mo.')),
    # Derive band center frequency and width from min/max edge columns
    mean_freq   = as.numeric((min_frequency_alpha + max_frequency_alpha) / 2),
    freq_range  = as.numeric(max_frequency_alpha - min_frequency_alpha)
  ) |>
  filter(inclusion == 1) |>
  dplyr::select(session_age, freq_group, mean_freq, freq_range) |>
  group_by(session_age, freq_group) |>
  summarise(
    mean_freq_avg  = mean(mean_freq,  na.rm = TRUE),
    freq_sd        = sd(mean_freq,    na.rm = TRUE),
    freq_range_avg = mean(freq_range, na.rm = TRUE),
    freq_range_sd  = sd(freq_range,   na.rm = TRUE),
    .groups        = "drop"
  ) |>
  mutate(
    mean_stats  = paste0(round(mean_freq_avg,  3), " (", round(freq_sd,       3), ")"),
    range_stats = paste0(round(freq_range_avg, 3), " (", round(freq_range_sd, 3), ")")
  ) |>
  dplyr::select(-c(mean_freq_avg, freq_sd, freq_range_avg, freq_range_sd)) |>
  pivot_wider(id_cols = session_age, names_from = freq_group,
              values_from = c(mean_stats, range_stats))

# Export Complementary Table
flextable(data_combined) |>
  set_header_labels(
    session_age    = "Visit",
    mean_stats_OB  = "Center Freq.\nOnly Burst",
    mean_stats_AC  = "Center Freq.\nAll Cycles",
    range_stats_OB = "Freq. Width\nOnly Burst",
    range_stats_AC = "Freq. Width\nAll Cycles"
  ) |>
  add_footer_lines("OB = Only-Burst alpha band (burst-cycle peak frequency); AC = All-Cycles alpha band. Values are Mean (SD) in Hz.") |>
  set_caption("Complementary Table: Alpha Frequency Band Parameters by Definition Strategy (OB vs. AC) and Visit") |>
  theme_booktabs() |> autofit() |>
  align(align = "left", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2tabs, "Table_SM1_Complementary_Table_AlphaLifespan_FrequencySelection_Ranges.docx"))


# =============================================================================
# SECTION 4: DATA LOADING & PREPARATION (for LMM analysis)
# =============================================================================

# --- EEG cleaning descriptives ---
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    # Remap session_age from acquisition labels to analysis ages (15→18, 40→42)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic/longitudinal age data ---
# age_months = exact chronological age at the visit (used as the GAMM smooth).
# SES covariates are z-scored; missing values are mean-imputed so no participant
# is lost from a model due to a single missing covariate.
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS, sujid %in% eeg_desc$sujid) |>
  mutate(
    # Remap session_age from acquisition labels to analysis ages (15→18, 40→42)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

desc_and_ages_wide = desc_and_ages |>
  dplyr::select(sujid, contains("mean"), GestationalAge_weeks) |>
  distinct() |>
  filter(sujid %in% eeg_desc$sujid) |>
  mutate(across(where(is.numeric), ~ scale(.x)[, 1]))

# --- Electrode map ---
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  # Remap region codes to full region labels (on-disk contract uses abbreviated codes)
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
    region == "T"  ~ "Temporal", region == "O" ~ "Occipital",
    TRUE         ~ "Central"
  )) |>
  filter(chinclu == 1)

# Combined descriptives frame (left anchor for all data merges)
descriptives = left_join(
  desc_and_ages |> dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
  eeg_desc
) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs))

# --- Lagged coherence — OB (Only-Burst) band definition ---
# Load LAcH data computed using OB (burst-based) alpha frequency band edges
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv")) |>
  merge(electrodes |> dplyr::select(label, chinclu, region), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  group_by(session_age, sujid, region, alpha_freqob) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop") |>
  rename(region = region) |>
  mutate(
    # Remap session_age from acquisition labels to analysis ages (15→18, 40→42)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Lagged coherence — AC (All-Cycles) band definition ---
# Load LAcH data computed using AC (all-cycles-based) alpha frequency band edges
hlcoh_data_allcycles = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long_allcycles.csv")) |>
  merge(electrodes |> dplyr::select(label, chinclu, region), by.x = "ch", by.y = "label") |>
  filter(chinclu == 1, epochs >= epochs_threshold) |>
  group_by(session_age, sujid, region, alpha_freqob) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop") |>
  rename(region = region) |>
  mutate(
    # Remap session_age from acquisition labels to analysis ages (15→18, 40→42)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# Combine OB and AC into one data frame for modelling and plotting
# Add freq_group label to track which band definition (OB vs. AC) was used for LAcH computation
hlcoh_data_combined = bind_rows(
  hlcoh_data           |> mutate(freq_group = "OB"),
  hlcoh_data_allcycles |> mutate(freq_group = "AC")
) |>
  left_join(
    descriptives |> dplyr::select(sujid, session_age, session_age, prop_epochs,
                                   Cohort, GestationalAge_weeks),
    by = c("sujid", "session_age")
  )

# =============================================================================
# SECTION 5: LAcH DESCRIPTIVES TABLE (OB vs. AC by Visit)
# =============================================================================

hlcoh_data_combined |>
  filter(session_age != -99)|>
  group_by(session_age, freq_group) |>
  summarise(alpha_LAcH_avg = mean(alpha_LAcH, na.rm = TRUE),
            alpha_LAcH_sd  = sd(alpha_LAcH,   na.rm = TRUE), .groups = "drop") |>
  mutate(LAcH_stats = paste0(round(alpha_LAcH_avg, 3), " (", round(alpha_LAcH_sd, 3), ")")) |>
  dplyr::select(-c(alpha_LAcH_avg, alpha_LAcH_sd)) |>
  pivot_wider(id_cols = session_age, names_from = freq_group, values_from = LAcH_stats) |>
  # Format session_age as factor with "mo." suffix for display
  mutate(session_age = factor(paste(session_age, "mo.")))|>
  flextable() |>
  set_header_labels(session_age = "Visit", OB = "Lifespan\nOnly Burst",
                    AC = "Lifespan\nAll Cycles") |>
  set_caption("Descriptive Summary: Alpha Lifespan (LAcH) by Frequency Band Definition and Visit") |>
  theme_booktabs() |> autofit() |> 
  set_table_properties(layout = "fixed") |>
  align(align = "left", part = "all") |>
  align(align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_docx(path = file.path(path2tabs, "Table_SM2_AlphaLifespan_Differences_FrequencySelection.docx"))


# =============================================================================
# SECTION 6: LMM — Does OB vs. AC Frequency Definition Affect LAcH and does this differ by visit?
# =============================================================================
# Model: alpha_LAcH ~ freq_group × session_age (i.e., visit) + prop_epochs + Cohort +
#                     GestationalAge_weeks + region + (1|sujid)
# Bootstrapped CIs for fixed-effect coefficients (1000 iterations).

fixed_effects = "alpha_LAcH ~ freq_group * session_age + prop_epochs + Cohort + GestationalAge_weeks + region + (1|sujid)"

m                 = lmer(as.formula(fixed_effects), data = hlcoh_data_combined, REML = TRUE)
performance_model = performance(m)
model_summary     = summary(m)
coefs             = as.data.frame(model_summary$coefficients)

# Extract relevant coefficient row names
burst_row = grep("^freq",        rownames(coefs), value = TRUE)[1]
age_row   = grep("^session_age",  rownames(coefs), value = TRUE)[1]
inter_row = grep(":",            rownames(coefs), value = TRUE)[1]

cat(sprintf("Running bootstrapping (%d simulations)...\n", n_bootstraps))
boot_ci = confint(m, parm = "beta_", method = "boot", nsim = n_bootstraps, oldNames = FALSE)

# Store main effects results: session_age (visit) and freq_group (OB vs. AC) fixed effects
# NOTE: Column names use standardized suffixes (_Estimate, _SE, _CI_Lower, _CI_Upper, _P_Value)
res_main = data.frame(
  metric         = "alpha_lifespan",
  Visit_Estimate   = coefs[age_row,   "Estimate"],
  Visit_SE         = coefs[age_row,   "Std. Error"],
  Visit_CI_Lower   = boot_ci[age_row,  "2.5 %"],
  Visit_CI_Upper   = boot_ci[age_row,  "97.5 %"],
  Visit_t          = coefs[age_row,   "t value"],
  Visit_df         = coefs[age_row,   "df"],
  Visit_P_Value    = coefs[age_row,   "Pr(>|t|)"],
  Burst_Estimate = coefs[burst_row, "Estimate"],
  Burst_SE       = coefs[burst_row, "Std. Error"],
  Burst_CI_Lower = boot_ci[burst_row, "2.5 %"],
  Burst_CI_Upper = boot_ci[burst_row, "97.5 %"],
  Burst_t        = coefs[burst_row, "t value"],
  Burst_df       = coefs[burst_row, "df"],
  Burst_P_Value  = coefs[burst_row, "Pr(>|t|)"]
)

inter_pval = coefs[inter_row, "Pr(>|t|)"]

res_inter = data.frame(
  metric          = "alpha_lifespan",
  Inter_Estimate  = coefs[inter_row, "Estimate"],
  Inter_SE        = coefs[inter_row, "Std. Error"],
  Inter_CI_Lower  = boot_ci[inter_row, "2.5 %"],
  Inter_CI_Upper  = boot_ci[inter_row, "97.5 %"],
  Inter_t         = coefs[inter_row, "t value"],
  Inter_df        = coefs[inter_row, "df"],
  Inter_P_Value   = inter_pval,
  Conditional_R2  = performance_model$R2_conditional,
  Marginal_R2     = performance_model$R2_marginal
)

# --- Numbered supplementary table SM2 ---
left_join(res_main, res_inter, by = "metric") |>
  nn_supp_table(id = "SM2", cols = c(
    metric         = "Metric",
    Burst_Estimate = "Beta (OB vs AC)",
    Burst_SE       = "s.e.",
    Burst_CI_Lower = "95% CI lower",
    Burst_CI_Upper = "95% CI upper",
    Burst_t        = "t",
    Burst_df       = "df",
    Burst_P_Value  = "P",
    Visit_Estimate = "Beta (visit)",
    Visit_SE       = "s.e. (visit)",
    Visit_t        = "t (visit)",
    Visit_df       = "df (visit)",
    Visit_P_Value  = "P (visit)",
    Inter_Estimate = "Beta (interaction)",
    Inter_SE       = "s.e. (interaction)",
    Inter_t        = "t (interaction)",
    Inter_df       = "df (interaction)",
    Inter_P_Value  = "P (interaction)",
    Marginal_R2    = "Marginal R2",
    Conditional_R2 = "Conditional R2"),
    note = paste(
      "Linear mixed model comparing alpha lifespan computed over an",
      "oscillation-based frequency range (OB) versus an all-cycles range (AC):",
      "Lifespan ~ Range definition x Visit + covariates + (1|child). Degrees of",
      "freedom from Satterthwaite's approximation. 95% CIs from parametric",
      "bootstrapping. OB, oscillation-based; AC, all cycles."))

left_join(res_main, res_inter, by = "metric") |>
  flextable() |>
  set_header_labels(
    metric = "Metric", Visit_SE = "Visit SE", Burst_SE = "OB vs. AC SE", Inter_SE = "Interaction SE",
    Visit_CI_Lower = "Visit CI Lower", Visit_CI_Upper = "Visit CI Upper",
    Burst_CI_Lower = "OB vs. AC CI Lower", Burst_CI_Upper = "OB vs. AC CI Upper"
  ) |>
  colformat_double(digits = 3) |>
  colformat_double(j = c("Visit_P_Value", "Burst_P_Value", "Inter_P_Value"), digits = 4) |>
  add_footer_lines("Bootstrapped 95% CIs (1000 iterations). OB vs. AC: OB = Only-Burst; AC = All-Cycles.") |>
  theme_booktabs() |> autofit() |>
  align(align = "left", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  save_as_html(path = file.path(path2tabs, "SupplementaryMaterial_Table_AlphaLifespan_FrequencySelection.html"))


# =============================================================================
# SECTION 7: PER-VISIT FOLLOW-UP LMMs (if interaction is significant)
# =============================================================================
# If freq_group × session_age interaction reaches p < .05, run a separate LMM
# per visit to identify which ages show an OB vs. AC difference.

emmeans_list = list()
if (!is.na(inter_pval) && inter_pval < 0.05) {

  cat("Interaction significant (p =", round(inter_pval, 4), ") — running per-visit follow-up LMMs.\n")

  sub_fixed_effects = "alpha_LAcH ~ freq_group + prop_epochs + Cohort + GestationalAge_weeks + (1|sujid)"

  for (age in unique(hlcoh_data_combined$session_age)) {
    sub_data    = hlcoh_data_combined |> filter(session_age == age)
    burst_descs = sub_data |>
      group_by(freq_group) |>
      summarise(m = mean(alpha_LAcH, na.rm = TRUE), s = sd(alpha_LAcH, na.rm = TRUE),
                .groups = "drop")

    sub_m             = lmer(as.formula(sub_fixed_effects), data = sub_data, REML = TRUE)
    performance_sub   = performance(sub_m)
    sub_coefs         = as.data.frame(summary(sub_m)$coefficients)
    burst_row_sub     = grep("^freq", rownames(sub_coefs), value = TRUE)[1]

    cat(sprintf("Running bootstrapping for visit %d months...\n", age))
    boot_ci_sub = confint(sub_m, parm = "beta_", method = "boot",
                          nsim = n_bootstraps, oldNames = FALSE)

    # Format OB and AC descriptive stats (mean ± SD) for output table
    desc_ob  = paste0(round(burst_descs$m[burst_descs$freq_group == "OB"], 3), " (",
                       round(burst_descs$s[burst_descs$freq_group == "OB"], 3), ")")
    desc_ac  = paste0(round(burst_descs$m[burst_descs$freq_group == "AC"], 3), " (",
                       round(burst_descs$s[burst_descs$freq_group == "AC"], 3), ")")

    emmeans_list[[as.character(age)]] = data.frame(
      metric         = "Lifespan",
      visit          = age,
      OB_MSD         = desc_ob,
      AC_MSD         = desc_ac,
      Estimate       = sub_coefs[burst_row_sub, "Estimate"],
      SE             = sub_coefs[burst_row_sub, "Std. Error"],
      t              = sub_coefs[burst_row_sub, "t value"],
      df             = sub_coefs[burst_row_sub, "df"],
      CI_Lower       = boot_ci_sub[burst_row_sub, "2.5 %"],
      CI_Upper       = boot_ci_sub[burst_row_sub, "97.5 %"],
      P_Value        = sub_coefs[burst_row_sub, "Pr(>|t|)"],
      Conditional_R2 = performance_sub$R2_conditional,
      Marginal_R2    = performance_sub$R2_marginal
    )
  }

  emmeans_results = bind_rows(emmeans_list) |>
    # Apply FDR correction across all age-groups, map p-values to significance stars
    mutate(p_fdr = p.adjust(P_Value, method = "fdr"),
           stars = case_when(p_fdr < 0.001 ~ "***", p_fdr < 0.01 ~ "**",
                              p_fdr < 0.05 ~ "*",   TRUE ~ "NS")) |>
    rename(session_age = visit)

  # --- Numbered supplementary table SM3 ---
  emmeans_results |>
    nn_supp_table(id = "SM3", cols = c(
      metric         = "Metric",
      session_age    = "Visit (months)",
      OB_MSD         = "OB, mean (s.d.)",
      AC_MSD         = "AC, mean (s.d.)",
      Estimate       = "Beta (OB vs AC)",
      SE             = "s.e.",
      t              = "t",
      df             = "df",
      CI_Lower       = "95% CI lower",
      CI_Upper       = "95% CI upper",
      P_Value        = "P",
      p_fdr          = "P (FDR)",
      Marginal_R2    = "Marginal R2",
      Conditional_R2 = "Conditional R2"),
      note = paste(
        "Visit-stratified linear mixed models comparing alpha lifespan computed over",
        "an oscillation-based frequency range (OB) versus an all-cycles range (AC).",
        "Degrees of freedom from Satterthwaite's approximation. 95% CIs from",
        "parametric bootstrapping. P values FDR-corrected (Benjamini-Hochberg)",
        "across visits. OB, oscillation-based; AC, all cycles."))

  emmeans_results |>
    rename(Visit = session_age) |>
    flextable() |>
    set_header_labels(
      metric = "Metric", visit = "Visit", OB_MSD = "OB Mean (SD)", AC_MSD = "AC Mean (SD)",
      Estimate = "\u03B2 (OB vs. AC)", SE = "SE", CI_Lower = "95% CI Lower",
      CI_Upper = "95% CI Upper", P_Value = "p", p_fdr = "p (FDR)"
    ) |>
    set_caption("Per-Age LMM Follow-Up: OB vs. AC Difference in Alpha Lifespan (LAcH)") |>
    colformat_double(digits = 3) |>
    colformat_double(j = c("P_Value", "p_fdr"), digits = 4) |>
    bold(j = "p_fdr", i = ~ p_fdr < 0.05, part = "body") |>
    add_footer_lines("FDR correction (Benjamini-Hochberg) applied across all ages.") |>
    theme_booktabs() |> autofit() |>
    align(align = "left", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    save_as_html(path = file.path(path2tabs, "SupplementaryMethods_Table_AlphaLifespan_FrequencySelection_by_Visit.html"))

} else {
  cat("Interaction not significant (p =", round(inter_pval, 4), ") — per-visit follow-up skipped.\n")
  emmeans_results = data.frame(session_age = integer(), stars = character())
}


# =============================================================================
# SECTION 8: BAR CHART — Alpha Lifespan (OB vs. AC) by Visit
# =============================================================================

hlcoh_plot_data = hlcoh_data_combined |>
  group_by(sujid, session_age, freq_group) |>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(freq_group))|>
  # Establish factor level order (OB first, AC second) for consistent plotting
  mutate(freq_group = factor(freq_group, levels = c('OB', 'AC')))

p_fig = ggplot(hlcoh_plot_data, aes(x = factor(session_age), y = alpha_LAcH, fill = freq_group)) +
  stat_summary(geom = "bar", fun = "mean", color = "black",
               position = position_dodge(0.8), alpha = 0.66, linewidth = 1) +
  geom_point(aes(color = freq_group),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
             alpha = 0.5) +
  stat_summary(geom = "errorbar", fun.data = median_se, color = "black",
               position = position_dodge(0.8), width = 0.2, linewidth = 0.8) +
  stat_summary(geom = "point", fun = "median", color = "black",
               position = position_dodge(0.8), size = 2) +
  stat_summary(geom = "line", fun = "median", aes(group = freq_group), color = "black",
               position = position_dodge(0.8), linewidth = 1) +
  THEME_BASE + THEME_TEXT +
  scale_fill_manual(values  = COLORS_MAIN) +
  scale_color_manual(values = COLORS_MAIN) +
  labs(x = "Visit", y = "Lifespan",
       color = "Frequency Band", fill = "Frequency Band")

# Overlay significance stars if follow-up LMMs were run
if (nrow(emmeans_results) > 0) {
  max_y = max(hlcoh_plot_data$alpha_LAcH, na.rm = TRUE)
  p_fig = p_fig +
    geom_text(data = emmeans_results,
              aes(x = factor(session_age), y = max_y - 0.1, label = stars),
              inherit.aes = FALSE, size = 5, vjust = 0)
}

ggsave(p_fig,
       filename = file.path(path2figs, "Fig_SM4_AlphaLifespan_FrequencySelection.jpeg"),
       width = 5, height = 4.5, dpi = 300)

