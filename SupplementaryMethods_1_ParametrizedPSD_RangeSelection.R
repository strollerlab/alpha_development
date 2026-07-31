# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 terms
#     r2value                 -> Model fit (R²)
#     mae                     -> Mean absolute error (MAE)
#     slope                   -> Aperiodic slope
#     offset                  -> Aperiodic offset
#     maxfreq                 -> Maximum fitting frequency cutoff (15, 30, 45 Hz)
#     epochs/epochs      -> Number of clean epochs per electrode
# -----------------------------------------------------------------------------
# Script: SupplementaryMethods_1_ParametrizedPSD_RangeSelection.R
# Purpose: Evaluates Specparam parametrization quality across different maximum
#          fitting frequency cutoffs (15 Hz, 30 Hz, 45 Hz). For each study visit:
#            (1) Computes per-subject R², MAE, and channel retention at each cutoff
#            (2) Exports an overall descriptives table (gtsummary)
#            (3) Runs Friedman test + pairwise Wilcoxon (FDR) across cutoffs
#            (4) Plots retained participants, R²/MAE boxplots, and PSD residual curves
#
#          The script helps justify the selection of 15 Hz 
# =============================================================================
# Inputs:
#   - Specparam parameter files (aperosc_parameters_*.csv) and PSD error files (aperosc_psds_*.csv) for each visit and max-frequency level. See README taxonomy for details.
# Outputs (all to path2save):
# - Fig_SM1_ParametrizedPSD_Range_Selection.jpeg
# - SupplementaryMethods_Complementary_Table_ParametrizedPSD_Range_Comparison_Analysis_**visit**_summary.html
# - SupplementaryMethods_Complementary_Table_ParametrizedPSD_Range_Comparison_Analysis_**visit**_pairwise_followup.html
# Dependencies: 00_Setup_PackageInstallation.R, 01_Utils_ProcFunctions.R
# =============================================================================

# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started
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

source(file.path(path2code, "00_Setup_PackageInstallation.R"))
source(file.path(path2code, "01_Utils_ProcFunctions.R"))  # Provides add_stat_pairwise_wilcox() and
                                     # add_stat_overall_friedman()

# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================

r2_thresh        = R2_THRESH         # Minimum Specparam R² threshold (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10; 'mae' column name)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)
final_prefit     = "no_prefit"       # Specparam fitting mode used in this dataset

# Analysis-level session visits (already remapped; no 15 or 40 in this list)
visits = c(1, 6, 12, 18, 30, 36, 42, 48)

# Colour palette for the three max-frequency levels (light → dark = 15 → 30 → 45 Hz)
my_gradual_palette = c("#B3E5FC", "#9FA8DA", "#7E57C2")

# Apply global flextable defaults (font size, padding, spacing)
flextable::set_flextable_defaults(font.size = 10, padding = 3, line_spacing = 1)

# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================

# path2sets comes from config_paths.R   # <<< CHECK: set this to the folder holding the per-visit input CSVs (see config_paths.R: path2sets)
path2data = path2sets    # EDIT — root data folder   # = Data/; this script prepends "" to each file name
# path2root comes from config_paths.R

# Supplementary Methods outputs (range-selection analysis). See README taxonomy.
path2tabs = file.path(path2root, "SupplementaryInformation",  "Methods", "Tables")  # range-comparison tables
path2figs = file.path(path2root, "SupplementaryInformation", "Methods", "Figures")  # Fig SM1

for (p in c(path2tabs, path2figs)) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# =============================================================================
# SECTION 3: DATA LOADING
# =============================================================================

# EEG cleaning descriptives (defines the denominator for retention calculations)
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic data ---
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
  mutate(across(where(is.numeric), ~ scale(.x)[, 1]))

# --- Electrode map ---
electrodes = read_csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  mutate(region = case_when(
    region == "Fr" ~ "Frontal", region == "P" ~ "Parietal",
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


# =============================================================================
# SECTION 4: PER-AGE PARAMETRIZATION LOOP
# =============================================================================
# For each study visit:
#   (1) Load all Specparam parameter and PSD files
#   (2) Apply quality filters and compute inclusion flags
#   (3) Build summary tables (overall + pairwise Friedman/Wilcoxon comparison)
#   (4) Plot participant retention, R²/MAE distributions, and PSD error curves

data_plot_fit = list()  # R²/MAE dual-axis boxplots (one per visit → combined at end)
data_psd_agg  = list()  # PSD residual error plots (one per visit → combined at end)

for (visit in visits) {
  visit_path = if      (visit == 18) c(15, 18)
             else if (visit == 42) 40
             else                visit

  # Load Specparam parameter files (all max-frequency levels combined)
  if (length(visit_path) > 1) {
    sets      = c(
      list.files(paste0(path2sets, "Age", visit_path[1]), pattern = "aperosc_parameters_*", recursive = TRUE, full.names = TRUE),
      list.files(paste0(path2sets, "Age", visit_path[2]), pattern = "aperosc_parameters_*", recursive = TRUE, full.names = TRUE)
    )
    sets_psds = c(
      list.files(paste0(path2sets, "Age", visit_path[1]), pattern = "aperosc_psds_*",       recursive = TRUE, full.names = TRUE),
      list.files(paste0(path2sets, "Age", visit_path[2]), pattern = "aperosc_psds_*",       recursive = TRUE, full.names = TRUE)
    )
  } else {
    sets      = list.files(paste0(path2sets, "Age", visit_path), pattern = "aperosc_parameters_*", recursive = TRUE, full.names = TRUE)
    sets_psds = list.files(paste0(path2sets, "Age", visit_path), pattern = "aperosc_psds_*",       recursive = TRUE, full.names = TRUE)
  }

  # Read and filter to analysis channels and eyes-closed resting state
  data     = lapply(sets,      read_csv) |> bind_rows() |>
    filter(ch %in% electrodes$label, is.na(block) | block == "ECrs")
  data_psd = lapply(sets_psds, read_csv) |> bind_rows() |>
    filter(ch %in% electrodes$label, is.na(block) | block == "ECrs")

  aper_voi = c("sujid", "ch", "epochs",
               "r2value", "mae", "slope", "maxfreq")

  # Standardize session age: remap raw labels (15, 40) to labels (18, 42) for consistency.
  data = data %>%
      dplyr::select(all_of(aper_voi)) |>
    mutate(session_age = visit) |>
    mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
    filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
           sujid %in% descriptives$sujid) |>
    left_join(descriptives, by = c("sujid", "session_age")) |>
    drop_na(slope) |>
    drop_na(any_of(c("age_months", "prop_epochs", "GestationalAge_weeks", "r2value"))) |>
    # Pre-filter to quality-passing electrodes before computing mean error curves
    filter(r2value >= r2_thresh, mae <= mae_thresh) |>
    left_join(electrodes, by = c("ch" = "label")) 

  # Compute quality-based inclusion flags for sensitivity sweep across maxfreq levels:
  # goodch           = count of channels passing R² + MAE + trial thresholds per subject-maxfreq
  # ch_inclusion     = binary flag: 1 if channel passes all three quality criteria
  # global_inclusion = binary flag: 1 if subject has ≥1 passing channel at this maxfreq level
  inclusion = data |>
    mutate(
      goodch        = sum(r2value >= r2_thresh & mae <= mae_thresh & epochs >= epochs_threshold),
      .by           = c(sujid, maxfreq)
    ) |>
    mutate(
      ch_inclusion  = if_else(r2value >= r2_thresh & mae <= mae_thresh & epochs >= epochs_threshold, 1, 0),
      .by           = c(sujid, maxfreq, ch)
    ) |>
    mutate(
      global_inclusion = if_else(goodch > 0, 1, 0),
      .by              = c(sujid, maxfreq)
    ) |>
    dplyr::select(sujid, maxfreq, ch, goodch, global_inclusion, ch_inclusion)

  data = left_join(data, inclusion, by = c("sujid", "maxfreq", "ch"))

  # --- Table: Overall summary ---
  data |>
    filter(global_inclusion == 1, ch_inclusion == 1, goodch >= CH_THRESHOLD) |>
    group_by(sujid, maxfreq) |>
    summarise(
      r2    = mean(r2value, na.rm = TRUE),
      mae   = mean(mae,     na.rm = TRUE),
      perch = mean(goodch) / 60,  
      .groups = "drop"
    ) |>
    tbl_summary(
      by        = maxfreq,
      include   = c(r2, mae, perch),
      missing   = "no",
      label     = list(r2 ~ "R\u00b2", mae ~ "MAE", perch ~ "Prop. of Electrodes\nRetained"),
      statistic = list(all_continuous() ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]")
    ) |>
    modify_header(label ~ "**Metric**") |>
    as_flex_table() |>
    flextable::autofit() |>
    flextable::colformat_double(digits = 3) |>
    theme_booktabs() |> autofit() |>
    align(align = "left", part = "all") |>
    align(j = 2:4, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
    flextable::save_as_html(path = file.path(path2tabs, paste0("SupplementaryMethods_Complementary_Table_ParametrizedPSD_Range_Comparison_Analysis_", visit, "_summary.html")))

  # Apply same quality filters to PSD data, then average error across channels per subject-freq-maxfreq
  data_psd = data_psd |>
    rename(ch = ch) |>
    left_join(inclusion, by = c("sujid", "maxfreq", "ch")) |>
    filter(global_inclusion == 1, ch_inclusion == 1,
           r2value >= r2_thresh, mae <= mae_thresh, epochs >= epochs_threshold, goodch >= CH_THRESHOLD) |>
    group_by(sujid, maxfreq, freq) |>
    summarise(error = mean(error, na.rm = TRUE), .groups = "drop")

  # Pairwise comparison: Friedman omnibus test + post-hoc Wilcoxon (FDR).
  # Restricts to fully-paired subjects with data at all three maxfreq levels (15, 30, 45 Hz).
  data |>
    filter(global_inclusion == 1, ch_inclusion == 1, r2value >= r2_thresh, goodch >= CH_THRESHOLD, mae <= mae_thresh) |>
    group_by(sujid, maxfreq) |>
    summarise(r2 = mean(r2value, na.rm = TRUE), mae = mean(mae, na.rm = TRUE),
              perch = mean(goodch) / 60, .groups = "drop") |>
    mutate(n_included = n(), .by = sujid) |>
    filter(n_included == length(unique(data$maxfreq[!is.na(data$maxfreq)]))) |>
    tbl_summary(
      by        = maxfreq,
      include   = c(r2, mae, perch),
      missing   = "no",
      label     = list(r2 ~ "R\u00b2", mae ~ "MAE", perch ~ "Prop. of Electrodes\nRetained"),
      statistic = list(all_continuous() ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]")
    ) |>
    add_stat(everything() ~ function(data, variable, by, ...) {
      add_stat_overall_friedman(data, variable, by, id_var = "sujid")
    }) |>
    add_stat(everything() ~ function(data, variable, by, ...) {
      add_stat_pairwise_wilcox(data, variable, by, id_var = "sujid", method = "fdr")
    }) |>
    modify_source_note(paste0(
      "Friedman test for repeated measures; chi-squared with df = k - 1 = ",
      dplyr::n_distinct(data$maxfreq) - 1,
      " (k = ", dplyr::n_distinct(data$maxfreq), " candidate upper frequency limits). ",
      "Pairwise comparisons are Wilcoxon signed-rank tests, FDR-corrected ",
      "(Benjamini-Hochberg); the V statistic is a rank sum and has no degrees of ",
      "freedom, so it is reported with n. Only subjects present at all frequency ",
      "levels are included.")) |>
    modify_header(label ~ "**Metric**") |>
    as_flex_table() |> flextable::autofit() |>
    flextable::autofit() |>
    flextable::colformat_double(digits = 3) |>
    theme_booktabs() |> autofit() |>
    align(align = "left", part = "all") |>
    align(j = 2:6, align = "center", part = "all") |>
    bold(part = "header")|>
    padding(padding = 1, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 9, part = "all") |>
    flextable::save_as_html(path = file.path(path2tabs, paste0("SupplementaryMethods_Complementary_Table_ParametrizedPSD_Range_Comparison_Analysis_", visit, "_pairwise_followup.html")))


  # Prepare per-visit R² and MAE summary for side-by-side boxplots (dual-axis visualization).
  summary_data_2 = data |>
    filter(global_inclusion == 1, ch_inclusion == 1, r2value >= r2_thresh, mae <= mae_thresh) |>
    group_by(sujid, maxfreq, region) |>
    summarise(r2 = mean(r2value, na.rm = TRUE), mae = mean(mae, na.rm = TRUE), .groups = "drop") |>
    group_by(sujid, maxfreq) |>
    summarise(r2 = mean(r2, na.rm = TRUE), mae = mean(mae, na.rm = TRUE), .groups = "drop")

  data_plot_fit[[paste0(visit, "mo")]] = summary_data_2 |>
    group_by(maxfreq) |>
    mutate(session_age = visit)

  # Prepare per-visit PSD error curves: aggregate residual error across subjects for each freq × maxfreq.
  # Compute mean error and SD envelope for visualization across frequency range.
  data_psd_agg[[paste0(visit, "mo")]] = data_psd |>
    group_by(maxfreq, freq) |>
    summarise(m_error  = mean(error, na.rm = TRUE),
              sd_error = sd(error,   na.rm = TRUE), .groups = "drop") |>
    mutate(sesion_age = visit)

  # Free per-visit objects before the next iteration
  rm(data_psd, data, inclusion)
}


scale_factor = 1
offset_val   = 0.8  # Shifts MAE boxes into the R² display range

fit_data = bind_rows(data_plot_fit) |>
  # Recode session_age to categorical labels for faceting (e.g., '1 mo.', '6 mo.', ..., '48 mo.')
  mutate(session_age = factor(session_age, labels = paste(visits, 'mo.'), levels = visits))

fit_plot = ggplot(fit_data, aes(x = factor(maxfreq))) +
  geom_boxplot(aes(y = r2),
               width = 0.35, position = position_nudge(x = -0.2),
               fill = "steelblue", color = "black", outlier.color = "steelblue", alpha = 0.5) +
  geom_boxplot(aes(y = (mae * scale_factor) + offset_val),
               width = 0.35, position = position_nudge(x = 0.2),
               fill = "indianred", color = "black", outlier.color = "indianred", alpha = 0.5) +
  
  # Removed the dot (.) from the facet_wrap formula
  facet_wrap(~session_age, nrow = 2) + 
  scale_y_continuous(
    name     = "R\u00b2",
    sec.axis = sec_axis(~ (. - offset_val) / scale_factor, name = "MAE") 
  ) +
  
  coord_cartesian(ylim = c(0.8, 1)) + 
  
  labs(x = "Max. Frequency (Hz)") +
  THEME_BASE + THEME_TEXT +
  
  #  Theme() calls for cleaner code
  theme(
    axis.title.y.left  = element_text(color = "steelblue", face = "bold", size = 12, margin = margin(r = 10)),
    axis.title.y.right = element_text(color = "indianred", face = "bold", size = 12, margin = margin(l = 10)),
    axis.line.y.right  = element_line(color = "black"),
    plot.margin        = margin(1.75, 0.1, 0.1, 0.1, "cm")
  )


psd_data = bind_rows(data_psd_agg)
psd_data = psd_data |>
  mutate(session_age = factor(sesion_age, labels = paste(visits, 'mo.'), levels = visits))

psds_plot = ggplot(psd_data,
                   aes(x = freq, y = m_error, color = factor(maxfreq),
                       fill = after_scale(color))) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = m_error + sd_error, ymax = m_error - sd_error), alpha = 0.20, linewidth = 0) +
  facet_wrap(.~session_age, nrow = 2) + 
  geom_hline(yintercept =  mae_thresh, linetype = "dashed", color = "darkred",  linewidth = 1.0) +
  geom_hline(yintercept = -mae_thresh, linetype = "dashed", color = "darkred",  linewidth = 1.0) +
  geom_hline(yintercept =  0,          linetype = "dashed", color = "black",    linewidth = 1.5) +
  THEME_BASE + THEME_TEXT +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Frequency (Hz)", y = "Absolute Error", color = "Max. Frequency (Hz)")

# =============================================================================
# SECTION 5: COMBINED MULTI-AGE FIGURES
# =============================================================================

combined_plot = ggpubr::ggarrange(
  fit_plot, psds_plot,
  ncol = 1, nrow = 2,
  labels        = c("a", "b"),
  hjust         = -0.1, vjust = 1.5)
ggsave(combined_plot,
       filename = file.path(path2figs, "Fig_SM1_ParametrizedPSD_Range_Selection.jpeg"),
       width = 30, height = 30, dpi = 300, unit = 'cm')


