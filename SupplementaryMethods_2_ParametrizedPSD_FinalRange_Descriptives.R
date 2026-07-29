# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term
#   Column identifiers below are an on-disk CSV / cross-script DATA CONTRACT and
#   are intentionally NOT renamed (would break the pipeline + OSF data). Display
#   labels and this map carry the manuscript terminology.
#     r2value                 -> Model fit (R²)
#     mae                     -> Mean squared error (MAE)
#     slope                   -> Aperiodic slope
#     offset                  -> Aperiodic offset
#     epochs                  -> Number of clean epochs per electrode
#     goodch                  -> Count of channels passing all quality thresholds
#     alpha_freq/alpha_ampl/alpha_osc/alpha_peak -> Oscillatory alpha band parameters
# -----------------------------------------------------------------------------
# Script: SupplementaryMethods_2_ParametrizedPSD_FinalRange_Descriptives.R
# Purpose: Characterisation of the Specparam model fit for the
#          SELECTED 15 Hz parametrization used in all
#          main analyses. Complements script SupplementaryMethods_1.
# =============================================================================
# Inputs:
#   - Data/Aperiodic_Oscillatory_ByCycle_Long.csv
#   - Data/psds_aperosc_long_region_*_15Hz*.csv
#   - Data/CrossVisit_EEG_CleaningDescriptives.csv
#   - Data/Sociodemographic_Descriptives_Long_Updated.csv
#   - Scripts/electrodes.csv
# Outputs (to path2figs, path2tabs, path2data)
#  - Fig_SM2_ParametrizedPSD_QualityMetrics_2to15Hz_range_extended.jpeg
#  - Table_S4_ParametrizedPSD_QualityMetrics_Descriptive.docx
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

r2_thresh        = R2_THRESH         # Minimum Specparam R² (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10)
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)

# Analysis ages (after session-label remapping)
ages               = c(1, 6, 12, 18, 30, 36, 42, 48)

# Age-indexed colour palette (oldest → youngest, sequential purples)
age_palette        = COLORS_GRADIENT
names(age_palette) = as.character(ages)

# Apply global flextable defaults
flextable::set_flextable_defaults(font.size = 10, padding = 3, line_spacing = 1)

# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================

path2data = path2sets    # In the anonymized data, path2sets = path2data as long as everything is in the same folder; otherwise, modify.
path2psd  = path2sets    # In the anonymized data, path2sets = path2psd as long as everything is in the same folder; otherwise, modify.

# path2root comes from config_paths.R

# Supplementary Methods outputs (final-range quality metrics). See README taxonomy.
path2tabs = file.path(path2root, "SupplementaryInformation",  "Methods", "Tables")  # range-comparison tables
path2figs = file.path(path2root, "SupplementaryInformation", "Methods", "Figures")  # Fig SM1

dir.create(path2figs, recursive = TRUE, showWarnings = FALSE)
dir.create(path2tabs, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 3: DATA LOADING & PREPARATION
# =============================================================================

# --- PSD files (15 Hz fitting, region-averaged) ---
# These were generated by script 02 and contain per-frequency oscillatory and
# aperiodic components for each subject × region × session age.
psd_files = list.files(path = path2psd, pattern = "psds_aperosc",
                        recursive = TRUE, full.names = TRUE)
psd_files = psd_files[!grepl("45Hz", psd_files) & grepl("region", psd_files)]

# --- EEG cleaning descriptives ---
eeg_desc = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv")) |>
  filter(clean_epochs_rest >= epochs_threshold, !sujid %in% EXCLUDED_SUBJECTS) |>
  mutate(inclusion_lmm = 1) |>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata, inclusion_lmm) |>
  mutate(
    prop_epochs  = scale(prop_epochs)[, 1],
    # Standardize session age labels: remap raw labels (15, 40) to labels (18, 42)
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  )

# --- Sociodemographic data ---
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv")) |>
  filter(dev_filter == 1, !sujid %in% EXCLUDED_SUBJECTS, sujid %in% eeg_desc$sujid) |>
  mutate(
    # Standardize session age labels to match EEG data nomenclature
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
  # Recode region abbreviations to full anatomical region names for display
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

# --- Aperiodic / oscillatory data ---
# Note: 'mae' is the raw Specparam MAE column name; displayed as 'mae' in tables/figures.
aper_voi = c("sujid", "session_age", "ch", "region", "chinclu", "epochs",
              "r2value", "mae", "goodch", "offset", "slope",
              "alpha_freq", "alpha_ampl", "alpha_osc", "alpha_peak",
              "inclusion_final_dummy")

aperiodic_raw = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv")) |>
  dplyr::select(all_of(aper_voi)) |>
  # Standardize session age (remap 15 -> 18 for consistency across data sources)
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% descriptives$sujid, goodch >= CH_THRESHOLD) |>
  left_join(descriptives, by = c("sujid", "session_age")) |>
  drop_na(slope) |>
  drop_na(any_of(c("age_months", "prop_epochs", "GestationalAge_weeks", "r2value")))

# Per-subject × visit good-channel count (used as a join key for PSD filtering)
goodchs = aperiodic_raw |>
  group_by(sujid, session_age) |>
  summarise(goodch = mean(goodch, na.rm = TRUE), .groups = "drop")

# --- PSD data ---
# Load all 15 Hz region-averaged PSD files; apply the same quality filters as
# the aperiodic dataset, then collapse to one error value per subject × visit × freq.
psd_data = lapply(psd_files, read_csv) |> bind_rows() |>
  mutate(
    session_age = if_else(session_age == 15, 18, session_age),
    session_age = if_else(session_age == 40, 42, session_age)
  ) |>
  left_join(goodchs, by = c("sujid", "session_age")) |>
  filter(r2value >= r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% eeg_desc$sujid, goodch >= CH_THRESHOLD) |>
  group_by(sujid, session_age, freq) |>
  summarise(abserror = mean(error, na.rm = TRUE), .groups = "drop")


# =============================================================================
# SECTION 4: RETAINED ELECTRODE BAR PLOT
# =============================================================================
# Stacked bar showing average number of retained electrodes per brain region
# at each visit. Each stack layer is one region; text labels show n per region.

electrodes_labels = aperiodic_raw |>
  filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% eeg_desc$sujid, chinclu == 1) |>
  group_by(sujid, session_age, region) |>
  summarise(n_elec = n(), .groups = "drop") |>
  group_by(session_age, region) |>
  summarise(n_elec_mean = mean(n_elec), .groups = "drop") |>
  # Format text labels for stacked bar chart: "n=" + rounded mean electrode count per region
  mutate(label = paste0("n=", round(n_elec_mean, 1)))

plot_elec = ggplot(electrodes_labels,
                   aes(x = factor(session_age), y = n_elec_mean, fill = region)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 1, color = "black", alpha = 0.33) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 8, size.unit = "pt") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette  = "Set1") +
  THEME_BASE + THEME_TEXT +
  theme(legend.position = "right", legend.direction = "vertical") +
  labs(x = "Age (months)", y = "Mean # Retained Electrodes", fill = "Region")


# =============================================================================
# SECTION 5: R² AND MAE DUAL-AXIS BOXPLOT
# =============================================================================
# Summarise R² and MAE to one row per subject × visit, then display as
# side-by-side boxplots (R² in blue, MAE in red) on a shared y-axis.
# scale_factor and offset are chosen so MAE boxes sit within [0.9, 1].

# Data collapse: electrode → region-subject-visit → subject-visit (for boxplot visualization)
aperiodic_region = aperiodic_raw |>
  filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% eeg_desc$sujid, chinclu == 1) |>
  group_by(sujid, session_age, region) |>
  summarise(r2_mean = mean(r2value, na.rm = TRUE), mae_mean = mean(mae, na.rm = TRUE),
            .groups = "drop") |>
  group_by(sujid, session_age) |>
  summarise(r2_mean = mean(r2_mean), mae_mean = mean(mae_mean), .groups = "drop")

scale_factor = 1
offset       = 0.9  # Shifts MAE boxes into the [0.9, 1] display range

fit_plot = ggplot(aperiodic_region, aes(x = factor(session_age))) +
  geom_boxplot(aes(y = r2_mean),
               width = 0.35, position = position_nudge(x = -0.2),
               fill = "steelblue", color = "black", outlier.colour = "steelblue", alpha = 0.5) +
  geom_boxplot(aes(y = (mae_mean * scale_factor) + offset),
               width = 0.35, position = position_nudge(x = 0.2),
               fill = "indianred", color = "black", outlier.colour = "indianred", alpha = 0.5) +
  scale_y_continuous(
    name     = "Model Fit (R\u00b2)",
    limits   = c(0.9, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Mean MAE")
  ) +
  THEME_BASE + THEME_TEXT +
  theme(
    axis.title.y.left  = element_text(color = "steelblue",  face = "bold", size = 12, margin = margin(r = 10)),
    axis.title.y.right = element_text(color = "indianred",  face = "bold", size = 12, margin = margin(l = 10)),
    axis.line.y.right  = element_line(color = "black")
  ) +
  labs(x = "Age (months)")


# =============================================================================
# SECTION 6: PSD RESIDUAL ERROR CURVES
# =============================================================================
# Mean absolute error per frequency per visit, with IQR ribbon.
# Dashed red lines mark ±MAE threshold; dashed black line marks zero.

# Compute mean and interquartile range of PSD residual error across subjects
psds_data_plot = psd_data |>
  group_by(session_age, freq) |>
  summarise(
    m_error     = mean(abserror,                       na.rm = TRUE),
    error_first = quantile(abserror, 0.25, na.rm = TRUE),
    error_third = quantile(abserror, 0.75, na.rm = TRUE),
    .groups     = "drop"
  )

psds_plot = ggplot(psds_data_plot,
                   aes(x = freq, y = m_error, color = factor(session_age),
                       fill = after_scale(color))) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = error_first, ymax = error_third), alpha = 0.20, linewidth = 0) +
  geom_hline(yintercept =  mae_thresh, linetype = "dashed", color = "darkred",  linewidth = 1.0) +
  geom_hline(yintercept = -mae_thresh, linetype = "dashed", color = "darkred",  linewidth = 1.0) +
  geom_hline(yintercept =  0,          linetype = "dashed", color = "black",    linewidth = 1.5) +
  scale_color_manual(values = age_palette) +
  THEME_BASE + THEME_TEXT +
  theme(legend.position = "right", legend.direction = "vertical") +
  labs(x = "Frequency (Hz)", y = "Absolute Error", color = "Age (months)")


# =============================================================================
# SECTION 7: HISTOGRAM — TOTAL RETAINED ELECTRODES BY VISIT
# =============================================================================
# Distribution of retained electrodes per subject at each visit.
# Useful for spotting under-representation of subjects at specific ages.

# Count distinct channels per subject-visit after all quality filters
electrode_check = aperiodic_raw |>
  filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% descriptives$sujid, chinclu == 1) |>
  group_by(sujid, session_age) |>
  summarise(n_elec = n_distinct(ch), .groups = "drop")

electrodes_plot_2 = ggplot(electrode_check, aes(x = n_elec)) +
  geom_bar(fill = COLORS_MAIN[1], color = "black", alpha = 0.7) +
  facet_wrap(~ session_age) +
  THEME_BASE + THEME_TEXT +
  labs(x = "Retained Electrodes per Subject", y = "Count")


# =============================================================================
# SECTION 8: ASSEMBLE AND SAVE FIGURE
# =============================================================================

fit_plot_top    = ggpubr::ggarrange(plot_elec, fit_plot,
                                     ncol = 2, nrow = 1, labels = c("a", "b"))
fit_plot_merged = ggpubr::ggarrange(fit_plot_top, psds_plot,
                                     ncol = 1, nrow = 2, heights = c(1, 1),
                                     labels = c("", "c"))

ggsave(fit_plot_merged,
       filename = file.path(path2figs, "Fig_SM2_ParametrizedPSD_QualityMetrics_2to15Hz_range_extended.jpeg"),
       width = 27, height = 18, units = "cm", dpi = 300)



# =============================================================================
# SECTION 9: DESCRIPTIVES TABLE (gtsummary → .docx)
# =============================================================================
# Reports mean R², MAE, and number of retained electrodes per visit.
# Saved as a formatted Word table with 2-decimal precision.

aperiodic_raw |>
  filter(r2value > r2_thresh, mae < mae_thresh, epochs >= epochs_threshold,
         sujid %in% descriptives$sujid, chinclu == 1) |>
  # Format session age for output labels (e.g., '1 mo.', '6 mo.', ...)
  mutate(session_age = paste0(session_age, " mo.")) |>
  group_by(sujid, session_age, region) |>
  summarise(
    r2_mean  = mean(r2value, na.rm = TRUE),
    mae_mean = mean(mae,     na.rm = TRUE),
    n_elec   = n_distinct(ch),
    .groups  = "drop"
  ) |>
  group_by(sujid, session_age) |>
  summarise(r2_mean = mean(r2_mean), mae_mean = mean(mae_mean),
            n_elec = sum(n_elec), .groups = "drop") |>
  # Order session_age factor levels chronologically for table display
  mutate(session_age = factor(session_age, levels = c('1 mo.', '6 mo.', '12 mo.', '18 mo.', '30 mo.', '36 mo.', '42 mo.', '48 mo.'))) |>
  tbl_summary(
    by        = session_age,
    include   = c(r2_mean, mae_mean, n_elec),
    missing   = "no",
    label     = list(r2_mean  ~ "R\u00b2",
                     mae_mean ~ "MAE",
                     n_elec   ~ " # Retained\nElectrodes"),
    statistic = list(all_continuous() ~ "{mean} ({sd})\n[{p2.5}, {p97.5}]")
  ) |>
  modify_header(label ~ "**Metric**") |>
  bold_labels() |>
  modify_caption("Parametrization model fit (R\u00b2, MAE) and retained electrode count by visit.") |>
  as_flex_table() |>
  flextable::autofit() |>
  flextable::colformat_double(digits = 3) |>
  theme_booktabs() |> autofit() |>
  align(align = "left", part = "all") |>
  align(j = 2:8, align = "center", part = "all") |>
  bold(part = "header")|>
  padding(padding = 1, part = "all") |>
  flextable::font(fontname = "Arial", part = "all") |> fontsize(size = 10, part = "all") |>
  flextable::save_as_docx(path = file.path(path2tabs, "Table_S4_ParametrizedPSD_QualityMetrics_Descriptive.docx"))


