# =============================================================================
# ALPHA BURST DEVELOPMENT PROJECT
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - code identifier -> Table 1 term
#   Naming convention (see NOMENCLATURE_CROSSWALK.md):
#     - proportions are `prop_*`  (never `per_*` / `percentage_*`)
#     - the unit of segmentation is the `epoch` (never "segment")
#     - the rhythmicity metric is `lifespan` (never "lifetime")
#
#   DATA columns -> DISPLAY LABELS (used in plot y-axis/legend labels):
#     offset                   -> "Offset"           (Aperiodic offset)
#     slope                    -> "Slope"            (Aperiodic exponent)
#     peak_freq                -> "Peak Freq."       (Peak Frequency, Hz)
#     peak_ampl                -> "Peak Amp."        (Peak Amplitude)
#     peak_per                 -> "Prop. of Peaks"   (Prop. of electrodes w/ oscillatory peak)
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
# Script: DatasetCreation_4_Topomaps.R
# Purpose: Creates channel-level averaged dataset for topographic heatmap visualization (aperiodic + coherence + burst).
# =============================================================================
# Dependencies: 00_Setup_PackageInstallation.R
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


# path2data comes from config_paths.R
# NOTE: This is a DATA-PREP script. Its only output is the intermediate CSV
# 'topological_heatmaps_data.csv' in path2data, which feeds the MATLAB topomap
# figures (Fig S3–S5, Fig SM3). It produces no Figures.
rsq_thresh       = R2_THRESH
epochs_incl      = EPOCHS_THRESHOLD
mae_thresh       = MAE_THRESH

# Load EEG quality descriptives; exclude pilots; z-score prop_epochs (% clean epochs) for uniform covariate scale
eeg_desc      = read_csv(file.path(path2data, "CrossVisit_EEG_CleaningDescriptives.csv"))|>
  filter(sujid != 'SUB-X2THXQ', sujid != 'SUB-YVEB6H')|>
  dplyr::select(sujid, session_age, prop_epochs, nsessions_withdata)|>
  mutate(prop_epochs = scale(prop_epochs)[,1])

# Load sociodemographic data; apply dev_filter for valid age records; center age_months to study baseline
desc_and_ages = read_csv(file.path(path2data, "Sociodemographic_Descriptives_Long_Updated.csv"))|>
  filter(dev_filter == 1, sujid != 'SUB-X2THXQ', sujid != 'SUB-YVEB6H')|>
  filter(sujid %in% eeg_desc$sujid)|>
  mutate(age_months = age_months - min(age_months))

# Extract one row per participant: demographic covariates (mean values) and gestational age; z-score all numeric columns
# ITN_mean missing values are imputed with the sample mean to preserve participants
desc_and_ages_wide = desc_and_ages|>
  dplyr::select(sujid, contains('mean'), GestationalAge_weeks)|>
  distinct()|>
  filter(sujid %in% eeg_desc$sujid)|>
  # Mean-impute ITN_mean (socioeconomic indicator) to avoid dropping participants
  # Mean-imputed only if present; ITN_mean is not a covariate in any model
  # here and is absent from the public release (see prepare_public_data.R).
  mutate(across(any_of("ITN_mean"), ~ if_else(is.na(.x), mean(.x, na.rm = TRUE), .x)))|>
  mutate(across(where(is.numeric), ~scale(.x)[,1]))

# Merge longitudinal visit info (per participant visit) with demographics; remove rows with missing EEG quality
descriptives = left_join(
  desc_and_ages|>dplyr::select(sujid, age_months, session_age, dev_filter, Cohort),
  eeg_desc) |>
  left_join(desc_and_ages_wide) |>
  filter(!is.na(prop_epochs))

# Aperiodic metrics: average across cycles within each electrode-visit cell, apply quality filters (R², MAE, trial count)
aper_voi       = c('sujid', 'session_age', 'ch', 'region', 'chinclu', 'mae', 'epochs', 'r2value', 'goodch', 'offset', 'slope', 'alpha_freq', 'alpha_ampl', 'alpha_osc', 'alpha_peak', 'inclusion_final_dummy', 'epochs')
aperiodic_data = read_csv(file.path(path2data, "Aperiodic_Oscillatory_ByCycle_Long.csv"))|>
  dplyr::select(all_of(aper_voi))|>
  # Harmonize session age coding
  mutate(session_age = if_else(session_age == 15, 18, session_age))|>
  group_by(session_age, ch)|>
  filter(r2value > rsq_thresh, mae < mae_thresh, sujid %in% eeg_desc$sujid)|>
  summarise(slope = mean(slope, na.rm = T),
            offset = mean(offset, na.rm = T),
            peak_ampl = mean(alpha_ampl, na.rm = T),
            peak_freq = mean(alpha_freq[alpha_peak ==1], na.rm = T),
            osc_ampl = mean(alpha_osc, na.rm = T),
            peak_prop = mean(alpha_peak, na.rm = T),
            epochs = mean(epochs, na.rm = T),
            r2value = mean(r2value, na.rm = T),
            mae = mean(mae, na.rm = T))

# Alpha lifespan (lagged coherence): average within each electrode-visit, harmonize age coding
hlcoh_data = read_csv(file.path(path2data, "LaggedCoh_Hilb_ByCycle_Long.csv"))|>
  group_by(session_age, ch)|>
  filter(sujid %in% eeg_desc$sujid, epochs >= epochs_incl)|>
  mutate(session_age = if_else(session_age == 15, 18, session_age))|>
  summarise(alpha_LAcH = mean(alpha_LAcH, na.rm = T))


# Burst properties: pivot Burst/NoBurst wide, compute amplitude correction ratios, rename for display
# Burst-exclusive metrics (prop_bursty_cycles_burst_Burst, prop_bursty_epochs_Burst) are renamed to user-friendly labels
voi = c('volt_amp', 'band_amp', 'frequency',
        'prop_bursty_cycles_burst','prop_bursty_epochs', 'avg_burst_duration')
burst_data = read_csv(paste0(path2data, '/BurstProperties_ByCycle_Long.csv'))|>
  filter(epochs >= epochs_incl, sujid %in% eeg_desc$sujid)|>
  dplyr::select(sujid, session_age, ch, is_burst, all_of(voi))|>
  # Harmonize session age coding
  mutate(session_age = if_else(session_age == 15, 18, session_age))|>
  group_by(session_age, is_burst, ch)|>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T)))|>
  # Pivot is_burst wide to enable Burst/NoBurst ratio computation
  pivot_wider(names_from = is_burst, values_from = c(volt_amp, band_amp, frequency, prop_bursty_cycles_burst, prop_bursty_epochs))|>
  rowwise()%>%
  # Compute amplitude correction factors (Burst / NoBurst ratios)
  mutate(volt_amp_corrected = volt_amp_Burst/volt_amp_NoBurst,
         band_amp_corrected = band_amp_Burst/band_amp_NoBurst)|>
  # Drop NoBurst proportions (not applicable to burst-exclusive metrics); strip the
  # "_Burst" suffix left by pivot_wider so the emitted topomap columns carry the same
  # canonical identifiers used everywhere else in the pipeline.
  # -> consumed by Figures_S3_S4_S5_topomaps_fieldtrip.m (groups(3).vars)
  dplyr::select(-c(prop_bursty_cycles_burst_NoBurst, prop_bursty_epochs_NoBurst))|>
  rename(prop_bursty_cycles_burst = prop_bursty_cycles_burst_Burst,   # -> "Prop. of Cycles w/ Burst"
         prop_bursty_epochs       = prop_bursty_epochs_Burst)         # -> "Prop. of Epochs w/ Burst"

# Merge all electrode-level metrics; create age label for topomap annotations (e.g., "12mo.")
topodata = left_join(aperiodic_data, hlcoh_data)|>
  left_join(burst_data)|>
  mutate(session_age_name = paste0(session_age, 'mo.'))


write.csv(topodata, paste0(path2data, '/topological_heatmaps_data.csv'))






  
