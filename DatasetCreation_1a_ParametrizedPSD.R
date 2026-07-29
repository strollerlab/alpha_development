# =============================================================================
# ALPHA BURST DEVELOPMENT PROJECT
# ----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term to columns
#
#   DATA CONTRACT columns to DISPLAY LABELS (used in plot y-axis/legend labels):
#     slope                   → "Slope" (Aperiodic slope)
#     offset                  → "Offset" (Aperiodic offset)
#     peak_freq               → "Peak Freq." (Alpha Peak Frequency, Hz)
#     peak_ampl               → "Peak Amp." (Alpha Peak Amplitude, μV²)
#     peak_prop                → "Prop. of Peaks" (Alpha Peak Proportion, 0–1)
#     osc_ampl                → "Band Power" (Oscillatory Alpha Band Power, μV²)
# -----------------------------------------------------------------------------
# Script: DatasetCreation_1a_ParametrizedPSD.R
# Purpose: For each study age, loads per-subject Specparam output CSVs, applies
#          EEG quality thresholds, defines individualized alpha frequency bands,
#          computes oscillatory band power from PSDs, and saves long-format
#          datasets (parameters + PSDs) per age group. These visit-stratified 
#          datasets are used in the merging code to generate the final datasets (as well as topomaps data). 
# =============================================================================
# Inputs:
#   - Data/Aperiodic/AgeNN/aperosc_parameters_*.csv   (per subject, per age)
#   - Data/Aperiodic/AgeNN/aperosc_psds_*.csv         (per subject, per age)
#   - electrodes.csv                                  (electrode map)
# Outputs:
#   - Data/aperiodic_oscillatory_long_allchannels_NNmo.csv  (one per age)
#   - Data/psds_aperosc_long_region_NNmo.csv                (one per age)
#   - Data/psds_aperosc_long_avg_NNmo.csv                   (one per age)
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

# =============================================================================
# SECTION 1: ANALYSIS PARAMETERS
# =============================================================================
# Local copies of global thresholds from 00_Setup_PackageInstallation.R.
# These allow the script to also be run as a standalone without sourcing 00_Setup.

r2_thresh        = R2_THRESH         # Minimum Specparam R-squared (0.900)
mae_thresh       = MAE_THRESH        # Maximum Specparam MAE (0.10)
                                     # NOTE: raw data column is called 'mae' 
epochs_threshold = EPOCHS_THRESHOLD  # Minimum clean epochs per electrode (5)
ch_threshold     = 0                 # Min. fraction of channels meeting inclusion
                                     # criteria (0 = retain all; track via goodch in next steps). 
                                     # This does not mean the criteria are not applied. Instead, it is filtered afterwards. 

# Specparam fitting settings
fiindx  = 15          # Specparam maximum frequency index: 15 Hz. See Supplementary Methods of the paper. 
fitmode = "no_prefit" # Specparam fitting mode used for this dataset. 

# Age-group-specific alpha band limits (Hz).
# Younger infants have slower alpha peaks; the window shifts upward with age.
# These defaults apply only when no individual peak frequency is available.
# visits 1  mo:    3–7 Hz
# visits 6  mo:    5–8 Hz
# visits 12–18 mo: 6–9 Hz
# visits 30+ mo:   7–10 Hz

visits = STUDY_VISITS  # All session visits (months) to process

# =============================================================================
# SECTION 2: PATH DEFINITIONS
# =============================================================================
# IMPORTANT: Update path2sets and path2data to match your local directory.
# path2sets: root folder for per-age Specparam output files (subfolders: Age1, Age6, …)
# path2data: destination folder for the merged CSV outputs

# path2data = "" Uncomment and add the path if a different path than those of config_paths.R is needed. 
# path2sets = ""

# Electrode map: links channel labels to brain region and inclusion flag.
# chinclu == 1 marks the 60 pre-selected analysis channels for consistency with other large studies. 
# region short codes (Fr, P, T, O) are expanded here for readability throughout the project.
electrodes = read.csv(file.path(path2data, "electrodes.csv")) |>
  dplyr::select(label, region, hemis, chinclu) |>
  mutate(region = case_when(
    region == "Fr" ~ "Frontal",
    region == "P"  ~ "Parietal",
    region == "T"  ~ "Temporal",
    region == "O"  ~ "Occipital",
    TRUE         ~ "Central"
  ))

# =============================================================================
# SECTION 3: MAIN PROCESSING LOOP (one iteration per age group)
# =============================================================================

for (visit in visits) {

  # --- 3a. Set age-specific alpha band limits ---
  # These defaults are used only for subjects without a detectable alpha peak.
  alpha = if      (visit == 1)        c(3, 7)
          else if (visit == 6)        c(5, 8)
          else if (visit < 30)        c(6, 9)
          else                        c(7, 10)

  # Remap session label 40 to 42 months
  session_age = if (visit != 40) visit else 42

  # --- 3b. Load Specparam parameter and PSD files for this age ---
  path2subsets = paste0(path2sets, "Age", visit, "/")

  sets      = list.files(path2subsets, pattern = "aperosc_parameters_*",
                          recursive = TRUE, full.names = TRUE)
  sets_psds = list.files(path2subsets, pattern = "aperosc_psds_*",
                          recursive = TRUE, full.names = TRUE)

  # Combine all subjects; keep only the 15 Hz Specparam fits and the chosen fitmode
  pow_and_aper = lapply(sets,      read.csv, header = TRUE) |> bind_rows() |>
    filter(maxfreq == fiindx, prefit == fitmode)

  psds = lapply(sets_psds, read.csv, header = TRUE) |> bind_rows() |>
    filter(maxfreq == fiindx, prefit == fitmode)

  # For visits > 36 months, retain only eyes-closed resting state ('ECrs') blocks
  # to match dim lights condition used at younger visits. Otherwise, we may have age-related differences due to variations partially explained by illumination levels. 
  if (visit > 36) {
    pow_and_aper = pow_and_aper |> filter(grepl("ECrs", block))
    psds         = psds         |> filter(grepl("ECrs", block))
  }

  pow_and_aper = merge(
    pow_and_aper |> mutate(ch = if_else(ch == "Cz", "E129", ch)),
    electrodes,
    by.x = "ch", by.y = "label"
  )

  # --- 3c. Compute per-subject channel-quality inclusion flags ---
  # goodch: number of channels passing both R² and MAE thresholds
  # ch_prop: proportion of included channels that pass both thresholds

  inclusion_fit = pow_and_aper |>
    filter(chinclu == 1) |>
    dplyr::select(sujid, ch, r2value, mae, epochs) |>
    mutate(
      goodch = sum(r2value > r2_thresh & mae < mae_thresh),
      ch_prop = mean(goodch) / n(),                          .by = sujid
    )

  # Electrode-level flags:
  #   inclusion_electrode_epochs: enough clean epochs
  #   inclusion_electrode_rsq    : R² above threshold & MAE
  #   inclusion_combined         : passes BOTH R² and epoch count
  inclusion_fit = inclusion_fit |>
    mutate(
      inclusion_electrode_epochs = epochs >= epochs_threshold,
      inclusion_electrode_rsq    = (r2value > r2_thresh & mae < mae_thresh),
      inclusion_combined         = as.integer(inclusion_electrode_rsq &
                                               epochs >= epochs_threshold),
      .by = sujid
    ) |>
    mutate(
      # Subject-level flag: enough channels pass combined criterion
      inclusion_final_nchan = sum(inclusion_combined),
      inclusion_final_dummy = inclusion_final_nchan >= round(ch_threshold * length(unique(ch))),
      .by = sujid
    ) |>
    dplyr::select(sujid, inclusion_final_nchan, inclusion_final_dummy, ch_prop, goodch) |>
    distinct()

  # Merge subject-level flags back and re-compute electrode-level flags
  pow_and_aper = left_join(pow_and_aper, inclusion_fit) |>
    mutate(
      inclusion_electrode_epochs = epochs >= epochs_threshold,
      inclusion_electrode_rsq    = r2value > r2_thresh,
      inclusion_combined         = as.integer(inclusion_electrode_rsq &
                                               epochs >= epochs_threshold),
      .by = sujid
    )

  # --- 3d. Detect alpha peaks ---
  # Create binary flags: 1 if Specparam detected a peak in that band (non-NA frequency), 0 otherwise
  pow_and_aper = pow_and_aper |>
    mutate(
      alpha_peak = as.integer(!is.na(alpha_freq)),
      theta_peak = as.integer(!is.na(theta_freq))
    )

  # --- 3e. Compute subject-level average alpha peak frequency and width ---
  # Used to define individualized (subject-specific) alpha band edges.
  # Inclusion criteria: high-quality R² fit, has detected alpha peak, pre-selected electrode set, width < 3 Hz.
  # Mean for frequency (balanced across channels); median for width (robust to outliers > 1.5 Hz).
  pow_and_aper = pow_and_aper |>
    mutate(
      alpha_freq_avg      = mean(  alpha_freq[r2value > r2_thresh & alpha_peak == 1 &
                                               chinclu == 1 & alpha_widt < 3], na.rm = TRUE),
      alpha_freq_widt_avg = median(alpha_widt[r2value > r2_thresh & alpha_peak == 1 &
                                               chinclu == 1 & alpha_widt < 3], na.rm = TRUE),
      .by = sujid
    )

  # --- 3f. Define individualised alpha band edges ---
  # Subjects with a detectable peak: band = [peak_freq - 0.5*width, peak_freq + 0.5*width]
  # Subjects without a detectable peak: use the age-group default alpha range
  pow_and_aper = pow_and_aper |>
    mutate(
      alpha_band_min = if_else(
        !is.na(alpha_freq_avg),
        alpha_freq_avg - 0.5 * alpha_freq_widt_avg,
        alpha[1]
      ),
      alpha_band_max = if_else(
        !is.na(alpha_freq_avg),
        alpha_freq_avg + 0.5 * alpha_freq_widt_avg,
        alpha[2]
      )
    ) |>
    # Round band edges to the nearest 0.5 Hz (Specparam frequency resolution)
    mutate(
      alpha_band_min = round(alpha_band_min * 2) / 2,
      alpha_band_max = round(alpha_band_max * 2) / 2
    )

  # --- 3g. Classify peak type per subject × channel ---
  # dummy_peak_ch: encoded peak classification per electrode (0=no peak, 1=theta only, 2=alpha only, 3=both)
  pow_and_aper = pow_and_aper |>
    mutate(
      dummy_peak_ch = case_when(
        theta_peak == 1 & alpha_peak == 1 ~ 3,
        theta_peak == 0 & alpha_peak == 0 ~ 0,
        theta_peak == 1                   ~ 1,
        TRUE                              ~ 2
      ),
      .by = c(sujid, ch)
    ) |>
    mutate(
      # Compute per-subject prevalence of each peak type across all channels
      prop_nopeak    = sum(dummy_peak_ch == 0) / length(unique(ch)),
      prop_bothpeaks = sum(dummy_peak_ch == 3) / length(unique(ch)),
      prop_thetapeak = sum(dummy_peak_ch == 1) / length(unique(ch)),
      prop_alphapeak = sum(dummy_peak_ch == 2) / length(unique(ch)),
      # type_of_peak: subject-level summary: 3=mixed (both peaks or theta+alpha coexist), else max dominant type
      type_of_peak      = if_else(
        prop_bothpeaks > 0 | (prop_thetapeak > 0 & prop_alphapeak > 0),
        3, max(dummy_peak_ch)
      ),
      .by = c(sujid, chinclu)
    )

  # Replace NA peak amplitudes with 0 (missing peak frequency to zero oscillatory amplitude by definition)
  pow_and_aper$theta_ampl[is.na(pow_and_aper$theta_ampl)] = 0
  pow_and_aper$alpha_ampl[is.na(pow_and_aper$alpha_ampl)] = 0

  # --- 3h. Compute alpha band power from PSDs ---
  # Extract individualised alpha band limits from parameter file; merge with PSD data to sum power
  # within each subject's personalised band. Integrates four PSD components: oscillatory, absolute, fooofed, aperiodic.
  frequencies = dplyr::select(pow_and_aper, sujid, alpha_band_min, alpha_band_max) |>
    distinct() |>
    mutate(session_age = session_age)

  band_powers = merge(psds, frequencies, by = "sujid") |>
    group_by(sujid, ch) |>
    summarise(
      alpha_osc = sum(oscillatory[freq >= alpha_band_min & freq <= alpha_band_max], na.rm = TRUE),
      alpha_abs = sum(absolute[   freq >= alpha_band_min & freq <= alpha_band_max], na.rm = TRUE),
      alpha_fit = sum(fooofed[    freq >= alpha_band_min & freq <= alpha_band_max], na.rm = TRUE),
      alpha_ape = sum(aperiodic[  freq >= alpha_band_min & freq <= alpha_band_max], na.rm = TRUE),
      .groups = "drop"
    )

  # Merge computed band powers back into the main parameter file (standardise ch to ch)
  pow_and_aper = merge(pow_and_aper, band_powers, by.x = c("sujid", "ch"), by.y = c("sujid", "ch"))

  # --- 3i. Save outputs ---
  # Full long-format channel-level file (all channels; one row per subject × channel)
  filename_long_all = file.path(path2data, "aperiodic_oscillatory_long_allchannels_", visit, "mo.csv")
  write.csv(pow_and_aper |> mutate(session_age = session_age), filename_long_all, row.names = FALSE)

  # PSD files for plotting: restrict to quality-passing, pre-selected channels only
  dummy_inclusion = dplyr::select(pow_and_aper, sujid, inclusion_final_dummy) |>
    group_by(sujid) |> distinct()

  psds = merge(psds, electrodes |> filter(chinclu == 1), by.x = "ch", by.y = "label") |>
    # Apply quality filters: R² criterion, epoch count, and fitting error (MAE)
    filter(r2value > r2_thresh, epochs >= epochs_threshold, mae <= mae_thresh)

  # Region-averaged PSDs (one row per subject × brain region × frequency) for plotting
  psds_region = psds |>
    group_by(sujid, region, freq) |>
    summarise(
      epochs      = mean(epochs,  na.rm = TRUE),
      r2value     = mean(r2value,      na.rm = TRUE),
      mae         = mean(mae,          na.rm = TRUE),
      absolute    = mean(absolute,     na.rm = TRUE),
      fooofed= mean(fooofed,      na.rm = TRUE),
      oscillatory = mean(oscillatory,  na.rm = TRUE),
      aperiodic   = mean(aperiodic,    na.rm = TRUE),
      error       = mean(error,        na.rm = TRUE),
      .groups     = "drop"
    )

  # Whole-brain average PSDs (one row per subject × frequency)
  psds_avg = psds |>
    group_by(sujid, freq) |>
    summarise(
      epochs      = mean(epochs,  na.rm = TRUE),
      r2value     = mean(r2value,      na.rm = TRUE),
      mae         = mean(mae,          na.rm = TRUE),
      absolute    = mean(absolute,     na.rm = TRUE),
      fooofed= mean(fooofed,      na.rm = TRUE),
      oscillatory = mean(oscillatory,  na.rm = TRUE),
      aperiodic   = mean(aperiodic,    na.rm = TRUE),
      error       = mean(error,        na.rm = TRUE),
      .groups     = "drop"
    )

  # Attach inclusion flag for easy filtering during plotting
  psds_avg  = merge(psds_avg,  dummy_inclusion, by = "sujid")
  psds_region = merge(psds_region, dummy_inclusion, by = "sujid")

  write.csv(psds_region |> mutate(session_age = session_age),
            file.path(path2data, "psds_aperosc_long_region_", visit, "mo.csv"), row.names = FALSE)

  write.csv(psds_avg  |> mutate(session_age = session_age),
            file.path(path2data, "psds_aperosc_long_avg_",  visit, "mo.csv"), row.names = FALSE)

}

