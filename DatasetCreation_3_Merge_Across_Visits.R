# =============================================================================
# ALPHA BURST DEVELOPMENT PROJECT
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - Table 1 term to columns
#
#   DATA CONTRACT columns to DISPLAY LABELS (used in plot y-axis/legend labels):
#     slope                  to"Slope" (Aperiodic slope)
#     offset                 to"Offset" (Aperiodic offset)
#     peak_freq              to"Peak Freq." (Alpha Peak Frequency, Hz)
#     peak_ampl              to"Peak Amp." (Alpha Peak Amplitude, μV²)
#     peak_prop               to"Prop. of Peaks" (Alpha Peak Proportion, 0–1)
#     osc_ampl               to"Band Power" (Oscillatory Alpha Band Power, μV²)
#     volt_amp / *_corrected to"Volt. Amp." (Peak-to-peak voltage, absolute/corrected)
#     band_amp / *_corrected to"Band Amp." (Burst band amplitude, absolute/corrected)
#     per_bursty_segments    to"Prop. of Epochs\nw/ Burst" (Proportion of epochs with alpha burst)
#     per_bursty_cycles_burstto"Prop. of Cycles\nw/ Burst" (Proportion of cycles with alpha burst)
#     avg_burst_duration     to"Burst Duration" (Consecutive cycles per burst)
#     alpha_LAcH             to"Lifespan" (Alpha lifespan in cycles; LAcH cumsum ≥90%)
#     is_burst / burst  to"Cycle Type" (categorical: "Burst" vs. "Non-Burst")
# -----------------------------------------------------------------------------
# Script: DatasetCreation_3_Merge_Across_Visits.R
# Purpose: Scans the Data folder, loads all per-age CSV files produced by
#          DtasetCreation_XX_R, and merges them into unified long-format datasets
#          ready for analysis.
# =============================================================================
# Inputs:
#   - Data/burst_properties_bycycle_long_*.csv     
#   - Data/lcohhilb_long_*.csv                     
#   - Data/aperiodic_oscillatory_long_allchannels_*.csv  
#   - Data/aperiodic_oscillatory_long_allch_*.csv       
#   - Data/hilb_coh_lifespan_long_burst*.csv
# Outputs:
#   - Data/BurstProperties_ByCycle_Long.csv
#   - Data/BurstProperties_ByCycle_Long_AdjEpoch.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long_allcycles.csv
#   - Data/Aperiodic_Oscillatory_ByCycle_Long.csv
#   - Data/Aperiodic_Oscillatory_ByCycle_Long_45Hz.csv
#   - Data/Aperiodic_Oscillatory_ByCycle_Long_Burst.csv
#   - Data/LaggedCoh_Hilb_ByCycle_Long_Burst.csv
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
# SECTION 1: PATH DEFINITIONS
# =============================================================================
# IMPORTANT:  Update path2sets to match your local data directory.
# path2sets : root folder containing all per-age Data CSV files
# path2save : destination folder for the merged output files

# path2sets comes from config_paths.R
path2save = paste0(path2sets, "")

# R² threshold applied at merge time only for the 45 Hz aperiodic dataset.
# All other merged files are filtered at analysis time (scripts 09–17).
r2_thresh = R2_THRESH  # 0.900 (from 00_Setup_PackageInstallation.R)
mae_thresh = MAE_THRESH  # 0.100 (from 00_Setup_PackageInstallation.R)

if (!dir.exists(path2save)) dir.create(path2save, recursive = TRUE)

# =============================================================================
# SECTION 2: BURST PROPERTIES — STANDARD AND EPOCH-ADJUSTED
# =============================================================================

# Discover all burst property CSV files
burst_sets = list.files(
  path       = path2sets,
  pattern    = "burst_properties_bycycle_long.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

# --- Epoch-adjusted burst file (sensitivity analysis) ---
# Contains burst metrics normalised for epoch-length differences across ages.
burst_sets_adj = burst_sets[grepl("adjepoch", burst_sets)]

burst_data_adj = lapply(burst_sets_adj, read.csv, header = TRUE) |>
  bind_rows() |>
  group_by(sujid) |>
  # Apply project naming convention: remap session label 15to18 months
  # (Python pipeline stores this visit as 15; project terminology uses 18)
  mutate(session_age = if_else(session_age == 15, 18, session_age)) |>
  ungroup()

# --- Standard burst file (main analyses) ---
# Filter to adapted frequency-band definition; exclude region-level, "Check", adjepoch
# and EOrs variants. NB: "area" here matches a legacy FILENAME on disk (pre-rename
# region-averaged outputs), not a column -- kept so stale files are still excluded.
burst_sets = burst_sets[
  !grepl("area",     burst_sets) &
   grepl("adapted",  burst_sets) &
  !grepl("Check",    burst_sets) &
  !grepl("adjepoch", burst_sets) &
  !grepl("EOrs",     burst_sets)
]

# NAMING LOGIC: Python burst pipeline stores 15-month visit as session_age=18.
# Remap 18to15 here; downstream scripts then apply standard remap 15to18 via PROJECT_NAMING_CONVENTIONS.
# This two-stage convention ensures files have consistent nomenclature across all R pipeline steps.
burst_data = lapply(burst_sets, read.csv, header = TRUE) |>
  bind_rows() |>
  group_by(sujid) |>
  # Invert Python pipeline's labeling: 18to15 (will be remapped 15to18 at analysis time)
  mutate(session_age = if_else(session_age == 18, 15, session_age)) |>
  ungroup()


# =============================================================================
# SECTION 3: LAGGED COHERENCE — ONLY-BURST AND ALL-CYCLES
# =============================================================================
lcoh_sets = list.files(
  path       = path2sets,
  pattern    = "lcohhilb_long_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)

# Remove region-averaged ("avg") and double-underscore variants (keep channel-level only)
lcoh_sets = lcoh_sets[!grepl("avg", lcoh_sets)]
lcoh_sets = lcoh_sets[!grepl("__",  lcoh_sets)]

# --- Only-burst coherence (adapted band; primary analysis) ---
# Filters to: burst-only cycles, adapted frequency bands, channel-level granularity (no "avg" or "__" variants).
lcoh_sets_ob = lcoh_sets[grepl("onlyburst", lcoh_sets) & grepl("adapted", lcoh_sets)]

lcoh_data = lapply(lcoh_sets_ob, read.csv, header = TRUE) |>
  bind_rows() |>
  group_by(sujid) |>
  # Apply naming remap (same convention as burst_data: invert Python pipeline's 18to15)
  mutate(session_age = if_else(session_age == 18, 15, session_age)) |>
  ungroup()

# --- All-cycles coherence (adapted band; used in supplement, script 17) ---
# Includes both burst and no-burst cycles; filtered to adapted frequency-band definition.
lcoh_sets_ac = lcoh_sets[grepl("allcycles", lcoh_sets) & grepl("adapted", lcoh_sets)]

lcoh_data_allcycles = lapply(lcoh_sets_ac, read.csv, header = TRUE) |>
  bind_rows() |>
  group_by(sujid) |>
  # Apply naming remap (same Python pipeline inversion: 18to15)
  mutate(session_age = if_else(session_age == 18, 15, session_age)) |>
  ungroup()


# =============================================================================
# SECTION 4: APERIODIC / OSCILLATORY — 15 Hz
# =============================================================================

aperosc_sets = list.files(
  path       = path2sets,
  pattern    = "aperiodic_oscillatory_long_allchannels_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)


# --- 15 Hz aperiodic/oscillatory dataset (used in select analyses requiring lower frequency ceiling) ---
aperosc_sets_15Hz = aperosc_sets[!grepl("adapted", aperosc_sets)]
aperosc_data = lapply(aperosc_sets_15Hz, read.csv, header = TRUE) |>
  bind_rows() |>
  dplyr::select(-any_of(c("epochsprop", "nepochs_group", "totalepochs", "interimr2",
                           "exclur2", "parameters_interim", "methodfooof",
                           "methodpsd", "banddef", "maxfreq", "block")))


# =============================================================================
# SECTION 5: APERIODIC — BURST-CONDITIONED VARIANT
# =============================================================================
# Separate FOOOF fits computed on burst vs. non-burst cycle segments.
# Applies same filtering and naming standardisation as Section 4.

aperosc_sets_burst = list.files(
  path       = path2sets,
  pattern    = "aperiodic_oscillatory_long_allch_.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)
aperosc_sets_burst = aperosc_sets_burst[!grepl("adapted", aperosc_sets_burst)]  # 15 Hz ceiling only

aperosc_data_burst = lapply(aperosc_sets_burst, read.csv, header = TRUE) |>
  bind_rows() |>
  # NB: the Python extractor already emits `epochs`; no rename needed (the old
  # `rename(trials = finaltrials)` step is obsolete.
  # Drop epoch-selection metadata: `epochsprop` (prop. of epochs kept by the interim
  # R2 selection) is not a model covariate; `epochs` (the count) is kept.
  dplyr::select(-any_of(c("epochsprop", "nepochs_group", "interimr2", "exclur2",
                           "parameters_interim", "methodfooof", "methodpsd",
                           "banddef", "maxfreq", "block")))

readr::write_csv(aperosc_data_burst, paste0(path2save, "Aperiodic_Oscillatory_ByCycle_Long_Burst.csv"))


# =============================================================================
# SECTION 6: LAGGED COHERENCE — BURST-CONDITIONED VARIANT
# =============================================================================
# Alpha Lifespan (LAcH) computed separately within burst and non-burst cycle segments.

lcoh_sets_burst = list.files(
  path       = path2sets,
  pattern    = "hilb_coh_lifespan_long.*\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)
lcoh_sets_burst = lcoh_sets_burst[grepl("_burst", lcoh_sets_burst)]  # Restrict to burst-conditioned files

lcoh_data_burst = lapply(lcoh_sets_burst, read.csv, header = TRUE) |>
  bind_rows() |>
  group_by(sujid) |>
  # Apply naming remap (same pipeline convention: 18to15)
  mutate(session_age = if_else(session_age == 18, 15, session_age)) |>
  ungroup() 

# =============================================================================
# SECTION 7: SAVE ALL MERGED FILES
# =============================================================================

readr::write_csv(burst_data,            paste0(path2save, "BurstProperties_ByCycle_Long.csv"))
readr::write_csv(burst_data_adj,        paste0(path2save, "BurstProperties_ByCycle_Long_AdjEpoch.csv"))
readr::write_csv(lcoh_data,             paste0(path2save, "LaggedCoh_Hilb_ByCycle_Long.csv"))
readr::write_csv(lcoh_data_allcycles,   paste0(path2save, "LaggedCoh_Hilb_ByCycle_Long_allcycles.csv"))
readr::write_csv(aperosc_data,          paste0(path2save, "Aperiodic_Oscillatory_ByCycle_Long.csv"))
readr::write_csv(lcoh_data_burst,       paste0(path2save, "LaggedCoh_Hilb_ByCycle_Long_Burst.csv"))
