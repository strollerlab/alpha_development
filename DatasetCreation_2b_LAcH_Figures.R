# =============================================================================
# ALPHA BURST DEVELOPMENT PROJECT
# ----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - code identifier -> Table 1 term
#     alpha_LAcH               -> "Lifespan"         (Alpha lifespan in cycles; LAcH cumsum >=90%)
# -----------------------------------------------------------------------------
# Script: DatasetCreation_2b_LAcH_Figures.R
# Purpose: Computes lagged coherence tile-plot, region-averaged, and cumulative LAcH datasets for figures. Supports burst and non-burst files.
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

# Load and remap electrode region abbreviations to full display labels (FrtoFrontal, PtoParietal, etc.)
electrodes = read.csv(file.path(path2data, "electrodes.csv"))|>
  filter(chinclu == 1)%>%dplyr::select(label, region)%>%
  mutate(region = case_when(
    region == "Fr" ~ "Frontal",
    region == "P" ~ "Parietal",
    region == "T" ~ "Temporal",
    region == "O" ~ "Occipital",
    TRUE ~ "Central"
  ))

# Alpha burst parameters to load the frequency in lagged coherence computation
# alpha_range: 'adapted' (visit-specific)
alpha_range = 'adapted'

# freqob = TRUE: use burst-only alpha range; FALSE: use all-cycles to create alpha range
freqob      = TRUE

# Lagged coherence parameters: Hilbert-based coherence (hilb) with coherence metric
lcoh_type   = 'hilb'
lcoh_metric = 'coh'
low_thresh  = FALSE

# Loop over burst vs. non-burst cycle types (FALSE=non-burst, TRUE=burst)
is_burst = c(FALSE, TRUE)

path2sets = path2sets   # <<< CHECK: set this to the folder holding the per-visit input CSVs (see config_paths.R: path2sets)
path2save = path2sets
if (!dir.exists(path2save)) dir.create(path2save, recursive = TRUE)

epochs_threshold   = EPOCHS_THRESHOLD
visits             = STUDY_VISITS

for (b in is_burst) {
  for (visit in visits) {

    # Harmonize session_age labels: raw naming convention uses 15 to 18mo and 40 to 42mo for display.
    # This maps internal visit codes to standardized session_age for output dataset naming.
    if (visit == 40) {
      session_age <- 42
    } else if (visit == 15) {
      session_age <- 18
    } else {
      session_age <- visit
    }

    # Construct age folder path and find lagged coherence CSV files
    age_name = paste0('Age', visit)
    sub_path2sets = paste0(path2sets, age_name)
    sets <- list.files(path = sub_path2sets,
                       pattern = "lagged_coh_py_hilb.*\\.csv$",
                       recursive = TRUE,
                       full.names = TRUE)


    if (b) {
      # BURST CYCLE PATH: filter for burst-specific files
      sets <- sets[grepl('burst', sets)]
      alpha_extension = paste0(alpha_range, '_alphafreq_variable_')
      suf = ""

      # Load alpha frequency bands for burst and all-cycle contexts (will use conditional logic below)
      if (session_age != 42) {frequencies_allc  <- read_csv(paste0(path2save, 'frequency_bands_allcycles_', alpha_range, '_', visit, 'mo.csv'))} else {frequencies_allc <- read_csv(paste0(path2save, 'frequency_bands_allcycles_', alpha_range, '_40mo.csv'))}
      if (session_age != 42) {frequencies_burst <- read_csv(paste0(path2save, 'frequency_bands_onlyburst_', alpha_range, '_', visit, 'mo.csv'))} else {frequencies_burst <- read_csv(paste0(path2save, 'frequency_bands_onlyburst_', alpha_range, '_40mo.csv'))}
      # Rename frequency columns to distinguish burst vs. all-cycle bands for later conditional filtering
      colnames(frequencies_allc)  <- gsub('alpha', 'alpha_allc', colnames(frequencies_allc))
      colnames(frequencies_burst) <- gsub('alpha', 'alpha_burst', colnames(frequencies_burst))
      frequencies_burst = frequencies_burst%>%dplyr::select(-inclusion)
      lcoh_extension = paste0(lcoh_type, '_', lcoh_metric, '_')

    } else {
      # NON-BURST CYCLE PATH: filter for non-burst files (exclude burst keyword)
      sets <- sets[!grepl('burst', sets)]
      sets <- sets[grepl(paste0(lcoh_type, '_', lcoh_metric), sets)]
      # Optional low-threshold variant (usually FALSE)
      if (low_thresh) {
        sets = sets[grepl('lowthresh', sets)]
        suf = "lowthresh"
        lcoh_extension = paste0(lcoh_type, '_', lcoh_metric, '_', suf, '')
      } else {
        sets <- sets[!grepl('lowthresh', sets)]
        lcoh_extension = paste0(lcoh_type, '_', lcoh_metric, '_')}

      # Select frequency band set: burst-only (if freqob=TRUE) or all-cycles (if FALSE)
      if (freqob) {
        alpha_extension = paste0('onlyburst_', alpha_range)
      } else {
        alpha_extension = paste0('allcycles_', alpha_range)
      }
      # Load single frequency band table for non-burst cycles
      if (session_age != 42) {frequencies  <- read_csv(paste0(path2save, 'frequency_bands_', alpha_extension, '_', session_age, 'mo.csv'))} else {frequencies <- read_csv(paste0(path2save, 'frequency_bands_',alpha_extension, '_42mo.csv'))}
    }

    # Initialize output lists for aggregation
    avg_data         = list()  # Full-brain averages (freq × cycle level)
    region_sdata       = list()  # Region-specific averages (freq × cycle × region level)
    ch_sujcum_lcoh   = list()  # Channel-level cumulative LAcH (alpha lifespan)
    region_sujcum_lcoh = list()  # Region-level cumulative LAcH (alpha lifespan)

    # Exclude problematic subjects at 48 months - These crashed during the computation
    if (visit == 48) {sets = sets[!grepl('SUB-STXSSR', sets)]}
    if (visit == 48) {sets = sets[!grepl('SUB-RXVYKG', sets)]}
    
    for (s in sets) {
      # Load lagged coherence CSV and filter to resting-state block; merge electrode metadata
      sdata <- readr::read_csv(s)|>filter(block == 'ECrs' | is.na(block))
      sdata <- filter(sdata, ch %in% electrodes$label)
      sdata <- merge(sdata, electrodes, by.x = 'ch', by.y = 'label')
      
      #if lcoh is present rename, otherwise

      if (!b) {
        # NON-BURST: Simple aggregations without burst_type stratification
        sdata <- sdata|>rename(LAcH = lcoh)
        # Average LAcH across channels and electrodes (full-brain level)
        avg_data[[s]] <- sdata%>%group_by(sujid, freq, cycle)%>%
          summarise(inclusion = mean(epochs >= epochs_threshold),
                    LAcH      = mean(LAcH))

        # Average LAcH by brain region (ROI level)
        region_sdata[[s]] <- sdata%>%group_by(sujid, region, freq, cycle)%>%
          summarise(inclusion = mean(epochs >= epochs_threshold),
                    LAcH      = mean(LAcH))

        # Compute LAcH lifespan (cumulative sum normalized by total): cycle where cumsum reaches ~90%
        sujcum_LAcH  <- sdata%>%filter(epochs >= epochs_threshold)%>%
          mutate(total = sum(LAcH), .by = c(sujid, freq, ch)) %>%
          mutate(clcoh = LAcH / total, .by = c(sujid, freq, ch)) %>%
          arrange(sujid, freq, ch, cycle) %>%
          group_by(sujid, freq, ch) %>%
          mutate(cumclcoh = cumsum(clcoh)) %>%
          ungroup()

        # Merge in alpha frequency bands for this subject
        sfreq       <- filter(frequencies, sujid == unique(sdata$sujid))
        sujcum_LAcH <- merge(sujcum_LAcH, sfreq, by = c("sujid"), all.x = TRUE)

        # Channel-level LAcH lifespan: average cumclcoh within alpha band, by cycle
        ch_sujcum_lcoh[[s]] <- sujcum_LAcH%>%
          group_by(sujid, cycle, ch)%>%
          summarise(cumLAcH = mean(cumclcoh[freq >= min_frequency_alpha & freq <= max_frequency_alpha], na.rm = T),
                    inclusion = mean(epochs>=epochs_threshold))

        # Region-level LAcH lifespan: average channel cumLAcH within alpha band, by region and cycle
        region_sujcum_lcoh[[s]] <- filter(sujcum_LAcH)%>%
          group_by(sujid, cycle, region)%>%
          summarise(cumLAcH = mean(cumclcoh[freq >= min_frequency_alpha & freq <= max_frequency_alpha], na.rm = T),
                    inclusion = mean(epochs>=epochs_threshold))

    } else {
      # BURST: Rename burst-specific column names to match non-burst schema (lagtocycle, n_epochstoepochs, etc.)
      sdata <- sdata%>%
        rename(cycle = lag, LAcH = lcoh_avg)

      # BURST: Aggregations stratified by burst_type (burst vs. non-burst cycle classification)

      # Average LAcH by burst type (full-brain level)
      avg_data[[s]] <- sdata%>%filter(epochs >= epochs_threshold)%>%group_by(sujid, burst_type, freq, cycle)%>%
        summarise(inclusion = mean(epochs >= epochs_threshold),
                  LAcH      = mean(LAcH))

      # Average LAcH by burst type and brain region (ROI level)
      region_sdata[[s]] <- sdata%>%filter(epochs >= epochs_threshold)%>%group_by(sujid, burst_type, region, freq, cycle)%>%
        summarise(inclusion = mean(epochs >= epochs_threshold),
                  LAcH      = mean(LAcH))

      # Load both all-cycles and burst-only frequency bands for conditional filtering
      sfreq             <- filter(frequencies_allc, sujid == unique(sdata$sujid))
      sfreq_burst       <- filter(frequencies_burst, sujid == unique(sdata$sujid))

      # Compute LAcH lifespan stratified by burst type
      sujcum_LAcH  <- sdata%>%filter(epochs >= epochs_threshold)%>%
        mutate(total = sum(LAcH), .by = c(sujid, burst_type, freq, ch)) %>%
        mutate(clcoh = LAcH / total, .by = c(sujid, burst_type, freq, ch)) %>%
        arrange(sujid, burst_type, freq, ch, cycle) %>%
        group_by(sujid, burst_type, freq, ch) %>%
        mutate(cumclcoh = cumsum(clcoh)) %>%
        ungroup()

      # Merge both frequency band tables to enable conditional filtering by burst type
      sujcum_LAcH <- merge(sujcum_LAcH, sfreq, by = c("sujid"), all.x = TRUE)
      sujcum_LAcH <- merge(sujcum_LAcH, sfreq_burst, by = c("sujid"), all.x = TRUE)

      # Channel-level LAcH lifespan: use conditional alpha bands (burst vs. allc) based on burst_type
      # Then collapse redundant rows if created by reframe
      ch_sujcum_lcoh[[s]] <-  sujcum_LAcH %>%
        group_by(sujid, burst_type, cycle, ch, region) %>%
        filter(epochs >= epochs_threshold) %>%
        reframe(
          # Burst cycles use burst-specific alpha bands; non-burst use all-cycles bands
          cumLAcH = if_else(
            burst_type == 'burst',
            mean(cumclcoh[freq >= min_frequency_alpha_burst & freq <= max_frequency_alpha_burst], na.rm = TRUE),
            mean(cumclcoh[freq >= min_frequency_alpha_allc & freq <= max_frequency_alpha_allc], na.rm = TRUE)
          ),
          inclusion = mean(epochs >= epochs_threshold),
          .groups = "drop")%>%
        group_by(sujid, burst_type, cycle, ch, region) %>%
        summarise(cumLAcH = mean(cumLAcH, na.rm = T),
                  inclusion = mean(inclusion, na.rm = T))

      # Region-level LAcH lifespan: average channel-level cumLAcH within each region and burst type
      region_sujcum_lcoh[[s]] <-  ch_sujcum_lcoh[[s]]%>%
        group_by(sujid, burst_type, region, cycle) %>%
        summarise(cumLAcH = mean(cumLAcH, na.rm = T),
                  inclusion = mean(inclusion, na.rm = T))

    }}

    # Construct output suffix for burst vs. non-burst datasets
    burst_name = ifelse(b, '_burst', '')

    # Combine across subjects and add session_age column; save full-brain tile-plot dataset
    avg_data   = dplyr::bind_rows(avg_data)%>%mutate(session_age = session_age)

    write.csv(avg_data,
              file = paste0(path2save, 'lcohhilb_tileplot_avg_', lcoh_extension, age_name, burst_name, '.csv'),
              row.names = FALSE)
    rm(avg_data)

    # Combine and save region-level tile-plot dataset (LAcH × freq × region × cycle)
    region_sdata       = dplyr::bind_rows(region_sdata)%>%mutate(session_age = session_age)

    write.csv(region_sdata,
              file = paste0(path2save, 'lcohhilb_tileplot_region_', lcoh_extension, age_name, burst_name, '.csv'),
              row.names = FALSE)
    rm(region_sdata)

    # Combine and save channel-level cumulative LAcH dataset (lifespan curves by channel)
    ch_sujcum_lcoh  = dplyr::bind_rows(ch_sujcum_lcoh)%>%mutate(session_age = session_age)
    write.csv(ch_sujcum_lcoh,
              file = paste0(path2save, 'lcohhilb_cumplot_allch_', lcoh_extension, age_name, burst_name, '.csv'),
              row.names = FALSE)
    rm(ch_sujcum_lcoh)

    # Combine and save region-level cumulative LAcH dataset (lifespan curves by region)
    region_sujcum_lcoh = dplyr::bind_rows(region_sujcum_lcoh)%>%mutate(session_age = session_age)
    write.csv(region_sujcum_lcoh,
              file = paste0(path2save, 'lcohhilb_cumplot_region_', lcoh_extension, age_name, burst_name, '.csv'),
              row.names = FALSE)
    rm(region_sujcum_lcoh) 
      
  }
}






