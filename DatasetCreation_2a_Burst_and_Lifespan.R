# =============================================================================
#  CODE Notes
# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE - code identifier -> Table 1 term
#   DATA columns -> DISPLAY LABELS (used in plot y-axis/legend labels):
#     volt_amp / *_corrected   -> "Volt. Amp."       (Voltage Amplitude, absolute/corrected)
#     band_amp / *_corrected   -> "Band Amp."        (Band Amplitude, absolute/corrected)
#     prop_bursty_epochs       -> "Prop. of Epochs w/ Burst"
#     prop_bursty_cycles_burst -> "Prop. of Cycles w/ Burst"
#     avg_burst_duration       -> "Burst Duration"   (Consecutive cycles per burst)
#     alpha_LAcH               -> "Lifespan"         (Alpha lifespan in cycles; LAcH cumsum >=90%)
#     is_burst / burst_type    -> "Cycle Type"       ("Burst" vs. "NoBurst")
#     prop_epochs              -> "Clean Epochs"     (EEG-quality covariate, z-scored)
#     mae                      -> reported as MAE (mean absolute error)
# -----------------------------------------------------------------------------
# Script: DatasetCreation_2a_Burst_and_Lifespan.R
# Purpose: Processes burst_properties_bycycle CSVs per age; computes burst duration, bursty cycles, corrected amplitudes, and frequency bands.
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

# --- Data Loading and Preparation ---------------------------------------
# Comment or uncomment depending on the datasets you want to process. The first path2sets and path2save are for the main burst properties dataset, while the second pair is for the epoch-adjusted dataset (see Methods)
path2sets = path2sets   # <<< CHECK: set this to the folder holding the per-visit input CSVs (see config_paths.R: path2sets)
path2save = path2sets
adapted_epoch = TRUE

if (adapted_epoch  == T) {
  path2save = paste0(path2save, 'adjepoch/')
}

if (!dir.exists(path2save)) dir.create(path2save, recursive = FALSE)

electrodes = readr::read_csv(file.path(path2data, "electrodes.csv"))%>%
  filter(chinclu == 1)%>%dplyr::select(label, region)%>%
  # Remap region codes to full region labels (on-disk contract uses abbreviated codes)
  mutate(region = case_when(
    region == "Fr" ~ "Frontal",
    region == "P" ~ "Parietal",
    region == "T" ~ "Temporal",
    region == "O" ~ "Occipital",
    TRUE ~ "Central"
  ))

mintime     = 15 # Minimum time in seconds of data to be included in the analysis. This is based on the minimum number of bursty cycles we can capture with a 3s epoch (see Methods). We apply this threshold at the subject level, so if a subject has at least one session with more than 15s of bursty data, they are included in the analysis.
# Notice this inclusion is different with respect to the other scripts in which we exclude based on the number of epochs. This is because of the adjusted epoch duration test. To unify the code in both non-normalized and length-normalized epochs, excluding based on time was a better approach than multiple conditionals in the code. 
ages        = STUDY_VISITS


bursti      = TRUE
alpha_range = 'adapted' # adapted (each age had one range)
freqob      = TRUE     # If true, the frequencies are only based on burst cycles for the lagged coherence lifespan. 

lcohi       = FALSE
lcoh_type   = 'hilb'
lcoh_metric = 'coh'
low_thresh  = FALSE


for (age in ages) {
  age_name = paste0('Age', age)
  # Session-age remapping: harmonize across data-acquisition vs. analysis naming.
  # Remap to analysis age for dataset consistency (15→18, 40→42).
  if (age != 40 & age != 18) {
    session_age = age
  } else if (age == 40) {
    session_age = 42
  } else {
    session_age = 15  # 18-month visit stored under nominal label '15' on disk
  }
  
  if (adapted_epoch) {
    if (age == 1) {epoch_t = 3} else if (age == 6) {epoch_t = 1.8} else if (age < 30) {epoch_t = 1.5} else {epoch_t = 1.3}
  } else {
    epoch_t = 3
  }
  
  sub_path2sets = paste0(path2sets, age_name)
  
  if (bursti) {
    sets = list.files(path = sub_path2sets,
                       pattern = "burst_properties_bycycle.*\\.csv$",
                       recursive = TRUE,
                       full.names = TRUE)
    
    if (alpha_range == "adapted") {sets = sets[(!grepl('stdalpha', sets) & !grepl('wider', sets))]} else if (alpha_range == 'stdalpha') {sets = sets[grepl('stdalpha', sets)]} else {sets = sets[grepl('wider', sets)]}
    
    for (s in sets) { # Loop through each set - do not apply lappy because RAM limitations in Age 1
      
      tdata = read.csv(s, header = TRUE)%>%filter((block != 'EOrs' | is.na(block)))
      
      burst_cycles_absolute = tdata|>
        # Convert is_burst string coding ("True"/"False") to boolean (TRUE/FALSE) for rle()
        mutate(burst_bool = if_else(is_burst == 'True' | is_burst == T, TRUE, FALSE))|>
        group_by(sujid, band, ch, epoch)%>%
        summarise(burst_lengths = list(with(rle(burst_bool), lengths[values == TRUE])),
                  num_bursts = length(unlist(burst_lengths)),
                  avg_burst_duration = ifelse(num_bursts > 0, mean(unlist(burst_lengths)), 0))|>
        filter(avg_burst_duration > 0)|>
        ungroup()%>%
        group_by(sujid, band, ch)%>%
        summarise(avg_burst_duration     = mean(avg_burst_duration, na.rm = TRUE))
      
      bursty_cycles = tdata%>%
        group_by(sujid, band, ch, epoch)%>%
        summarise(bursty_cycles      = sum(is_burst == "True" | is_burst == T, na.rm = TRUE),
                  prop_bursty_cycles = sum(is_burst == "True" | is_burst == T, na.rm = TRUE)/n())%>%
        ungroup()%>%
        group_by(sujid, band, ch)%>%
        # MANUSCRIPT Table 1 -> "Prop. of Cycles w/ Burst" (Prop. of cycles containing bursts
        # to total number of cycles). `prop_bursty_cycles_burst` conditions the mean on epochs
        # that contain >=1 bursty cycle; `prop_bursty_cycles` averages over all epochs.
        summarise(prop_bursty_cycles_burst = mean(prop_bursty_cycles[prop_bursty_cycles>0], na.rm = TRUE),
                  prop_bursty_cycles       = mean(prop_bursty_cycles, na.rm = TRUE))
      
      bursty_cycles = left_join(bursty_cycles, burst_cycles_absolute)
      
      # MANUSCRIPT Table 1 -> "Prop. of Epochs w/ Burst" (Prop. of epochs containing bursts
      # to total number of epochs). 
      bursty_epochs = tdata%>%
        group_by(sujid, band, ch, epoch)%>%
        summarise(bursty_epochs = if_else(sum(is_burst == "True" | is_burst == T, na.rm = TRUE)>0,1,0))%>%
        ungroup()%>%
        group_by(sujid, band, ch)%>%
        summarise(prop_bursty_epochs = mean(bursty_epochs, na.rm = TRUE))  # -> "Prop. of Epochs w/ Burst"
      
      tdata = tdata%>%
        group_by(sujid, band, ch, epoch)%>%
        mutate(volt_amp_corrected     = volt_amp/ mean(volt_amp[is_burst=="False" | is_burst == F], na.rm = TRUE))%>% # Legacy that is not congruent w/ other papers using this metric. Afterwards is recomputed.
        mutate(band_amp_corrected = band_amp/ mean(band_amp[is_burst=="False" | is_burst == F], na.rm = TRUE))%>%
        mutate(rise_decay_asym   = (time_rdsym / (time_rdsym + time_ptsym)) - .5,
               peak_trough_asym  = (time_peak / (time_peak + time_trough)) - .5)%>%
        dplyr::select(sujid, band, ch, epochs,
                      is_burst, volt_amp, band_amp, volt_amp_corrected, band_amp_corrected,
                      rise_decay_asym, peak_trough_asym, frequency)%>%
        group_by(sujid, band, ch, is_burst)%>%
        summarise(epochs = mean(epochs, na.rm = TRUE),
                  volt_amp = mean(volt_amp, na.rm = TRUE),
                  band_amp = mean(band_amp, na.rm = TRUE),
                  volt_amp_corrected     = mean(volt_amp_corrected, na.rm = TRUE),
                  band_amp_corrected = mean(band_amp_corrected, na.rm = TRUE),
                  rise_decay_asym   = mean(rise_decay_asym, na.rm = TRUE),
                  peak_trough_asym  = mean(peak_trough_asym, na.rm = TRUE),
                  sd_frequency      = sd(frequency, na.rm = TRUE),
                  frequency         = mean(frequency, na.rm = TRUE))%>%
        ungroup()%>%
        # Recode is_burst from string ("True"/"False") to display labels ("Burst"/"NoBurst")
        mutate(is_burst = if_else(is_burst == "True" | is_burst == T, "Burst", "NoBurst"))
      
      tdata = merge(tdata, bursty_cycles, by = c("sujid", "band", "ch"), all.x = TRUE)
      tdata = merge(tdata, bursty_epochs, by = c("sujid", "band", "ch"), all.x = TRUE)
      
      if (s == sets[1]) {
        data = tdata
        rm(tdata)
      } else {
        data = bind_rows(data, tdata)
        rm(tdata)
      }
      
    }
    
    inclusion = data%>%
        group_by(sujid)%>%summarise(inclusion = if_else(mean(epochs, na.rm = TRUE)*epoch_t >=mintime, 1, 0))
  
    data = merge(data%>%mutate(session_age = session_age,
                                alpha_range = alpha_range), inclusion, by = c("sujid"), all.x = TRUE)
  
    write.csv(data, 
              file = paste0(path2save, 'burst_properties_bycycle_long_', alpha_range, '_', age_name, '.csv'),
              row.names = FALSE)
    
    
    frequencies = data %>%
      group_by(sujid, inclusion, band, is_burst)%>%
      summarise(sd_frequency   = sd(frequency, na.rm = TRUE),
                mean_frequency = mean(frequency, na.rm = TRUE))%>%
      ungroup()%>%
      group_by(sujid, inclusion, band)%>%
      mutate(burst_sum  = sum(is_burst == 'Burst'))%>%
      # Select frequency band edges based on burst cycles (if present) or non-burst (if no burst observed).
      # Burst cycles occur at lower frequencies, so burst-based band is preferred for heterogeneous subjects.
      filter((burst_sum == 1 & is_burst == 'Burst') | (burst_sum == 0 & is_burst == 'NoBurst'))%>%
      dplyr::select(-is_burst)
    
    frequencies = frequencies %>%
      mutate(band          = 'alpha',
             min_frequency = round((mean_frequency - 1*sd_frequency)/0.5)*.5,
             max_frequency = round((mean_frequency + 1*sd_frequency)/0.5)*.5)%>%
      dplyr::select(sujid, inclusion, band, min_frequency, max_frequency)
    
    frequencies = reshape2::melt(frequencies, id.vars = c("sujid", "inclusion", "band"))
    frequencies = reshape2::dcast(frequencies, sujid + inclusion ~ variable + band,
                                    fill = NA)
    write_csv(frequencies, file = paste0(path2save, 'frequency_bands_onlyburst_', alpha_range, '_', age, 'mo.csv'))

    frequencies = data %>%
        group_by(sujid, inclusion, band)%>%
        summarise(sd_frequency   = sd(frequency, na.rm = TRUE),
                  mean_frequency = mean(frequency, na.rm = TRUE))%>%ungroup()
      
    frequencies = frequencies %>%
        mutate(band          = 'alpha',
               min_frequency = round((mean_frequency - 1*sd_frequency)/0.5)*.5,
               max_frequency = round((mean_frequency + 1*sd_frequency)/0.5)*.5)%>%
        dplyr::select(sujid, inclusion, band, min_frequency, max_frequency)
      
      frequencies = reshape2::melt(frequencies, id.vars = c("sujid", "inclusion", "band"))
      frequencies = reshape2::dcast(frequencies, sujid + inclusion ~ variable + band,
                                    fill = NA)
      write_csv(frequencies, file = paste0(path2save, 'frequency_bands_allcycles_', alpha_range, '_', age, 'mo.csv'))
  
    rm(list=c("data", "frequencies"))
    
    }
  
  if (lcohi) { 
    sets = list.files(path = sub_path2sets,
                         pattern = "lagged_coh_.*\\.csv$",
                         recursive = TRUE,
                         full.names = TRUE)
    
    sets = sets[!grepl('burst', sets)]
    sets = sets[grepl(paste0(lcoh_type, '_', lcoh_metric), sets)]
    
    if (low_thresh) {sets = sets[grepl('lowthresh', sets)]
    suf = "lowthresh"} else {sets = sets[!grepl('lowthresh', sets)] 
    suf = "" }
    
    lcoh_extension = paste(lcoh_type, lcoh_metric, suf, sep = '_')
    if (freqob ) {
       alpha_extension = paste0(alpha_range, '_onlyburst_') 
      } else { 
        alpha_extension = paste0(alpha_range, '_allcycles_') 
      }
    
      if (freqob) {
        if (session_age != 42) {frequencies = read_csv(paste0(path2save, 'frequency_bands_onlyburst_', alpha_range, '_', age, 'mo.csv'))} else {frequencies = read_csv(paste0(path2save, 'frequency_bands_onlyburst_', alpha_range, '_40mo.csv'))}
      } else { 
        if (session_age != 42) {frequencies = read_csv(paste0(path2save, 'frequency_bands_allcycles_', alpha_range, '_', age, 'mo.csv'))} else {frequencies = read_csv(paste0(path2save, 'frequency_bands_allcycles_', alpha_range, '_40mo.csv'))}
      }
    

    for (s in sets) {
      
      sdata = readr::read_csv(s)%>%filter((block != 'EOrs' | is.na(block)))
      
      # Harmonise the lagged-autocoherence column to the canonical `LAcH`.
      # Upstream emits `lcoh_avg` (burst-conditioned, jrpc.process_coherence_bursts)
      # or `lcoh` (all-cycles, EEG_metrics_rest_NN). Accept either; `any_of` makes
      # this a no-op if the file already ships `LAcH`.
      sdata = sdata %>% rename(any_of(c(LAcH = "lcoh_avg", LAcH = "lcoh")))
      stopifnot("LAcH" %in% colnames(sdata))
      
      sdata = sdata%>%
        # Normalize lagged coherence (cLAcH): LAcH / total LAcH per electrode-frequency
        mutate(total = sum(LAcH), .by = c(sujid, freq, ch)) %>%
        mutate(cLAcH = LAcH / total, .by = c(sujid, freq, ch)) %>%
        # Sort by cycle index to compute cumulative normalized coherence (LAcH cutoff logic)
        arrange(sujid, freq, ch, cycle) %>%
        # Cumulative sum of normalized coherence: identifies when alpha "lifespan" (90% cumulative coherence) is reached
        group_by(sujid, freq, ch) %>%
        mutate(cumclcoh = cumsum(cLAcH)) %>%
        ungroup()
      
      cycles_per_cycle = sdata%>%
        dplyr::select(sujid, freq, ch, cycle, cumclcoh)
      
      LAcH_lifespan = cycles_per_cycle %>%
        group_by(sujid, freq, ch) %>%
        # LAcH (Lifespan): cycle at which cumulative coherence first crosses 90% threshold
        filter(cumclcoh >= 0.90) %>%
        slice_head(n = 1) %>%
        ungroup() 
      
      LAcH_lifespan = merge(LAcH_lifespan, frequencies, by = c("sujid"), all.x = TRUE)
      
      LAcH_lifespan = LAcH_lifespan %>%
        group_by(sujid, inclusion, ch) %>%
        summarise(alpha_LAcH = mean(cycle[freq >= min_frequency_alpha & freq <= max_frequency_alpha]))
      
      
      sdata = merge(sdata, LAcH_lifespan, by = c("sujid", "ch"), all.x = TRUE)
      sdata = merge(sdata, frequencies, by = c("sujid", "inclusion"), all.x = TRUE)
      
      sdata = sdata %>%
        group_by(sujid, inclusion, ch)%>%
        summarise(epochs           = mean(epochs, na.rm = T), 
                  alpha_short_LAcH = mean(LAcH[freq >= min_frequency_alpha & freq <= max_frequency_alpha & cycle < .5], na.rm = TRUE),
                  alpha_long_LAcH  = mean(LAcH[freq >= min_frequency_alpha & freq <= max_frequency_alpha & cycle >= .5 & cycle <= 1.5], na.rm = TRUE),
                  epochs = mean(epochs, na.rm = TRUE),
                  alpha_LAcH = mean(alpha_LAcH, na.rm = TRUE))
      
      if (s == sets[1]) {
        data = sdata
      } else {
        data = bind_rows(data, sdata)
      }
      
      rm(sdata, cycles_per_cycle, LAcH_lifespan)
    }
    
  
    write.csv(data%>%mutate(session_age  = session_age,
                            alpha_range  = alpha_range,
                            alpha_freqob = if_else(freqob, 'AlphaOB', 'AlphaAllCycles'),
                            threshold    = if_else(low_thresh, '', '95%'),
                            lcoh_type    = lcoh_type,
                            lcoh_metric  = lcoh_metric), 
              file = paste0(path2save, 'lcohhilb_long_', lcoh_extension, alpha_extension, age_name, '.csv'),
              row.names = FALSE)
    
    rm(list=c("data"))}
  
}
