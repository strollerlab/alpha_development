"""
@author: J. Rico-Picó

SCRIPT ROLE (pipeline stage 1 of 2)
------------------------------------------------------------------------------
Primary per-participant feature-extraction script for the manuscript
manuscript. For every cleaned resting-state .set file (3 s epochs, MADE
preprocessing pipeline) this script optionally computes, per participant:
    1. Absolute/relative band power from the raw (non-burst-conditioned)
       power spectral density (PSD)                          [abspow]
    2. Aperiodic + oscillatory (specparam/FOOOF) parameterization of the PSD
       [aperosc]
    3. Cycle-by-cycle alpha burst detection and burst properties (bycycle)
       [bursbct]
    4. Hilbert-based lagged auto-coherence (LAcH), used downstream to derive
       "alpha lifespan" in R                                   [laggedcoh_hil]
Each measure is independently toggled via the 0/1 flags in the "Parameters
to extract" section below. Outputs are per-participant CSVs (and, for LAcH,
.npy matrices) written under `output_path/<sujid>/`.

The companion script is pipeline stage
2: it re-reads the burst-properties CSV produced here and re-splits the PSD /
aperiodic / LAcH computations by burst vs. no-burst epochs.
"""

# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE
#   specparam settings  -> "EEG Power-Spectrum Parameterization" (2-15 Hz, fixed
#                          aperiodic, peak_threshold=2, min_peak_height=0.01,
#                          max_n_peaks=7, peak_width_limits=[1,10]).
#   bycycle thresholds (all = 0.5 = median; "Alpha Burst Properties"):
#     amp_fraction        -> peak-to-peak voltage > median           (criterion 1)
#     period_consistency  -> phase/temporal consistency < median     (criterion 2)
#     amp_consistency     -> amplitude consistency < median          (criterion 2)
#     monotonicity        -> rise/fall symmetry < median             (criterion 3)
#   LAcH / lagged-coherence
# -----------------------------------------------------------------------------

"""
Interest Variables: determine the parameters and measures that you want to compute in this script
Measures:  Relative and Absolute Power (log or not logtrans), Aperiodic and Oscillatory, Burst, Lagged Coherence
"""

# ---------------------------------------------------------------------------
# CONFIGURATION — the only lines to edit (stage 1 (all epochs))
#   CODE_DIR : this folder (holds convenience_functions_jrp.py and
#              lagged_autocoherence.py)
#   SET_DIR  : EEGLAB .set files to process
#   DATA_DIR : where the per-subject CSVs are written. Must match the
#              path2sets set in config_paths.R, so the R stage finds them.
# ---------------------------------------------------------------------------
CODE_DIR = "~/AlphaBurstRhythm/Code/"
SET_DIR  = "~/AlphaBurstRhythm/Data/set/"
DATA_DIR = "~/AlphaBurstRhythm/Data/"

import os
CODE_DIR, SET_DIR, DATA_DIR = (os.path.expanduser(p) for p in (CODE_DIR, SET_DIR, DATA_DIR))

import subprocess
import sys

# Third-party dependencies required by this script (spectral / burst / MNE-Python stack).
packages = ['numpy', 'scipy', 'pandas', 'matplotlib', 'neurodsp', 'tqdm', 'pytest', 'mne', 'specparam', 'neurodsp']

# Ensure all packages are installed
# Dependency check: on a fresh machine/environment this pip-installs
# any missing package before the real imports below run, so the script is runnable
for package in packages:
    try:
        __import__(package)
        print(f"{package} is already installed.")
    except ImportError:
        print(f"{package} is not installed. Installing now...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

import os
import numpy as np
import mne
import pandas as pd
import bycycle as bc  # cycle-by-cycle burst detection (Cole & Voytek, 2019)

# Local helper modules (not on PyPI): PSD/FOOOF convenience wrappers and the
# lagged auto-coherence (LAcH) implementation used for "alpha lifespan".
sys.path.append(CODE_DIR)   # convenience_functions_jrp.py / lagged_autocoherence.py
import convenience_functions_jrp as jrpc
import lagged_autocoherence as la


###############################################################################
# Path and basic parameters ###################################################
###############################################################################

# Specify folder containing the .set files
# (For the manuscript EEGLAB-format epoched, cleaned resting-state recordings; 3 s epochs, MADE (Raw data needs to be in .set format)
#  preprocessing pipeline output). Set SET_DIR in the CONFIGURATION block at the
#  top of this file. Run the script once per cohort/visit, pointing SET_DIR at
#  that visit's folder and DATA_DIR at the matching output folder.

folder_path = SET_DIR

# Where do you want to save the files
# Destination root; a per-participant subfolder (named after `sujname`, see
# renaming block below) is created under this path for every processed file.
output_path  = os.path.join(DATA_DIR)

os.makedirs(output_path, exist_ok = True)

# Set files in the folder, modify how the file ends
# Only files ending in "rest.set" (i.e., the resting-state condition) are queued for processing.
set_files = [f for f in os.listdir(folder_path) if f.endswith('rest.set')]

# ---- Run-control / batching flags -------------------------------------------------
skipi      = 0   # 1 = skip participants whose expected output files already exist (resumable batch runs); 0 = always reprocess
upsidedown = 0   # 1 = reverse the file processing order (useful for running two machines from opposite ends of the list concurrently)
multicomp  = 0   # 1 = split set_files across parallel jobs/machines using the quartile logic below; 0 = process the full list
stepm      = 2   # Which half/quarter of the split this job should take when multicomp == 1 (1 or 2; see block below)
rewritei   = 1   # Passed through to jrpc.save_dataframe: 1 = overwrite existing output CSVs, 0 = do not overwrite and generates csv adding the version (version = previous version found in the subject folder + 1)
numchar    = 8   # Number of leading characters of the raw filename used as sujid when rename == 0 (e.g., P123456_Resting_Paradigm.set numchar = 7)
nversion   = 0   # Which is the version you are running? If set to 0 no version added, otherwise vZ to skip files



# ---- Participant-ID renaming (maps raw filenames -> manuscript subject IDs) --------
rename     = 1        # 1 = build sujid as prename + substring(namepos1:namepos2) + posname; 0 = use the first `numchar` characters of the filename as-is
namepos1   = 1
namepos2   = 4
prename    = 'BUDDY'
posname    = ''

# Reverse processing order (see `upsidedown` above).
if upsidedown == 1:
   set_files = set_files[::-1]

# Multi-parallel batching: carve the file list into a subset for this run.
# stepm == 2 takes (roughly) the second half of the list; stepm == 1 takes the first
# quarter (plus a remainder correction when `upsidedown == 0`). This is a manual,
# non-deterministic load-splitting convenience for running the pipeline on >1 machine
# it does not affect which participants are *eventually* processed
if multicomp == 1 :
        half_participants = int(round(len(set_files)/4))
        if upsidedown == 0 and stepm == 1 and half_participants*4 < len(set_files):
            half_participants = half_participants + (len(set_files)-half_participants*4)

        if stepm == 2:
            startpoint = int(round(len(set_files)/2))
            if upsidedown == 0 :
                startpoint = startpoint + 1
        elif stepm == 1:
            startpoint = 0


        set_files = set_files[startpoint:startpoint+half_participants-1]


##############################################################################
# Parameters to extract: 0 = No compute, 1 = Compute (for aperosc, abspow is necessary)
###############################################################################
abspow         = 1   # Absolute/relative band power from the raw PSD (Methods: "Power Spectrum Density")
aperosc        = 1   # Aperiodic + periodic (specparam/FOOOF) parameterization (Methods: "EEG Power-Spectrum Parameterization")
bursbct        = 1   # Cycle-by-cycle alpha burst detection via bycycle (Methods: "Alpha Burst Properties")
laggedcoh_hil  = 1   # Hilbert-based lagged auto-coherence (Methods: LAcH / "alpha lifespan" precursor)

###############################################################################
# Signal per blocks? - MODIFY
################################################################################
blocks = 1   # 1 = split epochs by condition block (e.g., eyes-closed/eyes-open) using block_events below; 0 = treat all epochs as one pool
block_events  = {'ECrs',
                 'EOrs'} #Events in the EEG
block_names   = {'ECrs': 'EyesClosed', #Equivalence for writing
                 'EOrs': 'EyesOpen'}

filter_event    = 1        # 1 = drop epochs whose event name matches `event_to_filter` before any computation
event_to_filter = 'SOCV'   # Substring of event names to exclude (e.g., non-resting-state / task markers)

################################################################################
# Freq Range, Notice min should be at least twice higher than the s in the epoch
################################################################################

# Only computed/needed if abspow or aperosc is requested (bursbct/laggedcoh_hil use
# their own frequency ranges defined further below, independent of this block).
if abspow == 1 or aperosc == 1: #We save this information compute this information if necessary
    fmin, fmax     = 2, 45     # PSD frequency window retained for band-power/FOOOF (Hz)
    nranges        = 4         # Number of canonical bands defined below (theta/alpha/beta/gamma)
    range_theta    = [3, 6]
    range_alpha    = [7, 10]
    range_beta     = [11, 20]
    range_gamma    = [21, 45]

    # Integer-index lookup (used when iterating band_names/bands_ranges by position)...
    band_names = {
        0: 'theta',
        1: 'alpha',
        2: 'beta',
        3: 'gamma'}

    # ...and the matching frequency-range lookup, kept in the same index order as band_names.
    bands_ranges = {
        0: range_theta,
        1: range_alpha,
        2: range_beta,
        3: range_gamma}

    # Name-keyed version of the same band definitions, persisted verbatim into
    # output CSVs (column `banddef`) as the analysis "receipt" for reproducibility.
    band_definition = {
        'theta': range_theta,
        'alpha': range_alpha,
        'beta':  range_beta,
        'gamma': range_gamma}

##############################################################################
# Power spectrum computation parameters - MODIFY
##############################################################################

    methodpsd       = 2     #Hanning window = 1, Welch = 2, Spectopo like = 3, Multitaper Auto = 4, Multitaper Manual = 5
    badnwidth_res   = 'NA'  # Frequency resolution when multitapers used, recommended 1 or 2
    n_tapers        = 'NA'  # Number tapers in manual - recomended 3
    n_perseg        = 1000  #Welch and spectopo number of timepoints to run the fft (e.g., 1000 with 1000 srate = 1s)
    n_overlap       = 500   #Welch and spectopo points of overlapping (e.g., in 1000 timewindow 500 will be 50% overlap)

##############################################################################
# Define the Parameters and Methods to Compare - DON'T MODIFY
###############################################################################
# Derived/bookkeeping dicts built from the settings above. These are what gets
# passed to jrpc.process_psd() and what gets serialized into the `psdmethod`
# column of the output CSVs. This way, R scripts can confirm which PSD
# estimator produced a given row.

    extra_params = {
        'nperseg': n_perseg,      # For Welch and Spectopo
        'noverlap': n_overlap,        # For Welch and Spectopo
        'time_bandwidth': badnwidth_res,    # For multitaper methods
        'num_tapers': n_tapers,        # For manual multitaper
        'fmin': fmin,              # For built-in multitaper
        'fmax': fmax              # For built-in multitaper
    }

    method_names = {
        1: "FFT/Hanning",
        2: "Welch",
        3: "Spectopo",
        4: "Multitaper (MNE)",
        5: "Multitaper (Manual)"
    }

    # Build the method-specific parameter dict actually used for PSD computation
    # and for the human-readable `psdmethod` provenance string in the outputs.
    if method_names[methodpsd] == "FFT/Hanning":
        psd_parameters = {
            'method': 'FFT/Hanning',
            'overlap': 'noover'}
    elif method_names[methodpsd] == "Welch":
        psd_parameters = {
            'method': 'Welch',
            'nperseg': n_perseg,
            'overlap': n_overlap/n_perseg}
    elif method_names[methodpsd] == "Spectopo":
        psd_parameters = {
            'method': 'Spectopo',
            'nperseg': n_perseg,
            'overlap': n_overlap/n_perseg}
    elif method_names[methodpsd] == "Multitaper (MNE)":
        psd_parameters = {
            'method': 'MultitaperMNE',
            'bw': badnwidth_res}
    else: 
        psd_parameters = {
            'method': 'MultitaperManuel',
            'bw': badnwidth_res,
            'ntapers': n_tapers}
    
 ##############################################################################
 # Aperiodic computation parameters - MODIFY
 ##############################################################################
 # These map directly onto the "EEG Power-Spectrum Parameterization" settings
 # in the Methods/nomenclature box at the top of this file (specparam/FOOOF).
if aperosc == 1:
    psd_fmin          = 2       # Lower fit bound (Hz)
    peak_width_limits = [2, 12] # Allowed oscillatory-peak bandwidth (Hz)
    max_n_peaks       = 5       # Max number of periodic peaks specparam may fit
    min_peak_height   = 0.05    # Minimum peak height (above the aperiodic fit) to be retained
    peak_threshold     = 2       # Peak selection threshold (in SD of the flattened spectrum)
    rsquare_threshold = 0.90    # Final-fit R^2 acceptance threshold (goodness of aperiodic+periodic fit)
    aperiodic_mode    = 'fixed' # Aperiodic component form: 'fixed' (no spectral knee). It is not adapted to the 'knee' yet!

    ## If you are performing a reduction by excluding the bad fit epochs
    # epoch-selection / pre-fitting strategy used to guard against noisy single-trial
    # spectra biasing the aperiodic fit (see jrpc.interim_fooof_epochselection):
    fit_pregroups             = 0    # 0 = fit the mean PSD only; 1 = pre-fit groups of epochs and discard poorly-fit ones (used here); 2 = fit every trial individually
    rsquare_threshold_interim = .90  # R^2 cutoff applied during the interim (pre-group) fit used to flag/discard bad epochs
    perelect                  = 0.80 # Minimum fraction of electrodes that must pass rsquare_threshold_interim for a epoch-group to be retained
    epochs_multi              = 1    # 1 = iterate over every group size in `nepochs_group_multi`; 0 = use the single `nepochs_group_prec` value instead
    nepochs_group_multi       = [9999, 10]  # epoch-group sizes to test when epochs_multi == 1 (9999 is a sentinel meaning "no grouping / fit all epochs together" in our studies. If you foresee more than 9999 epochs in your study, just increase it to an unreasonable number)
    nepochs_group_prec        = 5    # Fallback single group size when epochs_multi == 0
    psd_fmaxs                 = [15, 20, 40]  # Upper fit-frequency bounds compared across runs (noise check on fmax)

 ##############################################################################
 # Burst Properties - MODIFY
 ##############################################################################
 # bycycle (Cole & Voytek) cycle-by-cycle burst-detection settings
if bursbct == 1:
    #theta_burst = [4, 7]
    burst_alpha = [7, 10]   # Alpha band used for burst detection ONLY (independent of range_alpha)
    #burst_beta  = [10.6, 20.90]
    #burst_gamma = [21, 45]
    alpha = burst_alpha     # kept for backward-compat with references below; sourced from burst_alpha only


    # Burst detection parameters (bycycle threshold_kwargs; all thresholds are
    # fractions of the corresponding feature's distribution, 0.5 = median split):
    amp_fraction       = .5  #Amplitude to be considered a burst (.5 = median)     -> criterion 1
    amp_consistency    = .5  # amplitude consistency threshold                     -> criterion 2
    period_consistency = .5  # phase/temporal consistency threshold                -> criterion 2
    monotonicity       = .5  # rise/fall symmetry threshold (median = .5); Methods "amplitude consistency"  -> criterion 3
    mincycles          =  3  # Minimum consecutive cycles required to accept a burst run
    treshold_features  =  1  # if 1 keep the values of the treshold parameters in each cycle, ej retains amp_fraction
    only_bursts        =  0  # if set to 1, export data only of the bursts

    # Cycle-level feature columns retained from bc.features.compute_features()
    # output (see bycycle docs); persisted verbatim as CSV columns (data contract).
    features_keep = [
             'is_burst', 'period', 'time_decay', 'time_rise',
             'time_peak', 'time_trough', 'volt_decay', 'volt_rise',
             'volt_amp', 'volt_peak', 'volt_trough', 'time_rdsym',
             'time_ptsym', 'band_amp']

    # NOTE: this locally re-defines `band_names` (shadowing the theta/alpha/beta/gamma
    # dict built above under `if abspow == 1 or aperosc == 1:`) to alpha-only, since
    # burst detection in this manuscript is restricted to the alpha band.
    band_names_burst = {
        0: 'alpha'}

    bands_ranges_burst = {
        0: alpha,}

    band_definition_burst = {
        'alpha': alpha}

##############################################################################
# Lagged coherence parameters - MODIFY
##############################################################################
# Hilbert-based lagged auto-coherence (LAcH) which combined with the R pipeline computes "alpha lifespan" (cycle count at 90% cumulative LAcH).
if laggedcoh_hil == 1 :
    lcohtype   = 'hilb'  # 'hilb' = Hilbert-transform based LAcH (vs. 'fft' alternative supported by the loop below)
    plv_or_coh = 'coh'   # Metric flavor: 'coh' = coherence, otherwise phase-locking value (PLV)
    fminch     = 3       # Min carrier frequency for the LAcH frequency sweep (Hz)
    fmaxch     = 12      # Max carrier frequency for the LAcH frequency sweep (Hz)
    freqfacch  = .5      # Frequency step of the sweep (Hz)

    cminh      = 0.05    # Minimum lag, in oscillation cycles
    cmaxh      = 4       # Maximum lag, in oscillation cycles
    cfach      = .05     # Lag step, in oscillation cycles

    thresh_prctile_hilb  =  None# During the subrogation, which is the threshold (None = no shuffle)

    freq_rangeh = np.arange(fminch, fmaxch, freqfacch)  # Frequency list used
    lag_rangeh  = np.arange(cminh, cmaxh, cfach)        # Lag (cycle) sweep values actually used

    lcohhil_parameters = {
        'ncyclesmin': cminh,
        'ncyclesmax':cmaxh,
        'cyclestep': cfach,
        'freqs': [fminch, fmaxch]}

    # Base output filename stem; suffix (`suf`, set below) is appended to encode
    # the surrogate-threshold choice actually used for a given run.
    basefile_lcoh = "lagged_coh_py_" + lcohtype + '_' + plv_or_coh

    if lcohtype == 'hilb':
        if thresh_prctile_hilb is None:
            niter = 1
            suf = '_notresh'
        elif thresh_prctile_hilb  < 95:
            niter = 1000 # Modify
            suf   = "_lowtresh"
        else:
            niter = 1000 # Modify
            suf   = ""
    else:
        suf = ""
else:
    basefile_lcoh = ""
    suf           = ""


def files2process(setfile, output_path, abspow, aperosc, bursbct, laggedcoh_hil,
                  basefile_lcoh,
                  extension = '.csv', numversion = 0):
            """
            Check, for one participant, which of the requested output CSVs
            already exist on disk under `output_path/<sujname>/`.

            Used to implement resumable batch runs (`skipi == 1`, see main
            loop): if all requested measures already have an output file, the
            participant is skipped instead of being reprocessed.

            NOTE: relies on the module-level `sujname` variable set in the
            enclosing loop rather than the `setfile` argument (pre-existing
            behavior, left unchanged).

            Parameters
            ----------
            setfile : str
                Unused positionally (see note above) -- kept for call-site
                compatibility.
            output_path : str
                Root output directory; `output_path/sujname/` is checked.
            abspow, aperosc, bursbct, laggedcoh_hil : int (0/1)
                Which measures are requested for this run (mirrors the
                module-level flags of the same name).
            basefile_lcoh : str
                Fully-resolved LAcH output filename stem (including any
                threshold suffix) to check for.
            extension : str, default '.csv'
                File extension to check.
            numversion : int, default 0
                If > 0, checks for a versioned filename ('..._v<numversion>')
                instead of the default (unversioned) filename.

            Returns
            -------
            (totaldone, totalfiles) : tuple[int, int]
                Count of expected output files that already exist vs. the
                total number of output files expected for the requested
                measures. Main loop skips the participant when these are equal.
            """
            subfolder_path  = os.path.join(output_path, sujname)
            if numversion == 0:
                version = ''
            else:
                version = '_v' + str(numversion)

            print(f"Checking what files the participant {sujname} has")
            totalfiles = abspow + aperosc + bursbct + laggedcoh_hil
            totaldone  =  0
            if abspow  == 1:
                absfile    = 'psds' + version + extension
                abspowpath = os.path.join(subfolder_path, absfile)
                abs_dummy  = os.path.exists(abspowpath)
                totaldone  = totaldone + abs_dummy
                print(f"Power Spectrum Computation: {abs_dummy}")
            if bursbct == 1:
                burstfile    = 'burst_properties_bycycle' + version + extension
                burstpath = os.path.join(subfolder_path, burstfile)
                burst_dummy = os.path.exists(burstpath)
                totaldone = totaldone + burst_dummy
                print(f"Burst Computation: {burst_dummy}")
            if aperosc == 1:
                aperfile = 'aperosc_parameters' + version + extension
                ape_dummy = os.path.exists(os.path.join(subfolder_path, aperfile))
                print(f"Power Spectrum Parametrization: {ape_dummy}")
                totaldone = totaldone + ape_dummy
            if laggedcoh_hil == 1 :
                lagfilehil = basefile_lcoh + version + extension
                lcohpath  = os.path.join(subfolder_path, lagfilehil)
                lco_dummy = os.path.exists(lcohpath)
                print(f"Lagged Coherence: {lco_dummy}")
                totaldone = totaldone + lco_dummy
            return totaldone, totalfiles
        
 ## Starts the loops
 # Main per-participant processing loop: resolve the subject ID, decide
 # whether to skip (already-processed) participants, then read/clean the
 # EEGLAB epochs and dispatch to each requested measure's computation block.

errori     = 0
error_list = []
for s, set_file in enumerate(set_files):

            # Resolve the manuscript-facing participant ID (sujid) from the raw filename.
            if rename == 0:
                   sujname = set_file[0:numchar]
            else:
                   sujname = prename + set_file[namepos1:namepos2] + posname

            # Resumability check: skip participants who already have every
            # requested output file on disk (see files2process docstring above).
            if skipi == 1:
                basefile_lcoh_full = basefile_lcoh + suf
                filesdone, totalfiles = files2process(sujname, output_path, abspow=abspow, 
                                     aperosc=aperosc, bursbct=bursbct, 
                                     laggedcoh_hil=laggedcoh_hil,
                                     basefile_lcoh = basefile_lcoh_full,
                                     numversion = nversion, 
                                     extension = '.csv')
            else:
                totalfiles = 1
                filesdone  = 0
                
            if totalfiles-filesdone == 0 and skipi == 1:
                print(f"Participant {s} - {sujname} - was already processed: Skipping")
            else:
                print(f"Participant {s} - {sujname}: Starting processing")



                try:
                    filepath = os.path.join(folder_path, set_file)
                    subfolder_path = []
                    subfolder_path = os.path.join(output_path, sujname)
                    os.makedirs(subfolder_path, exist_ok=True)

                    # Load the cleaned, epoched EEGLAB (.set) recording for this participant.
                    raw = mne.io.read_epochs_eeglab(filepath)
                    raw.load_data()

                    # Optionally drop epochs whose event label matches `event_to_filter`
                    # (For this project non-resting-state task markers concatenated).
                    if filter_event:
                        events_to_keep = [key for key in raw.event_id.keys() if event_to_filter not in key]
                        raw            = raw[events_to_keep]


                    ch_names       = raw.ch_names
                    electrodes     = len(raw.ch_names)
                    electrode_names=raw.ch_names
                    num_epochs     = len(raw)
                    raw._data     *= 1000000  # Convert MNE's native Volts to microvolts (uV), the unit used throughout Methods/figures.
                    data           = raw._data
                    fs             = raw.info['sfreq'] # Extract the sampling rate

                    print(f"Starting the processing of {sujname}")
                    if bursbct == 1:
                         # bycycle expects (channels, epochs, time); MNE stores epochs as
                         # (epochs, channels, time), so swap the first two axes here.
                         data_copy  = np.copy(data)
                         burst_data = np.transpose(data_copy, axes=(1, 0, 2))

                    # Derive a per-epoch condition-block label ('EyesClosed'/'EyesOpen', etc.)
                    # so below measures can be computed overall and/or split by block.
                    if blocks == 1:
                        event_ids   = raw.events[:, 2]  # Third column contains the event IDs
                        event_names = [key for event_id in event_ids for key, value in raw.event_id.items() if value == event_id]
                        new_event_names = []
                        event_id_selection = {}
                        for name in event_names:
                                normalized_name = name  # Default to the original event name.
                                for block in block_events:
                                    if block in name:  # If the block string is a substring of the event name...
                                        normalized_name = block  # ...set it to that block.
                                        break  # Stop checking once a match is found.
                                new_event_names.append(normalized_name)

                    else:
                        # No block splitting requested: still need a per-epoch label array
                        # (used later for counting epochs/epoch_block indexing), so fall back
                        # to the raw (un-normalized) event name for every epoch.
                        event_ids   = raw.events[:, 2]  # Third column contains the event IDs
                        new_event_names = [key for event_id in event_ids for key, value in raw.event_id.items() if value == event_id]

                    event_array = np.array(new_event_names)

                    # ------------------------------------------------------------------
                    # MEASURE 1: Absolute / relative band power from the raw PSD
                    # (Methods: "Power Spectrum Density"; independent of burst status).
                    # ------------------------------------------------------------------
                    if abspow == 1 or aperosc == 1:
                        powspcpath = os.path.join(subfolder_path, "psds.csv")
                        print(f"Starting Power Spectrum computation of {sujname}")

                        # Compute the PSD once here; both the raw band-power measure (abspow)
                        # and the aperiodic/oscillatory fit (aperosc, further below) reuse it.
                        freqs, psds = jrpc.process_psd(raw, methodpsd, raw.info['sfreq'], **extra_params) # Returns PSD as epochs, epochs, freqs
                        freq_mask = (freqs >= fmin) & (freqs <= fmax)
                        foi = np.where(freq_mask)[0]
                        psds = psds[:, :, foi]  # Restrict to the analysis frequency window [fmin, fmax].

                        freqs = freqs[foi]

                        # Only (re)compute the trial/band-level breakdown if the output doesn't
                        # already exist, or a forced rewrite (skipi == 0) was requested.
                        if not os.path.exists(powspcpath) or skipi == 0:
                                    # Relative power = each frequency bin's power divided by the
                                    # total power summed across the full [fmin, fmax] window.
                                    psds_sum = np.nansum(psds[:,:,:], axis = 2, keepdims = True)
                                    relpsds  = psds/psds_sum #Here we compute the power-spectrum relative first correcting each freq. bin by the sum of all freqs. (i.e., f1 / (f1 + .... fN))
                                    column_names = ['ch', 'band', 'sujid', 'psdmethod', 'freq', 'Abspow', 'Relpow', 'nepochs']
                                    df = pd.DataFrame(columns=column_names)

                                    # Collapse the PSD into per-band absolute/relative power, one
                                    # row per (trial x channel x band), for every canonical band.
                                    for bi, band in enumerate(band_names):
                                        band_range = bands_ranges[bi]

                                        psd_band                     = {}
                                        psd_band_rel                 = {}
                                        freq_mask                    = (freqs >= band_range[0]) & (freqs <= band_range[1])
                                        foi = np.where(freq_mask)[0]  # Get indices for frequencies of interest
                                        # Band power = sum (not integral) of PSD across the band's frequency bins.
                                        psd_band[band_names[bi]]     = np.nansum(psds[:, :, foi], axis=2)
                                        psd_band_rel[band_names[bi]] = np.nansum(relpsds[:, :, foi], axis=2)

                                        if bi == 0:
                                            psds_bands_data = []

                                        # NOTE: `band` here is the dict *key* yielded by `enumerate(band_names)`
                                        # (0,1,2,3 for theta/alpha/beta/gamma as currently defined)
                                        for ch_idx, ch_name in enumerate(ch_names):

                                            for epoch_idx in range(np.size(psds,0)):


                                                # One output row per trial x channel x band: absolute + relative power,
                                                # tagged with the block label and full PSD-method/band-definition provenance.
                                                _blocklabel = new_event_names[epoch_idx]
                                                _nepochs_block = new_event_names.count(_blocklabel)
                                                rowpow = {'Epoch': epoch_idx + 1, 'ch': ch_name, 'sujid': sujname, 'psdmethod': str(psd_parameters), 'banddef': str(band_definition),'pow': 'abs', 'block':new_event_names[epoch_idx], 'band': band_names[band], 'nepochs':_nepochs_block}
                                                rowpow['Abspow']  = psd_band[band_names[band]][epoch_idx, ch_idx]
                                                rowpow['Relpow']  = psd_band_rel[band_names[band]][epoch_idx, ch_idx]

                                                psds_bands_data.append(rowpow)

                                    df_abs = pd.DataFrame(psds_bands_data)
                                    jrpc.save_dataframe(df_abs, subfolder_path, base_filename="trial_psdpow", rewrite=rewritei, ext=".csv")


                                    # Frequency-resolved (per-Hz, not band-collapsed) average PSD,
                                    # used for spectral plots (e.g., group-average spectra figures).
                                    if blocks == 0:

                                        avg_df_abs = df_abs.pivot_table(index=['ch', 'band', 'sujid', 'psdmethod', 'banddef'],
                                            values=['Abspow', 'Relpow', 'epochs']).reset_index()

                                        avg_psd_burst_abs_file_path = os.path.join(subfolder_path, "avg_pow_band.csv")
                                        avg_df_abs.to_csv(avg_psd_burst_abs_file_path, index=False)

                                        # Assuming avgpsd and avgrelpsd are 2D arrays with dimensions (n_channels, n_freq)
                                        n_channels = len(ch_names)
                                        n_freq = len(freqs)

                                        # No block splitting: average PSD/relative-PSD across all epochs.
                                        avgpsd = np.average(psds[:,:,:], axis = 0)
                                        avgrelpsd = np.average(relpsds[:,:,:], axis = 0)

                                        # Create the columns for channels and frequency by repeating and tiling the original arrays.
                                        data = {
                                            'ch': np.repeat(ch_names, n_freq),          # Repeat each channel name n_freq times
                                            'freq': np.tile(freqs, n_channels),         # Tile frequency values for each channel
                                            'sujid': sujname,                           # Constant subject ID for all rows
                                            'psdmethod': method_names[methodpsd],       # Constant PSD method for all rows
                                            'abs': avgpsd.flatten(),                    # Flatten the 2D array into 1D
                                            'rel': avgrelpsd.flatten(),                 # Flatten the 2D array into 1D
                                            'block': None,
                                            'nepochs':psds.shape[0]
                                        }

                                        # Create the DataFrame from the dictionary.
                                        df_psddata = pd.DataFrame(data)

                                    else:
                                        # Block splitting requested: repeat the same frequency-resolved
                                        # averaging, but separately for each condition block (e.g.,
                                        # eyes-closed vs. eyes-open), using only that block's epochs.
                                        avg_df_abs = df_abs.pivot_table(index=['ch', 'block', 'band', 'sujid', 'psdmethod', 'banddef'],
                                            values=['Abspow', 'Relpow', 'epochs']).reset_index()

                                        df_list = []
                                        for block in block_events:
                                            blockpos = np.where(event_array == block)[0]
                                            if blockpos.size == 0 :
                                                continue  # This participant has no epochs for this block; skip it.

                                            n_channels = len(ch_names)
                                            n_freq     = len(freqs)

                                            avgpsd_block    = np.nanmean(psds[blockpos, :, :], axis=0)
                                            avgrelpsd_block = np.nanmean(relpsds[blockpos, :, :], axis=0)

                                            df_block = pd.DataFrame({
                                                'ch': np.repeat(ch_names, n_freq),
                                                'freq': np.tile(freqs, n_channels),
                                                'sujid': sujname,  # constant subject ID
                                                'psdmethod': method_names[methodpsd],  # constant PSD method
                                                'block': block,
                                                'abs': avgpsd_block.flatten(),
                                                'rel': avgrelpsd_block.flatten(),
                                                'epochs':len(blockpos)
                                            })

                                            df_list.append(df_block)

                                            # Combine the DataFrames for all blocks into one
                                            df_psddata = pd.concat(df_list, ignore_index=True)

                                    jrpc.save_dataframe(df_psddata, subfolder_path, base_filename="psds", rewrite=rewritei, ext=".csv")
                                    jrpc.save_dataframe(avg_df_abs, subfolder_path, base_filename="avg_pow_band", rewrite=rewritei, ext=".csv")
        
                    # ------------------------------------------------------------------
                    # MEASURE 2: Aperiodic + oscillatory parameterization (specparam/FOOOF)
                    # (Methods: "EEG Power-Spectrum Parameterization").
                    # Iterates over every combination of epoch-group size
                    # (epochs_group_multi) x upper fit-frequency bound (psd_fmaxs) x block,
                    # so the resulting CSV lets the R pipeline pick/compare fit settings
                    # rather than baking a single choice in at extraction time.
                    # ------------------------------------------------------------------
                    if aperosc == 1:
                        powspcpath   = os.path.join(subfolder_path, "psds.csv")
                        aperosc_path = os.path.join(subfolder_path, "aperosc_parameters.csv")
                        if not os.path.exists(aperosc_path) or skipi == 0:
                            print(f"Starting Power Spectrum Parametrization of {sujname}")
                            l = 1
                            aperosc_results_df = []
                            aperosc_psd_df = []

                            # Choose which epoch-group sizes to sweep over (see the
                            # "Aperiodic computation parameters" block for definitions).
                            if epochs_multi == 1:
                                nepochs_group_multi = nepochs_group_multi
                            else:
                                nepochs_group_multi = nepochs_group_prec

                            for nepochs_group in nepochs_group_multi:
                                for f_idx, fidx in enumerate(psd_fmaxs):

                                    jrpc.update_progress(l, (len(psd_fmaxs)*len(nepochs_group_multi)), bar_length=(len(psd_fmaxs)), message=f"{sujname} PSD Parametrization Progress: ")
                                    l = l + 1

                                    # specparam/FOOOF call signature for this (nepochs_group, fidx) combination.
                                    parameters = {
                                        "fmin": fmin,
                                        "fmax": fidx,
                                        "peak_width_limits": peak_width_limits,
                                        "max_n_peaks": max_n_peaks,
                                        "min_peak_height": min_peak_height,
                                        "peak_threshold": peak_threshold,
                                        "aperiodic_mode": aperiodic_mode}

                                    if blocks == 0:
                                        block_events = {'NA'}


                                    for blocki in block_events:
                                        if blocks == 0:
                                            tpsds = psds
                                        else:
                                            blockpos = np.where(event_array == blocki)[0]
                                            if blockpos.size == 0 :
                                                continue  # No epochs for this block for this participant; skip.
                                            else:
                                                tpsds    = psds[blockpos, :, :]

                                        # Core specparam/FOOOF fit, with optional pre-fit epoch-group
                                        # screening (fit_pregroups) to discard poorly-fit trial subsets
                                        # before the final fit -- see jrpc.interim_fooof_epochselection.
                                        aperosc_results, fooofpsd = jrpc.interim_fooof_epochselection(tpsds, freqs, ch_names, fit_pregroups, fmin, fidx,
                                                                           peak_width_limits, max_n_peaks, min_peak_height, peak_threshold,
                                                                           rsquare_threshold_interim, perelect, bands_ranges, aperiodic_mode=aperiodic_mode, nepochs_group=nepochs_group)

                                        # Tag the fit-result rows with full provenance (method + band
                                        # definitions + which sweep setting produced them).
                                        aperosc_results['methodfooof'] = str(parameters)
                                        aperosc_results['methodpsd']   = str(psd_parameters)
                                        aperosc_results['banddef']     = str(band_definition)
                                        aperosc_results['maxfreq']     = fidx
                                        aperosc_results['sujid']       = sujname
                                        aperosc_results['block']       = blocki
                                        
                                        # `nepochs_group` sentinel decoding (see nepochs_group_multi definition above)
                                        # 1 = epoch-by-trial pre-fit; any other value = pre-fit in groups of
                                        # that size. Recorded verbatim so the R pipeline can select/compare.
                                        if nepochs_group == 9999:
                                            aperosc_results['prefit']        = 'no_prefit'
                                            aperosc_results['nepochs_group'] = 'NA'
                                        elif nepochs_group == 1:
                                            aperosc_results['prefit']        = 'prefit_by_trial'
                                            aperosc_results['nepochs_group'] = 1
                                        else:
                                            aperosc_results['prefit']        = 'prefit_by_group'
                                            aperosc_results['nepochs_group'] = nepochs_group



                                        # Mask frequencies based on range
                                        freq_mask = (freqs >= fmin) & (freqs <= fidx)
                                        foi = freqs[freq_mask]  # Frequencies of interest (masked)

                                            # Initialize an empty list to store row
                                        # Extract the fitted (model + aperiodic-only) spectra at foi, alongside
                                        # the same fit-result metadata, for spectral-plot reconstruction downstream.
                                        aperosc_psd   = jrpc.fooofpsd_extraction(fooofpsd, aperosc_results, foi, ch_names, parameters=parameters)
                                        aperosc_psd   = pd.DataFrame(aperosc_psd)

                                        aperosc_psd['sujid']       = sujname
                                        aperosc_psd['methodpsd']   = str(psd_parameters)
                                        aperosc_psd['methodfooof'] = str(parameters)
                                        aperosc_psd['block']       = blocki
                                        aperosc_psd['parameters']  = parameters
                                        aperosc_psd['maxfreq']     = fidx

                                        if nepochs_group == 9999:
                                            aperosc_psd['prefit']        = 'no_prefit'
                                            aperosc_psd['nepochs_group'] = 'NA'
                                        elif nepochs_group == 1:
                                            aperosc_psd['prefit']        = 'prefit_by_trial'
                                            aperosc_psd['nepochs_group'] = 1
                                        else:
                                            aperosc_psd['prefit']        = 'prefit_by_group'
                                            aperosc_psd['nepochs_group'] = nepochs_group



                                        aperosc_results_df.append(aperosc_results)
                                        aperosc_psd_df.append(aperosc_psd)
                            # Stack every (nepochs_group x fidx x block) combination into two long-format
                            # tables: per-channel fit parameters, and the reconstructed fitted spectra.
                            df_aperosc = pd.concat(aperosc_results_df)
                            df_aperosc_psd = pd.concat(aperosc_psd_df)

                            jrpc.save_dataframe(df_aperosc, subfolder_path, base_filename="aperosc_parameters", rewrite=rewritei, ext=".csv")
                            jrpc.save_dataframe(df_aperosc_psd, subfolder_path, base_filename="aperosc_psds", rewrite=rewritei, ext=".csv")

                    # ------------------------------------------------------------------
                    # MEASURE 3: Cycle-by-cycle alpha burst detection (bycycle)
                    # (Methods: "Alpha Burst Properties"); this is the burst-detection
                    # step whose output (`burst_properties_bycycle.csv`) is consumed by
                    # the companion script EEG_metrics_rest_NN_basedonburst.py to split
                    # PSD/aperiodic/coherence measures by burst vs. no-burst epochs.
                    # ------------------------------------------------------------------
                    if bursbct == 1:
                        burstpath = os.path.join(subfolder_path, "burst_properties_bycycle.csv")
                        if not os.path.exists(burstpath) or skipi == 0:
                            # bycycle threshold_kwargs dict (see "Burst Properties" parameter
                            # block above for the meaning of each threshold).
                            burst_parameters = {
                                'amp_fraction_threshold': amp_fraction,
                                'amp_consistency_threshold': amp_consistency,
                                'period_consistency_threshold': period_consistency,
                                'monotonicity_threshold': monotonicity,
                                'min_n_cycles': mincycles}

                            print(f"Starting Bycycle Burst detection of {sujname}")
                            bm_list = []
                            # Burst detection is run independently per band (currently alpha-only,
                            # see bands_ranges_burst), per channel, per epoch/trial -- bycycle
                            # operates on one continuous 1-D signal at a time.
                            for band_idx in range(len(bands_ranges_burst)):
                                    print(f"{sujname} {band_names_burst[band_idx]} Detection in Progress: ")
                                    for ch, ch_name in enumerate(ch_names):
                                        prev_epochs = []
                                        for epoch in range(num_epochs):
                                            signal = raw._data[epoch,ch,:]
                                            bw = bands_ranges_burst[band_idx]
                                            # Core bycycle cycle-by-cycle feature extraction + burst
                                            # classification for this trial/channel/band using the
                                            # amp/period/monotonicity thresholds defined above.
                                            bm = bc.features.compute_features(signal, fs, f_range = bw, threshold_kwargs = burst_parameters)
                                            if treshold_features == 0:
                                                bm = pd.DataFrame(bm)[features_keep]
                                            else:
                                                bm = pd.DataFrame(bm)
                                            bm['epoch'] = epoch
                                            bm['ch']    = ch_name
                                            bm['chnum'] = ch
                                            bm['band']  = band_names_burst[band_idx]
                                            bm['block'] = new_event_names[epoch]
                                            # `epoch_block`: the epoch's ordinal position *within its own
                                            # condition block* (e.g., 3rd EyesClosed epoch), as opposed to
                                            # `epoch`, which is the ordinal position within the whole
                                            # recording. Needed to align
                                            # burst epochs to the block-restricted PSD arrays. Built
                                            # incrementally via a running tally (`prev_epochs`) of block
                                            # labels seen so far for this channel/trial loop.
                                            if (blocks == 0) or (epoch == 0):
                                                bm['epoch_block'] = epoch
                                                prev_epochs.append(new_event_names[epoch])
                                            else:
                                                prev_epochs.append(new_event_names[epoch])
                                                bm['epoch_block'] = prev_epochs.count(new_event_names[epoch])-1 # We need the position, and python starts with 0. Given that count gives you the actual number of epochs and not their position, removing one we obtain the position

                                            if blocks == 0:
                                                bm['epochs'] = num_epochs
                                            else:
                                                bm['epochs'] = new_event_names.count(new_event_names[epoch])

                                            bm_list.append(bm)

                            # Stack every (band x channel x trial) cycle table into one long-format
                            # burst-properties table and tag it with method/band provenance.
                            dfs_burst = pd.concat(bm_list, ignore_index=True)
                            dfs_burst['sujid'] = sujname
                            dfs_burst['methodburst'] = str(burst_parameters)
                            dfs_burst['banddef'] = str(band_definition_burst)


                            # bycycle's native 'period' column is in samples; rename to 'frequency'
                            # and convert samples -> Hz (cycle frequency = sampling rate / period).
                            dfs_burst.rename(columns={'period': 'frequency'}, inplace=True)
                            dfs_burst['frequency'] = dfs_burst['frequency'].apply(lambda x: fs / x)

                            if only_bursts:
                                dfs_burst = dfs_burst[dfs_burst['is_burst'] == True]

                            jrpc.save_dataframe(dfs_burst, subfolder_path, base_filename="burst_properties_bycycle", rewrite=rewritei, ext=".csv")

                    # ------------------------------------------------------------------
                    # MEASURE 4: Hilbert-based lagged auto-coherence (LAcH)
                    # (Methods: precursor to "alpha lifespan", cycle at 90% cumulative LAcH).
                    # ------------------------------------------------------------------
                    if laggedcoh_hil == 1:

                        # Resolve the surrogate-threshold suffix again here (mirrors the
                        # module-level block above) in case thresh_prctile_hilb/suf need to be
                        # re-derived at this point in participant-loop scope.
                        if lcohtype == 'hilb':
                            if thresh_prctile_hilb is None:
                                thresh_value = None
                            else:
                                thresh_value = thresh_prctile_hilb 

                            
                        basefile = basefile_lcoh + suf

                        lcohhilpath = os.path.join(subfolder_path, basefile + ".csv")
                        if not os.path.exists(lcohhilpath) or skipi == 0 or rewritei == 1:
                            print(f"Starting Lagged Coherence  Hilbert Based Computation of {sujname}. Take a seat, this could take a while")
                            df_lcoh = []
                            l = 1
                            if blocks == 0:
                                    block_events = {'NA'}

                            for blocki in block_events:
                                if blocks == 0:
                                    total_lcoh_data = raw._data
                                else:
                                    blockpos  = np.where(event_array == blocki)[0]
                                    if blockpos.size == 0 :
                                        continue  # No epochs for this block; skip.
                                    else:
                                        total_lcoh_data = raw._data[blockpos, :, :]

                                srate = raw.info['sfreq']

                                # lcoh_ar_total: per-trial LAcH values (epochs x channels x carrier-freqs x lags).
                                # lcoh_ar_mat:   epoch-averaged LAcH (channels x carrier-freqs x lags), what
                                #                gets exported as the long-format CSV below.
                                lcoh_ar_total = np.full((total_lcoh_data.shape[0], total_lcoh_data.shape[1], len(freq_rangeh), len(lag_rangeh)), np.nan) #epochs, channels, frequencies, lags
                                lcoh_ar_mat   = np.full((len(ch_names), len(freq_rangeh), len(lag_rangeh)), np.nan)

                                # Compute LAcH per channel, per trial (Hilbert-based autocoherence is
                                # inherently single-channel/single-trial), sweeping the full
                                # carrier-frequency x lag grid in one call per trial via `la.*`.
                                for ch_idx, ch in enumerate(ch_names):
                                        lcoh_data = total_lcoh_data[:,ch_idx,:]
                                        lcoh_data = np.squeeze(lcoh_data)
                                        epoch_data      = np.zeros((1,lcoh_data.shape[1]))
                                        for epoch in range(lcoh_data.shape[0]):
                                            epoch_data[0,:] = lcoh_data[epoch,:]

                                            if lcohtype == 'hilb':
                                                lcoh_ar_total[epoch,ch_idx,:,:]  = np.squeeze(la.lagged_hilbert_autocoherence(epoch_data, freq_rangeh, lag_rangeh, srate, type = plv_or_coh, thresh_prctile=thresh_value, n_shuffles=niter))
                                            elif lcohtype == 'fft':
                                                lcoh_ar_total[epoch,ch_idx,:,:]  = np.squeeze(la.lagged_fourier_autocoherence(epoch_data, freq_rangeh, lag_rangeh, srate, type = plv_or_coh))


                                            jrpc.update_progress(l, len(ch_names)*raw._data.shape[0], bar_length=20, message=f"{sujname} Lagged Coherence Progress: ")
                                            l = l + 1

                                        # epoch-average this channel's LAcH surface for the exported summary matrix.
                                        lcoh_ar_mat[ch_idx, :, :] = np.nanmean(np.squeeze(lcoh_ar_total[:,ch_idx,:,:]), axis = 0)


                                # Build a long-format (one row per electrode x freq x lag) table from the
                                # epoch-averaged 3-D matrix by taking every (I,J,K) index triple and mapping
                                # each axis back to its real-world label (channel name / Hz / lag-in-cycles).
                                I, J, K = np.indices(lcoh_ar_mat.shape)  # I: electrodes, J: frequency, K: lag

                                # Output filename encodes: method (hilb/fft) x metric (coh/plv) x block
                                # (if >1 block) x surrogate-threshold suffix -- so every parameter
                                # combination that could be run gets its own traceable output file.
                                if len(block_events) > 1:
                                    if lcohtype == 'hilb':
                                        if plv_or_coh == 'coh':
                                            lcohpath_mat_name = "lagged_coh_py_mat_hilb_" + blocki + suf + ".npy"
                                        else:
                                            lcohpath_mat_name = "lagged_coh_py_mat_hilb_plv_" + blocki +  suf + ".npy"
                                    else:
                                        if plv_or_coh == 'coh':
                                            lcohpath_mat_name = "lagged_coh_py_mat_fft_" + blocki + suf + ".npy"
                                        else:
                                            lcohpath_mat_name = "lagged_coh_py_mat_fft_plv_" + blocki + suf + ".npy"
                                else:
                                    if lcohtype == 'hilb':
                                        if plv_or_coh == 'coh':
                                            lcohpath_mat_name = "lagged_coh_py_mat_hilb" + suf + ".npy"
                                        else:
                                            lcohpath_mat_name = "lagged_coh_py_mat_hilb_plv" + suf + ".npy"
                                    else:
                                        if plv_or_coh == 'coh':
                                            lcohpath_mat_name = "lagged_coh_py_mat_fft" + suf + ".npy"
                                        else:
                                            lcohpath_mat_name = "lagged_coh_py_mat_fft_plv" + suf + ".npy"
                                    
                                lcohpath_mat = os.path.join(subfolder_path, lcohpath_mat_name)
                                # Persist the full per-trial LAcH tensor (not just the trial average)
                                # as a .npy matrix -- consumed directly by the companion
                                # basedonburst_v3 script to re-split LAcH by burst vs. no-burst epochs.
                                np.save(lcohpath_mat, lcoh_ar_total)
                                # Flatten the index arrays and map them to the actual labels/values:
                                df_long = pd.DataFrame({
                                                        'ch': np.array(raw.ch_names)[I.flatten()],
                                                        'freq': freq_rangeh[J.flatten()],
                                                        'cycle': lag_rangeh[K.flatten()],
                                                        'lcoh': lcoh_ar_mat.flatten()
                                                    })

                                df_long["sujid"]      = sujname
                                df_long['parameters'] = lcohhil_parameters
                                df_long['block']      = blocki
                                if blocks == 0:
                                    df_long['epochs']    = lcoh_data.shape[0]
                                else:
                                    df_long['epochs']    = new_event_names.count(blocki)

                                df_long['method']     = lcohtype
                                df_long['metric']     = plv_or_coh

                                if lcohtype == 'hilb':
                                    df_long['threshold'] = thresh_prctile_hilb
                                else:
                                    df_long['threshold'] = 'NA'


                                df_lcoh.append(df_long)

                            # Combine all blocks (or the single pseudo-block 'NA') into one CSV per participant.
                            df_lcoh_t = pd.concat(df_lcoh)
                            jrpc.save_dataframe(df_lcoh_t, subfolder_path, base_filename=basefile, rewrite=rewritei, ext=".csv")


                # Any ValueError raised anywhere in the try block above (e.g., malformed/
                # unreadable .set file, incompatible event structure) is caught here so one
                # participant's failure does not halt the whole batch; the error is logged
                # to `error_list` for post-hoc review instead of being silently swallowed.
                except ValueError as e:
                    print(f"Skipping {sujname} : {e}")
                    error_message = f"Participant {sujname} crashed because: {e}"
                    error_list.append(error_message)
                
           
