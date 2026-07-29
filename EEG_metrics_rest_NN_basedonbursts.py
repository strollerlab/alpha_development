"""
@author: J. Rico-Picó

SCRIPT ROLE (pipeline stage 1 of 2)
------------------------------------------------------------------------------
Burst-conditioned feature-extraction script for the Alpha Burst manuscript.
This script re-reads the same cleaned resting-state .set files as
`EEG_metrics_rest_v5.py`, but re-derives PSD-based measures split by burst
status (burst vs. no-burst epochs, per channel), using the epoch-level
`burst_properties_bycycle.csv` and the LAcH `.npy` matrix produced by that
first-stage script as inputs:
    1. Absolute/relative band power, computed separately for burst vs.
       no-burst epochs                                        [abspow]
    2. Aperiodic + oscillatory (specparam/FOOOF) parameterization, fit
       separately on burst vs. no-burst epochs                [aperosc]
    3. Lagged auto-coherence (LAcH), re-aggregated by burst vs. no-burst
       epoch (via jrpc.process_coherence_bursts)               [laggedcoh_hil]
There is no burst-*detection* step here (`bursbct` from stage 1 is not
repeated) -- burst status is read in from stage 1's output CSV
(`burst_filename = 'burst_properties_bycycle'`) and used purely as a
classifier to split epochs before computing the measures above. 
"""

# -----------------------------------------------------------------------------
# MANUSCRIPT NOMENCLATURE
#   specparam settings  -> "EEG Power-Spectrum Parameterization" (2-15 Hz, fixed
#                          aperiodic, peak_threshold=2, min_peak_height=0.01,
#                          max_n_peaks=7, peak_width_limits=[1,10]).
#   LAcH / lagged-coherence
# -----------------------------------------------------------------------------

"""
Interest Variables: determine the parameters and measures that you want to compute in this script
Measures:  Relative and Absolute Power (log or not logtrans), Aperiodic and Oscillatory, Burst, Lagged Coherence
Notice that in absolute and relative power it is necessary to compute previously the absolute and relative power
"""

# ---------------------------------------------------------------------------
# CONFIGURATION — the only lines to edit (stage 2 (burst-conditioned))
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

# Third-party dependencies required by this script (spectral / MNE-Python stack;
packages = ['numpy', 'scipy', 'pandas', 'matplotlib', 'neurodsp', 'tqdm', 'pytest', 'mne', 'specparam']

# Ensure all packages are installed
# Dependency check: pip-installs any missing package before
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

# Local helper module (not on PyPI): PSD/FOOOF convenience wrappers, including
# the channel-based FOOOF fitting and burst/coherence-splitting helpers used
# only by this stage-2 script (interim_fooof_epochselection_chbased,
# fooofpsd_extraction_chbased, process_coherence_bursts).
sys.path.append(CODE_DIR)   # convenience_functions_jrp.py / lagged_autocoherence.py
import convenience_functions_jrp as jrpc

###############################################################################
# Path and basic parameters ###################################################
###############################################################################

# Specify folder containing the .set files
# (EEGLAB-format epoched, cleaned resting-state recordings; must be the SAME
#  cohort/visit path used for the corresponding stage-1 run, since this script
#  re-reads that run's `burst_properties_bycycle.csv` / LAcH .npy outputs.)

folder_path = SET_DIR   # must be the SAME visit folder used for the stage-1 run

# Where do you want to save the files
# Destination root; must match the stage-1 output_path for this
# cohort/visit so the burst-properties CSV and LAcH .npy files are found.
output_path  = os.path.join(DATA_DIR)

# Set files in the folder, modify how the file ends
# Only files ending in "rest.set" (i.e., the resting-state condition) are queued for processing.
set_files = [f for f in os.listdir(folder_path) if f.endswith('rest.set')] # MODIFY

# ---- Run-control / batching flags -------------------------------------------------
skipi      = 0   # 1 = skip participants whose expected output files already exist (resumable batch runs); 0 = always reprocess
upsidedown = 1   # 1 = reverse the file processing order (useful for running two machines from opposite ends of the list concurrently)
multicomp  = 0   # 1 = split set_files across parallel jobs/machines using the quartile logic below; 0 = process the full list
stepm      = 1   # Which half/quarter of the split this job should take when multicomp == 1 (1 or 2; see block below)
rewritei   = 1   # Passed through to jrpc.save_dataframe: 1 = overwrite existing output CSVs, 0 = do not overwrite
numchar    = 8   # Number of leading characters of the raw filename used as sujid when rename == 0
nversion   = 0   # Which is the version you are running? If set to 0 no version added, otherwise vZ to skip files

# ---- Participant-ID renaming (maps raw filenames manuscript subject IDs; MUST
#      match the mapping used in the stage-1 run for the same cohort/visit) --------
rename     = 0        # 1 = build sujid as prename + substring(namepos1:namepos2) + posname; 0 = use the first `numchar` characters of the filename as-is
namepos1   = 1
namepos2   = 4
prename    = 'BUDDY'
posname    = ''

# Reverse processing order (see `upsidedown` above).
if upsidedown == 1 :
   set_files = set_files[::-1]

# Multi-machine/parallel batching: carve the file list into a subset for this run
# (see EEG_metrics_rest_v5.py for the full explanation of this splitting logic;
# it does not affect which participants are eventually processed, only which
# machine handles which subset in a given run).
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
# Parameters to extract: 0 = No compute, 1 = Compute (for aperosc, abspow is necessary because we need epoch x electrode psds)
###############################################################################
abspow         = 1   # Burst-conditioned absolute/relative band power (Methods: "Power Spectrum Density", split by burst status)
aperosc        = 1   # Burst-conditioned aperiodic + periodic (specparam/FOOOF) parameterization
laggedcoh_hil  = 1   # Burst-conditioned Hilbert-based lagged auto-coherence (LAcH)
###############################################################################
# Signal per blocks?
################################################################################
blocks = 0   # 1 = split epochs by condition block (e.g., eyes-closed/eyes-open) using block_events below; 0 = treat all epochs as one pool
block_events  = {'ECrs', 'EOrs'} #Events in the EEG
block_names   = {'ECrs': 'EyesClosed', #Equivalence for writing
                 'EOrs': 'EyesOpen'}

filter_event    = 0        # 1 = drop epochs whose event name matches `event_to_filter` before any computation
event_to_filter = 'SOCV'   # Substring of event names to exclude (e.g., non-resting-state / task markers)

################################################################################
# Freq Range, Notice min should be at least twice higher than the s in the epoch
################################################################################
theta = [3, 5]
alpha = [6, 9]
beta  = [12, 20]
gamma = [21, 45]

# Burst-detection band actually used in stage-1 (documented here so stage-2 can
# check independence). Set this to whatever stage-1's `burst_alpha` was.
burst_detect_band = [7, 10]
if list(burst_detect_band) == list(alpha):
    print("NOTE [PB02]: burst-detection band {} equals the stage-2 power alpha "
          "band {}. Allowed, but the burst-conditioned power contrast is then "
          "not independent of the detection band.".format(burst_detect_band, alpha))

# Only computed/needed if abspow or aperosc is requested - MODIFY
if abspow == 1 or aperosc == 1: #We save this information compute this information if necessary
    fmin, fmax     = 2, 45     # PSD frequency window retained for band-power/FOOOF (Hz)
    nranges        = 4         # Number of canonical bands defined below (theta/alpha/beta/gamma)
    range_theta    = theta
    range_alpha    = alpha
    range_beta     = beta
    range_gamma    = gamma

    # Integer-index lookup (used when iterating band_names/bands_ranges by position)...
    band_names = {
        0: 'theta',
        1: 'alpha',
        2: 'beta',
        3: 'gamma'}

    # The burst-detection band(s) to condition on -- must match (a subset of)
    # the bands present in the upstream burst_properties_bycycle.csv `band`
    # column, since that column is filtered against these names below.
    band_names_burst = {
        0: "alpha"}

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
    n_perseg        = 2000  #Welch and spectopo number of timepoints to run the fft (e.g., 1000 with 1000 srate = 1s)
    n_overlap       = 1000  #Welch and spectopo points of overlapping (e.g., in 2000 timewindow 1000 will be 50% overlap)

##############################################################################
# Define the Parameters and Methods to Compare - DON'T MODIFY
###############################################################################

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
    # NOTE: must match the PSD method actually used in the stage-1 run
    # (EEG_metrics_rest_v5.py) if burst-conditioned and unconditioned PSD
    # values are to be compared downstream.
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
 # Aperiodic computation parameters
 ##############################################################################
 # specparam/FOOOF settings for the burst-conditioned re-fit ("EEG Power-Spectrum
 # Parameterization" applied separately to burst vs. no-burst epoch subsets).
 # NOTE: several values intentionally differ from the stage-1 unconditioned
 # fit below, since burst/no-burst subsets contain fewer epochs per participant

if aperosc == 1:
    psd_fmin          = 2       # Lower fit bound (Hz)
    peak_width_limits = [1, 12] # Allowed oscillatory-peak bandwidth (Hz); narrower floor than stage-1's [2, 12]
    max_n_peaks       = 5       # Max number of periodic peaks specparam may fit
    min_peak_height   = 0.01    # Minimum peak height (above the aperiodic fit) to be retained
    peak_threshold     = 2      # Peak selection threshold (in SD of the flattened spectrum)
    rsquare_threshold = 0.90    # Final-fit R^2 acceptance threshold (goodness of aperiodic+periodic fit)
    aperiodic_mode    = 'fixed' # Aperiodic component form: 'fixed' (no spectral knee) per Methods


    ## If you are performing a reduction by excluding the bad fit epochs
    # epoch-selection / pre-fitting strategy used to guard against noisy single-epoch
    # spectra biasing the aperiodic fit (see jrpc.interim_fooof_epochselection_chbased,
    # the channel-based variant used in this burst-conditioned script):
    fit_pregroups             = 1    # 0 = fit the mean PSD only; 1 = pre-fit groups of epochs and discard poorly-fit ones (used here); 2 = fit every epoch individually
    rsquare_threshold_interim = .90  # R^2 cutoff applied during the interim (pre-group) fit used to flag/discard bad epochs
    perelect                  = 0.75 # Minimum fraction of electrodes that must pass rsquare_threshold_interim for a epoch-group to be retained
    epochs_multi              = 1    # 1 = iterate over every group size in `nepochs_group_multi`; 0 = use the single `nepochs_group_prec` value instead
    nepochs_group_multi       = [9999]  # epoch-group size(s) tested; 9999 is a value meaning "no grouping / fit all (burst or no-burst) epochs together"
    nepochs_group_prec        = 5      # Fallback single group size when epochs_multi == 0
    psd_fmaxs                 = [15, 45]  # Upper fit-frequency bounds compared across runs

 ##############################################################################
 # Burst Properties
 ##############################################################################


burst_filename = 'burst_properties_bycycle'  # Stem of stage-1's per-participant bycycle output CSV (must match jrpc.save_dataframe's base_filename there)

###############################################################################
# Lagged hilbert autocoherence
###############################################################################
# Hilbert-based lagged auto-coherence (LAcH), re-split by burst vs. no-burst
# epoch status. Sweep settings (frequency/lag ranges) mirror stage-1's; what
# differs here is that the *input* being coherence-analyzed is the stage-1
# LAcH matrix subset by burst status (see main loop / `matfile_lcoh` below),
# not a Hilbert computation.
if laggedcoh_hil == 1 :
    lcohtype   = 'hilb'  # 'hilb' = Hilbert-transform based LAcH (must match stage-1's lcohtype for the on-disk matrix to be readable)
    plv_or_coh = 'coh'   # Metric flavor: 'coh' = coherence, otherwise phase-locking value (PLV)
    fminch     = 3       # Min carrier frequency for the LAcH frequency sweep (Hz)
    fmaxch     = 12      # Max carrier frequency for the LAcH frequency sweep (Hz)
    freqfacch  = .5      # Frequency step of the sweep (Hz)

    cminh      = 0.05    # Minimum lag, in oscillation cycles
    cmaxh      = 4       # Maximum lag, in oscillation cycles
    cfach      = .05     # Lag step, in oscillation cycles

    thresh_prctile_hilb  =  95# During the subrogation, which is the threshold (None = no subrrogate)

    fft_rangeh = np.arange(fminch, fmaxch, freqfacch)  # Carrier-frequency sweep values actually used
    lag_rangeh = np.arange(cminh, cmaxh, cfach)         # Lag (cycle) sweep values actually used

    lcohhil_parameters = {
        'ncyclesmin': cminh,
        'ncyclesmax':cmaxh,
        'cyclestep': cfach,
        'freqs': [fminch, fmaxch]}

    # Two output-name stems are built here (stage-1 only needed one):
    #   basefile_lcoh -> stem for the burst-conditioned long-format LAcH CSV
    #                    written by THIS script (see jrpc.save_dataframe call
    #                    in main loop); gets the '_burst' suffix appended below
    #                    so it is never confused with stage-1's unconditioned CSV.
    #   matfile_lcoh  -> stem used to LOCATE stage-1's raw per-block LAcH
    #                    `.npy` matrices on disk (input to this script, read
    #                    via np.load in the main loop) -- must exactly match
    #                    the filename stage-1 wrote them under.
    basefile_lcoh = "lagged_coh_py_" + lcohtype + '_' + plv_or_coh
    matfile_lcoh  = "lagged_coh_py_mat_" + lcohtype + '_' + plv_or_coh

    # ==============================================================================
    if lcohtype == 'hilb':
        if thresh_prctile_hilb is None:
            niter = 1
            suf   = '_nothresh'
        elif thresh_prctile_hilb < 95:
            niter = 1000
            suf   = "_lowthresh"
        else:
            niter = 1000
            suf   = ""
    else:
        suf = ""

    # Append the surrogate-threshold suffix to both stems; `basefile_lcoh`
    # additionally gets '_burst' so the output CSV name self-documents that
    # it holds burst-conditioned (not raw) LAcH values, while `matfile_lcoh`
    # deliberately does NOT get '_burst' since it must still resolve to
    # stage-1's original (unconditioned) .npy filename on disk.
    basefile_lcoh = basefile_lcoh + suf + '_burst'
    matfile_lcoh  = matfile_lcoh + suf
else:
    basefile_lcoh = ""

#########################################
# Function methods and deffinition Methods Function Definitions
#########################################

def files2process(setfile, output_path, abspow, aperosc, laggedcoh_hil,
                  laggedcoh_name,
                  extension = '.csv', numversion = 0):
            """
            Check, for one participant, which of the requested burst-conditioned
            output CSVs already exist on disk under `output_path/<sujname>/`.

            Used to implement resumable batch runs (`skipi == 1`, see main
            loop): if all requested measures already have an output file, the
            participant is skipped entirely.

            NOTE: like stage-1's files2process(), this relies on the
            module-level `sujname` variable set in the enclosing loop rather
            than the `setfile` argument (pre-existing behavior, left
            unchanged). Signature also differs from stage-1: no `bursbct`
            flag (this script performs no burst detection) and a
            `laggedcoh_name` parameter (expects the burst-suffixed
            `basefile_lcoh` stem) in place of stage-1's `basefile_lcoh`.

            Parameters
            ----------
            setfile : str
                Path to the participant's .set file (unused directly; see
                NOTE above).
            output_path : str
                Root output directory containing one subfolder per participant.
            abspow, aperosc, laggedcoh_hil : int (0 or 1)
                Which burst-conditioned measures were requested for this run.
            laggedcoh_name : str
                Base filename stem for the burst-conditioned LAcH CSV
                (i.e. `basefile_lcoh`, already including the '_burst' suffix).
            extension : str, default '.csv'
                Output file extension to check for.
            numversion : int, default 0
                Optional version suffix ('_v{numversion}') appended to output
                filenames; 0 means no suffix.

            Returns
            -------
            (totaldone, totalfiles) : tuple[int, int]
                Count of requested outputs already present vs. total requested;
                `totaldone == totalfiles` signals this participant can be skipped.
            """
            subfolder_path  = os.path.join(output_path, sujname)
            if numversion == 0:
                version = ''
            else:
                version = '_v' + str(numversion)
                 
            print(f"Checking what files the participant {sujname} has")
            totalfiles = abspow + aperosc + laggedcoh_hil
            totaldone  =  0
            if abspow  == 1:
                absfile    = 'psds_burst' + version + extension
                abspowpath = os.path.join(subfolder_path, absfile)
                abs_dummy  = os.path.exists(abspowpath)
                totaldone  = totaldone + abs_dummy
                print(f"Power Spectrum Computation: {abs_dummy}")
            if aperosc == 1:
                aperfile = 'aperosc_parameters_burst' + version + extension
                ape_dummy = os.path.exists(os.path.join(subfolder_path, aperfile))
                print(f"Power Spectrum Parametrization: {ape_dummy}")
                totaldone = totaldone + ape_dummy
            if laggedcoh_hil == 1 :
                lagfilehil = laggedcoh_name + version + extension
                lcohpath  = os.path.join(subfolder_path, lagfilehil)
                lco_dummy = os.path.exists(lcohpath)
                print(f"Lagged Coherence: {lco_dummy}")
                totaldone = totaldone + lco_dummy 
            return totaldone, totalfiles
        
 ## Starts the loops
 # Main per-participant processing loop (stage 2): resolve the subject ID,
 # decide whether to skip (already-processed) participants, re-load the same
 # cleaned EEGLAB epochs used in stage 1, load stage-1's burst-classification
 # CSV, then dispatch to each requested measure's burst-conditioned
 # computation block (PSD/abspow, aperiodic/FOOOF, LAcH).

errori     = 0
error_list = []
for s, set_file in enumerate(set_files):

            # Resolve the manuscript-facing participant ID (sujid) from the raw filename.
            if rename == 0:
                   sujname = set_file[0:numchar]
            else:
                   sujname = prename + set_file[namepos1:namepos2] + posname

            # Resumability check: skip participants who already have every
            # requested burst-conditioned output file on disk (see
            # files2process docstring above).
            filesdone, totalfiles = files2process(sujname, output_path, abspow=abspow,
                                     aperosc=aperosc,
                                     laggedcoh_hil = laggedcoh_hil,
                                     laggedcoh_name = basefile_lcoh,
                                     numversion = nversion,
                                     extension = '.csv')


            if totalfiles-filesdone == 0 and skipi == 1:
                print(f"Participant {s} - {sujname} - was already processed: Skipping")
            else:
                print(f"Participant {s} - {sujname}: Starting processing")

                    # Function to read raw EEG data
                try:
                    filepath = os.path.join(folder_path, set_file)
                    subfolder_path = []
                    subfolder_path = os.path.join(output_path, sujname)
                    os.makedirs(subfolder_path, exist_ok=True)

                    # Corrected function to read raw EEG data
                    # Load the same cleaned, epoched EEGLAB (.set) recording used in stage 1
                    # (epoch order/count must match so the burst-status CSV rows align).
                    raw = mne.io.read_epochs_eeglab(filepath)
                    raw.load_data()

                    # Optionally drop epochs whose event label matches `event_to_filter`
                    # (e.g., non-resting-state task markers accidentally retained upstream).
                    if filter_event:
                        events_to_keep = [key for key in raw.event_id.keys() if event_to_filter not in key]
                        raw            = raw[events_to_keep]

                    ch_names       = raw.ch_names
                    electrodes     = len(raw.ch_names)
                    electrode_names=raw.ch_names
                    num_epochs     = len(raw)
                    raw._data     *= 1000000  # Convert MNE's native Volts to microvolts (uV), the unit used throughout Methods/figures.
                    data           = raw._data
                    fs             = raw.info['sfreq']

                    print(f"Starting the processing of {sujname}")

                    # Load stage-1's per-epoch, per-band bycycle burst classification
                    # (columns include epoch/epoch index, `band`, and burst-run summary
                    # stats); used below purely as a boolean epoch classifier -- no
                    # burst detection is (re-)run here.
                    subburstpath    = os.path.join(subfolder_path, burst_filename + '.csv')
                    burstdata_total = pd.read_csv(subburstpath)

                    # Derive a per-epoch condition-block label ('EyesClosed'/'EyesOpen', etc.)
                    # so downstream measures can be computed overall and/or split by block.
                    if blocks == 1:
                        event_ids   = raw.events[:, 2]  # Third column contains the event IDs
                        event_names = [key for event_id in event_ids for key, value in raw.event_id.items() if value == event_id]
                        new_event_names = []
                        for name in event_names:
                                normalized_name = name  # Default to the original event name.
                                for block in block_events:
                                    if block in name:  # If the block string is a substring of the event name...
                                        normalized_name = block  # ...set it to that block.
                                        break  # Stop checking once a match is found.
                                new_event_names.append(normalized_name)
                    else: 
                        event_ids   = raw.events[:, 2]  # Third column contains the event IDs
                        new_event_names = [key for event_id in event_ids for key, value in raw.event_id.items() if value == event_id]
                             
                    if abspow == 1 or aperosc == 1:
                        # abspow and aperosc share the same underlying PSD computation
                        # (psds/relpsds below); only the downstream aggregation differs
                        # (band-summed power here vs. specparam fit further down).
                        powspcpath   = os.path.join(subfolder_path, "psds_burst.csv")
                        aperosc_path = os.path.join(subfolder_path, "aperosc_parameters_burst.csv")
                        if not os.path.exists(powspcpath) or skipi == 0 or aperosc == 1:
                            print(f"Starting Absolute/Relative Power computation of {sujname}")

                            # Process PSDs: calculate power spectra (psds) and frequencies.
                            freqs, psds = jrpc.process_psd(raw, methodpsd, raw.info['sfreq'], **extra_params)  # psds shape: epochs x channels x freqs
                            freq_mask = (freqs >= fmin) & (freqs <= fmax)
                            foi = np.where(freq_mask)[0]
                            psds = psds[:, :, foi]  # restrict frequencies of interest
                            psds_sum = np.nansum(psds, axis=2, keepdims=True)
                            relpsds  = psds / psds_sum
                            freqs = freqs[foi]
                            if abspow == 1:
                                # Prepare a list to store epoch-level output.
                                psds_bands_data = []
                                # Outer loop: the burst-DETECTION band(s) whose classification
                                # we're conditioning on (`band_names_burst`, alpha-only in this
                                # manuscript -- see "Burst Properties" parameter block above).
                                for bi_b, band_b in enumerate(band_names_burst):
                                    # Restrict stage-1's burst CSV to rows for this detection band,
                                    # then further restrict to only the (channel, epoch) pairs that
                                    # bycycle actually flagged as a burst (`is_burst == True`); the
                                    # remaining rows define the "burst" epoch/channel set below.
                                    burst_band_df = burstdata_total[burstdata_total['band'] == band_names_burst[bi_b]]
                                    burst_band_df = burst_band_df[burst_band_df['is_burst'] == True]

                                    # Inner loop: the POWER band(s) we're computing PSD/aperiodic
                                    # measures for (`band_names`/`bands_ranges`, e.g. theta/alpha/beta/
                                    # gamma); independent of, and possibly different from, the burst-
                                    # detection band(s) looped above.
                                    for bi, band in enumerate(band_names):
                                        band_range = bands_ranges[bi]
                                        # Filter the burst data for the current band.


                                        # Group by channel and remove duplicate epochs.
                                        #  per-channel dict of the *set* of epoch indices classified
                                        # as "burst" for the current detection band (used as the
                                        # burst/no-burst mask below).
                                        burst_band = (
                                            burst_band_df.groupby('ch')['epoch']
                                            .apply(lambda x: list(set(x)))
                                            .to_dict()
                                        )

                                        # Next, sum PSD values over the frequencies of interest for this band.
                                        freq_mask_band = (freqs >= band_range[0]) & (freqs <= band_range[1])
                                        foi_band       = np.where(freq_mask_band)[0]
                                        # For absolute power.
                                        psd_band       = {band_names[bi]: np.nansum(psds[:, :, foi_band], axis=2)}
                                        # For relative power.
                                        psd_band_rel  = {band_names[bi]: np.nansum(relpsds[:, :, foi_band], axis=2)}

                                        # Loop over each channel.
                                        for ch_idx, ch_name in enumerate(ch_names):
                                            # Initialize an array of zeros for the burst status of each epoch.
                                            ch_burst = np.zeros(psds.shape[0])

                                            try:
                                                # 1. Retrieve the indices and ensure they are a numpy array
                                                indices = np.array(burst_band[ch_name])

                                                # 2. Filter out indices that are out of bounds (larger than the array size)
                                                valid_indices = indices[(indices >= 0) & (indices < len(ch_burst))]

                                                # 3. Assign 1 only to the valid indices
                                                ch_burst[valid_indices] = 1

                                            except KeyError:
                                                # If the channel is not found in burst data, leave as zeros
                                                pass
                                            except Exception as e:
                                                # (Optional) Catch any other unexpected errors to prevent crash
                                                print(f"Warning: Could not process burst for {ch_name}: {e}")

                                            # Compute counts per condition for the current channel:
                                            ntr_burst   = np.sum(ch_burst == 1)
                                            ntr_noburst = np.sum(ch_burst == 0)

                                            ntr_noburst = ntr_noburst - sum(np.isnan(psds[:,ch_idx,0])) #Designed for low-density montage correction. In regular EEG with signal in all electrodes x epoch this does not change a thing

                                            # Loop over each epoch (epoch).
                                            for epoch_idx in range(psds.shape[0]):
                                                if ch_burst[epoch_idx] == 1:
                                                    condition = 'burst'
                                                    count = ntr_burst  # total burst epochs for this electrode
                                                else:
                                                    condition = 'noburst'
                                                    count = ntr_noburst  # total non-burst epochs for this electrode

                                                if psd_band[band_names[bi]][epoch_idx, ch_idx] < 0.00001 :
                                                    abspowval = np.nan
                                                    relpowval = np.nan
                                                else:
                                                    abspowval = psd_band[band_names[bi]][epoch_idx, ch_idx]
                                                    relpowval = psd_band_rel[band_names[bi]][epoch_idx, ch_idx]

                                                # Create a dictionary with the epoch-level information.
                                                # Column names are the data contract consumed by the
                                                # downstream R pipeline; `burst`/`burstband` are the two
                                                # columns added relative to stage-1's (unconditioned)
                                                # epoch_psdpow.csv equivalent.
                                                rowpow = {
                                                    'Epoch': epoch_idx + 1,
                                                    'ch': ch_name,
                                                    'sujid': sujname,
                                                    'psdmethod': str(psd_parameters),
                                                    'banddef': str(band_definition),
                                                    'pow': 'abs',
                                                    'block': new_event_names[epoch_idx],
                                                    'band': band_names[bi],
                                                    'burst': condition,           # 'burst' or 'noburst' classification for this (channel, epoch)
                                                    'burstband': band_names_burst[bi_b],  # which detection band this classification came from
                                                    'Abspow': abspowval,
                                                    'Relpow': relpowval,
                                                    'epochs': count  # number of valid epochs for that electrode & condition
                                                }
                                                psds_bands_data.append(rowpow)
                                
                                # Convert the list of dictionaries into a DataFrame.
                                df_abs = pd.DataFrame(psds_bands_data)
                                if abspow == 1 or aperosc == 1:
                                    jrpc.save_dataframe(df_abs, subfolder_path, base_filename="epoch_psdpow_burst", rewrite=rewritei, ext=".csv")
                        
                                # If no blocks, compute additional, frequency-resolved averages.
                                if blocks == 0:
                                    avg_df_abs = df_abs.pivot_table(
                                        index=['ch', 'band', 'sujid', 'psdmethod', 'banddef', 'burst', 'burstband'],
                                        values=['Abspow', 'Relpow', 'epochs']
                                    ).reset_index()

                                    # Calculate frequency-resolved averages per band, channel, and condition.
                                    # This block re-derives the same burst/no-burst epoch classification as
                                    # above (independently, from `epoch_block` rather than `epoch`) in order
                                    # to average the FULL frequency-resolved PSD/relPSD spectrum (not just the
                                    # band-summed scalar), for spectral-shape plots/comparisons in the manuscript.
                                    df_list = []

                                    for bi_b, band_b in enumerate(band_names_burst):
                                        burst_band = burstdata_total[burstdata_total['band'] == band_names_burst[bi_b]]
                                        burst_band = burst_band[burst_band['is_burst'] == True]

                                        # BUGFIX (PB01): index by 'epoch', not 'epoch_block'.
                                        # This branch is the `blocks == 0` case: `psds`/`relpsds` below are
                                        # the FULL-recording arrays and are NOT sliced to a block, so the
                                        # mask must be built from the epoch's ordinal position within the
                                        # whole recording ('epoch'). 'epoch_block' is the position WITHIN a
                                        # condition block and restarts at 0 for every block, so using it
                                        # against a full-recording array silently mislabels burst/no-burst
                                        # epochs whenever more than one block is present.
                                        # NB: the two sibling branches are already correct and differ from
                                        # this one exactly here -- the block-sliced PSD branch (`blocks != 0`)
                                        # legitimately uses 'epoch_block' *after* filtering on ['block'],
                                        # and the aperosc branch guards with `if blocks > 0: ...['block']`.
                                        # OLD: burst_band.groupby('ch')['epoch_block']
                                        burst_band = (
                                            burst_band.groupby('ch')['epoch']
                                            .apply(lambda x: list(set(x)))  # Remove duplicates using `set`
                                            .to_dict())

                                        # ======================================================================
                                        burst_names = ['burst', 'noburst']

                                        for  bi, burst in enumerate(burst_names):
                                            for ch_idx, ch in enumerate(ch_names):
                                                            ch_data = psds[:, ch_idx, :]

                                                            try:
                                                                # Create a boolean mask indicating burst epochs for the current channel.
                                                                mask = np.array([epoch in burst_band[ch] for epoch in range(ch_data.shape[0])])

                                                                # For the non-burst condition, invert the mask so that "True" indicates epochs without bursts.
                                                                if burst != 'burst':
                                                                    mask = ~mask
                                                                    # If no epoch remains (i.e. all epochs are bursts), warn and skip this channel.
                                                                    if np.sum(mask) == 0:
                                                                        print(f"Warning: Participant {sujname} did not have any epochs without bursts, skipping electrode {ch}")
                                                                        continue
                                                            except KeyError:
                                                                # If the channel is not even in burstdata, handle accordingly.
                                                                print(f"Warning: Participant {sujname} did not have any bursts in electrode {ch}: all epochs selected for no burst condition")
                                                                if burst == 'burst':
                                                                    # For burst condition, skip the channel if no burst info exists.
                                                                    continue
                                                                else:
                                                                    mask = np.ones(ch_data.shape[0], dtype=bool)  # otherwise, select all epochs

                                                            # Create a boolean mask for epochs that have valid (non-NaN) entries in the first column.
                                                            valid_mask = ~np.isnan(ch_data[:, 0])

                                                            # Combine the burst mask with the valid data mask.
                                                            final_mask = mask & valid_mask

                                                            # Average the frequency-resolved PSD (abs/rel) across all epochs
                                                            # satisfying `final_mask` (burst-or-noburst AND non-NaN) for this channel.
                                                            ch_data_burst     = np.nanmean(ch_data[final_mask,:], axis = 0)
                                                            ch_data_burst_rel = np.nanmean(relpsds[final_mask, ch_idx, :], axis = 0)
                                                            
                                                            data_temp = pd.DataFrame({'ch': [ch] * len(freqs),
                                                                    'freq': freqs,
                                                                    'sujid': [sujname] * len(freqs),
                                                                    'psdmethod': [method_names[methodpsd]] * len(freqs),
                                                                    'abs': ch_data_burst,
                                                                    'rel': ch_data_burst_rel,
                                                                    'block': [None] * len(freqs),
                                                                    'burst': [burst] * len(freqs),
                                                                    'burstband': [band_b] * len(freqs)})
                                                            df_list.append(data_temp)
    
                                    # After all loops are finished, concatenate the list of DataFrames into one.
                                    if df_list:
                                        df_psddata = pd.concat(df_list, ignore_index=True)
                                    else:
                                        df_psddata = pd.DataFrame()
                                else:
                                    # Blocks case: also average separately by block. Identical
                                    # burst/no-burst masking logic to the `blocks == 0` branch
                                    # above, just additionally sliced to the epochs of the
                                    # current `blocki` (e.g. 'EyesClosed') before masking.
                                    event_array = np.array(new_event_names)
                                    avg_df_abs = df_abs.pivot_table(
                                        index=['ch', 'block', 'band', 'sujid', 'psdmethod', 'banddef', 'burst', 'burstband'],
                                        values=['Abspow', 'Relpow', 'epochs']
                                    ).reset_index()


                                    df_list   = []

                                    for blocki in block_events:
                                        blockpos = np.where(event_array == blocki)[0]

                                        if blockpos.size == 0 :
                                            continue
                                        else:
                                            tpsds     = psds[blockpos, :, :]
                                            tpsds_rel = relpsds[blockpos, :, :]

                                            # PB04 (see first occurrence): use a list, not a set, for deterministic order.
                                            # OLD: burst_names = {'burst', 'noburst'}
                                            burst_names = ['burst', 'noburst']
                                            for bi_b, band_b in enumerate(band_names_burst):
                                                burst_band = burstdata_total[burstdata_total['band'] == band_names_burst[bi_b]]
                                                burst_band = burst_band[burst_band['is_burst'] == True]
                                                burst_band = burst_band[burst_band['block']==blocki]  # further restrict to epochs from this block
                                                            
                                                burst_band = (
                                                    burst_band.groupby('ch')['epoch_block']
                                                    .apply(lambda x: list(set(x)))  # Remove duplicates using `set`
                                                    .to_dict())
                                            
                                            
                                                for  bi, burst in enumerate(burst_names):                       
                                                    for ch_idx, ch in enumerate(ch_names):
                                                                    ch_data = tpsds[:, ch_idx, :]
                                                                    try:
                                                                        # Create a boolean mask indicating burst epochs for the current channel.
                                                                        mask = np.array([epoch in burst_band[ch] for epoch in range(ch_data.shape[0])])
                                                                        
                                                                        # For the non-burst condition, invert the mask so that "True" indicates epochs without bursts.
                                                                        if burst != 'burst':
                                                                            mask = ~mask
                                                                            # If no epoch remains (i.e. all epochs are bursts), warn and skip this channel.
                                                                            if np.sum(mask) == 0:
                                                                                print(f"Warning: Participant {sujname} did not have any epochs without bursts and block {blocki}, skipping electrode {ch}")
                                                                                continue
                                                                    except KeyError:
                                                                        # If the channel is not even in burstdata, handle accordingly.
                                                                        print(f"Warning: Participant {sujname} did not have any bursts in electrode {ch} and block {blocki}: all epochs selected for no burst condition")
                                                                        if burst == 'burst':
                                                                            # For burst condition, skip the channel if no burst info exists.
                                                                            continue
                                                                        else:
                                                                            mask = np.ones(ch_data.shape[0], dtype=bool)  # otherwise, select all epochs
                                                                    
                                                                    # Create a boolean mask for epochs that have valid (non-NaN) entries in the first column.
                                                                    valid_mask = ~np.isnan(ch_data[:, 0])
                                                                    
                                                                    # Combine the burst mask with the valid data mask.
                                                                    final_mask = mask & valid_mask    
                                                                
                                                                    ch_data_burst     = np.nanmean(ch_data[final_mask,:], axis = 0)
                                                                    ch_data_burst_rel = np.nanmean(tpsds_rel[final_mask, ch_idx, :], axis = 0)
                                                                    
                                                                    data_temp = pd.DataFrame({'ch': [ch] * len(freqs),
                                                                            'freq': freqs,
                                                                            'sujid': [sujname] * len(freqs),
                                                                            'psdmethod': [method_names[methodpsd]] * len(freqs),
                                                                            'abs': ch_data_burst,
                                                                            'rel': ch_data_burst_rel,
                                                                            'block': blocki * len(freqs),
                                                                            'burst': [burst] * len(freqs),
                                                                            'burstband': [band_b] * len(freqs)})
                                                                    df_list.append(data_temp)
                                                                    
                                    if df_list:
                                        df_psddata = pd.concat(df_list, ignore_index=True)
                                    else:
                                        df_psddata = pd.DataFrame()
    
                                # Finally, save the new PSD DataFrame and the averaged epoch-based (pivot) DataFrame.
                                if abspow ==    1:
                                    jrpc.save_dataframe(df_psddata, subfolder_path, base_filename="psds_burst", rewrite=rewritei, ext=".csv")
                                    jrpc.save_dataframe(avg_df_abs, subfolder_path, base_filename="avg_pow_band_burst", rewrite=rewritei, ext=".csv")

                    if aperosc == 1:
                        aperosc_path = os.path.join(subfolder_path, "aperosc_parameters_burst.csv")

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
                                        tpsds = psds

                                        if blocks != 0:
                                            blockpos = np.where(event_array == blocki)[0]
                                            if blockpos.size == 0 :
                                                continue  # No epochs for this block for this participant; skip.
                                            else:
                                                tpsds    = psds[blockpos, :, :]

                                        # Burst/no-burst epoch classification for this detection band
                                        # (and, if blocked, this block) -- same pattern as the abspow
                                        # block above, but here the mask is applied per-channel BEFORE
                                        # the FOOOF fit (via the channel-based `_chbased` jrpc variants),
                                        # rather than after computing a channel-wide PSD summary.
                                        # PB04 (see first occurrence): use a list, not a set, for deterministic order.
                                        # OLD: burst_names = {'burst', 'noburst'}
                                        burst_names = ['burst', 'noburst']
                                        for band in range(len(band_names_burst)):
                                                burst_band = burstdata_total[burstdata_total['band'] == band_names_burst[band]]
                                                burst_band = burst_band[burst_band['is_burst'] == True]

                                                if blocks > 0:
                                                    burst_band = burst_band[burst_band['block'] == blocki]

                                                burst_band = (
                                                            burst_band.groupby('ch')['epoch_block']
                                                            .apply(lambda x: list(set(x)))  # Remove duplicates using `set`
                                                            .to_dict())

                                                for  bi, burst in enumerate(burst_names):
                                                    for ch_idx, ch in enumerate(ch_names):
                                                        ch_data = tpsds[:, ch_idx, :]

                                                        try:
                                                            # Create a boolean mask indicating burst epochs for the current channel.
                                                            mask = np.array([epoch in burst_band[ch] for epoch in range(ch_data.shape[0])])

                                                            # For the non-burst condition, invert the mask so that "True" indicates epochs without bursts.
                                                            if burst != 'burst':
                                                                mask = ~mask
                                                                # If no epoch remains (i.e. all epochs are bursts), warn and skip this channel.
                                                                if np.sum(mask) == 0:
                                                                    print(f"Warning: Participant {sujname} did not have any epochs without bursts, skipping electrode {ch}")
                                                                    continue
                                                        except KeyError:
                                                            # If the channel is not even in burstdata, handle accordingly.
                                                            print(f"Warning: Participant {sujname} did not have any bursts in electrode {ch}: all epochs selected for no burst condition")
                                                            if burst == 'burst':
                                                                # For burst condition, skip the channel if no burst info exists.
                                                                continue
                                                            else:
                                                                mask = np.ones(ch_data.shape[0], dtype=bool)  # otherwise, select all epochs

                                                        # Create a boolean mask for epochs that have valid (non-NaN) entries in the first column.
                                                        valid_mask = ~np.isnan(ch_data[:, 0])

                                                        # Combine the burst mask with the valid data mask.
                                                        final_mask = mask & valid_mask

                                                        # Warn and skip the channel if no valid epoch remains.
                                                        if np.sum(final_mask) == 0:
                                                            print(f"Warning: Participant {sujname} did not have any valid epochs for electrode {ch}, skipping")
                                                            continue

                                                        # Finally, apply the final mask to select only valid (burst-or-noburst,
                                                        # non-NaN) epochs before the FOOOF fit for this channel.
                                                        ch_data = ch_data[final_mask, :]

                                                        # Core specparam/FOOOF fit for this single channel's burst-conditioned
                                                        # epoch subset, with optional pre-fit epoch-group screening
                                                        # (fit_pregroups) -- see jrpc.interim_fooof_epochselection_chbased
                                                        # (the channel-based variant used only in this stage-2 script,
                                                        # since epochs here are pre-filtered per-channel by burst status).
                                                        aperosc_results, fooofpsd = jrpc.interim_fooof_epochselection_chbased(ch_data, freqs, ch, fit_pregroups, fmin, fidx,
                                                                                           peak_width_limits, max_n_peaks, min_peak_height, peak_threshold,
                                                                                           rsquare_threshold_interim, bands_ranges, aperiodic_mode=aperiodic_mode, nepochs_group=nepochs_group)

                                                        # Tag the fit-result rows with full provenance (method + band
                                                        # definitions + which sweep setting + burst condition produced them).
                                                        aperosc_results['methodfooof'] = str(parameters)
                                                        aperosc_results['methodpsd']   = str(psd_parameters)
                                                        aperosc_results['banddef']     = str(band_definition)
                                                        aperosc_results['maxfreq']     = fidx
                                                        aperosc_results['sujid']       = sujname
                                                        aperosc_results['block']       = blocki
                                                        aperosc_results['burst']       = burst
                                                        aperosc_results['burstband']   = band_names_burst[band]

                                                        # `nepochs_group` sentinel decoding (see nepochs_group_multi definition
                                                        # above): 9999 = fit performed with no epoch-subgroup pre-screening;
                                                        # 1 = epoch-by-epoch pre-fit; any other value = pre-fit in groups of
                                                        # that size. Recorded verbatim so the R pipeline can select/compare.
                                                        if nepochs_group == 9999:
                                                            aperosc_results['prefit']        = 'no_prefit'
                                                            aperosc_results['nepochs_group'] = 'NA'
                                                        elif nepochs_group == 1:
                                                            aperosc_results['prefit']        = 'prefit_by_epoch'
                                                            aperosc_results['nepochs_group'] = 1
                                                        else:
                                                            aperosc_results['prefit']        = 'prefit_by_group'
                                                            aperosc_results['nepochs_group'] = nepochs_group

                                                        # Mask frequencies based on range
                                                        freq_mask = (freqs >= fmin) & (freqs <= fidx)
                                                        foi = freqs[freq_mask]  # Frequencies of interest (masked)

                                                        # Extract the fitted (model + aperiodic-only) spectrum at foi for
                                                        # this single channel, alongside the same fit-result metadata, for
                                                        # spectral-plot reconstruction downstream (channel-based variant,
                                                        # since the fit itself was per-channel above).
                                                        aperosc_psd   = jrpc.fooofpsd_extraction_chbased(fooofpsd, aperosc_results, foi, ch, parameters)
                                                        aperosc_psd   = pd.DataFrame(aperosc_psd)

                                                        aperosc_psd['sujid']           = sujname
                                                        aperosc_psd['methodpsd']       = str(psd_parameters)
                                                        aperosc_psd['methodfooof']     = str(parameters)
                                                        aperosc_psd['block']           = blocki
                                                        aperosc_psd['parameters']      = parameters
                                                        aperosc_psd['maxfreq']         = fidx
                                                        aperosc_psd['burst']           = burst
                                                        aperosc_psd['burstband']       = band_names_burst[band]

                                                        # Same `nepochs_group` sentinel decoding as `aperosc_results` above.
                                                        if nepochs_group == 9999:
                                                            aperosc_psd['prefit']        = 'no_prefit'
                                                            aperosc_psd['nepochs_group'] = 'NA'
                                                        elif nepochs_group == 1:
                                                            aperosc_psd['prefit']        = 'prefit_by_epoch'
                                                            aperosc_psd['nepochs_group'] = 1
                                                        else:
                                                            aperosc_psd['prefit']        = 'prefit_by_group'
                                                            aperosc_psd['nepochs_group'] = nepochs_group

                                                        # ===========================================================
                                                        aperosc_results_df.append(aperosc_results)
                                                        aperosc_psd_df.append(aperosc_psd)

                            # Concatenate every (nepochs_group x fmax x block x burst-condition x
                            # channel) fit result into the two long-format output tables.
                            df_aperosc = pd.concat(aperosc_results_df)
                            df_aperosc_psd = pd.concat(aperosc_psd_df)

                            jrpc.save_dataframe(df_aperosc, subfolder_path, base_filename="aperosc_parameters_burst", rewrite=rewritei, ext=".csv")
                            jrpc.save_dataframe(df_aperosc_psd, subfolder_path, base_filename="aperosc_psds_burst", rewrite=rewritei, ext=".csv")

                    if laggedcoh_hil == 1:
                        # Re-split stage-1's raw per-epoch LAcH matrix (`.npy`, loaded below)
                        # into burst vs. no-burst subsets using stage-1's burst classification
                        # CSV (`burstdata_total`), via jrpc.process_coherence_bursts(). No LAcH
                        # is (re-)computed here
                        fbandsh  = ['alpha']
                        if blocks == 0:
                            if set_file == set_files[0]:
                                # One-time on the first participant: finalize the matrix filename
                                # stem (`matfile_lcoh`) with its .npy extension for the no-blocks case.
                                matfile_lcoh = matfile_lcoh + '.npy'

                            lcohdata = np.load(os.path.join(subfolder_path, matfile_lcoh))
                            lcoh_df_total = jrpc.process_coherence_bursts(burstdata_total=burstdata_total,
                                                     lcohdata=lcohdata,
                                                     ch_names=ch_names,
                                                     fft_range = fft_rangeh,
                                                     lag_range= lag_rangeh,
                                                     bands= fbandsh,
                                                     burst_names=['burst', 'noburst'])

                            lcoh_df_total['sujid'] = sujname
                            lcoh_df_total['block'] = None
                            lcoh_df_total['method']=lcohtype
                            lcoh_df_total['metric']= plv_or_coh
                        else:
                            # Blocked case: stage-1 wrote one .npy matrix per block
                            # (`matfile_lcoh + "_" + block + ".npy"`); load and re-split
                            # each block's matrix separately, skipping any block the
                            # participant lacks epochs for.
                            lcoh_df_total = []
                            for blocki, block in enumerate(block_names):
                                try:
                                    lcoh_hilb_filename = matfile_lcoh + "_" + block + ".npy"
                                    # This is where the crash happens if the file isn't there
                                    lcohdata = np.load(os.path.join(subfolder_path, lcoh_hilb_filename))

                                except FileNotFoundError:
                                    print(f"Warning: The file for block {blocki} ('{lcoh_hilb_filename}') was not found. Participant {sujname} did not have epochs in this condition.")
                                    continue

                                burstdata_band = burstdata_total[burstdata_total['block']==block]
                                blockpos  = np.where(event_array == block)[0]
                                
                                if blockpos.size > 0:
                                    lcoh_df = jrpc.process_coherence_bursts(burstdata_total=burstdata_band,
                                                             lcohdata=lcohdata,
                                                             ch_names=ch_names,
                                                             fft_range = fft_rangeh,
                                                             lag_range= lag_rangeh,
                                                             bands= fbandsh,
                                                             burst_names=['burst', 'noburst'])
                                    
                                    lcoh_df['sujid'] = sujname
                                    lcoh_df['block'] = block
                                    lcoh_df['method']=lcohtype
                                    lcoh_df['metric']= plv_or_coh
                                    lcoh_df_total.append(lcoh_df)
                            lcoh_df_total = pd.concat(lcoh_df_total, ignore_index=True)
                            
                        jrpc.save_dataframe(lcoh_df_total, subfolder_path, base_filename=basefile_lcoh, rewrite=rewritei, ext=".csv")



                # Any ValueError raised anywhere in the try block above (e.g., malformed/
                # unreadable .set file, missing stage-1 burst CSV/LAcH matrix, incompatible
                # event structure) is caught here so one participant's failure does not
                # halt the whole batch; the error is logged to `error_list` for post-hoc
                # review instead of being silently swallowed.
                except ValueError as e:
                    print(f"Skipping {sujname} : {e}")
                    error_message = f"Participant {sujname} crashed because: {e}"
                    error_list.append(error_message)
                
           
    
    
                
            
    
       