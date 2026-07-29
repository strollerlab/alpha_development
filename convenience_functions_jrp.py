#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 14:38:19 2025

@author: J. Rico-Picó
"""
import numpy as np 
import pandas as pd 
import warnings
import sys
import re
import specparam
from   scipy.signal import welch 
import os
from   scipy.signal.windows import dpss
from   mne.time_frequency   import psd_array_multitaper
from   scipy.signal import fftconvolve
import mne
from scipy.signal.windows import hann

  
######## LAGGED COHERENCE FUNCTIONS FOR FFT - MTMCONVOL PLUS CLASSIC LAGGED COHERENCE (FRANSEN ET AL. 2015) AND CONVENIENCE FUNCTIONS           
def ft_specest_mtmconvol(dat, time, toi, t_ftimwin, foi, pad, polyorder=0, verbose=False):
    """
    Compute a single-frequency time–frequency transform analogous to FieldTrip's ft_specest_mtmconvol.
    
    Parameters
    ----------
    dat : ndarray, shape (n_epochs, n_channels, n_times)
        Input data (e.g., from an MNE object) for each epoch.
    time : 1D ndarray
        Time vector (in seconds) corresponding to the n_times samples.
    toi : 1D array_like
        Time points of interest (in seconds) at which to extract the spectrum.
    t_ftimwin : float or scalar
        Duration (in seconds) of the taper (time–window). A single value is expected.
    foi : float or single-element array_like
        Frequency of interest (Hz). Only a single frequency is allowed.
    pad : float
        Padding duration (in seconds). Must be at least as long as the epoch duration.
    polyorder : int, optional
        The polynomial order for detrending (if ≥0, demeaning is performed).
    verbose : bool, optional
        If True, prints extra information.
        
    Returns
    -------
    spectrum : ndarray, shape (n_epochs, n_channels, 1, n_toi)
        Complex time–frequency transform at the given time–points of interest.
    ntaper : int
        Number of tapers used (here always 1).
    foi_out : ndarray, shape (1,)
        The frequency of interest (returned as a one–element array).
    toi_used : ndarray
        The final set of time–points-of–interest (subset of input 'toi' that were valid).
    """
    # --- 1. Preprocessing ---
    if dat.dtype not in [np.float32, np.float64]:
        dat = dat.astype(np.float64)
    
    # Remove polynomial trend (here simply demean along time)
    if polyorder >= 0:
        dat = dat - np.mean(dat, axis=2, keepdims=True)
        
    n_epochs, n_channels, n_times = dat.shape
    fsample = 1.0 / np.mean(np.diff(time))
    epoch_duration = n_times / fsample
    if verbose:
        print(f"n_epochs: {n_epochs}, n_channels: {n_channels}, n_times: {n_times}")
        print(f"Sampling frequency: {fsample:.3f} Hz, epoch duration: {epoch_duration:.3f} s")
        
    # --- 2. Padding check ---
    if round(pad * fsample) < n_times:
        raise ValueError("The padding specified is shorter than the data.")
    # (Note: here we use fftconvolve with mode='same', so no explicit zero-padding is required.)
    
    # --- 3. Frequency of interest ---
    if np.isscalar(foi):
        foi_val = foi
        foi_out = np.array([foi_val])
    else:
        foi = np.asarray(foi)
        if foi.size != 1:
            raise ValueError("Expected a single frequency of interest for 4D output [epochs, electrodes, 1, toi].")
        foi_val = foi.item()
        foi_out = np.array([foi_val])
        
    # --- 4. Time points of interest ---
    toi = np.asarray(toi)
    # Select those time points within the original time range
    valid_mask = (toi >= time[0]) & (toi <= time[-1])
    toi_used = toi[valid_mask]
    if toi_used.size == 0:
        raise ValueError("No time points of interest fall within the data time range.")
    # For each requested toi, find the closest index in the original time vector
    toi_indices = np.array([np.argmin(np.abs(time - t)) for t in toi_used])
    n_toi = len(toi_indices)
    if verbose:
        print(f"Using {n_toi} time–of–interest samples.")
    
    # --- 5. Taper and Wavelet construction ---
    if np.isscalar(t_ftimwin):
        t_win = t_ftimwin
    else:
        t_ftimwin = np.asarray(t_ftimwin)
        if t_ftimwin.size != 1:
            raise ValueError("Expected a single t_ftimwin for this analysis.")
        t_win = t_ftimwin.item()
        
    win_samples = int(round(t_win * fsample))
    if win_samples < 1:
        raise ValueError("Time window is too short (in samples).")
    
    # Construct a symmetric Hanning taper and normalize it
    taper  = np.hanning(win_samples)
    taper  = taper / np.linalg.norm(taper)
    ntaper = 1  # only one taper is used here

    # Create a time vector (in seconds) for the wavelet. Note that we center around zero.
    half_range = (win_samples - 1) / 2.0
    t_vec = np.linspace(-half_range, half_range, win_samples) / fsample
    angle = 2 * np.pi * foi_val * t_vec
    # Construct the complex wavelet: multiply taper with complex exponential
    wavelet = taper * (np.cos(angle) + 1j * np.sin(angle))
    
    # Scaling (as in MATLAB: sqrt(2 / t_ftimwinsample))
    scaling = np.sqrt(2.0 / win_samples)
    
    if verbose:
        print(f"Wavelet built with t_ftimwin: {t_win:.3f}s ({win_samples} samples), Foi: {foi_val} Hz")
    
    # --- 6. Convolution ---
    # Preallocate the spectrum array: [n_epochs, n_channels, 1, n_toi]
    spectrum = np.empty((n_epochs, n_channels, 1, n_toi), dtype=np.complex64)
    
    for epoch in range(n_epochs):
        for ch in range(n_channels):
            # Convolve data with the wavelet using FFT-based convolution (mode='same' returns original length)
            conv_result = fftconvolve(dat[epoch, ch, :], wavelet, mode='same')
            conv_result *= scaling
            # Extract only the time–points of interest
            spectrum[epoch, ch, 0, :] = conv_result[toi_indices]
            
    return spectrum, ntaper, foi_out, toi_used

def compute_epoch_based_lagged_coherence(ft_coefs, delta=1):
    """
    Compute epoch-based lagged coherence from 4D time-frequency coefficients.
    
    The input array is assumed to be four-dimensional with shape:
         (n_epochs, n_channels, 1, n_timepoints)
    This function first reduces the data by removing the singleton frequency dimension.
    Then, it computes lagged coherence using a one-window shift (delta steps).
    
    The calculation per channel proceeds as:
       1. Let X1 be the coefficients at time indices [0, ..., n_timepoints-delta-1]
          and X2 be the coefficients at time indices [delta, ..., n_timepoints-1].
       2. For each epoch and channel, compute:
                num_epoch = sum( X1 * conj(X2) )       over time
                den1_epoch = sum( |X1|^2 )               over time
                den2_epoch = sum( |X2|^2 )               over time
       3. Sum these quantities across epochs to yield a total numerator and power sums.
       4. For each channel, if the product of powers is positive, then:
                coh = |total_num| / sqrt(total_den1 * total_den2)
           Otherwise, coherence is NaN.
    
    Parameters
    ----------
    ft_coefs : np.ndarray
         Four-dimensional complex array with shape (n_epochs, n_channels, 1, n_timepoints)
         representing time-frequency coefficients.
    delta : int, optional
         The number of time-window shifts to apply (default is 1).
    
    Returns
    -------
    coh : np.ndarray, shape (n_channels,)
         Lagged coherence computed per channel after summing over epochs.
    """
    # If the input is 4D but the third dimension is singleton, squeeze it.
    if ft_coefs.ndim == 4:
        if ft_coefs.shape[2] == 1:
            ft_coefs = np.squeeze(ft_coefs, axis=2)
        else:
            # Alternatively, if more than one frequency is present,
            # you might average across the frequency dimension:
            ft_coefs = np.mean(ft_coefs, axis=2)
    
    n_epochs, n_channels, n_timepoints = ft_coefs.shape
    if n_timepoints <= delta:
        return np.full(n_channels, np.nan)
    
    # Split the coefficients into two time-shifted sets.
    X1 = ft_coefs[:, :, :-delta]  # shape: (n_epochs, n_channels, n_timepoints-delta)
    X2 = ft_coefs[:, :, delta:]   # shape: (n_epochs, n_channels, n_timepoints-delta)
    
    # Compute numerator and denominators for each epoch and channel.
    # Here, summing over the time window axis.
    num   = np.sum(X1 * np.conjugate(X2), axis=2)    # shape: (n_epochs, n_channels)
    den1  = np.sum(np.abs(X1) ** 2, axis=2)          # shape: (n_epochs, n_channels)
    den2  = np.sum(np.abs(X2) ** 2, axis=2)          # shape: (n_epochs, n_channels)
    
    # Sum over epochs:
    total_num  = np.sum(num, axis=0)   # shape: (n_channels,)
    total_den1 = np.sum(den1, axis=0)
    total_den2 = np.sum(den2, axis=0)
    
    # Compute lagged coherence per channel.
    coh = np.full(n_channels, np.nan, dtype=np.float64)
    valid = (total_den1 * total_den2) > 0
    coh[valid] = np.abs(total_num[valid]) / np.sqrt(total_den1[valid] * total_den2[valid])
    
    return coh


def process_coherence_bursts(burstdata_total: pd.DataFrame,
                             lcohdata: np.ndarray,
                             ch_names: list[str],
                             fft_range: np.ndarray,
                             lag_range: np.ndarray,
                             bands: list[str],
                             burst_names: list[str] = None) -> pd.DataFrame:
    """
    Processes lagged coherence data by averaging based on burst occurrences
    in different frequency bands and reshapes the results into a long DataFrame.
    """
    if burst_names is None:
        burst_names = ['burst', 'noburst']

    # --- Derive Metadata from Data ---
    n_epochs, n_channels, n_freqs, n_lags = lcohdata.shape

    # --- Validation ---
    if n_channels != len(ch_names):
        raise ValueError(f"Mismatch between lcohdata channels ({n_channels}) and ch_names length ({len(ch_names)})")
    if n_freqs != len(fft_range):
        raise ValueError(f"Mismatch between lcohdata frequencies ({n_freqs}) and fft_range length ({len(fft_range)})")
    if n_lags != len(lag_range):
        raise ValueError(f"Mismatch between lcohdata lags ({n_lags}) and lag_range length ({len(lag_range)})")
    if not all(band in burstdata_total['band'].unique() for band in bands):
        print("Warning: Some requested bands are not found in burstdata.")

    # --- Data Processing Initialization ---
    lcoh_burst_avg = np.full((n_channels, n_freqs, n_lags, len(burst_names), len(bands)), np.nan)
    lcoh_burst_counts = np.zeros((n_channels, len(burst_names), len(bands)), dtype=int)

    print(f"Processing {n_channels} channels, {len(bands)} bands...")

    # --- Main Loop ---
    for bi, band in enumerate(bands):
        # CORRECTED: Use `band` directly, not the undefined `band_names_burst`
        burst_band_df = burstdata_total[burstdata_total['band'] == band]
        burst_band_df = burst_band_df[burst_band_df['is_burst'] == True]

        burst_epochs_by_ch = (burst_band_df.groupby('ch')['epoch_block']
                              .apply(lambda x: list(set(x)))
                              .to_dict())

        for bursti, burst_condition in enumerate(burst_names):
            for ch_idx, ch in enumerate(ch_names):
                mask = np.zeros(n_epochs, dtype=bool) # Start with an empty mask

                try:
                    # Get burst epochs for the current channel
                    burst_epochs = burst_epochs_by_ch[ch]
                    # CORRECTED: Use `n_epochs`, not `ch_data.shape[0]`
                    mask = np.array([epoch in burst_epochs for epoch in range(n_epochs)])

                    if burst_condition != 'burst':
                        mask = ~mask

                    if np.sum(mask) == 0:
                        print(f"Warning: No epochs for condition '{burst_condition}' on electrode {ch} in band {band}.")
                        continue
                except KeyError:
                    # CORRECTED: Removed undefined `sujname`
                    print(f"Warning: No burst data for electrode {ch} in band {band}.")
                    if burst_condition == 'burst':
                        continue 
                    else:

                        mask = np.ones(n_epochs, dtype=bool)

                # --- Create valid mask and combine ---
                valid_mask = ~np.isnan(lcohdata[:, ch_idx, 0, 0])
                
                # CORRECTED LOGIC: Use a combined `final_mask` for all subsequent operations
                final_mask = mask & valid_mask

                current_n_epochs                      = np.sum(final_mask)
                lcoh_burst_counts[ch_idx, bursti, bi] = current_n_epochs

                # Check if any valid epochs remain
                if current_n_epochs > 0:
                    # CORRECTED LOGIC: Slice data using the `final_mask`
                    data_slice = lcohdata[final_mask, ch_idx, :, :]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        lcoh_burst_avg[ch_idx, :, :, bursti, bi] = np.nanmean(data_slice, axis=0)
                # If no epochs, the values remain NaN and count remains 0, so no `else` is needed

    print("\nProcessing finished.")
    print("Creating long DataFrame...")

    # --- Reshape to Long DataFrame ---
    I, J, K, L, M = np.indices(lcoh_burst_avg.shape)

    df_long = pd.DataFrame({
        'ch':         np.array(ch_names)[I.flatten()],
        'freq':       np.array(fft_range)[J.flatten()],
        'lag':        np.array(lag_range)[K.flatten()], # Renamed 'cycle' to 'lag' for consistency
        'burst':      np.array(burst_names)[L.flatten()], # Renamed for clarity
        'burst_band': np.array(bands)[M.flatten()], # Renamed for clarity
        'lcoh':       lcoh_burst_avg.flatten(),
        'n_epochs':   lcoh_burst_counts[I.flatten(), L.flatten(), M.flatten()]
    })

    print("DataFrame created successfully.")
    return df_long

def compute_epoch_based_lagged_coherence_chbased(ft_coefs, delta=1):
    """
    Compute epoch-based lagged coherence from time-frequency coefficients for a single channel.
    
    The function assumes the input is either:
      - A 4D complex array with shape (n_epochs, 1, 1, n_timepoints), or
      - A 2D complex array with shape (n_epochs, n_timepoints).
    
    Steps:
      1. Remove singleton dimensions if present.
      2. Define two time-shifted arrays:
           X1 corresponds to time indices [0, ..., n_timepoints-delta-1]
           X2 corresponds to time indices [delta, ..., n_timepoints-1]
      3. For each epoch, compute:
              num_epoch = sum( X1 * conj(X2) ) over time
              den1_epoch = sum( |X1|^2 ) over time
              den2_epoch = sum( |X2|^2 ) over time
      4. Sum these quantities across epochs.
      5. If (total_den1 * total_den2) > 0, compute coherence as:
              coh = |total_num| / sqrt(total_den1 * total_den2)
         Otherwise, return NaN.
    
    Parameters
    ----------
    ft_coefs : np.ndarray 
        Complex array with shape either (n_epochs, 1, 1, n_timepoints) or (n_epochs, n_timepoints).
    delta : int, optional
        The number of time-window shifts to apply (default is 1).
    
    Returns
    -------
    coh : float
        Lagged coherence computed across epochs.
    """
    # If a 4D array is provided, remove singleton dimensions.
    if ft_coefs.ndim == 4:
        ft_coefs = np.squeeze(ft_coefs, axis=(1, 2))
    
    # After squeezing, we expect the shape to be (n_epochs, n_timepoints).
    if ft_coefs.ndim != 2:
        raise ValueError("Input array must be 2D after squeezing, got shape: {}".format(ft_coefs.shape))
    
    n_epochs, n_timepoints = ft_coefs.shape
    
    # If the time dimension is too short for the specified shift, return NaN.
    if n_timepoints <= delta:
        return np.nan
    
    # Create time-shifted versions of the coefficients.
    X1 = ft_coefs[:, : -delta]  # shape: (n_epochs, n_timepoints - delta)
    X2 = ft_coefs[:, delta:]    # shape: (n_epochs, n_timepoints - delta)
    
    # Compute quantities for each epoch by summing over the time axis.
    num   = np.sum(X1 * np.conjugate(X2), axis=1)  # Complex-valued sum per epoch.
    den1  = np.sum(np.abs(X1)**2, axis=1)           # Power in X1 per epoch.
    den2  = np.sum(np.abs(X2)**2, axis=1)           # Power in X2 per epoch.
    
    # Sum these values across all epochs.
    total_num  = np.sum(num)
    total_den1 = np.sum(den1)
    total_den2 = np.sum(den2)
    
    # Compute lagged coherence if the denominators are valid.
    if total_den1 * total_den2 > 0:
        coh = np.abs(total_num) / np.sqrt(total_den1 * total_den2)
    else:
        coh = np.nan
    
    return coh

def create_noise(data):
    noise_arr = np.random.normal(size=data._data.shape)
    nchans = noise_arr.shape[0]

    ch_types = ['seeg']*nchans
    ch_names = ['noise #{}'.format(i) for i in range(nchans)]

    info = mne.create_info(ch_names=ch_names, sfreq=data.info['sfreq'], ch_types=ch_types)
    noise = mne.io.RawArray(noise_arr, info, verbose=False)

    noise.notch_filter(np.arange(50, data.info['sfreq']//2, 50), trans_bandwidth=0.1, verbose=False)

    return noise

def update_progress(current, total, bar_length=10, message = "Progress"):
    """Prints a progress bar to the console.

    Args:
        current (int): The current iteration count.
        total (int): The total number of iterations.
        bar_length (int): The length of the progress bar in characters.
    """
    percent = current / total * 100
    filled_length = int(bar_length * current // total)
    bar = '#' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write(f'\r{message}: |{bar}| {percent:5.1f}% Complete')
    sys.stdout.flush()

#### POWER SPECTRUM FUNCTIONS ################################################
def compute_fft_hanning(epochs, sfreq):
    """
    Compute PSD using an FFT after applying a Hanning window on each epoch.
    Returns PSD estimates with units V²/Hz for a one-sided spectrum.

    Parameters:
        epochs : MNE Epochs object.
        sfreq  : Sampling frequency.
    Returns:
        freqs : Frequency bins.
        psds  : PSD for each epoch as a 3D array (n_epochs, n_channels, n_freq_bins).
    """
    data = epochs.get_data()  # (n_epochs, n_channels, n_times)
    n_epochs, n_channels, n_times = data.shape

    # Create Hanning window and calculate its power
    window = np.hanning(n_times)
    U = np.sum(window**2)  # Window power

    # Frequency bins for one-sided spectrum
    freqs = np.fft.rfftfreq(n_times, 1 / sfreq)

    # Apply window and compute PSD using FFT
    psds = np.empty((n_epochs, n_channels, len(freqs)))  # Pre-allocate array for efficiency
    for i, epoch in enumerate(data):
        epoch_windowed = epoch * window[np.newaxis, :]  # Apply window to each channel
        fft_vals = np.fft.rfft(epoch_windowed, axis=-1)
        

        psd = np.abs(fft_vals)**2 / (sfreq * U)

        if n_times % 2 == 0:
            # If even, Nyquist is at the end (do not double Nyquist)
            psd[:, 1:-1] *= 2
        else:
            # If odd, Nyquist is not hit exactly (double everything after DC)
            psd[:, 1:] *= 2

        psds[i] = psd  # Store PSD for each epoch
    
    return freqs, psds

def compute_fft_hanning_normlength(epochs, sfreq):
    """
    Computes the power spectrum for multi-epoch, multi-channel data based on
    the specific scaling from the original MATLAB script.

    Parameters:
    -----------
    epochs : MNE Epochs object or numpy.ndarray
        If an MNE Epochs object, `get_data()` will be called.
        If a NumPy array, it should have the shape (n_epochs, n_channels, n_times).
    sfreq : int or float
        The sampling frequency of the data in Hz.

    Returns:
    --------
    freqs : numpy.ndarray
        A 1D array of frequency bins (excluding DC and Nyquist).
    amp_squared : numpy.ndarray
        A 3D array of SQUARED, LENGTH-NORMALISED AMPLITUDE for each epoch,
        channel and frequency bin, shape (n_epochs, n_channels, n_freqs).

        *** Q12 FIX (2026-07-15) -- THIS IS NOT A PSD. ***
        There is no division by sfreq and no window-power (U) normalisation, so
        the units are amplitude^2 (arbitrary), NOT V^2/Hz. The function
        reproduces the original MATLAB script's scaling and is kept for parity
        with it, so the NUMBERS ARE DELIBERATELY UNCHANGED -- only the name and
        this documentation were corrected. Do not mix these values with
        compute_welch() output in the same analysis. For a PSD in V^2/Hz, use
        compute_welch().
    """
    # Check if the input is an MNE-Python Epochs object and get data
    if hasattr(epochs, 'get_data'):
        data = epochs.get_data()
    else:
        data = epochs # Assume it's already a NumPy array

    if data.ndim != 3:
        raise ValueError("Input data must be a 3D array (epochs, channels, times).")

    n_epochs, n_channels, n_times = data.shape

    # Create a single hANNING window and apply it to all epochs and channels via broadcasting.
    window = hann(n_times)
    data_windowed = data * window

    # --- 2. Fourier Transform ---
    # Use rfft for real-valued input for efficiency.
    fft_vals = np.fft.rfft(data_windowed, axis=-1)

    # --- 3. Get Frequencies ---
    # Calculate the frequency bins once.
    all_freqs = np.fft.rfftfreq(n_times, 1 / sfreq)

    # --- 4. Get Amplitude and Normalize ---
    # Normalize amplitude by the length of the data (n_times).
    amp_vals = np.abs(fft_vals) / n_times

    # --- 5. Select Frequencies and Scale Amplitude ---
    # Select only the frequencies between DC (0 Hz) and Nyquist.
    freqs = all_freqs[1:-1]
    amp_vals_scaled = 2 * amp_vals[..., 1:-1]

    amp_squared = amp_vals_scaled**2

    return freqs, amp_squared

def compute_welch(epochs, sfreq, nperepoch=None, noverlap=None):
    """
    Compute PSD using SciPy's Welch method.
    
    Parameters:
        epochs  : MNE Epochs object.
        sfreq   : Sampling frequency.
        nperepoch : Epoch length (if None, full epoch is used).
        noverlap: Overlap between Epochs (if None, defaults to nperepoch//2).
    
    Returns:
        freqs : Frequency bins.
        psds  : PSD for each epoch as a 3D array.
    """
    data = epochs.get_data()  # (n_epochs, n_channels, n_times)
    n_epochs, n_channels, n_times = data.shape

    if nperepoch is None:
        nperepoch = n_times
    if noverlap is None:
        noverlap = nperepoch // 2

    all_psd = []
    freqs = None
    for epoch in data:
        psd_channels = []
        for channel in epoch:
            f, p = welch(channel, fs=sfreq, window='hann',
                         nperepoch=nperepoch, noverlap=noverlap)
            psd_channels.append(p)
            if freqs is None:
                freqs = f
        all_psd.append(np.array(psd_channels))
    all_psd = np.array(all_psd)
    return freqs, all_psd

def compute_spectopo(epochs, sfreq, nperepoch=None, noverlap=None):
    """
    Compute a Spectopo-like PSD by using Welch and converting its output to dB.
    Note: This method now returns its PSD in dB scale.
    Parameters:
        epochs  : MNE Epochs object.
        sfreq   : Sampling frequency.
        nperepoch : Epoch length.
        noverlap: Overlap between Epochs.
    
    Returns:
        freqs  : Frequency bins.
        psd_db : PSD for each epoch (in dB) as a 3D array.
    """
    freqs, psd = compute_welch(epochs, sfreq, nperepoch, noverlap)
    psd_db = 10 * np.log10(psd)
    return freqs, psd_db

def compute_multitaper(epochs, sfreq, time_bandwidth=4.0, fmin=0, fmax=np.inf):
    """
    Compute PSD using MNE's built-in multitaper method.
    Parameters:
        epochs        : MNE Epochs object.
        sfreq         : Sampling frequency.
        time_bandwidth: Time–bandwidth product (NW).
        fmin          : Minimum frequency for analysis.
        fmax          : Maximum frequency for analysis.
    
    Returns:
        freqs : Frequency bins.
        psds  : PSD for each epoch as a 3D array.
    """
    data = epochs.get_data()
    all_psd = []
    freqs = None
    for epoch in data:
        psd, freqs = psd_array_multitaper(epoch, sfreq=sfreq,
                                          fmin=fmin, fmax=fmax,
                                          bandwidth=time_bandwidth,
                                          adaptive=False, normalization='full',
                                          verbose=False)
        all_psd.append(psd)
    all_psd = np.array(all_psd)
    return freqs, all_psd

def compute_manual_multitaper(epochs, sfreq, time_bandwidth=4.0, num_tapers=None):
    """
    Compute PSD manually using a multitaper approach with DPSS tapers.
    Parameters:
        epochs        : MNE Epochs object.
        sfreq         : Sampling frequency.
        time_bandwidth: Time–bandwidth product (NW).
        num_tapers    : Number of tapers (default: 2*NW - 1).
    Returns:
        freqs : Frequency bins.
        psds  : PSD for each epoch as a 3D array.
    """
    
    data = epochs.get_data()  # (n_epochs, n_channels, n_times)
    n_epochs, n_channels, n_times = data.shape
    if num_tapers is None:
        num_tapers = int(2 * time_bandwidth) - 1
    freqs = np.fft.rfftfreq(n_times, 1/sfreq)

    all_psd = []
    for epoch in data:
        psd_channels = []
        for channel in epoch:
            tapers, _ = dpss(n_times, NW=time_bandwidth, Kmax=num_tapers, return_ratios=True)
            psd_tapers = []
            for taper in tapers:
                tapered_signal = channel * taper
                fft_vals = np.fft.rfft(tapered_signal)
                U_t = np.sum(taper**2)
                psd_val = np.abs(fft_vals)**2 / (sfreq * U_t)
                if n_times % 2 == 0:
                    psd_val[1:-1] *= 2
                else:
                    psd_val[1:] *= 2
                psd_tapers.append(psd_val)
            psd_channel = np.mean(psd_tapers, axis=0)
            psd_channels.append(psd_channel)
        all_psd.append(np.array(psd_channels))
    all_psd = np.array(all_psd)
    return freqs, all_psd


##### FOOOF related functions #################################################

def extract_band_params(peaks, band_range):
    """
    Extracts parameters of the first peak within a specified frequency band range.

    Parameters:
        peaks : list or np.ndarray
            List of peak parameters, where each peak is an iterable [CF, PW, BW].
        band_range : tuple
            Frequency band range as (min_freq, max_freq).

    Returns:
        list
            Parameters of the first peak within the range [CF, PW, BW], or [np.nan, np.nan, np.nan] if no peaks are found.
    """
    # Ensure peaks is iterable; if it's a scalar, return NaN for all parameters
    if isinstance(peaks, list) or isinstance(peaks, np.ndarray):
        band_peaks = [peak for peak in peaks if band_range[0] <= peak[0] <= band_range[1]]
        if band_peaks:
            return band_peaks[0]  # Return the first peak in the range
        else:
            return [np.nan, np.nan, np.nan]  # No peak found in the range
    else:
        return [np.nan, np.nan, np.nan]  # If peaks is not iterable, return NaN
    
def interim_fooof_epochselection(psds, freqs, ch_names, fit_pregroups, fmin, fmax,
                                 peak_width_limits, max_n_peaks, min_peak_height, peak_threshold,
                                 rsquare_threshold_interim, perelect, bands_ranges, aperiodic_mode, nepochs_group=None):
    """
    Compute included epoch absolute power and select epochs (or groups) based on FOOOF fits.

    Parameters
    ----------
    psds : np.ndarray
        PSD absolute power array with shape (epochs, electrodes, freqs).
    freqs : np.ndarray
        Frequency vector corresponding to the PSD.
    fit_pregroups : int
        If 2, process each epoch individually; if 1, process epochs grouped by nepochs_group;
        otherwise return mean_abspower.
    fmin : float
        Lower frequency bound for FOOOF fitting.
    fmax : float
        Upper frequency bound for FOOOF fitting.
    peak_width_limits : iterable
        Parameter for FOOOF specifying minimum and maximum allowable peak widths.
    max_n_peaks : int
        Maximum number of peaks allowed in specparam.
    min_peak_height : float
        Minimum peak height for specparam.
    peak_threshold : float
        Peak threshold for specparam.
    rsquare_threshold_interim : float
        R² threshold for considering a epoch or group as valid.
    perelect : float
        Proportion of electrodes required to pass the r² threshold.
    nepochs_group : int, optional
        Number of epochs to average when processing groups (fit_pregroups==1).
    mean_abspower : np.ndarray, optional
        Mean absolute power to return if fit_pregroups is not 1 or 2.

    Returns
    -------
    included_psds : np.ndarray
        Averaged PSD of selected epochs (shape: electrodes x freqs).
    epochs_included : int or None
        Number of epochs (or groups) that were included.
    epochs_percentage : float or None
        Proportion of total epochs that were included.
    r2_included : float or None
        Average r² of the included (selected) fits.
    r2_discard : float or None
        Average r² of the discarded (non-selected) fits.
    """
    electrodes = psds.shape[1]
    num_epochs = psds.shape[0]
    
    # ------------------------
    # Method 2: Process each epoch individually.
    
    lenpsd         = len(np.arange(fmin,fmax,freqs[1]-freqs[0])) + 1
        
    
    if fit_pregroups   == 0 or num_epochs <= nepochs_group:
        included_psds  = np.average(psds[:, :, :], axis=0)

    elif fit_pregroups == 2:
        r2_value = np.zeros((electrodes, num_epochs))
        
        # Loop over each epoch and electrode to perform FOOOF fits.
        for epoch in range(num_epochs):
            for elec in range(electrodes):
                psd_data = psds[epoch, elec, :]
                fm = specparam.SpectralModel(
                    peak_width_limits=peak_width_limits,
                    max_n_peaks=max_n_peaks,
                    min_peak_height=min_peak_height,
                    peak_threshold=peak_threshold, 
                    aperiodic_mode=aperiodic_mode)
                # Use the provided frequency bounds.
                fm.fit(freqs=freqs, power_spectrum=psd_data, freq_range=[fmin, fmax])
                r2_value[elec, epoch] = fm.r_squared_
                
        # Determine which epochs meet the r² criteria:
        r2_dummy = (r2_value >= rsquare_threshold_interim)
        epoch_counts = np.sum(r2_dummy, axis=0)
        valid_epochs = epoch_counts >= (perelect * electrodes)
        
        r2_included = np.average(r2_value[:, valid_epochs]) if np.any(valid_epochs) else np.nan
        r2_discard  = np.average(r2_value[:, ~valid_epochs]) if np.any(~valid_epochs) else np.nan
        
        indices = np.where(valid_epochs)[0]
        epochs_included   = len(indices)
        epochs_percentage = epochs_included / num_epochs

        if epochs_percentage > 0:
            included_psds = np.average(psds[indices, :, :], axis=0)
        else:
            included_psds = np.average(psds[:, :, :], axis=0)
            
    # ------------------------
    # Method 1: Process epochs grouped by nepochs_group.
    elif fit_pregroups == 1:
        if nepochs_group is None:
            raise ValueError("nepochs_group must be provided when fit_pregroups==1")
        
        # Determine the number of groups.
        if num_epochs <= nepochs_group:
            ngroups =1 
        else:
            ngroups = int(np.round(num_epochs / nepochs_group))
            if (ngroups * nepochs_group) > num_epochs or (num_epochs - ngroups * nepochs_group) < nepochs_group:
                ngroups -= 1
        
        epoch1_loop = np.zeros(ngroups, dtype=int)
        epoch2_loop = np.zeros(ngroups, dtype=int)
        r2_value = np.zeros((electrodes, ngroups))
        
        # Loop over groups.
        for grp in range(ngroups):
            epoch1 = nepochs_group * grp
            epoch2 = nepochs_group * grp + nepochs_group
            # Adjust the second index if necessary to avoid overflow.
            if (epoch2 + nepochs_group) >= num_epochs:
                epoch2 = num_epochs
            
            epoch1_loop[grp] = epoch1
            epoch2_loop[grp] = epoch2
            
            # For each electrode, average within the group and perform FOOOF fitting.
            for elec in range(electrodes):
                psd_data_group = np.average(psds[epoch1:epoch2, elec, :], axis=0)
                fm = specparam.SpectralModel(
                    peak_width_limits=peak_width_limits,
                    max_n_peaks=max_n_peaks,
                    min_peak_height=min_peak_height,
                    peak_threshold=peak_threshold, aperiodic_mode=aperiodic_mode)
                fm.fit(freqs=freqs, power_spectrum=psd_data_group, freq_range=[fmin, fmax])
                r2_value[elec, grp] = fm.r_squared_
                
        # Check group passes criterion:
        r2_dummy = (r2_value >= rsquare_threshold_interim)
        group_counts = np.sum(r2_dummy, axis=0)
        valid_groups = group_counts >= (perelect * electrodes)
        
   
        r2_included = np.average(r2_value[:, valid_groups]) if np.any(valid_groups) else np.nan
        r2_discard  = np.average(r2_value[:, ~valid_groups]) if np.any(~valid_groups) else np.nan
        
        # Mark individual epochs for exclusion based on group-level analysis.
        selected_epochs = np.ones(num_epochs, dtype=bool)
        for grp in range(ngroups):
            if not valid_groups[grp]:
                pos1 = epoch1_loop[grp]
                pos2 = epoch2_loop[grp]
                selected_epochs[pos1:pos2] = False
        
        indices = np.where(selected_epochs)[0]
        
        epochs_included   = len(indices)
        epochs_percentage = epochs_included / num_epochs
        
        if epochs_percentage > 0:
            included_psds = np.average(psds[indices, :, :], axis=0)
        else:
            included_psds = np.average(psds[:, :, :], axis=0)

    aperi_psd = np.zeros((electrodes, lenpsd))
    absol_psd = np.zeros((electrodes, lenpsd))
    fooof_psd = np.zeros((electrodes, lenpsd))
    oscil_psd = np.zeros((electrodes, lenpsd))
    error_psd = np.zeros((electrodes, lenpsd))
        
    # Initialize arrays for storing FOOOF results
    r2value = np.zeros(electrodes)
    mae = np.zeros(electrodes)
    offset = np.zeros(electrodes)
    slope = np.zeros(electrodes)
        
    # Initialize arrays for storing peak parameters (default to NaN)
    alpha_freq = np.full(electrodes, np.nan)
    alpha_ampl = np.full(electrodes, np.nan)
    alpha_widt = np.full(electrodes, np.nan)
        
    theta_freq = np.full(electrodes, np.nan)
    theta_ampl = np.full(electrodes, np.nan)
    theta_widt = np.full(electrodes, np.nan)
        
    beta_freq = np.full(electrodes, np.nan)
    beta_ampl = np.full(electrodes, np.nan)
    beta_widt = np.full(electrodes, np.nan)
        
    gamma_freq = np.full(electrodes, np.nan)
    gamma_ampl = np.full(electrodes, np.nan)
    gamma_widt = np.full(electrodes, np.nan)
        
     # Loop over electrodes to compute PSDs and extract FOOOF parameters
    for elec_idx in range(electrodes):
            # Extract PSDs for the current electrode
        psds_elec = included_psds[elec_idx, :]  # Use correct indexing
            
            # Initialize and configure FOOOF
        fm = specparam.SpectralModel(
                peak_width_limits=peak_width_limits,
                max_n_peaks=max_n_peaks,
                min_peak_height=min_peak_height,
                peak_threshold=peak_threshold, 
                aperiodic_mode=aperiodic_mode)
            
            # Fit the FOOOF model on the electrode's PSD data
        fm.fit(freqs=freqs, power_spectrum=psds_elec, freq_range=[fmin, fmax])
        
            # Store overall results
        r2value[elec_idx] = fm.r_squared_
        # `fm.error_` = specparam's fit error. specparam.compute_error() defaults to
        # error_metric='mae' (mean ABSOLUTE error) and this pipeline never overrides it
        mae[elec_idx] = fm.error_
        offset[elec_idx] = fm.aperiodic_params_[0]
        slope[elec_idx] = fm.aperiodic_params_[1]
            
            # Extract components from the FOOOF model
        aperi_psd[elec_idx, :] = fm._ap_fit
        absol_psd[elec_idx, :] = fm.power_spectrum
        fooof_psd[elec_idx, :] = fm.modeled_spectrum_
        oscil_psd[elec_idx, :] = fm.power_spectrum - fm._ap_fit
        error_psd[elec_idx, :] = fm.power_spectrum - fm.modeled_spectrum_
            
    
        peaks = fm.get_params('peak_params')

        peaks_array = np.atleast_2d(np.asarray(peaks, dtype=float))

        if peaks_array.size:
            peaks_valid = peaks_array[~np.isnan(peaks_array).any(axis=1)]
        else:
            peaks_valid = peaks_array
        if peaks_valid.size:
                # Alpha band
                alpha_params = extract_band_params(peaks_valid, bands_ranges[1])
                alpha_freq[elec_idx], alpha_ampl[elec_idx], alpha_widt[elec_idx] = alpha_params
            
                # Theta band
                theta_params = extract_band_params(peaks_valid, bands_ranges[0])
                theta_freq[elec_idx], theta_ampl[elec_idx], theta_widt[elec_idx] = theta_params
            
                # Beta band
                beta_params = extract_band_params(peaks_valid, bands_ranges[2])
                beta_freq[elec_idx], beta_ampl[elec_idx], beta_widt[elec_idx] = beta_params
            
                # Gamma band
                gamma_params = extract_band_params(peaks_valid, bands_ranges[3])
                gamma_freq[elec_idx], gamma_ampl[elec_idx], gamma_widt[elec_idx] = gamma_params    
        
    psd_dict = {
        'aperiodic': aperi_psd,
        'absolute': absol_psd,
        'fooofed': fooof_psd,
        'oscillatory': oscil_psd,
        'error': error_psd
    }

    if fit_pregroups != 0 and num_epochs > nepochs_group:
        totalepochs    = num_epochs
        epochs    = epochs_included
        epochsprop      = epochs_percentage
        interimr2      =r2_included
        interimr2exclu = r2_discard
        interimr2tresh = rsquare_threshold_interim
        perelectresh   = perelect 
        
        if fit_pregroups == 1 :
            nepochs_group = nepochs_group
            prefit         = 'prefit_bygroup'
        else :
            nepochs_group = 1    
            prefit        = 'prefit_byepoch'
    else:
        totalepochs    = np.size(psds, axis = 0)
        epochs         = np.size(psds, axis = 0)
        epochsprop      = 'NA'
        interimr2      = 'NA'
        interimr2exclu = 'NA'
        perelectresh   = perelect
        interimr2tresh = rsquare_threshold_interim
        nepochs_group  = nepochs_group
        prefit         = 'noprefit'
    
    parameters = {"interim_r2treshold": interimr2tresh,
           "nepochs_group": nepochs_group,
           "perelectresh": perelectresh}
       
    aperosc_results = {'freq': [fmax]*electrodes,
                       'ch': ch_names,
                         'prefit_specific': [prefit]*electrodes,
                         'nepochs_group': [nepochs_group]*electrodes,
                         'totalepochs': [totalepochs]*electrodes,
                         'epochs': [epochs]*electrodes,
                         'epochsprop': [epochsprop]*electrodes,
                         'interimr2': [interimr2]*electrodes,
                         'exclur2': [interimr2exclu]*electrodes,
                         'r2value': r2value, 
                         'mae': mae,
                         'offset': offset, 
                         'slope': slope, 
                         'theta_freq': theta_freq,
                         'theta_ampl': theta_ampl,
                         'theta_widt': theta_widt,
                         'alpha_freq': alpha_freq,
                         'alpha_ampl': alpha_ampl,
                         'alpha_widt': alpha_widt,
                         'beta_freq':  beta_freq,
                         'beta_ampl':  beta_ampl,
                         'beta_widt':  beta_widt,
                         'gamma_freq': gamma_freq,
                         'gamma_ampl': gamma_ampl,
                         'gamma_widt': gamma_widt,
                         'parameters_interim': [str(parameters)]*len(ch_names)}
                
    results_df = pd.DataFrame(aperosc_results)
    
    return results_df, psd_dict

def fooofpsd_extraction(fooofpsd, aperosc_results, foi, ch_names, parameters):
    aperosc_psd = []
    # Iterate over each channel and frequency
    for ch, ch_name in enumerate(ch_names):
          for f_idx, freq in enumerate(foi):  # Iterate over frequencies of interest
         # Construct a dictionary to hold values for this combination of channel and frequency
             row = {'ch': ch_name,  # Channel name
                    'freq': freq,  # Current frequency value
                    'epochs': aperosc_results['epochs'][ch],  # Final epochs for this channel
                    'r2value': aperosc_results['r2value'][ch],  # R² value for this channel
                    'mae': aperosc_results['mae'][ch]}
             
             row['absolute']    = fooofpsd['absolute'][ch, f_idx]
             row['aperiodic']   = fooofpsd['aperiodic'][ch, f_idx]
             row['fooofed']     = fooofpsd['fooofed'][ch, f_idx]
             row['oscillatory'] = fooofpsd['oscillatory'][ch, f_idx]
             row['error']       = fooofpsd['error'][ch, f_idx]
             row['parameters']  = str(parameters)
             
             aperosc_psd.append(row)      
    return aperosc_psd      

def interim_fooof_epochselection_chbased(psds, freqs, ch, fit_pregroups, fmin, fmax,
                                 peak_width_limits, max_n_peaks, min_peak_height, peak_threshold,
                                 rsquare_threshold_interim, bands_ranges, aperiodic_mode, nepochs_group=None):
    """
    Compute included epoch absolute power and select epochs (or groups) based on FOOOF fits.

    Parameters
    ----------
    psds : np.ndarray
        PSD absolute power array with shape (epochs, freqs), i.e. a SINGLE
        channel's spectra. This variant is per-channel; there is no electrode
        axis.

    freqs : np.ndarray
        Frequency vector corresponding to the PSD.
    fit_pregroups : int
        If 2, process each epoch individually; if 1, process epochs grouped by nepochs_group;
        otherwise return mean_abspower.
    fmin : float
        Lower frequency bound for FOOOF fitting.
    fmax : float
        Upper frequency bound for FOOOF fitting.
    peak_width_limits : iterable
        Parameter for FOOOF specifying minimum and maximum allowable peak widths.
    max_n_peaks : int
        Maximum number of peaks allowed in specparam.
    min_peak_height : float
        Minimum peak height for specparam.
    peak_threshold : float
        Peak threshold for specparam.
    rsquare_threshold_interim : float
        R² threshold for considering a epoch or group as valid.
    perelect : float
        Proportion of electrodes required to pass the r² threshold.
    nepochs_group : int, optional
        Number of epochs to average when processing groups (fit_pregroups==1).
    mean_abspower : np.ndarray, optional
        Mean absolute power to return if fit_pregroups is not 1 or 2.

    Returns
    -------
    included_psds : np.ndarray
        Averaged PSD of selected epochs (shape: electrodes x freqs).
    epochs_included : int or None
        Number of epochs (or groups) that were included.
    epochs_percentage : float or None
        Proportion of total epochs that were included.
    r2_included : float or None
        Average r² of the included (selected) fits.
    r2_discard : float or None
        Average r² of the discarded (non-selected) fits.
    """

    # Note: We assume psds.shape is (epochs, freqs)
    num_epochs = psds.shape[0]
    
    # ------------------------
    # Method 2: Process each epoch individually.
    
    lenpsd         = len(np.arange(fmin,fmax,freqs[1]-freqs[0])) + 1
        
    
    if fit_pregroups   == 0 or num_epochs <= nepochs_group:
        included_psds  = np.average(psds, axis=0)

    elif fit_pregroups == 2:
        r2_value = np.zeros(num_epochs)
        
        # Loop over each epoch and electrode to perform FOOOF fits.
        for epoch in range(num_epochs):
                psd_data = psds[epoch, :]
                fm = specparam.SpectralModel(
                    peak_width_limits=peak_width_limits,
                    max_n_peaks=max_n_peaks,
                    min_peak_height=min_peak_height,
                    peak_threshold=peak_threshold, aperiodic_mode=aperiodic_mode
                )
                # Use the provided frequency bounds.
                fm.fit(freqs=freqs, power_spectrum=psd_data, freq_range=[fmin, fmax])
                r2_value[epoch] = fm.r_squared_
                
        # Determine which epochs meet the r² criteria:
        r2_dummy     = (r2_value >= rsquare_threshold_interim)
        valid_epochs = r2_dummy == 1
        
        r2_included = np.average(r2_value[valid_epochs]) if np.any(valid_epochs) else np.nan
        r2_discard  = np.average(r2_value[~valid_epochs]) if np.any(~valid_epochs) else np.nan
        
        indices = np.where(valid_epochs)[0]
        epochs_included   = len(indices)
        epochs_percentage = epochs_included / num_epochs

        if epochs_percentage > 0:
            included_psds = np.average(psds[indices, :], axis=0)
        else:
            included_psds = np.average(psds[:, :], axis=0)
   
    # ------------------------
    # Method 1: Process epochs grouped by nepochs_group.
    elif fit_pregroups == 1:
        if nepochs_group is None:
            raise ValueError("nepochs_group must be provided when fit_pregroups==1")
        
        # Determine the number of groups.
        if num_epochs <= nepochs_group:
            ngroups =1 
        else:
            ngroups = int(np.round(num_epochs / nepochs_group))
            if (ngroups * nepochs_group) > num_epochs or (num_epochs - ngroups * nepochs_group) < nepochs_group:
                ngroups -= 1
        
        epoch1_loop = np.zeros(ngroups, dtype=int)
        epoch2_loop = np.zeros(ngroups, dtype=int)
        r2_value = np.zeros((ngroups))
        
        # Loop over groups.
        for grp in range(ngroups):
            epoch1 = nepochs_group * grp
            epoch2 = nepochs_group * grp + nepochs_group
            # Adjust the second index if necessary to avoid overflow.
            if (epoch2 + nepochs_group) >= num_epochs:
                epoch2 = num_epochs
            
            epoch1_loop[grp] = epoch1
            epoch2_loop[grp] = epoch2
            
            # For each electrode, average within the group and perform FOOOF fitting.
            psd_data_group = np.average(psds[epoch1:epoch2, :], axis=0)
            fm = specparam.SpectralModel(
                    peak_width_limits=peak_width_limits,
                    max_n_peaks=max_n_peaks,
                    min_peak_height=min_peak_height,
                    peak_threshold=peak_threshold, aperiodic_mode=aperiodic_mode)
            fm.fit(freqs=freqs, power_spectrum=psd_data_group, freq_range=[fmin, fmax])
            r2_value[grp] = fm.r_squared_
                
        # Check group passes criterion:
        r2_dummy     = (r2_value >= rsquare_threshold_interim)
        valid_groups =  r2_dummy == 1
        
   
        r2_included = np.average(r2_value[valid_groups]) if np.any(valid_groups) else np.nan
        r2_discard  = np.average(r2_value[~valid_groups]) if np.any(~valid_groups) else np.nan
        
        # Mark individual epochs for exclusion based on group-level analysis.
        selected_epochs = np.ones(num_epochs, dtype=bool)
        for grp in range(ngroups):
            if not valid_groups[grp]:
                pos1 = epoch1_loop[grp]
                pos2 = epoch2_loop[grp]
                selected_epochs[pos1:pos2] = False
        
        indices = np.where(selected_epochs)[0]
        
        epochs_included   = len(indices)
        epochs_percentage = epochs_included / num_epochs
        
        if epochs_percentage > 0:
            included_psds = np.average(psds[indices, :], axis=0)
        else:
            included_psds = np.average(psds[:, :], axis=0)

    aperi_psd = np.zeros((lenpsd))
    absol_psd = np.zeros((lenpsd))
    fooof_psd = np.zeros((lenpsd))
    oscil_psd = np.zeros((lenpsd))
    error_psd = np.zeros((lenpsd))
        
    # Initialize arrays for storing FOOOF results
    r2value = 0
    mae    = 0
    offset = 0
    slope  = 0
        
    # Initialize arrays for storing peak parameters (default to NaN)
    alpha_freq = np.full(1, np.nan)
    alpha_ampl = np.full(1, np.nan)
    alpha_widt = np.full(1, np.nan)
        
    theta_freq = np.full(1, np.nan)
    theta_ampl = np.full(1, np.nan)
    theta_widt = np.full(1, np.nan)
        
    beta_freq = np.full(1, np.nan)
    beta_ampl = np.full(1, np.nan)
    beta_widt = np.full(1, np.nan)
        
    gamma_freq = np.full(1, np.nan)
    gamma_ampl = np.full(1, np.nan)
    gamma_widt = np.full(1, np.nan)
        
     # Loop over electrodes to compute PSDs and extract FOOOF parameters
            # Extract PSDs for the current electrode
    psds_elec = included_psds[:]  # Use correct indexing
            
            # Initialize and configure FOOOF
    fm = specparam.SpectralModel(
                peak_width_limits=peak_width_limits,
                max_n_peaks=max_n_peaks,
                min_peak_height=min_peak_height,
                peak_threshold=peak_threshold
            )
            
            # Fit the FOOOF model on the electrode's PSD data
    fm.fit(freqs=freqs, power_spectrum=psds_elec, freq_range=[fmin, fmax])
        
            # Store overall results
    r2value = fm.r_squared_
    # See note in fooofpsd_extraction: fm.error_ is MAE (specparam default error_metric='mae').
    mae     = fm.error_
    offset  = fm.aperiodic_params_[0]
    slope  = fm.aperiodic_params_[1]
            
            # Extract components from the FOOOF model
    aperi_psd[:] = fm._ap_fit
    absol_psd[:] = fm.power_spectrum
    fooof_psd[:] = fm.modeled_spectrum_
    oscil_psd[:] = fm.power_spectrum - fm._ap_fit
    error_psd[:] = fm.power_spectrum - fm.modeled_spectrum_
            
        
            # Extract peak parameters (if peaks are present, otherwise assign NaN)
    peaks = fm.get_params('peak_params') 
    peaks_array = np.array(peaks)# Returns peaks as [CF, PW, BW] or an empty array
    if not np.any(np.isnan(peaks_array)):
                # Alpha band
                alpha_params = extract_band_params(peaks, bands_ranges[1])
                alpha_freq, alpha_ampl, alpha_widt = alpha_params
            
                # Theta band
                theta_params = extract_band_params(peaks, bands_ranges[0])
                theta_freq, theta_ampl, theta_widt = theta_params
            
                # Beta band
                beta_params = extract_band_params(peaks, bands_ranges[2])
                beta_freq, beta_ampl, beta_widt = beta_params
            
                # Gamma band
                gamma_params = extract_band_params(peaks, bands_ranges[3])
                gamma_freq, gamma_ampl, gamma_widt = gamma_params    
                
    psd_dict = {
        'aperiodic': aperi_psd,
        'absolute': absol_psd,
        'fooofed': fooof_psd,
        'oscillatory': oscil_psd,
        'error': error_psd
    }

    if fit_pregroups != 0 and num_epochs > nepochs_group:
        totalepochs    = num_epochs
        epochs    = epochs_included
        epochsprop      = epochs_percentage
        interimr2      =r2_included
        interimr2exclu = r2_discard
        interimr2tresh = rsquare_threshold_interim
        
        if fit_pregroups == 1 :
            nepochs_group = nepochs_group
            prefit         = 'prefit_bygroup'
        else :
            nepochs_group = 1    
            prefit        = 'prefit_byepoch'
    else:
        totalepochs    = np.size(psds, axis = 0)
        epochs    = np.size(psds, axis = 0)
        epochsprop      = 'NA'
        interimr2      = 'NA'
        interimr2exclu = 'NA'
        interimr2tresh = rsquare_threshold_interim
        nepochs_group  = nepochs_group
        prefit         = 'noprefit'
    
    parameters = {"interim_r2treshold": interimr2tresh,
           "nepochs_group": nepochs_group}
       
    aperosc_results = {'freq': [fmax],
                       'ch': ch,
                         'prefit': [prefit],
                         'totalepochs': [totalepochs],
                         'epochs': [epochs],
                         'epochsprop': [epochsprop],
                         'interimr2': [interimr2],
                         'exclur2': [interimr2exclu],
                         'r2value': r2value, 
                         'mae': mae,
                         'offset': offset, 
                         'slope': slope, 
                         'theta_freq': theta_freq,
                         'theta_ampl': theta_ampl,
                         'theta_widt': theta_widt,
                         'alpha_freq': alpha_freq,
                         'alpha_ampl': alpha_ampl,
                         'alpha_widt': alpha_widt,
                         'beta_freq':  beta_freq,
                         'beta_ampl':  beta_ampl,
                         'beta_widt':  beta_widt,
                         'gamma_freq': gamma_freq,
                         'gamma_ampl': gamma_ampl,
                         'gamma_widt': gamma_widt,
                         'parameters_interim': [str(parameters)]}
                
    results_df = pd.DataFrame(aperosc_results)
    
    return results_df, psd_dict

def fooofpsd_extraction_chbased(fooofpsd, aperosc_results, foi, ch_names, parameters):
        aperosc_psd = []
        # Iterate over each channel and frequency
        for f_idx, freq in enumerate(foi):  # Iterate over frequencies of interest
             # Construct a dictionary to hold values for this combination of channel and frequency
                 row = {'ch': ch_names,  # Channel name
                        'freq': freq,  # Current frequency value
                        'epochs': aperosc_results['epochs'][0],  # Final epochs for this channel
                        'r2value': aperosc_results['r2value'][0],  # R² value for this channel
                        'mae': aperosc_results['mae'][0]}
                 
                 row['absolute']    = fooofpsd['absolute'][f_idx]
                 row['aperiodic']   = fooofpsd['aperiodic'][f_idx]
                 row['fooofed']     = fooofpsd['fooofed'][f_idx]
                 row['oscillatory'] = fooofpsd['oscillatory'][f_idx]
                 row['error']       = fooofpsd['error'][f_idx]
                 row['parameters']  = str(parameters)
                 
                 aperosc_psd.append(row)      
        return aperosc_psd
    

#########################################
# Function to call PSD functions
#########################################

def process_psd(epochs, method, sfreq, **kwargs):
    """
    Runs one of the PSD computation methods.
    
    Parameters:
        epochs : MNE Epochs object.
        method : Integer selection:
                 1 = FFT with Hanning,
                 2 = Welch,
                 3 = Spectopo-like,
                 4 = Built-in multitaper,
                 5 = Manual multitaper.
        sfreq  : Sampling frequency.
        kwargs : Additional parameters for each method.
    
    Returns:
        freqs : Frequency bins.
        psds  : PSD for each epoch.
                *For methods 1,2,4,5 the output is in linear units (V²/Hz).
                *For method 3 (Spectopo), the output is in dB.
    """
    if method == 1:
        return compute_fft_hanning(epochs, sfreq)
    elif method == 2:
        return compute_welch(epochs, sfreq,
                             nperepoch=kwargs.get('nperepoch', None),
                             noverlap=kwargs.get('noverlap', None))
    elif method == 3:
        return compute_spectopo(epochs, sfreq,
                                nperepoch=kwargs.get('nperepoch', None),
                                noverlap=kwargs.get('noverlap', None))
    elif method == 4:
        return compute_multitaper(epochs, sfreq,
                                  time_bandwidth=kwargs.get('time_bandwidth', 4.0),
                                  fmin=kwargs.get('fmin', 0),
                                  fmax=kwargs.get('fmax', np.inf))
    elif method == 5:
        return compute_manual_multitaper(epochs, sfreq,
                                         time_bandwidth=kwargs.get('time_bandwidth', 4.0),
                                         num_tapers=kwargs.get('num_tapers', None))
    elif method == 6 :
        return compute_fft_hanning_normlength(epochs,sfreq)
        
    else:
        raise ValueError("Invalid method specified.")

def save_dataframe(df, subfolder_path, base_filename="lagged_coh_py", rewrite=1, ext=".csv"):
    """
    Save a pandas DataFrame to a CSV file, optionally versioning the filename.
    
    Parameters:
      df             : pandas DataFrame to save.
      subfolder_path : Path to the folder where the file will be saved.
      base_filename  : The base name for the file (default "lagged_coh_py").
      rewrite        : If 1, the file is saved as <base_filename><ext> (overwriting any existing file).
                       If 0, the function searches for existing files that start with <base_filename>_v<number><ext>
                       and saves the file with the next version (e.g., _v3) if _v1 and _v2 exist.
      ext            : The file extension (default ".csv").
    
    Returns:
      full_path      : The full file path where the DataFrame was saved.
    """
    
    if not os.path.exists(subfolder_path):
        raise FileNotFoundError(f"The specified folder {subfolder_path} does not exist.")
        
    if rewrite == 1:
        file_name = base_filename + ext
    else:
        # List all files in the folder
        files = os.listdir(subfolder_path)
        # Pattern to match files like "<base_filename>_v<number><ext>"
        pattern = re.compile(r"^" + re.escape(base_filename) + r"_v(\d+)" + re.escape(ext) + r"$")
        
        # Extract version numbers from matching files
        version_numbers = []
        for f in files:
            m = pattern.match(f)
            if m:
                version_numbers.append(int(m.group(1)))
        
        # Get the next version number (if none exist, it starts with 1)
        next_version = max(version_numbers, default=0) + 1
        file_name = f"{base_filename}_v{next_version}{ext}"

    full_path = os.path.join(subfolder_path, file_name)
    df.to_csv(full_path, index=False)
    
    return full_path
