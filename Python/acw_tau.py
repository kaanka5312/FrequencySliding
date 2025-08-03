import pickle
import numpy as np
from statsmodels.tsa.stattools import acf
from scipy.integrate import simpson
from scipy import optimize

def autocorr_decay(l, A, tau, B):
    return A * (np.exp(-(l / tau)) + B)

def calc_int(ts, fs):
    n = len(ts)
    i_acf = acf(ts, nlags=n - 1)
    lags = np.arange(0, n)
    acw_50_samples = np.argmax(i_acf <= 0.5)
    acw_0_samples = np.argmax(i_acf <= 0)
    acw_1_over_e_samples = np.argmax(i_acf <= 1 / np.e)

    acw_50 = acw_50_samples / fs
    acw_0 = acw_0_samples / fs
    acw_1_over_e = acw_1_over_e_samples / fs

    acf_before_0 = i_acf[:acw_0_samples]
    acw_integral = simpson(acf_before_0, dx= 1 / fs)

    acf_first_100 = i_acf[:101]  # As in Ito et al paper
    try:
        A, tau, B = optimize.curve_fit(
            autocorr_decay,
            lags[1:101],
            acf_first_100[1:],
            p0=[0, np.random.rand(1)[0] + 0.01, 0],
            bounds=(([0, 0, -np.inf], [np.inf, np.inf, np.inf])),
            method="trf",
        )[0]
    except:
        tau = np.nan
    
    return acw_50, acw_0, acw_integral, acw_1_over_e, tau

def calc_int_array(ts, fs):
    """
    Assumes ts is times x channels
    """
    n_times, n_channels = ts.shape
    int_array = np.zeros((5, n_channels))
    for i_channel in range(n_channels):
        int_array[:, i_channel] = calc_int(ts[:, i_channel], fs)
    return int_array

import numpy as np

def acw_slidingwindow(eeg, fs, window, overlap=50, lag=None):
    """
    Python equivalent of the MATLAB ACW_kaan function.

    Parameters
    ----------
    eeg : 1d array-like
        EEG time course for a single channel.
    fs : float
        Sampling frequency in Hz.
    window : float
        Window length in seconds. If zero, uses the full time series as one block.
    overlap : float, optional
        Percentage overlap between consecutive windows (default 50).
    lag : float or None, optional
        Maximum lag in seconds to compute autocorrelation. If None, defaults to (N-1)/fs.

    Returns
    -------
    ACW0 : float
        Width until first zero crossing (as implemented in original).
    ACW50 : float
        Full-width-at-half-maximum of the main lobe of the autocorrelation.
    ACF_mean : 1d numpy array
        Averaged autocorrelation function across windows (lags from -lag to +lag).
    time : 1d numpy array
        Time axis corresponding to ACF_mean (in seconds).
    """
    eeg = np.asarray(eeg, dtype=float)
    if window is None:
        raise ValueError("window must be provided (in seconds)")
    if overlap is None:
        overlap = 50.0
    if lag is None:
        lag = (len(eeg) - 1) / fs

    # samples
    window_samps = int(round(window * fs)) if window != 0 else len(eeg)
    if window_samps <= 0:
        raise ValueError("Window length in samples must be positive.")
    lag_samps = int(np.floor(lag * fs))

    step = int(window_samps * overlap / 100.0)
    if step <= 0:
        raise ValueError("Overlap too small or window leads to non-positive step size.")

    acf_list = []
    start = 0
    while True:
        end = start + window_samps
        if end >= len(eeg):
            break
        segment = eeg[start:end]
        x = segment - np.mean(segment)
        denom = np.dot(x, x)
        if denom == 0:
            acf = np.zeros(2 * lag_samps + 1)
        else:
            full = np.correlate(x, x, mode="full") / denom  # normalized so zero-lag is 1
            center = len(full) // 2
            acf = full[center - lag_samps : center + lag_samps + 1]
        acf_list.append(acf)
        start += step

    if len(acf_list) == 0:
        raise ValueError("No windows were generated: check signal length vs window/overlap.")

    ACF = np.vstack(acf_list)
    ACF_mean = ACF.mean(axis=0)

    # ACW50 (full-width at half max of main lobe)
    peak_idx = np.argmax(ACF_mean)
    half_max = ACF_mean[peak_idx] / 2.0
    above_half = ACF_mean >= half_max
    right_side = above_half[peak_idx:]
    under_half_idxs = np.where(~right_side)[0]
    if under_half_idxs.size == 0:
        first_under_half = len(right_side) + 1  # mimic absence handling
    else:
        first_under_half = under_half_idxs[0] + 1  # +1 to match MATLAB's 1-based logic

    acw50_samps = 2 * (first_under_half - 1) - 1
    ACW50 = acw50_samps / fs

    # ACW0 (first zero crossing on right of peak)
    negative = ACF_mean < 0
    right_neg = negative[peak_idx:]
    neg_idxs = np.where(right_neg)[0]
    if neg_idxs.size == 0:
        acw0_index = window_samps + 1  # as in original: if not found
    else:
        acw0_index = neg_idxs[0] + 1

    acw0_samps = 2 * (acw0_index - 1) - 1
    ACW0 = acw0_samps / fs

    time = np.linspace(-lag_samps / fs, lag_samps / fs, 2 * lag_samps + 1)

    return ACW0, ACW50, ACF_mean, time

def acw_slidingwindow_array(ts, fs, window):
    """
    ts: times x channels
    Returns a dict with per-channel ACW0, ACW50, and the full ACF_mean/time.
    """
    n_times, n_channels = ts.shape
    acw0 = np.zeros(n_channels)
    acw50 = np.zeros(n_channels)
    acf_means = []
    time_vector = None

    for i in range(n_channels):
        acw0[i], acw50[i], acf_mean, time = acw_slidingwindow(ts[:, i], fs, window)
        acf_means.append(acf_mean)
        if time_vector is None:
            time_vector = time  # same for all channels

    acf_means = np.stack(acf_means, axis=1)  # shape: (lags, channels)

    return {
        "ACW0": acw0,                # (n_channels,)
        "ACW50": acw50,              # (n_channels,)
        "ACF_mean": acf_means,       # (2*lag+1, n_channels)
        "time": time_vector          # (2*lag+1,)
    }
