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