function filtered_ts = bandpass_cheby1(ts, f_low, f_high, srate)
% BANDPASS_CHEBY1 applies a zero-phase Chebyshev Type I bandpass filter
% Inputs:
%   ts      - time series (1D or 2D matrix: T x N, column vector for 1)
%   f_low   - lower cutoff frequency in Hz (e.g., 0.01)
%   f_high  - upper cutoff frequency in Hz (e.g., 0.073)
%   srate   - sampling rate in Hz (e.g., 1/3 for TR = 3s)
% Output:
%   filtered_ts - bandpass filtered time series

    order = 4;            % Filter order
    rp = 0.1;             % Passband ripple (in dB)

    % Normalize frequencies to Nyquist (fs/2)
    Wn = [f_low f_high] / (srate / 2);  % normalized cutoff frequencies

    % Design Chebyshev Type I filter using second-order sections
    [sos, g] = cheby1(order, rp, Wn, 'bandpass');
    
    % Apply zero-phase filtering (forward & backward)
    filtered_ts = filtfilt(sos, g, ts);
end
