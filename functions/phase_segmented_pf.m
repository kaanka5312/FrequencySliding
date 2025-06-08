function PF_struct = phase_segmented_pf(x_phase, peak_freq)
% PHASE_SEGMENTED_PF handles matrix input where each column is one region.
%
% Inputs:
%   x_phase    - T x N matrix of phase (unwrapped Hilbert phase)
%   peak_freq  - (T-1) x N matrix of instantaneous frequency
%
% Output:
%   PF_struct  - struct with fields, each 1 x N (region-wise averages):
%       PF_peak, PF_trough, PF_rise, PF_fall, PF_whole

    [~, N] = size(peak_freq);
    
    % Preallocate output
    PF_struct.PF_peak   = nan(1, N);
    PF_struct.PF_trough = nan(1, N);
    PF_struct.PF_rise   = nan(1, N);
    PF_struct.PF_fall   = nan(1, N);
    PF_struct.PF_whole  = nan(1, N);

    for r = 1:N
        x = x_phase(:, r);
        pf = peak_freq(:, r);

        % Phase bins (shift by -1 for PF indexing)
        idx_peak = find(x >= -pi/12 & x <= pi/12) - 1;
        idx_trough = find((x >= 11*pi/12 & x <= pi) | (x >= -pi & x <= -11*pi/12)) - 1;
        idx_rise = find(x >= -7*pi/12 & x <= -5*pi/12) - 1;
        idx_fall = find(x >= 5*pi/12 & x <= 7*pi/12) - 1;

        idx_peak(idx_peak < 1) = [];
        idx_trough(idx_trough < 1) = [];
        idx_rise(idx_rise < 1) = [];
        idx_fall(idx_fall < 1) = [];

        % Store per-region PF means
        PF_struct.PF_peak(r)   = mean(pf(idx_peak));
        PF_struct.PF_trough(r) = mean(pf(idx_trough));
        PF_struct.PF_rise(r)   = mean(pf(idx_rise));
        PF_struct.PF_fall(r)   = mean(pf(idx_fall));
        PF_struct.PF_whole(r)  = mean(pf);
    end
end
