function MF = mf(x, fs, freqrange)
%% Median Frequency
% The frequency that divides the power spectrum to two equal halves. 
% Input:
%   x: time series (1D vector)
%   fs: sampling frequency (scalar)
%   freqrange: frequency of interest (1x2 vector, like [1 50] for EEG).
% Output:
%   MF: median frequency

% The beginning is the same as the PLE function:
% pwelch uses welch method that is useful when you have many sampling
% points (like EEG / MEG data) since it uses a sliding window approach when
% calculating PSD (take PSDs of chunks and average them, see the 
% documentation of pwelch for details). If your data has low number of time 
% points, you might want to switch to periodogram. A useful rule of thumb:
% your data should cover 2-3 full cycles of the lowest frequency you are
% interested in (for example, if your lowest frequency is 1 Hz, your chunks
% should be 2-3 seconds). 

[psd, freq] = pwelch(x, [], [], [], fs); 
% [psd, freq] = periodogram(x, [], [], fs);

% Next step is taking the frequencies you are interested in. The Nyquist
% for EEG data can go up to 500 Hz, but neuronally those frequencies are
% not thought to be very informative. 

foi = (freq > freqrange(1)) & (freq < freqrange(2));
freq2 = freq(foi);
psd2 = psd(foi);
    
cum_sum = cumsum(psd2); 
[~, where] = max(cum_sum >= (cum_sum(end) / 2));
MF = freq2(where);
end