function [peak_freq,x_phase] = pf_cohen(x,srate,varargin)

%   This function returns the peak frequency/frequency sliding of the
%   filtered data X. If X is a matrix, then pf_cohen operates along the columns
%   of X.

%   pf_cohen(x,n): specifies the order of the median filter. Default: 10
%   srate : Sampling rate (0.5 Hz for TR=2)
%   pf_cohen(x,[],n): Remove the first and last n time points of outputs to reduce the
%   inaccuracies from median filter, especially for fMRI.
%   Default:0
  
%   Reference:Cohen, M. X. (2014). Fluctuations in oscillation frequency control spike timing and coordinate neural networks. Journal of Neuroscience, 34(27), 8988-8998. doi:10.1523/JNEUROSCI.0261-14.2014

%   Time: 2023.07.11
%   Author: Yujia Ao

% get filter order
filt_order = 10;
if ~isempty(varargin)
  filt_order = varargin{1};
end

x_complex = hilbert(x);
x_phase = angle(x_complex); % phase in pi-border for output and visualizing
x_phase_cont = unwrap(x_phase); % contineous phase over pi-border for calculating derivative

% cohen method
dx_phase = diff(x_phase_cont); % calculating derivative of phase
freqslide = srate*dx_phase/(2*pi); % transfer to Instantaneous Frequency(Hz)

peak_freq = medfilt1(freqslide,filt_order,'omitnan','truncate'); % median filter %'zeropad'


% get time points to be removed
rm_time = [];
if numel(varargin)>1
  rm_time = varargin{2};
  peak_freq(1:1+(rm_time-1),:) = [];
  peak_freq(end-(rm_time-1):end,:) = [];
  x_phase(1:1+(rm_time-1),:) = [];
  x_phase(end-(rm_time-1):end,:) = [];
end


