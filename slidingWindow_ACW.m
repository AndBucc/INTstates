function [result, result_mean] = slidingWindow_ACW(EEG,fs,window,overlap)
%
%   Input:
%       -EEG:       time course of the EEG in a single channel
%       -FS:        sample frequency or sample rate(in Hz)
%       -WINDOW:    size (in seconds) of the windows or block to compute
%                   the AC function. Zero means that the block is all the
%                   time serie. Zero is the default.
%       -OVERLAP:   % of overlapping. 50% by default.
%
%   Output:
%       -RESULT: Vector with the values of each window
%       -RESULT_MEAN: Average across windows (depends on the measure)
%
%
% Version: 1.0
%
% Date: September 12, 2017
% Last update: September 12, 2017
%


%% Set initial defaults (user-specified)
if nargin<3 || isempty(window),
    error('Not enough input arguments.');
end
if nargin<4 || isempty(overlap),
    overlap=50;
end


% Time to samples
window=window*fs;


% Sliding window for computing what you want
ii=1;   % Windows counter
while true
    % Begining and ending of the current window
    SWindow=[1+(ii-1)*window*(100-overlap)/100, window+(ii-1)*window*(100-overlap)/100];

    % Check if index exceeds vector dimensions. If so, break!
    if SWindow(2)>=length(EEG), break; end

    % Measure computation into the window
    [result(ii), ~, ~, ~] = yasir_acw(EEG(SWindow(1):SWindow(2)),fs);


    % Next window
    ii=ii+1;
end

% for convenience, we also compute the mean of all the values in each window
result_mean=nanmean(result);


