function [SEM] = SEM_time_series(Signal,per_autocorr_decrease)
% SEM of time series (which has coorealtions)

% Inputs:
% Signal - the time series to analyze
% per_autocorr_decrease - drop value of autocorrelation to estimate the actual variance

% Outputs:
% SEM - corresponding standart error of the mean

%% The SEM
Signal = Signal-mean(Signal); % subtracting average

[AC,l] = xcorr(Signal,Signal,'coeff'); % autocorrelation normalized
AC = AC*var(Signal,1); % multiply autocorrelation by the variance of the signal

z0 = find(l==0); % 0 auto-correlation index = variance
Z = z0; % will count index which we drop beloew per_autocorr_decrease

Thresh_per = AC(z0); % will contian sum of norms autocorrealation (starts its variance)
VS = Thresh_per; % starting value

while VS>per_autocorr_decrease*Thresh_per % checking at which point we decrease per_autocorr_decrease %
    Z = Z+1;
    VS = AC(Z);
end

VAR = AC(z0) + 2*sum(AC(z0+1:Z-1)); % actual variance
SEM = sqrt(VAR/numel(Signal)); % SEM

end