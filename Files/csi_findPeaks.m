function [peak_position, peak_values] = csi_findPeaks(spec, nPeaks)
%%%% Description:                               Find peaks in a spectrum
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2018-07
%%% --------------------------------------------------------------------
%%%
%%% Uses the real part of the spectrum to find local maxima e.g. peaks
%%% in the spectrum.
%%% 
%%% nPeaks = number of peaks, default is set to 10;


% Process input
if nargin == 1, nPeaks = 10; end

% Real part of data
rp = real(spec);

%Plot it
% figure(); hold on; plot(rp); hold on;

% Prominence minimum
mrp = abs(std(diff(rp))); 


% % Find peaks e.g. local maxima
% [loc_max_bool ]= islocalmax(rp,'MaxNumExtrema',nPeaks,...
%      'MinSeparation',0);
 [loc_max_bool ]= islocalmax(rp,'MaxNumExtrema',nPeaks,'MinProminence',mrp,...
     'MinSeparation',0);
% Convert boolean to linear indexes
[peak_position,~] = find(loc_max_bool == 1); 
% Values at position
peak_values = rp(peak_position);

%Plot found maxima
%plot(peak_position, peak_values,'or');
