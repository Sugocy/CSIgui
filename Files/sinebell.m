function sb_window = sinebell(L, shift)
% Returns the sinebell window for filtering of NMR spectra.
% Causes narrow lines but induces strong oscillations and 
% disappearing integrals.
%
% Input L:  window size e.g. sample size of signal.
%   shift:  if given, the sinebell is shifted, default shift = 0;
%           Shift is defined as numbers of data samples to shift.
%
% Source information:
% sb = @(t, tmax) sin(pi.*t ./ tmax);
% https://www.chemie.uni-hamburg.de/nmr/insensitive/tutorial/
% en.lproj/window_function.html
%
% Shifted SB (shb): See sine-bell, but added shift over time-dimensions.
%                   sin(pi*(t+b)/(AQ+b)). b = shift, AQ = end point of bell
%                   on time axis.
%
% Contact: qhoutum2@umcutrecht.nl

if nargin == 1, shift = 0; end

% Formula SineBell 1.
% sb = @(t, tmax) sin(pi.*t ./ tmax);
% Formula SineBell 2.
sb = @(t, tmax, b) sin(pi.*(t+b)./(tmax+b));

% Window size
t = 0:L-1;
% Calc window
sb_window = sb(t, max(t), shift);

