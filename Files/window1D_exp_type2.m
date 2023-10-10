function window = window1D_exp_type2(L, EL, ES)
% Create expontential apodization window given window size L and
% exponential length (EL) and strength (EL) or stretch.
%
% If only window size is given, exponential length is set to 4 and
% strength, e.g. stretched exponential is set to 1.


% Returng if L is zero or negative
if L <= 0, return; end

if nargin == 1, EL = 4; ES = 1; end

% Define x-axis
x = 0:(EL/L):EL;

% Calculate exp values
tmp = exp(-1.*x.^ES);
tmp = tmp./(max(tmp));
    
% Plot window
% figure(); plot(tmp); hold on;

% Output.
window = tmp;