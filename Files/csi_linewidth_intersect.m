function linewidth = csi_linewidth_intersect(data, ax, peak_range)
%%% Calculate linewidth of a peak at a given frequency peak_range using an
%%% intersect FWHM method.


% Baseline | Samples
N = size(data,1); 
% Baseline
bv = mean(real([data(1:round(N/20)); data(N-round(N/20):N)]));  % Value
% bvline = cat(2,ax',repmat(bv,N,1));                           % Line

% Maximum
data_peak_ranged = real(data(peak_range(1):peak_range(2)));
[mv, mi] = max(data_peak_ranged);

% Set maximum as peak-centre
peak_pos = mi+peak_range(1)-1; % Data index
xv = ax(peak_pos);             % In unit (ppm/unitless etc.)
ax_peak_ranged = ax(peak_range(1):peak_range(2));

% FWHM
fwhm = mv - ((mv-bv)/2);
fwhmline = cat(2,ax',repmat(fwhm,N,1));


% Intersect FWHM line and DATA (peak_range) (Estimate)
[xest,yest] = polyxpoly(ax_peak_ranged',data_peak_ranged, ...
    fwhmline(:,1)', fwhmline(:,2)');

% The answer in PPM or unit-less: ALWAYS use the first two - if more
% intersects are found, data is wrong anyways... 
if size(xest,1) == 1 % The peak is not within the peak_range probabaly
    linewidth = NaN; 
else
    xest = [xest(1) xest(end)]; % yest = [yest(1) yest(end)];
    linewidth = diff(xest); 
end