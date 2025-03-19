function [area, peak_est] = csi_AUC(data, poi, bv_mask, strict, doDisp)
% Estimate the area under the curve of a peak definied by peak-range
% (unit-less) and the peak-maximum plus intersection with a base-line
% defined by baseline-value mask size (integer).
%
% Input:    data            full spectrum [1 x N];
%           poi             range of peak in data [1 x 2], unitless; 
%           bv_mask         baseline value mask-size: the number of points
%                           to determine the base-line value using the two 
%                           tails of the spectrum. If empty, 35% of the
%                           data-vector its length.
%           strict          Allow the peak-range to increase (0, default) 
%                           if no intersection with the base-value within 
%                           poi is found. Max peak-range increase is 2x
%                           given peak-range.
%           doDisp          display spectrum, found peak and included
%                           points for both methods.
%           
% Output:
%           area            the absolute-cumalitive summated area under
%                           the curve
%           peak_est        calculated peak-location, related to the
%                           data vector.
%
% Method:
% First, the baseline value (bv) is calculated using the bv_mask size. The
% x-axis location of intersection with bv is estimated using polyxpoly. If
% no or only a single intersect value is found - the peak range will
% iteratively be increased until two intersect values are found; unless
% strict is on (1). 
% A second peak-estimate is performed by checking if the change of the
% slope of the peak-values is significantly decreased e.g. almost a
% straight line - without crossing/touching the bv. The point in this area
% closest to the bv will be used as the peak-range estimate.
%
%
%
% Q. van Houtum, PhD; version 1.0 10/2024;
% quincyvanhoutum@gmail.com

if nargin < 5, doDisp = 0; end
if nargin < 4, strict = 0; end
if nargin < 3 || isempty(bv_mask), bv_mask = round(size(data,1) * 0.35); 
end




% Peak-estimation using given peak-of-interest range.
if ~strict
    [peak_est, bv] = csi_PeakEstimate(data, poi, bv_mask, strict, doDisp);
else
    peak_est = poi;
end

%% --- Area-under-the-curve
% Interpolate the data within the peak range by 10% and calculate the area
% under the curve via cumptrapz.

% Interpolate between (found) peak-range
N = size(data,1); ax = 1:N;
peak_ind = floor(peak_est(1)):ceil(peak_est(2)); 
dx = mean(diff(ax))./10; ax_interp = peak_ind(1):dx:peak_ind(end);
data_interp = interp1(ax(peak_ind), real(data(peak_ind)), ax_interp);

% x-axis peak-estimates closest to interpolated x-axis
[~, peak_est_interp(1)] = min(abs( peak_est(1) - ax_interp ));
[~, peak_est_interp(2)] = min(abs( peak_est(2) - ax_interp ));

% Interpolated data
xarea = ax_interp(peak_est_interp(1):peak_est_interp(2));
parea = data_interp(peak_est_interp(1):peak_est_interp(2));

% AUC: Absolute
auc{1} = cumtrapz(xarea, abs(parea));
% AUC: bv-corrected
auc{2} = cumtrapz(xarea(parea >= bv(1)), abs( parea(parea >= bv(1)) ));

% Safety for output generation
for kk = 1:2, if isempty(auc{kk}), auc{kk} = 0; end; end

% Set output.
try
    area = cell2mat(cellfun(@(x) x(end), auc, 'UniformOutput', 0));
catch
    area = nan;
end
   

% --- Display
if doDisp
    % Plot data and baseline
    figure();  pd = plot(data); hold on; pbv = plot(bv, '-.k');
    
    % Plot found peak x-axis
    mx = max(data(peak_ind));
    pra1 = plot(xarea, parea, 'og', 'markersize', 5); 
    pra2 = plot(xarea(parea <= bv(1)), parea(parea < bv(1)),...
        'or', 'markersize', 6);      
    peak_range_y = bv(1) + 0.15.* [-mx mx];    
    
    % Plot peak-points included in AUC markers 
    pra0 = plot(repmat(peak_est(1),1,2),  peak_range_y, '--r');
    plot(repmat(peak_est(2),1,2), peak_range_y , '--r');   
    % set(pd.Parent,'XLim', poi);
    
    % Legend.
    legend_objects = {pd, pbv, pra0, pra1, pra2};
    legend_objects_empty = cellfun(@isempty, legend_objects);
    legend_descriptions = {'Data', 'Base-value', 'peak-range',...
        'Included (default)', 'Excluded (corrected)'};
    legend_descriptions(legend_objects_empty) = [];
    legend_objects(legend_objects_empty) = [];
    legend([legend_objects{:}], legend_descriptions{:}, 'Location', 'NW');

    % Title
    title(sprintf('AUC default | corrected: %2.2e | %2.2e', area));
end