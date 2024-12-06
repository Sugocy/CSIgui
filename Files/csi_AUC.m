function [area, peak_est] = csi_AUC(data, poi, bv_mask, doDisp)
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
%           doDisp          display spectrum, found peak and included
%                           points for both methods.
%           
% Output:
%           area            (1) the absolute-cumalitive summated area under
%                           the curve and (2) the area under the curve
%                           excluding values in the peak range below the
%                           baseline-value.
%           peak_est        calculated peak-location, related to the
%                           data vector.
%
% Q. van Houtum, PhD; version 1.0 10/2024;
% quincyvanhoutum@gmail.com

if nargin < 4, doDisp = 0; end
if nargin < 3 || isempty(bv_mask)
        bv_mask = round(size(data,1) * 0.35); 
end

% --- Peak estimation

% Baseline value
N = size(data,1); ax = 1:N; mask = round(bv_mask./2);
bv = mean(real(cat(1,data(1:mask),data(N-mask:N)))); 
bv = repmat(bv,N,1); 

% Peak-range-data (sub)
data_sub = real(data(poi(1):poi(2))); ax_sub = ax(poi(1):poi(2));
bvline_sub = bv(poi(1):poi(2));

% Peak maximum in range-data
[mx, mx_ind] = max(data_sub);

% Intersect bv-line and peak-range (Estimate)
[xest, ~] = polyxpoly(ax_sub, data_sub, ax_sub, bvline_sub);

% Intersection closest to peak-max location
xest_low = xest(xest < ax_sub(mx_ind));
[~, I] = sort(abs(xest_low  - ax_sub(mx_ind)));
if isempty(I)
    
    % First point to reach bv larger than given peak range
%     data_sub_mx = data(poi(1):mx_ind);
%     full_area_bool = data_sub_mx <= bv(1);
%     ind = find(full_area_bool == 1);
%     peak_est(1) = poi(1) + ind(1); 
%     poi(1) = peak_est(1);
    
    % Point after maximum closest to bv
    [~, closest] = min(abs(data_sub(poi(1):mx_ind) - bv(1)));
    if isempty(closest), closest = 0; end
    peak_est(1) = poi(1) + closest;
    
    
else
    peak_est(1) = xest_low(I(1))-1;  % Peak location estimate below
end


xest_high = xest(xest > ax_sub(mx_ind));
[~, I] = sort(abs(xest_high  - ax_sub(mx_ind)));
if isempty(I)
    
    % First point to reach bv larger than given peak range
%     data_sub_mx = data(mx_ind+poi(1):end);
%     full_area_bool = data_sub_mx <= bv(1);
%     ind = find(full_area_bool == 1);
%     peak_est(2) = mx_ind + poi(1) + ind(1) - 1; 
%     poi(2) = peak_est(2);

    % Point after maximum closest to bv
    [~, closest] = min(abs(data_sub(mx_ind:end) - bv(1)));
    peak_est(2) = poi(1) + mx_ind+closest;
    
    
else
    peak_est(2) = xest_high(I(1)); % Peak location estimate above
end

% --- Area-under-the-curve

% Peak location
peak_est = round(peak_est); 
xarea = ax(peak_est(1) : peak_est(2)); parea = real(data(xarea));

% Default
auc = cumtrapz(abs(parea));
if isempty(auc), area(1) = NaN; else, area(1) = auc(end); end

% bv-corrected
auc = cumtrapz(abs(parea(parea >= bv(1))));
if isempty(auc), area(2) = NaN; else, area(2) = auc(end); end

% --- Display
if doDisp
    % Plot data and baseline
    figure();  pd = plot(real(data)); hold on; pbv = plot(bv, '-.k');
    
    % Plot found peak x-axis
    pra1 = plot(xarea, parea, 'og', 'markersize', 5); 
    pra2 = plot(xarea(parea <= bv(1)), parea(parea < bv(1)),...
        'or', 'markersize', 6);      
    peak_range_y = bv(1) + 0.15.* [-mx mx];    
    
    % Plot peak-points included in AUC markers (default/corrected)
    pra0 = plot(repmat(peak_est(1),1,2),  peak_range_y, '--r');
    plot(repmat(peak_est(2),1,2), peak_range_y , '--r');   
    set(pd.Parent,'XLim', poi);
    
    % Legend
    legend([pd, pbv, pra0, pra1, pra2],...
        'Data', 'Base-value', 'peak-range',...
        'Included (default)', 'Excluded (corrected)', 'Location',...
        'NW');
    title(sprintf('AUC default | corrected: %2.2e | %2.2e', area));
end