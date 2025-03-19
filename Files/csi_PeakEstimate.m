function [peak_est, bv] = csi_PeakEstimate(data, poi, bv_mask, strict, doDisp)
% Find an estimate of the peak in data in the range poi - with freedom to
% increase the poi-range for improved peak-estimation.
%
%
% Im an idiot... csi_findPeaks already existed.
% Even bigger idiot: findpeaks is a Matlab function designed like this
% script. 
%
% HOWEVER, it doesnt give proper x-locations of the peak...
%
% Im not an idiot at all, actually quite smart what I did :P Was not as
% fast - now faster (up to 4-times and more..)

if nargin < 5, doDisp = 0; end
if nargin < 4, strict = 0; end
if nargin < 3 || isempty(bv_mask), bv_mask = round(size(data,1) * 0.35); 
end

%% --- Baseline

% Baseline value
N = size(data,1); ax = 1:N; mask = round(bv_mask./2);
bv = mean((cat(1,data(1:mask),data(N-mask:N)))); bv = repmat(bv,N,1); 

% POI rang vector.
poir = poi(1):poi(2);

% Maximum in peak-range - this is the "peak" of interest
[mx, mx_loc] = max(data(poir)); mx_loc = mx_loc + poi(1) - 1 ; 

%% --- Peak estimation I: intersect with base value
% Check for intersections with base value within given peak range. If
% strict calculations are off - increase peak range until an intersection
% is found, limiting growth on either side of the maximum by half the
% original peak range.

% This has been compared to an all-data-at-once approach, but the latter
% requires up to 4-times as much computation time!
xest = []; poi_tmp = poi; lStep = -1; rStep = 1; n = 0; dpoi = diff(poi)./2;
while isempty(xest) || sum(xest == 0) >= 1 || isscalar(xest)
    
    % Cut data
    poir = poi_tmp(1):poi_tmp(2);
    data_sub = data(poir); ax_sub = ax(poir); bvline_sub = bv(poir);

    % Find intersection with base-value
    xest_tmp = polyxpoly(ax_sub, data_sub, ax_sub, bvline_sub);
    
    % Check if xest is found on either side of the peak maximum.
    % Set step-size of side with an intersect to 0 if found.
    if ~isempty(xest_tmp)
        % Find intersection before and after maximum in original range
        Il = find((xest_tmp < mx_loc), 1, 'last');          
        Ir = find((xest_tmp > mx_loc), 1, 'first');      
        % Save and set step value to 0.
        if ~isempty(Il) ~= 0, xest(1) = xest_tmp(Il); lStep = 0; end           
        if ~isempty(Ir) ~= 0, xest(2) = xest_tmp(Ir); rStep = 0; end
    end    

    % Break if new maximum in range.
    [mx_tmp, mx_ind_tmp] = max(data(poir)); 
    mx_ind_tmp = mx_ind_tmp + poi_tmp(1) - 1;
    if ~isequal(mx_ind_tmp, mx_loc) 
        % Only if the difference in maxima is more than 5%
        % and 4 data points away - we "label" it as a true new maximum        
        if (mx_ind_tmp - mx_loc) <= 3 && ((mx_tmp - mx)./mx_tmp) <= 0.05
            mx_loc = mx_ind_tmp; mx = mx_tmp;
        else
            % If it is a new maximum - stop.
            if mx_ind_tmp < mx_loc
                poi_tmp(1) = mx_ind_tmp + 1; 
                if numel(xest) >= 1, xest(1) = poi_tmp(1); end
            elseif mx_ind_tmp > mx_loc
                poi_tmp(2) = mx_ind_tmp - 1; 
                if numel(xest) == 2, xest(2) = poi_tmp(2); end
            end
        end
    end

    % Limit increase of peak range on either side of peak by half original
    % peak range size. Set grown poi-range to x-estimate; If strict
    % calculations is on - stop after a single iteration.
    if n > dpoi || strict
        if lStep ~= 0, xest(1) = poi_tmp(1); end
        if rStep ~= 0, xest(2) = poi_tmp(2); end
        break;     
    end

    % Increase poi-range
    poi_tmp = poi_tmp + [lStep rStep];

    % Count iterations
    n = n + 1;
end

if doDisp
    fprintf('\n@AUC - original peak-range: %i:%i.\n', poi);
    fprintf('@AUC - updated peak-range (intersect): %2.2f:%2.2f.\n', xest);
end      


%% --- Peak estimation II: local minima within range
% check local minima within peak range - if prominence is large enough -
% this could mean a "split" peal - and the minima should be used as the
% peak-estimate.
%
% Peak prominence set to >20% of maximum prominence possible within range.

% From the previous peak-estimation, xest is always created.
poi_user = poi; poi = [floor(xest(1)) ceil(xest(2))];

% figure(); plot(data); hold on; plot(xest, data(poi), 'or')

% Local minima within (new) poi-range
poir = poi(1):poi(2); data_sub = data(poir); 
[~, mx_ind] = max(data_sub); % Peak maximum in range-data

% Include prominence to find throughs/valleys with significant depth 
% within range
minProm = 0.2 .* (max(data_sub) - min(data_sub));
[loc_min, prom] = islocalmin(data_sub, 'MinProminence', minProm);
loc_min = find(loc_min == 1); prom = prom(loc_min);   

% Before and after maximum
prom_oi = {prom(loc_min < mx_ind), prom(loc_min > mx_ind)};

for kk = 1:numel(prom_oi)
    if ~isempty(prom_oi{kk})
        % Local minima of interest
        loc_min_oi = loc_min(max(prom_oi{kk}) == prom) + (poi(1) - 1);

        % Current peak-range estimate vs local-minima:         
        ind_opt = [loc_min_oi poi(kk)];

        % Unused: closest to BV
        % [~, ind] = min(abs( [data(ind_opt)] - bv(1) )); % closest to bv

        % Within user-given peak range
        if kk == 1,     withinUserRange = ind_opt >= poi_user(kk);
        elseif kk == 2, withinUserRange = ind_opt <= poi_user(kk);
        end
            
        % If local minimum not within user-range - use user-value
        if ~withinUserRange(1)            
            xest(kk) = poi(kk); % Use user value
        else            
            xest(kk) = loc_min_oi; % Use local minimum. 
        end

    end
end       

% Data and maximum-location within (new) poi-range
poi = [floor(xest(1)) ceil(xest(2))]; poir = poi(1):poi(2);
data_sub = data(poir); ax_sub = ax(poir); [~, mx_ind] = max(data_sub); 

if doDisp
fprintf('@AUC - updated peak-range (local minima):  %3.2f:%3.2f.\n', xest);
end



% Maximum possible prominence in peak-only-data
prom_max = max(data_sub) - min(data_sub);

% Local min/max
[locmin, prommn] = islocalmin(data_sub); locmin_ind = find(locmin == 1);
[locmax, prommx] = islocalmax(data_sub); locmax_ind = find(locmax == 1);

% Visualize
% figure();plot(bv(ax_sub), '--k'); hold on; plot(data_sub); 
% plot(locmin_ind, data_sub(locmin_ind), 'or')
% plot(locmax_ind, data_sub(locmax_ind), 'og')

% If more local maxima than minima - we have "mini" peaks.

if numel(locmax_ind) > numel(locmin_ind)
    minind = (locmin_ind + poi(1) - 1);
    maxind = (locmax_ind + poi(1) - 1);
    
    % Closest to maximum
    [mx, mx_ind] = max(data_sub);

    % Ease of indexing
    isBeforeMx = (locmin_ind < mx_ind); isAfterMx = (locmin_ind > mx_ind);

    % Near base value
    abs(data_sub(locmin_ind(isBeforeMx)) - bv(1));
    abs(data_sub(locmin_ind(isAfterMx)) - bv(1));
    
end





%% --- Set Output

% Set Output
peak_est = xest;

