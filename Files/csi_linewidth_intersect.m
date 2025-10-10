function [linewidth, xest, yest] = csi_linewidth_intersect(data, ax, peak_range, doPlot)
%%% Calculate linewidth of a peak at a given frequency peak_range using an
%%% intersect FWHM method.

if nargin < 4, doPlot = 0; end

try

% Baseline | Samples
N = size(data,1); 
% Baseline
bv = mean(real([data(1:round(N/20)); data(N-round(N/20):N)]));  % Value
bvline = cat(2,ax',repmat(bv,N,1));                          % Line


% -------- Range correction via peak-maximum

% Center the peak-range around maximum inside the given peak-range.
[~, mi] = max(real(data(peak_range(1):peak_range(2))));
% Set max as peak position.
peak_pos = mi+peak_range(1);

% Peak width: keep same width as given range
peak_width = round(diff(peak_range)./2);
if peak_width < 3, peak_width = 3; end

% Set range by adding and subtracting peak width from peak position.
peak_range = peak_pos + [-peak_width peak_width];

% -----------------------------------------

% Recalculate maximum and location in updated range
data_peak_ranged = real(data(peak_range(1):peak_range(2)));
[mv, mi] = max(data_peak_ranged);

% Set maximum as peak-centre
% peak_pos = mi+peak_range(1)-1; % Data index
ax_peak_ranged = ax(peak_range(1):peak_range(2));
xv = ax_peak_ranged(mi);  % In unit (ppm/unitless etc.)

% FWHM
fwhm = mv - ((mv-bv)/2);
fwhmline = cat(2,ax',repmat(fwhm,N,1));

% Intersect FWHM line and DATA (peak_range) (Estimate)
[xest,yest] = polyxpoly(ax_peak_ranged',data_peak_ranged, ...
    fwhmline(:,1)', fwhmline(:,2)');

% The answer in PPM or unit-less: ALWAYS use the first two - if more
% intersects are found, data is wrong anyways... 
if size(xest,1) == 1 || size(xest,1) == 0 % The peak is not within the peak_range probabaly
    linewidth = NaN; xest = NaN; yest = NaN;
else
    xest = [xest(1) xest(end)]; % yest = [yest(1) yest(end)];
    linewidth = diff(xest); 
end

% Plot results
if doPlot
   
    % Create new figure for user input
    fh = figure(); 
    % Plot the spectrum (reverse xaxis)
    spObj = plot(ax, real(data)); spObj.Parent.XDir = 'reverse'; hold on;
    % Axes
    axObj = spObj.Parent;
    % PLOT BASELINE
    baseObj = plot(axObj, bvline(:,1),bvline(:,2),'--k');
    % datapoint
    cursorMode = datacursormode(fh); hdtip = cursorMode.createDatatip(baseObj);
    hdtip.Position = [bvline(1,1) bv 0]; updateDataCursors(cursorMode);

    % Plot Maximum
    plot(axObj, xv,  mv, 'or');

    % Plot FWHM line and estimate value
    fwhmObj = plot(axObj, fwhmline(:,1), fwhmline(:,2),'--c');
    % Datapoint
    cursorMode = datacursormode(fh); 
    hdtip = cursorMode.createDatatip(fwhmObj);
    hdtip.Position = [fwhmline(1,1) fwhm 0]; 
    updateDataCursors(cursorMode);

    % Plot FWHM intersect points
    plot(axObj, xest, yest, '-om');

end


catch % catches errors due to full nan-values data
    linewidth = NaN;
    xest = NaN; yest = NaN;
end
