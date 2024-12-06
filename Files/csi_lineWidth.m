function [fwhm, fwhm_val, fwhm_pos] = csi_lineWidth(data, ppm, poi, show_plot) 
% Returns the FWHM linewidth, value and position by interpolating the data
% before calculating the FWHM at the given peak of interest.
% 
% Input:
% Data      - Is the spectrum, complex or real.
% Xaxis     - The ppm range of the given spectrum.
% Poi       - The peak of interest given as a numerical index range
%
% Output:
% FWHM      - Linewidth at the half maximum of given range (1x1)
% Value     - Half maximum value (1x2)
% Position  - Position of the half value (1x2)
%
% Q. van Houtum

% Default no graphs displayed
if nargin == 3, show_plot = 0; end
    
interp_factor = 10^-3;              % Interpolation factor
doi = data(poi,:); xoi = ppm(poi)'; % Data at POI

% Maximum value of peak % -------------------------------------------- %
% Keeping in mind negative values with the max-function.
mi = min(real(doi)); 
if mi < 0
    doi_real = real(doi) + abs(mi); 
    doi_imag = imag(doi);
    doi = complex(doi_real,doi_imag);
end

% Interpolate DATA and X-AXIS
xoii = interp1(xoi,1:interp_factor:size(xoi,1),'Pchip'); % New xaxis
doii = interp1(xoi,doi,xoii,'linear');                   % New data array


% Find maximum position and half-maximum
[~, max_ind] = max(real(doii)); 
max_half = abs((doii(max_ind))/2);

% Continue with real part.
doii = real(doii);

% LEFT SIDE (to pos ppm) % ------------------------------------------ %
% Find position closest to half-max / left side

close_to_HM = abs( (doii(max_ind:end)) - (max_half) );
xoii_left   = xoii(max_ind:end);
bool = islocalmin(close_to_HM); ind = find(bool == 1);
if isempty(ind)
    fwhm_indl = 1;      no_local_minimum = 1;
else
    fwhm_indl = ind(1); no_local_minimum = 0;
end

if max_ind == size(doii,2)
    fwhm_val(2) = doii(max_ind); fwhm_pos(2) = xoii(max_ind);
else
    fwhm_val(2) = doii(fwhm_indl + (max_ind -1) ); 
    fwhm_pos(2) = xoii(fwhm_indl + (max_ind -1) );
end
% RIGHT SIDE (negative) % ------------------------------------------ %
% Find position closest to half-max / right side

% Calculate difference between HM and each point.
% Flip, so we begin at the maximum and ride the peak down 
close_to_HM = flip(abs( doii(1:max_ind) - max_half),2);
xoii_right  = xoii(1:max_ind);

% Find the lowest  e.g. almost zero! (Almost no difference)
bool = islocalmin(close_to_HM); ind = find(bool == 1);
% Take the first minimum e.g. first index to be half maximum from maximum
if isempty(ind)
    fwhm_indr = 1;      no_local_minimum = 1;
else
    fwhm_indr = ind(1); no_local_minimum = 0;
end

% Correct for the flip we applied
fwhm_indr = abs( size(doii(1:max_ind),2) - fwhm_indr ) + 1;

if max_ind == size(doii,2)
    fwhm_val(1) = doii(max_ind); 
    fwhm_pos(1) = xoii(max_ind);
else
  
    fwhm_val(1) = doii(fwhm_indr); 
    fwhm_pos(1) = xoii(fwhm_indr);
end

% Full Width at Half Maximum
if mi < 0, fwhm_val = fwhm_val-abs(mi); end
fwhm = abs(fwhm_pos(1) - fwhm_pos(2));

% NB: Corrections for inaccurate measurements
if no_local_minimum, fwhm = NaN;
end

% ----------------------------------------------------------------------- %

if show_plot 
    figure(); 
    
    % Ylimit value
    ylim_val = [-2.*std(real(data)) max(real(data))+(2.*std(real(data)))];
    
    % Show original
    ax1 = subplot(3,1,1); 
    plot(ax1, ppm, real(data)); ylim(ylim_val);     
    if no_local_minimum, title('Original - LineWidth not accurate');
    else,                title('Original');
    end
    % Plot fwhm in original
    hold(ax1,'on'); plot(ax1, fwhm_pos, fwhm_val,'or');
    ax1.XDir = 'Reverse';
    
    % Plot peak of interest
    ax2 = subplot(3,1,2); 
    plot(ax2, xoi, real(data(poi,:)),'s-'); 
    title('Peak of Interest');
    hold(ax2,'on'); plot(ax2, fwhm_pos, fwhm_val,'-or', 'LineWidth', 1.5);
    ax2.XDir = 'Reverse';
    
    % Show interpolated
    ax3 = subplot(3,1,3); 
    plot(ax3, xoii, interp1(xoi,real(data(poi,:)),xoii,'linear'),'-'); 
    title('Peak of Interest (Interpolated)');
    % Plot fwhm in interpolated
    hold(ax3,'on'); plot(ax3, fwhm_pos, fwhm_val,'-or', 'LineWidth', 1.5);
    ax3.XDir = 'Reverse';
    
    % Show LW in last plot
    lw_str = sprintf('FWHM: %2.3f',fwhm);
    text(ax3, ax3.XLim(1)+(diff(ax3.XLim).*0.75),...
              ax3.YLim(1)+(diff(ax3.YLim).*0.75),lw_str);    
end