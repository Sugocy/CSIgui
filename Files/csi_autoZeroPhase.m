function spec = csi_autoZeroPhase(spec, poi, method, disp_on)
% Quick and dirty method to apply zeroth order phase correcton.
% The peak of interest can be described by a range or a single point. 
%
% input: spectrum, range, method, display results
%
% Method 1: maximize the real component of the spectrum (default)
% Method 2: match the real component to the spectrums maximum magnitude
%
% qhoutum2@umcutrecht.nl

% Default no visualization of data
if nargin == 2, method = 1; disp_on = 0; end
if nargin == 3, disp_on = 0; end

if method == 1

% Split complex
magn = abs(spec); pha = angle(spec); % Magnitude and phase of complex spec

% Calculate all phases at this peak
pha_step  = (2*pi/360); 
pha_range = pha(poi,:)+((1:360).*pha_step);

% Real of each phase for all poi given
real_corr = NaN(size(poi,2),size(pha_range,2));
for pii = 1:size(poi,2)
    real_corr(pii,:) = magn(poi(pii),:).*cos(pha_range(pii,:));
end

% Max phase position for each poi in range.
[maxReal_inRange, maxPhaStep_inRange] = max(real_corr,[],2);
% Take max of max; e.g. take in poi range one with highest maximum.
[~, maxReal_poi_pos] = max(maxReal_inRange);
% Calculate corrected phase for each point in the spectrum by adding 
% the found phase change which maximizes real part of the complex 
% at the points in the given range. 
pha_corrected = repmat(pha,[1 size(poi,2)]) + (pha_step.*maxPhaStep_inRange)';
spec_all = complex(magn.*cos(pha_corrected), magn.*sin(pha_corrected));



% Show if requested.
if disp_on == 1
    figure();

    subplot(3,1,1); 
    plot(real(spec)); title('Original');
    
    subplot(3,1,2); 
    plot(real_corr'); hold on;
    plot(maxPhaStep_inRange,maxReal_inRange,'or'); 
    title('Real @ all phases');
    
    subplot(3,1,3); 
    plot(real(spec_all)); title('Corrected'); hold on;
    ci = plot(real(spec),'--m','LineWidth', 1.2); 
    legend(ci,'Selected spectrum');
    
    drawnow;
end

% Set output as spectrum with phase corrected resulting in the highest 
% real value at that point in the given range.
spec = spec_all(:,maxReal_poi_pos);

elseif method == 2
    
    
    % Get temporary spectrum: phase and angle & cut proper range
    tmp = spec; tmp = tmp(poi,:);
    
    % Get phase and angle
    ph = angle(tmp); mg = abs(tmp);

    % Absolute and its max amplitude
    spec_abs = abs(tmp); target_realcomp = max(spec_abs);
    
    % Calculate all phase changes to add to original phase
    dphi = (2.*pi./360);           % Phase step
    dphi_range = (1:360).*dphi;    % Phase range

    % Add the increasing phase change to a repeated column of the phases of 
    % each point in the spectrum.
    ph_new = repmat(ph,1,360) + dphi_range;
    % Calculate the new maximum in for the newly phased real components
    rl_new = mg.*cos(ph_new); mx_new = max(rl_new,[],1);
    % Get difference between new real components maximum and the 
    % abs maximum of the original; the target
    df_new = abs(target_realcomp - mx_new);
    % Minimum distance real part to the absolute maximum.
    [~, ind] = min(df_new);
    
    % Show if requested.
    if disp_on == 1
        figure(); 
        
        subplot(2,1,1); plot(rl_new,'b'); hold on;
        plot(rl_new(:,ind),'r','LineWidth',1.2);
        title('Rephased real component');
        
        subplot(2,1,2); plot(df_new); 
        title('Difference with target absolute spectrum maximum');
    end
    
    % Zeroth order phase correction of interest
    poi = dphi_range(:,ind);
    
    % New complex spectrum
    mag = abs(spec); pha = angle(spec);
    spec = complex( mag.*cos(pha+poi),...
                    mag.*sin(pha+poi) );

end
