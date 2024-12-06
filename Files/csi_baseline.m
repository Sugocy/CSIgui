function [data, bline, blevel, best] = csi_baseline(data, method, preview)
% Apply baseline correction using user input for base line curve detection.
%
% Input: data       - Spectrum
%        method     - (1) Use clicked locations (0) use data locations
%        preview    - (1) View outcome
%
% Output: [data, bline,     baseline correction
%                blevel,    baseline default level 
%                best       baseline correction estimate from user
%
% Contact: qhoutum2@umcutrecht.nl

% Process # of input
if     nargin == 1, method = 1; preview = 0;
elseif nargin == 2, preview = 0;
end

% Process data

% Base level from begin and end of spectrum
blevel = mean(real(data([1:10 end-10:end]))); 

x4fit = 1:size(data,1);

%% Userinput

% Plot original
figh = figure(); plot(x4fit, real(data)); title('Spectrum'); hold on;
% Plot: BASELEVEL
plot(figh.Children, x4fit, repmat(blevel,size(data,1)),'--k');
% figh.Children.XDir = 'Reverse';
% Get user input
title(figh.Children,'Click 6 points to define baseline offset');

pos = ginputQ(5,figh); if isnan(pos),return;end

x = pos(:,1); y = pos(:,2); x = round(x); 

title('Spectrum'); 

if method % Method 1: Use clicked locations
    % Base line estimate user.
    best = cat(2, [1; x(1)-10; x; x(end)+10; size(data,1)],...
                  [blevel; blevel; y; blevel; blevel]);
else
    % Base line estimate user.
    tmp = real(data(x,:));
    best = cat(2, [1; x(1)-10; x; x(end)+10; size(data,1)],...
                  [blevel; blevel; tmp; blevel; blevel]);
end

% Create baseline 
bline = interp1(best(:,1), best(:,2), x4fit, 'pchip');
bline = -1.* (bline');

% Create corrected data
tmp = real(data) + bline ;
data = complex(tmp, imag(data));

if preview 
    % Plot: ORIGINAL
    plot(figh.Children, real(data),'-b'); 
    title('Baseline correction'); hold on;
    % Plot: BASELINE
    plot(figh.Children, x4fit, -1.*bline,'-g');
    % Plot: CORRECTED
    plot(figh.Children, x4fit, tmp, '-r')
    % Legend
    legend('Data','Baselevel','Baseline','Corrected');
else
    close(figh);
end
