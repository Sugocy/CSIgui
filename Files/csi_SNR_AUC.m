function [SNR, peak_est] = csi_SNR_AUC(data, poi, noise, bv_mask, unit, strict, doDisp)
% Calculate SNR using the area under the curve of the peak within the given
% peak-of-interest (poi) range.
%
% If noise is a single value, it is a noise-mask size, otherwise it is
% assumed noise-array with equal size to data was given.
%
% Input:
%   data        data-array with spectra with index 1 as sample index.
%   poi         unitless index of peak of interest range in data.
%   noise       noise-mask size or noise per voxel as an array.
%   bv_mask     noise-mask for baseline calculation required for AUC.
%   unit        signal unit magnitude (0), real (1, default) or
%               imaginary (2).
%   strict      Allow the peak-range to increase (0, default) 
%               if no intersection with the base-value within 
%               poi is found. See csi_AUC for more info
%   doDisp      display AUC per voxel - can be graphically intensive!
%
% Q. van Houtum, PhD; version 1.0 01/2025
% quincyvanhoutum@gmail.com

if nargin < 7, doDisp = 0; end
if nargin < 6, strict = 0; end
if nargin < 5, unit = 1; end 

% --- STRUCTURE
% Exclude first index dimensions
% Create vector of ones equal in size to each remaining dimension.
sz = size(data); 

% Data coverted to unit of interest: magn = 0; real = 1; imaginary = 2;
if     unit == 1, data = real(data); 
elseif unit == 2, data = imag(data);  
else,             data = abs(data); 
end

% Creat cell-array of each spectrum
cell_layout = arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),...
    'UniformOutput',0);
specs = mat2cell(data, sz(1), cell_layout{:});


% --- NOISE

% Noise = noise component(s) itself or
% Noise = size of noise mask (go into if-statement)
if sum((noise(:) - floor(noise(:))) == 0) > 0 && isscalar(noise) 

    % Calculate voxel-noise from data-area
    mask_size_double = [round(noise./2) noise-round(noise./2)];
    mask = [1:mask_size_double(1) (sz(1) - mask_size_double(2) + 1):sz(1)];

    % Precalculate the noise - already in the correct unit!
    noise = cellfun(@voxel_noise,...
        specs, repmat({mask}, size(specs)),'UniformOutput',0);
    noise = cell2mat(noise);

end

% --- AUC
[auc, peak_est] = cellfun(@csi_AUC, specs, ...
         repmat({poi}, size(specs)), ...
         repmat({bv_mask}, size(specs)),...
         repmat({strict}, size(specs)),...
         repmat({doDisp}, size(specs)),'UniformOutput',0);    
auc = cellfun(@(x) x', auc, 'UniformOutput', 0);
auc = cell2mat(auc);

% --- SNR
ind = arrayfun(@(x) 1:x, sz(2:end),'UniformOutput', 0);
SNR = auc(1, ind{:}) ./ noise;



function noise = voxel_noise(spectrum, mask)
% Calculate noise
noise = abs( std( spectrum(mask,:),[],1 ) );