function SNR = csi_SNR_AUC(data, poi, noise, bv_mask, doDisp)
% Calculate SNR using the area under the curve of the peak within the given
% peak-of-interest (poi) range.
%
% If noise is a single value, it is a noise-mask size, otherwise it is
% assumed noise-array with equal size to data was given.
%
% Input:
%   data        data-array with spectra with index 1 as sample index.
%   poi         unitless index of peak of interest range in data.
%   noise       noise-mask size or noise per voxel array.
%   bv_mask     noise-mask for baseline calculation for the AUC.
%   doDisp      display AUC per voxel - can be intensive!
%
% Q. van Houtum, PhD; version 1.0 10/2024
% quincyvanhoutum@gmail.com

if nargin < 5, doDisp = 0; end

% --- STRUCTURE
% Exclude first index dimensions
% Create vector of ones equal in size to each remaining dimension.
sz = size(data);
cell_layout = arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),...
    'UniformOutput',0);
% Creat cell of data +  Time gone from 1st dim. Inside cells now.
specs = mat2cell(data, sz(1), cell_layout{:});


% --- NOISE
if numel(noise) == 1 % Calculate voxel-noise from data-area
    mask_size_double = [round(noise./2) noise-round(noise./2)];
    mask = [1:mask_size_double(1) (sz(1) - mask_size_double(2) + 1):sz(1)];

    % Precalculate the noise
    noise = cellfun(@voxel_noise,...
        specs, repmat({mask}, size(specs)),'UniformOutput',0);
    noise = cell2mat(noise);
end

% --- AUC
auc = cellfun(@csi_AUC, specs, ...
         repmat({poi}, size(specs)), ...
         repmat({bv_mask}, size(specs)),...
         repmat({doDisp}, size(specs)),'UniformOutput',0);    
auc = cellfun(@(x) x', auc, 'UniformOutput', 0);
auc = cell2mat(auc);


% --- SNR
ind = arrayfun(@(x) 1:x, sz(2:end),'UniformOutput', 0);
SNR = auc(1,ind{:}) ./ noise;



function noise = voxel_noise(spectrum, mask)
% Calculate noise
noise = abs( std( spectrum(mask,:),[],1 ) );