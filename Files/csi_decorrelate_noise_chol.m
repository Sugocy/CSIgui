function [spec, nCov_Chol, scaleMatrix] = csi_decorrelate_noise_chol(spec, chan_ind, nCov)
% Decorrelate signal according to Cholesky decomposition for decorrelation
% CG-SENSE, Pruessman et al MRM 2001.
%
% Adapted from Ayhan Gursan, 2021, UMC Utrecht. a.gursan@umcutrecht.nl.
%
% Input: spectra,  index for channels, noise_covariance matrix,
%   If nCov is not an input argument, it will be calculated using
%   noise in the data.
%
% Quincy van Houtum, PhD; 10/2023
% quincyvanhoutum@gmail.com

if nargin < 2, nCov = csi_noisecov_usingdata(spec, chan_ind); end
dim = size(spec);


% Generate virtual coils using Chlesky-Decomp, resulting in unit-noise.
nCov_Chol = cellfun(@chol, nCov, 'UniformOutput', 0);

% Reshape to {nS x nChan} x nVoxels
[spec, permv, szr] = csi_combine_reshape(spec, chan_ind);
nCov_Chol = reshape(nCov_Chol, [numel(nCov_Chol) 1]);

% Decompose signal
[spec, scaleMatrix] = cellfun(@decorrelate_noise_equation, ...
                              spec, nCov_Chol, 'UniformOutput', false);

% Reshape back to matrix with same dimensions
if numel(szr) > 2 && (numel(dim)-2 > 2)
    spec = reshape(spec, szr(3:end)); vec_sz = size(szr(3:end),2);
    spec = permute(spec,[vec_sz + [1 2] 1:vec_sz]);
    spec = cell2mat(spec);
else
    spec = cat(numel(size(spec))+1, spec{:});
end

% Permute to correct order as before
restore_permv = permv;
for kk = 1:size(permv,2), restore_permv(kk) = find(kk == permv); end
spec = permute(spec, restore_permv); 

% Reshape nCovariance Cholensky decomposition matrix.
if numel(szr) >= 3 && (numel(dim)-2 > 2)
    nCov_Chol = reshape(nCov_Chol,szr(3:end));
end

% --- Decorrelate noise equation
function [dSignal, scaleMatrix] = decorrelate_noise_equation(signal, nCov_Chol)
% Decorrelated-signal = pinv(Cholesky-decomposed Cov-Matrix) * Signal;
% i.e., Pseudo-inverse of cholesky decomposed nCov matrix-times the 
% signal.
scaleMatrix = pinv(nCov_Chol'); dSignal = ( scaleMatrix * signal')';

% --- Noise Covariance matrix using data
function noise_cov = csi_noisecov_usingdata(spec, chan_ind, noise_mask)
% Returns a noise covariance matrix with respect to each channel for 
% every voxel using the data itself.
%
% output noise_cov = {nChan x nChan} x Spatial Dimensions ...
%
% !! If no separate noise-data is available, use this function
% !! If no noise_mask given - it will be calculated as 2x(1/6)% of the
%                              start and end of the spectrum.

dim = size(spec);
if nargin < 3
    nS = dim(1); half_nm_size = round(nS./6);
    noise_mask = [1:half_nm_size (nS - half_nm_size + 1):nS];
end

% Index for cutting data
cut_ind = arrayfun(@(x) 1:x, dim,'UniformOutput',0);
cut_ind{1} = noise_mask;
noise_data = spec(cut_ind{:});

% Reshape channel index: {nS x nChan} x spatial dimensions...
if isempty(chan_ind), chan_ind = numel(dim)+1; end

% Cell-ify and reshape data to {nS x nChan} x nVox
[noise_data, ~, szr] = csi_combine_reshape(noise_data, chan_ind);

% Noise covariance Matrix
noise_cov = cellfun(@cov, noise_data,'UniformOutput',0);
noise_cov = cellfun(@diag, noise_cov,'UniformOutput',0);
noise_cov = cellfun(@diag, noise_cov,'UniformOutput',0);

% Reshape to {noise_cov} x Spatial Dimensions
if numel(szr) >= 3 && (numel(dim)-2 > 2)
    noise_cov = reshape(noise_cov, szr(3:end)); 
end