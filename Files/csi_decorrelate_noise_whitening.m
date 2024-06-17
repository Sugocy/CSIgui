function [data, scaleMatrix] = csi_decorrelate_noise_whitening(data, chan_ind, nCov)
% Apply whitening/decorrelation to the data array using the noise 
% covariance matrix. The whitening transformation is equal to the WSVD
% method. The scaleMatrix is required for WSVD and therefor returned.
% However, the scaleMatrix can be recalculated using the same
% noiseCovariance matrix - but without applying it to the spectra.
%
% Input
%   data      = array of spectra [M x N x ...]
%   chan_ind  = channel index in data-array.
%   nCov      = noise covariance matrix for every voxel of the data-array
%               thus with size [M x N x ...] equal to data.
%
% If no noise-covariance matrix is given, it will be calculated via a 
% noise-mask of 2x (1/6) of the #samples from the tails of the spectrum. 
%
% Quincy van Houtum, PhD; 05/2024
% quincyvanhoutum@gmail.com


% Calculate noise-cov using data if necessary
if nargin < 3, nCov = csi_noisecov_usingdata(data, chan_ind); end

% Reshape to {nS x nChan} x nVoxels % -------------------------- %
[data, permv, szr] = csi_combine_reshape(data, chan_ind);


% Apply whitening % -------------------------------------------- %

% Eigenvalue and vector of the noise-cov matrix
% e.g. nCov * nVec = nVec * nVal
[nVec, nVal]= cellfun(@eig, nCov, 'UniformOutput', 0);

% Calculate the scale matrix
scaleMatrix = cellfun(@csi_combine_ZCA_scaleMatrix, nVec, nVal, ...
    'UniformOutput',0);

% Apply scaling
data = cellfun(@(x,y) x * y, data, reshape(scaleMatrix, size(data)),...
    'UniformOutput',0);
        

% Reshape to orignal index order % ----------------------------- %

% Calculate permute vector to set non-cell indexes at correct
% dimension.
perm_index = [numel(size(data)) + [1 2] 1:numel(size(data))];
data = permute(data, perm_index);

data = cell2mat(data); % Convert to matrix
data = reshape(data, szr); % Reshape to original size

% Create undo-permute-permute-vector
nindex = numel(permv); restore_permv = NaN(1,nindex);
for kk = 1:nindex, restore_permv(kk) = find(kk == permv); end

% Permute channel index back to orignal position
data = permute(data, restore_permv);

% --- Noise Covariance matrix using data
function nCov = csi_noisecov_usingdata(spec, chan_ind, noise_mask)
% Returns a noise covariance matrix with respect to each channel for 
% every voxel using the data itself.
%
% output nCov = {nChan x nChan} x Spatial Dimensions ...
%
% !! If no separate noise-data is available, use this function
% !! If no noise_mask given - it will be calculated as 2x(1/6)% of the
%                              start and end of the spectrum.

dim = size(spec);
if nargin < 3
    nS = dim(1); half_nm_size = round(nS./6);
    noise_mask = [1:half_nm_size nS - half_nm_size];
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
nCov = cellfun(@cov, noise_data,'UniformOutput',0);
nCov = cellfun(@diag, nCov,'UniformOutput',0);
nCov = cellfun(@diag, nCov,'UniformOutput',0);

% Reshape to {nCov} x Spatial Dimensions
if numel(szr) >= 3 && (numel(dim)-2 > 2)
    nCov = reshape(nCov, szr(3:end)); 
end

function scaleMatrix = csi_combine_ZCA_scaleMatrix(nVec,nVal)
% Calculate the scale matrix for ZCA whitening
scaleMatrix = nVec * diag(sqrt(0.5) ./ sqrt(diag(nVal)));
