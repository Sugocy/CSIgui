function [spec, Q, W, A] = WSVD2(spec, noiseCov)
% Whitened singular value decomposition to calculate weights for combining
% multiple coil/channels in MR-data.
%
% This algorithm is based on the WSVD algorithm described in:
% Rodgers et al., http://dx.doi.org/10.1002/mrm.22230
%
% Expected input data format: 
%           spec     = [ nSamples  x nChannels ] x 1
%           noiseCov = [ nChannels x nChannels ] x 1
%
% Output:
%           spec = combined spectra
%           Q    = quality of svd
%           W    = weights per channel
%           A    = amplitudes per channel
%
% Q. van Houtum, PhD.
% quincyvanhoutum@gmail.com

% calculate noise-covariance matrix using data
if nargin < 2, noiseCov = csi_noisecov_usingdata(spec,2); end

% Preperation % -------------------------------------------------------- %
nCoils = size(spec,2);
refChannel = 1; % Used for scaling, could be any of the channels.

% Whitening % ---------------------------------------------------------- %

% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V 
% whose columns are the corresponding right eigenvectors, so that A*V = V*D.
[noiseVec, noiseVal] = eig(noiseCov);

% Scale spectra using eigen-values
scaleMatrix = noiseVec * diag(sqrt(0.5)./sqrt(diag(noiseVal)));
spec = spec*scaleMatrix; 


% SVD % ---------------------------------------------------------------- %
[u,s,v] = svd(spec,'econ');

% SVD Quality
Q(1) = ((s(1,1)/norm(diag(s)))*sqrt(nCoils)-1)/(sqrt(nCoils)-1);

% Coil amplitudes
A = v(:,1)'/scaleMatrix; % (updated code)


% Arbitrary scaling such that the first coil weight is real and positive
svdRescale = norm(A) * normalise(A(refChannel));

% Correct for scaling
A = A / svdRescale;

% Calculate applied weights.
W = 0.5 * (noiseCov \ A') * conj(svdRescale) * svdRescale; % (updated code)

% Calulate combination % ----------------------------------------------- %
spec = u(:,1) * s(1,1) * svdRescale;



function [val] = normalise(val)
% Normalise a vector.
val = val / norm(val);


function noise_cov = csi_noisecov_usingdata(spec, chan_ind, noise_mask)
% Returns a noise covariance matrix with respect to each channel for 
% every voxel using the data itself.
%
% output noise_cov = {nChan x nChan} x Spatial Dimensions ...
%
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