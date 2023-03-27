% Find maximum likelihood combined spectrum from a receive array using:
%
% (i) the original WSVD algorithm;
%
% USAGE:
% [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, noiseMask, varargin)
%
% rawSpectra: must be an N x M matrix (N = frequency points, M = receive array elements)
% noiseMask: must be an N x 1 logical matrix (N = frequency points).
%              true --> this point is baseline noise.
%              false --> this point contains signal.
%            or it can be [] and a full noise covariance matrix supplied
%            via the name/value options in varargin.
%
% 
% (ii) the extended WSVD+Apod algorithm;
%
% USAGE:
% [wsvdApodCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApod(rawSpectra, noiseMask, apod, timeAxis, varargin)
%
% rawSpectra: must be an N x M matrix (N = frequency points, M = receive array elements)
% noiseMask: must be an N x 1 logical matrix (N = frequency points).
%              true --> this point is baseline noise.
%              false --> this point contains signal.
%            or it can be [] and a full noise covariance matrix supplied
%            via the name/value options in varargin.
% apod: Apodization factor (in Hz).
% timeAxis: Times corresponding to each FID point (in s). (Must be computed
%           from scan parameters.)
%
%
% (iii) the extended WSVD+Apod+Blur algorithm.
%
% USAGE:
% [wsvdApodBlurCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApodBlur(allData, noiseMask, apod, timeAxis, blurRadius, blurExtent, varargin)
%
% allData: must be CHEM_SHIFT_PTS x COIL x NCOL x NROW
%          i.e. a whole CSI slice must be supplied instead of a single
%          voxel at a time
% noiseMask: must be an N x 1 logical matrix (N = frequency points).
%              true --> this point is baseline noise.
%              false --> this point contains signal.
%            or it can be [] and a full noise covariance matrix supplied
%            via the name/value options in varargin.
% apod: Apodization factor (in Hz).
% timeAxis: Times corresponding to each FID point (in s). (Must be computed
%           from scan parameters.)
% blurRadius : Blur radius / voxels.
% blurExtent : 0 -> 1x1, 1 -> 3x3, 2 -> 5x5, etc.
%
%
% REFERENCE:
% 1. Rodgers and Robson, Magnetic Resonance in Medicine, 2010.
%    (http://dx.doi.org/10.1002/mrm.22230).
% 2. Rodgers and Robson, Magnetic Resonance in Medicine, YYYY.
%    (http://dx.doi.org/INSERT FINAL DOI HERE).
% 
% Copyright Chris Rodgers, University of Oxford, 2014.

function [wsvdCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvd(rawSpectra, noiseMask, varargin)

options = processVarargin(varargin{:});

if ~isfield(options,'debug')
    options.debug = 0;
end

if ~isfield(options,'phaseRefChannel')
    options.phaseRefChannel = 1;
end

if numel(size(rawSpectra))>2
    error('rawSpectra must be an N x M matrix (N = frequency points, M = receive array elements)')
end

nCoils = size(rawSpectra,2);

%% Measure noise statistics from data (unless instructed otherwise)
if isfield(options,'noiseCov') && ~ischar(options.noiseCov)
  % Noise covariance matrix has been supplied
  noiseCov=options.noiseCov;
elseif isfield(options,'noiseCov') && strcmp(options.noiseCov,'disable')
    % Disable noise prewhitening
    noiseCov=eye(nCoils)*0.5;
elseif isfield(options,'noiseCov') && strcmp(options.noiseCov,'diag')
    % Use only the noise variances (not off diagonal elements)...
    noiseCov=diag(diag(cov(rawSpectra(noiseMask,:))));
elseif ~isfield(options,'noiseCov') || strcmp(options.noiseCov,'normal')
    % Estimate from data (DEFAULT)
    noiseCov=cov(rawSpectra(noiseMask,:));
else
    error('Option "noiseCov" has an unknown value.')
end

[noiseVec, noiseVal] = eig(noiseCov);
% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V 
% whose columns are the corresponding right eigenvectors, so that A*V = V*D.

scaleMatrixFft = noiseVec*diag(sqrt(0.5)./sqrt(diag(noiseVal)));
invScaleMatrixFft = inv(scaleMatrixFft);

scaledSpectra=rawSpectra*scaleMatrixFft;

if options.debug
    disp('Noise covariance matrix eigenvalues:')
    disp(noiseVal)
end

%% Compute optimal reconstruction using SVD
[u,s,v]=svd(scaledSpectra,'econ');
% SVD quality indicator
wsvdQuality(1) = ((s(1,1)/norm(diag(s)))*sqrt(nCoils)-1)/(sqrt(nCoils)-1);

% Coil amplitudes
wsvdCoilAmplitudes=v(:,1)'*invScaleMatrixFft;
% There's an arbitrary scaling here such that the first coil weight is
% real and positive

svdRescale = norm(wsvdCoilAmplitudes)*normalise(wsvdCoilAmplitudes(options.phaseRefChannel));

wsvdCoilAmplitudes=wsvdCoilAmplitudes / svdRescale;

if options.debug
    fprintf('DEBUG: SVD quality indicator = %g\n',wsvdQuality);
    fprintf('SVD coil amplitudes:\n');
    disp(wsvdCoilAmplitudes)
end

wsvdCombination = u(:,1)*s(1,1)*svdRescale;

wsvdWeights = 0.5*inv(noiseCov) * wsvdCoilAmplitudes' * conj(svdRescale) * svdRescale;

end


function [wsvdApodCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApod(rawSpectra, noiseMask, apod, timeAxis, varargin)
% Two-step variant of WSVD. (1) Compute weights from apodized data.
% (2) Apply those weights to the original single-element spectra.
%
% Note that this should be IDENTICAL to "WSVD" combination if
% apod = 0.

if ~exist('apod','var')
    error('You must specify the apodization factor.')
end

if ~exist('timeAxis','var')
    error('You must supply the timeAxis.')
end

if numel(size(rawSpectra))>2
    error('rawSpectra must be an N x M matrix (N = frequency points, M = receive array elements)')
end

% Step 1: Compute weights from apodized data.
[~, ...
    wsvdQuality, ...
    wsvdCoilAmplitudes, ...
    wsvdWeights] ...
    = wsvd(specApodize(timeAxis, rawSpectra, options.apod), [], 'noiseCov', noiseCov);

% Step 2: Apply weights to the raw single-element spectra.
wsvdApodCombination = rawSpectra * wsvdCoilWeights;
end


function [wsvdApodBlurCombination, wsvdQuality, wsvdCoilAmplitudes, wsvdWeights] = wsvdApodBlur(allData, noiseMask, apod, timeAxis, blurRadius, blurExtent, varargin)
% Two-step variant of WSVD. (1) Compute weights from apodized and
% blurred data. (2) Apply those weights to the original
% single-element spectra.
%
% apod : Apodization constant / Hz. (30 Hz).
% blurRadius : Blur radius / voxels. (1)
% blurExtent : 0 -> 1x1, 1 -> 3x3, 2 -> 5x5, etc. (1)
%
% This time a whole CSI slice must be supplied instead of a single voxel at
% a time. So rawSpectra must be: CHEM_SHIFT_PTS x COIL x NCOL x NROW

if ~exist('apod','var')
    error('You must specify the apodization factor.')
end
        
if ~exist('timeAxis','var')
    error('You must supply the timeAxis.')
end

if ~exist('blurRadius','var')
    error('You must specify the blur radius (in units of 1 voxel dimension).')
end

if ~exist('blurExtent','var')
    error('You must specify the blur matrix extent (in voxels).')
end

if numel(size(allData))~=4
    error('allData must be an #frequency points x #receive array elements x #columns (in CSI) x #rows (in CSI) matrix')
end

% Set up the blur matrix
blurExtentVec = -blurExtent:blurExtent;
[blur_x,blur_y] = ndgrid(blurExtentVec,blurExtentVec);
blur_r = sqrt(blur_x.^2+blur_y.^2);
blurMatrix = exp(-(blur_r/options.blurRadius).^2);

% Set up output matrix
numSamples = size(allData,1);
numCoils = size(allData,2);
numColumns = size(allData,3);
numRows = size(allData,4);
numSlices = 1;

outputSize = [size(allData,3) size(allData,4)];
wsvdApodBlurCombination = zeros([ size(allData,1) outputSize ]);

voxNums = 1:prod(outputSize);
for voxDx = 1:prod(outputSize)
    [vox_crs{voxDx}(1), vox_crs{voxDx}(2), vox_crs{voxDx}(3)] = ind2sub([outputSize 1],voxNums(voxDx));
end

% Loop over each voxel in the CSI slice
for voxDx=1:prod(outputSize)
    % Step 1: Apodize and concatenate ("blur") the raw data.
    
    % Find the voxel numbers needed for this blur matrix
    [blurVox_C, blurVox_R, blurVox_S] = ndgrid(vox_crs{voxDx}(1) + blurExtentVec, vox_crs{voxDx}(2) + blurExtentVec, vox_crs{voxDx}(3));
    
    % Ignore blur that would exceed data dimensions
    blurVox_mask = (blurVox_C < 1 | blurVox_C > numColumns | ...
        blurVox_R < 1 | blurVox_R > numRows | ...
        blurVox_S < 1 | blurVox_S > numSlices);
    
    blurVox_C(blurVox_mask) = NaN;
    blurVox_R(blurVox_mask) = NaN;
    blurVox_S(blurVox_mask) = NaN;
    
    % Remember allData has dimensions:
    % allData = zeros(numSamples,numCoils,numColumns,numRows);
    
    data_blurred = zeros(numSamples,...
        numel(blurMatrix),...
        numCoils);
    
    for blurDx=1:numel(blurVox_C)
        if isnan(blurVox_C(blurDx)), continue, end
        
        for coilDx=1:numCoils
            data_blurred(:,blurDx,coilDx) ...
                = blurMatrix(blurDx) * ...
                specApodize(timeAxis,...
                allData(:, coilDx, blurVox_C(blurDx), blurVox_R(blurDx)),...
                apod);
        end
    end
    
    % Step 2: Compute weights from apodized and blurred (concatenated) data.
    [~, ...
        wsvdQuality(:,voxDx), ...
        wsvdCoilAmplitudes(:,voxDx), ...
        wsvdWeights(:,voxDx)] ...
        = wsvd(...
        reshape(data_blurred,...
        [numSamples*numel(blurMatrix),...
        numCoils]), ...
        noiseMask, varargin{:}); %#ok<*AGROW>
    
    % Step 3: Apply weights to the raw single-element spectra.
    wsvdApodBlurCombination(:,voxDx) = allData(:,:,voxDx) * wsvdWeights(:,voxDx);
end

end

function [outSpec] = specApodize(timevals,spec,amount)
% Apodize spectra.
%
% [outSpec] = specApodize(timevals,spec,amount)

% Check dimensions

% Prepare the window function
outSpec = bsxfun(@(oneSpec,windowFunc) specFft(specInvFft(oneSpec).*windowFunc),spec,exp(-amount*timevals));
end

function [spec] = specFft(fid,dim)
% FFT for spectroscopy. Converts fid --> spectrum.
%
% Accounts correctly for the 0.5x weighting for the t=0 FID point arising
% from the one-sided definition of the Fourier Transform relations in
% spectroscopy.
%
% Optional parameter: dim is the dimension to be FFT'd. Default is 1.
%
% EXAMPLE:
% zz=randn(2048,1)+1i*randn(2048,1);maxdiff(zz,specFft(specInvFft(zz)))
% OR
% zz=randn(2048,128)+1i*randn(2048,128);maxdiff(zz,specFft(specInvFft(zz)));maxdiff(zz,specFft(specInvFft(zz,2),2))

if nargin<2
    dim = 1;
end

perm = [dim 1:(dim-1) (dim+1):numel(size(fid))];

% Re-use variable name "spec" to economise on RAM.
spec = permute(fid,perm);

% t=0 point is treated differently by convention for a FID
spec(1,:) = spec(1,:) * 0.5;

spec = ipermute(fftshift(fft(spec,[],1),1)/size(spec,1),perm);
end

function [fid] = specInvFft(spec,dim)
% Inverse FFT for spectroscopy. Converts spectrum --> fid.
%
% Accounts correctly for the 0.5x weighting for the t=0 FID point arising
% from the one-sided definition of the Fourier Transform relations in
% spectroscopy.
%
% Optional parameter: dim is the dimension to be FFT'd. Default is 1.
% 
% EXAMPLE:
% zz=randn(2048,1)+1i*randn(2048,1);maxdiff(zz,specFft(specInvFft(zz)))
% OR
% zz=randn(2048,128)+1i*randn(2048,128);maxdiff(zz,specFft(specInvFft(zz)));maxdiff(zz,specFft(specInvFft(zz,2),2))

if nargin<2
    dim = 1;
end

fid = ifft(fftshift(spec,dim),[],dim)*size(spec,dim);

perm = [dim 1:(dim-1) (dim+1):numel(size(fid))];

% Re-use variable name "fid" to economise on RAM.
fid = permute(fid,perm);

% t=0 point is treated differently by convention for a FID
fid(1,:) = 2*fid(1,:);

fid = ipermute(fid,perm);

end

% Allow varargin to hold either a struct or a list of field/value pairs
% or a mixture of structs and field/value pairs.
%
% If parameters are specified more than once, the last value is used.
%
% [options] = processVarargin(varargin)
%
% Example 1:
%
% To use this in a function, add code as follows:
%
% function [ ... ] = myfunction(param1, param2, ... , varargin)
% ...
% % Read options from varargin
% options = processVarargin(varargin{:});
%
%
% Example 2:
%
% inputStruct = struct('abc',1);
% inputStruct.def = {23};
% options = processVarargin(inputStruct);
% options2 = processVarargin('abc',1,'def',{23});
% isequal(options,options2)
%
% After running this code, inputStruct, options and options2 will all be
% the same.

% Copyright Chris Rodgers, University of Oxford, 2008-13.
% $Id: processVarargin.m 6515 2013-05-14 10:57:39Z crodgers $

function [opt] = processVarargin(varargin)

full_f = cell(0,1);
full_v = cell(0,1);

idx = 1;
while idx <= numel(varargin)
    if isstruct(varargin{idx}) && numel(varargin{idx}) == 1
        f = fieldnames(varargin{idx});
        v = struct2cell(varargin{idx});
        idx = idx+1;
    elseif ischar(varargin{idx}) && numel(varargin) > idx
        f = varargin(idx);
        v = varargin(idx+1);
        idx = idx+2;
    else
        error('Bad input format! Expected a struct or field/value pair.')
    end

    % Now append to the main list
    full_f(end+1:end+numel(f),1) = f;
    full_v(end+1:end+numel(f),1) = v;
end

% Sort fields and catch duplicate field names
[full_f, fnDx] = myUnique(full_f);
full_v = full_v(fnDx);

opt = cell2struct(full_v,full_f);
end

function [out,dx] = myUnique(in)
% unique is changing in R2012a so roll our own

[tmpSort, tmpSortDx] = sort(in);
mask = true(size(tmpSort));
for idx=numel(tmpSort)-1:-1:1
    if isequal(tmpSort(idx+1),tmpSort(idx))
        mask(idx) = false;
    end
end

out = tmpSort(mask);
dx = tmpSortDx(mask);

end
