function SNR = csi_SNR(spectrum, mask_size, method, range)
% Description:                 Calculate SNR of multidimensional MRSI data
% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-09
% --------------------------------------------------------------------
% Calculate SNR for given spectra using the formula shown below input
% description.
%
% Input:
%       spectrum;   Array of all spectra of interest.(dHz x N x M x etc.)
%       mask_size;  Value for the mask size used to calculate noise. If no 
%                   mask_size is given a defeault mask size of 50 is used.
%       method;     Either using the absolute (0) or real(1) part of the
%                   spectrum
%       range;      Range of spectrum to calculate SNR from. Noise is still
%                   calculated using the noise mask. range = [start end]
% 
% Formula: SNR = mean( max ( abs ( spec ))) / abs(std(spec))) --> (1)
% Formula: SNR = mean( max ( real( spec ))) / abs(std(spec))) --> (0)
%
% Contact: qhoutum2@umcutrecht.nl

% Process #input arguments
if     nargin == 1, mask_size = 50; method = 0; 
elseif nargin == 2, method = 0;
end
if nargin <= 3, range = [1 size(spectrum,1)]; end


% Noise mask
sz = size(spectrum); mask = (sz(1)-mask_size):sz(1)-1;

% Exclude first index dimensions
% Create vector of ones equal in size to each remaining dimension.
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

% Creat cell of data +  Time gone from 1st dim. Inside cells now.
specs = (mat2cell(spectrum, sz(1), cell_layout{:}));


% Precalculate the noise
noise = cellfun(@spectrum_noise,...
    specs, repmat({mask}, size(specs)),'UniformOutput',0);

% SNR of each spectrum.
if method == 0 
    SNR = cellfun(@SNRfunc,...
        specs, noise, repmat({range},size(specs)),...
        'UniformOutput', 0); 
elseif method == 1
    SNR = cellfun(@SNRfunc_real,...
        specs, noise, repmat({range},size(specs)),...
        'UniformOutput', 0); 
end

% Convert to matrix.
SNR = cell2mat(SNR);


function SNR = SNRfunc(spectrum, noise, range)
% Calculate SNR
SNR = mean(max( abs(spectrum(range(1):range(2))),[],1) ./ noise);

function SNR = SNRfunc_real(spectrum, noise, range)
% Calculate SNR
SNR = mean( max(real(spectrum(range(1):range(2))),[],1) ./ noise);

function noise = spectrum_noise(spectrum, mask)
% Calculate noise
noise = abs( std( spectrum(mask,:),[],1 ) );