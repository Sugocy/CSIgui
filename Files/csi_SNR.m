function SNR = csi_SNR(spectrum, mask_sz, snr_unit, range, mask_side)
% Description:                 Calculate SNR of multidimensional MRSI data
% Creator: Dr. Q. van Houtum       Version: 1.9          Date: 2025-01
% --------------------------------------------------------------------
% Calculate SNR for given spectra using the formula shown below input
% description.
%
% Input:
%       spectrum;   Array of all spectra of interest.(dHz x N x M x etc.)
%       mask_size;  Value for the mask size used to calculate noise. If no 
%                   mask_size is given a defeault mask size of 50 is used.
%                   OR the noise given per  voxel.
%       snr_unit;   Either using the absolute (0), real(1, default) or
%                   imaginary part of spectrum
%                   spectrum
%       range;      Range of spectrum to calculate SNR from. Noise is still
%                   calculated using the noise mask. range = [start end]
%       mask_side;  Calculate the noise mask using (0) both sides, 
%                   (1) start or (2) end of the spectrum. Both-sides is 
%                   the default. 
%
% Formula: SNR = max ( abs ( spec ))) ./ abs(std( abs(  spec)))) --> (0)
% Formula: SNR = max ( real( spec ))) ./ abs(std( real( spec)))) --> (1)
% Formula: SNR = max ( imag( spec ))) ./ abs(std( imag( spec)))) --> (2)
%
% Q. van Houtum, PhD; version 2.3 01/2025
% Contact: quincyvanhoutum@gmail.com

% Process #input arguments
if     nargin == 1, mask_sz = 50; snr_unit = 0; 
elseif nargin == 2, snr_unit = 0; 
end
if nargin <= 3, range = [1 size(spectrum,1)]; end
if nargin <= 4, mask_side = 0; end

% Data size
sz = size(spectrum); 

% NaN position
nan_ind = isnan(spectrum);

% Noise mask
doMask = 0;
if sum((mask_sz(:) - floor(mask_sz(:))) == 0) > 0 && numel(mask_sz) == 1 

    if mask_side == 0
        mask_size_double = [round(mask_sz./2) mask_sz-round(mask_sz./2)];
        mask = [1:mask_size_double(1) (sz(1) - mask_size_double(2) + 1):sz(1)];
    elseif mask_side == 1
        mask = 1:mask_sz;
    elseif mask_side == 2
        mask = (sz(1)-mask_sz + 1):sz(1);
    end
    
    doMask = 1;
end

% Create cell-array of data to use cellfun
% Exclude first index dimensions
% Create vector of ones equal in size to each remaining dimension.
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
specs = (mat2cell(spectrum, sz(1), cell_layout{:}));

if doMask
    % Pre-calculate the noise
    noise = cellfun(@spectrum_noise,...
        specs, repmat({mask}, size(specs)), ...
        repmat({snr_unit}, size(specs)),'UniformOutput',0);
else    
    % Use given noise-value(s)
    if ~iscell(mask_sz)
        noise = (mat2cell(mask_sz, 1, cell_layout{:}));
    end
end

% SNR of each spectrum.
SNR = cellfun(@SNRfunc,...
    specs, noise, repmat({range},size(specs)),...
    repmat({snr_unit}, size(specs)), 'UniformOutput', 0); 

% Convert to matrix.
SNR = cell2mat(SNR); SNR(nan_ind) = NaN;


function SNR = SNRfunc(spectrum, noise, range, unit)
% Calculate SNR
if unit == 1
    SNR =  max(real(spectrum(range(1):range(2))),[],1) ./ noise;
elseif unit == 2
    SNR =  max(imag(spectrum(range(1):range(2))),[],1) ./ noise;
else
    SNR = max( abs(spectrum(range(1):range(2))),[],1) ./ noise;
end 


function noise = spectrum_noise(spectrum, mask, unit)
% Calculate noise
if unit == 1
    noise = abs( std( real(spectrum(mask,:)),[], 1) );
elseif unit == 2
    noise = abs( std( imag(spectrum(mask,:)),[], 1) );
else
    noise = abs( std( abs(spectrum(mask,:)),[], 1) );
end