function spec = csi_fft(fid, corr_N, corr_onesided, double_shift)	
%%%% Description:                    Forward Fourier of FIDs in CSI volume.
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%% Backward fourier of spectra ordered in space (1-ND). Expects the
%%% frequency dimension/FID samples to be on index 1.
%%%
%%% Input:      spec - Array with each FID on the first dimension.
%%%                    (fres x M x N x P x ...) 
%%%
%%% Used algorithm: fftshift( fft( FID ));
%%%                 Calculate the fft of the FID and shift the zero-freq to
%%%                 the center of the spectrum. Each FID will be converted 
%%%                 to a cell-array with size (M x N x P ...) whereafter 
%%%                 the fft-method is applied on each cell. No loops, fast.
%%%
%%% corr_N:         Correct for the 1/#samples (N) in the FFT. Off by
%%%                 default.
%%% double_shift:   Shift the FID twice for handling (symmetric) echoes. 
%%%                 Disables onesided_cor if enabled.
%%% corr_onesided:  Correct for onesided FFT in spectroscopy. On by
%%%                 default. Sample at t = 0 is halfed before fft.
%%%
%%% Contact: quincyvanhoutum@gmail.com
%%% See also: csi_ifft(spec);

if     nargin == 1,     corr_N = 0; corr_onesided = 1; double_shift = 0; 
elseif nargin == 2,                 corr_onesided = 1; double_shift = 0;
elseif nargin == 3,                                    double_shift = 0;
end

%% Convert FID to cell if array

if ~iscell(fid)
    was_cell = 0;
    
    % Exclude first-e.g. time-dimension.
    % Create vector of ones for each dimension. Cell/FID.
    sz = size(fid); % [ Dim1 Dim2 ... DimN ]; Dim1 = FID-signal.
    cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

    % Creat cell of data.
    fid = mat2cell(fid, sz(1), cell_layout{:}); 
else
    was_cell = 1;
end

%% FFT over all FIDS.
% 1. Check for pre-fft shift
% 2. Check for onesided FFT correction
% 3. FFT
% 4. Shift 
% 5. Check for N-factor correction

% Apply pre-shift (applicable to echoes)
if double_shift
    fid = cellfun(@fftshift, fid, 'UniformOutput', 0); 
    corr_onesided = 0;
end

% Correct for t = 0 because of non-symmetrical e.g. one-sided fft. 
if corr_onesided
    fid = cellfun(@(x) x .* [0.5 ones(1,size(x,1)-1)]', fid, 'Uniform', 0);
end

% FFT
spec = cellfun(@fft,      fid,  'UniformOutput', 0);                        % Fast Forward Fourier

% FFT Shift
spec = cellfun(@fftshift, spec, 'UniformOutput', 0);                        % Shift zero-freq to centrum of spectrum

% Correct N-factor 
if corr_N
    spec = cellfun(@times, ...
        spec, repmat({1/sz(1)},size(spec)),'UniformOutput',0);              % Swap the 1/N factor convention in the inverse Fourier
end


%% Convert SPEC to array
% ONLY IF REQUIRED

if ~was_cell
    % Undo squeeze of first-dimension.
    % spec = permute(spec, [numel(size(spec))+1 1:numel(size(spec))]);
    % Convert.
    spec = cell2mat(spec);
    if numel(sz) <= 2, spec = squeeze(spec); end
end