function spec = csi_fft(fid, correct_N, double_shift)	
%%%% Description:                    Forward Fourier of FIDs in CSI volume.
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%% Backward fourier of spectra ordered in space (1-ND). Expects the
%%% frequence dimension to be on index 1.
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
%%% Contact: qhoutum2@umcutrecht.nl    
%%% See also: csi_ifft(spec);

if nargin == 1,     correct_N = 0; double_shift = 0;    
elseif nargin == 2,                double_shift = 0;
    
end

%% Convert FID to cell
% ONLY IF REQUIRED

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
% 2. FFT
% 3. Shift 
% 4. Check for N-factor correction

% Apply pre-shift (applicable to echoes)
if double_shift, fid = cellfun(@fftshift, fid, 'UniformOutput', 0); end

% FFT
spec = cellfun(@fft,      fid,  'UniformOutput', 0);                        % Forward Fast Fourier

% FFT Shift
spec = cellfun(@fftshift, spec, 'UniformOutput', 0);                        % Shift zero-freq to centrum of spectrum

% Correct N-factor 
if correct_N
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