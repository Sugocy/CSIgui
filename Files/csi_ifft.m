function fid = csi_ifft(spec, corr_onesided)	
%%%% Description:                            Inverse Fourier of CSI data.
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%% Inverse or backward fourier of spectra ordered in space (1-ND). 
%%% Expects the frequence to be on index-dimension one!
%%%
%%% Input:      spec - Array with each spectrum on the first dimension.
%%%                    (fres x M x N x P x ...)
%%%
%%% Used algorithm: ifft( ifftshift( spectra ));
%%%                 Shift zero-freq to center before calculating the ifft.
%%%                 Each spectrum will be converted to a cell-array with
%%%                 size (M x N x P ...) whereafter the ifft-method is 
%%%                 applied on each cell. No loops, fast.
%%%
%%% corr_onesided:  Correct for onesided FFT in spectroscopy. On by
%%%                 default. Sample at t = 0 is double after ifft.
%%% 
%%% Contact: quincyvanhoutum@gmail.com                 
%%% See also: csi_fft(fid);               

if nargin == 1, corr_onesided = 1; end

%% Convert SPEC to cell
% ONLY IF REQUIRED

if ~iscell(spec)
    was_cell = 0;
    % Exclude first-e.g. time-dimension.
    % Create vector of ones for each dimension. Cell/FID.
    sz = size(spec); 
    cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

    % Creat cell of data.
    spec = mat2cell(spec, sz(1), cell_layout{:}); 
    % spec = squeeze(spec); % Time gone.
else
    was_cell = 1;
end

%% FFT over all FIDS.
% First inverse the data shift around 0 in csi_fft prior to ifft.
fid = cellfun(@ifftshift, spec,  'UniformOutput', 0);  % Inverse shift.   
fid = cellfun(@ifft,      fid,   'UniformOutput', 0);  % Inverse FFT                     

% Correct for t = 0 because of non-symmetrical e.g. one-sided fft. 
if corr_onesided
    fid = cellfun(@(x) x .* [2 ones(1,size(x,1)-1)]', fid, 'Uniform', 0);
end

%% Convert FID to array
% ONLY IF REQUIRED

if ~was_cell
    % Undo squeeze
    % fid = permute(fid, [numel(size(fid))+1 1:numel(size(fid))]);
    % Convert.
    fid = cell2mat(fid);

    % Undo if 1D data.
    if numel(sz) <= 2, fid = squeeze(fid); end
end