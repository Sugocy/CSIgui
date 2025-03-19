function fid_zf = csi_zeroFill(fid, N, append)
%%%% Description:                                Zero fill spectral data.
%%% Creator: Ir. Q. van Houtum       Version: 1.3          Date: 2017-07
%%% --------------------------------------------------------------------
%%% Zero fill data to increase sensitivity for lower SNR peaks. Data
%%% expected to be in timedomain e.g. FID/Echo.
%%%
%%% fid (dt x M x N ...);  time domain data (FID)
%%% N   (1xN)           ;  requested size of data/fid after zero filling.
%%% append              ;  'pre'  append zeroes to end (default)
%%%                        'post' append zeroes to begin
%%%                        'both' append on both sides
%%%                        If append set to both, N/2 is added at the start
%%%                        and end of the FID.
%%% 
%%%
%%% Contact: qhoutum2@umcutrecht.nl

if nargin == 2, append = 'post'; end
append = lower(append);

% Prepare data and options % ------------------------------------------- %

% Data size.
sz = size(fid); % [ Dim1 Dim2 ... DimN ]; Dim1 = FID-signal.

% Calculate nr of zeroes to add
Nadd = N-sz(1); 
if strcmp(append,'both'), Nadd = Nadd/2; end

% Go from array to cell/FID
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
% Creat cell of data.
fid = mat2cell(fid, sz(1), cell_layout{:}); fid = squeeze(fid); % Time gone.


% Pad the FID array in each cell % ------------------------------------- %
fid_zf = cellfun(@padarray, fid, ...
    repmat({Nadd},size(fid)),repmat({append},size(fid)),...
    'UniformOutput', 0);  


% Recreate array from cell % ------------------------------------------- %

% Undo squeeze of first-dimension.
fid_zf = permute(fid_zf, [numel(size(fid_zf))+1 1:numel(size(fid_zf))]);
% Convert to array
fid_zf = cell2mat(fid_zf);
% 1D Correction
if numel(sz) <= 2, fid_zf = squeeze(fid_zf); end