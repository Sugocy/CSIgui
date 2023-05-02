function [data, win] = csi_apodization3D(data, index, type)
%%%% Description:                            Apply 3D filter to MRSI data.
%%% Creator: Ir. Q. van Houtum       Version: 1.3          Date: 2017-07
%%% --------------------------------------------------------------------
%%%
%%% [data, window] = csi_apodization3D(data, spatial indexes, window type);
%%%
%%% Available types: Hamming, Hann, Blackman or Flattop.
%%%
%%% Calculate the correct 3D filter window from 3 individual 1D windows
%%% with sizes equal to data sizes @ index and apply the filter to every 
%%% dimension.
%%%
%%% Not be used for 1D spectra. The 3D filter for apodization is applied to
%%% spatial MRS data in k-space prior to spatial FFT. 
%%%
%%% SOURCES: 
%%%
%%% Correct sinc-like PSF with wide side-lobes resulting from FFT on 
%%% truncated data.
%%% Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2235821/
%%%
%%% 3D Window mathematics
%%% https://dsp.stackexchange.com/questions/19519/extending-1d-window- ...
%%% functions-to-3d-or-higher

if nargin == 2, type = 'hamming'; end

% 1D/2D/3D
dim = numel(index); 

% Loop each dimension and create 1D window
win_vec = cell(1,size(index,2));
for kk = 1:size(index,2) 
    win_vec{kk} = window1D(size(data,index(kk)),type, 'Periodic');
end

% Create 2D or 3D window.
if     dim == 2, win = window2D(win_vec{:});    
elseif dim == 3, win = window3D(win_vec{:});
end

% Permutation of window array to specified index in data.
sz = size(data); index_bool = zeros(1,numel(sz)); index_bool(index) = 1;
n = numel(sz); win_perm = 1:n; win_perm(index) = 1:numel(index); 
win_perm(~index_bool) = numel(index)+1:numel(sz);

% Permute the window to "index" dimensions
win = permute(win,win_perm);

% Repmat vector to repeat it for all values on index(1); 
% For a spectrum this means the sample index.
sz = size(data);
index_bool = zeros(1,numel(sz)); index_bool(index) = 1;
rep_vec = ones(1,numel(sz)); rep_vec(~index) = 0;
rep_vec(~index_bool) = sz(index_bool == 0 );
% Apply repetition vector to window
win = repmat(win, rep_vec);

% Apply hamming
data = data .* win;

