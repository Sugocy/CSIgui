function win2D = window2D(wx, wy)
%%%% Description:                     Create a 2D filter window.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-06
%%% --------------------------------------------------------------------
%%%
%%% Input: WX - 1D window in row/x-direction
%%%        WY - 1D window in col/y-direction
%%%
%%% WY, according to convention, should be transposed 1D window.
%%%
%%% If no WY or WZ is given, the transposed 1D window WX is used for WY.
%%% This method of passing input allows for differently sized windows. 
%%%
%%%

%% Analyse input.
if nargin == 1, wy = wx'; end

%% Outerproduct method
% Create 2D window using x and y e.g. r and c.
win2D =  bsxfun(@times, wx, wy);

% No transpose necessary if method:
% [mr , mc] = meshgrid(wx, wy);
% win2D = mr .* mc;


