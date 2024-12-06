function win3D = window3D(wx, wy, wz)
%%%% Description:                     Create a 3D filter window.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-06
%%% --------------------------------------------------------------------
%%%
%%% Input: WX - 1D window in row/x-direction
%%%        WY - 1D window in col/y-direction
%%%        WZ - 1D window in pag/z-direction
%%%
%%% If no WY is given, the transposed 1D window WX is used for WY and
%%% 3rd dimension permuted WX is used for WZ.
%%%
%%% This allows for differently sized windows.
%%%
%%% !!!!! Be aware: WY and WZ need permutations! This is done during
%%% calculations.
%%%
%%% https://dsp.stackexchange.com/questions/19519/extending-1d-window-functions-to-3d-or-higher

%% Analyse input.
if nargin == 1, wy = wx; wz = wx; end

%% Permute data
wy = wy'; 
wz = permute(wz,[3 2 1]);

%% Outerproduct method

% Create 2D window using x and y e.g. r and c.
win2D = bsxfun(@times, wx, wy);                   % Transpose WY for outerpr
win3D = bsxfun(@times, win2D, wz);                % Transpose WZ for... 
