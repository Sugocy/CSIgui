function win3D = window3D(wx, wy, wz)
% Return a 3D  filter window given the 1D filter windows for every 
% dimension.
%
% Input: 
%   wx, wy, wz:     A 1D window for each dimension.
%
% Output:
%   win3D	        A 3D filter window constructed using given 1D
%                   window(s).
%
% If wy and  wz are NOT given, the transposed 1D window wx is used for wy 
% and 3rd-dimension permuted wx is used for wz. This allows for differently 
% sized windows in every dimension.
%
% NB. Windows for wy and wz need permutations which are applied during 
% calculations.
%
% Source:
% Stackexchange - extending-1d-window-functions-to-3d-or-higher.
%
% Quincy van Houtum. v06.2017
% quincyvanhoutum@gmail.com

% Analyse input.
if nargin == 1, wy = wx; wz = wx; end

% Permute data
wy = wy'; 
wz = permute(wz,[3 2 1]);

% Outerproduct method

% Create 2D window using x and y e.g. r and c.
win2D = bsxfun(@times, wx, wy);                   
win3D = bsxfun(@times, win2D, wz);                
