function [data, phaseroll] = csi_voxelshift(data, shift, ind, isCentered)
% Apply a voxelshift in the spatial domain by adding a linear phaseroll to 
% the data in the frequency domain (k-space).
%
% Input
%   data:       multidimensional complex array with spatial dimension at  
%               indexes described in ind. Assumes the k-space data is 
%               centered i.e. zero-frequency component in the center, 
%               otherwise "isCentered" should be set to 0.
%   shift:      Integer #voxels of to shift "data" for each spatial 
%               dimension in ind as a list [1 x numel(ind)].
%               
%               Shift direction, axes origin at top left of volume.
%               Voxel left/right:   [-1, +1]
%               Voxel up/down:      [+1, -1]
%               Voxel forw/backw:   [-1, +1] 
%                       forward = toward the last "slice". Slice #1 will
%                       become #Last.
%
%   dim:        [1 x N] list of indexes of the spatial dimensions in data.
%   isCentered  Optional input to shift the zero-frequency component to the
%               center. Default, 1. Set to 0 if a shift is required.
%   
% Output
%   data:       input data multiplied by a linear phaseroll
%   
%
% Method
% Applying a phase-roll in k-space is equivalent to a linear shift 
% (translation) in the image domain. In MRI, a phase-roll can be applied 
% to k-space data by multiplying the complex k-space matrix by a linear 
% phaseroll. Multiply k-space voxel(s) by a phaseroll defined by:
%
%  phaseroll(n,m,l) = exp((-/+)1i .* 2pi .* (kx[n]dx + ky[m]dy + kz[l]dz))
%
% where n,m,l are the k-space domain index per dimension and dx, dy, dz
% the shift in space. (-/+) defines the direction of the shift. Rewriting 
% it more descrete for a per voxel shift:
%
%  phaseroll(n,...) = exp((-/+)1i .* 2pi .* ( [n - Nx/2 / Nx] + ...)
% 
% where n is the k-space domain index (per dim) and Nx is the matrix size
% for that dimension.
%
%
% Quincy van Houtum, v01.2026
% Contact: quincyvanhoutum@gmail.com


if nargin < 3, isCentered = 1; end

% Possible safety implementation
% isInteger = mod(shift,1) == 0;

% Move the "ind" dimensions to the first indexes.
if ~isequal(ind, 1:numel(ind))
    % Reorder the data array
    sz = size(data); pvec = 1:numel(sz); pvec(ismember(pvec,ind)) = [];
    pvec = cat(2,ind, pvec); data = permute(data, pvec);
end

% Sizes --> What if 2D etc!
sz = size(data);    % Kx, Ky, Kz, ...
sz_spat = sz(1:3);  % Nx, Ny, Nz - kspace matrix size

% Create range vectors for k-space
% Depends on centered k-space i.e. if the zero frequency comp at center.
if isCentered
    range = arrayfun(@(x) (-floor(x./2):(ceil(x./2) - 1)) ./ x, ...
        sz_spat, 'UniformOutput', 0);
else
    range = arrayfun(@(x) (0:x - 1)./x, sz_spat, 'UniformOutput', 0);
end

% Grid for each k-space voxel
[mx, my, mz] = ndgrid(range{:}); % NB. mesh vs. ndgrid differences!  

% K-space grid and #voxels linear shift
phaseroll = (shift(1) .* mx + shift(2) .* my + shift(3) .* mz);
phaseroll = exp(1i .* 2.*pi .* phaseroll);

% Multiply every FID per voxel by the phase roll.
% tmp = permute(tmp, [2:numel(sz) 1]); % already away.
data = data .* phaseroll;

% Revert to original dimension order
if ~isequal(ind, 1:numel(ind))
    data = ipermute(data, pvec); phaseroll = ipermute(phaseroll, pvec);
end


% // --- Previous checksum 2025
%
% Direction of shift: 
% minus = move right +x, plus = move left  -x      // checked
% minus = move down  -y, plus = move up    +y      // checked
% minus = move forw  +z, plus = move forw  -z      // checked
% backwards = move volume towards the first "slice" (#1 --> #End, etc,)
% forwards  = move volume towards the last "slice"  (#End --> #1, etc,)
%
% shiftd = +1; % Shift direction from sign of shiftv            
% Original: phaseroll = exp(shiftd .* 1i .* 2.*pi .* phaseroll);