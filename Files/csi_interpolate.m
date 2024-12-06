function data = csi_interpolate(data, new_dim, spat_dim, interp_method)
%%% Interpolate 1D, 2D or 3D data to a new resolution specified with 
%%% new_dim. The first index must contain the FID/spectrum samples, the 
%%% consequent indexes are the spatial dimensions for either 1D/2D/3D space. 
%%% Dimensionality is found by the #elements in new_dim. If this condition 
%%% is not met, the spat_dim input must be used to set the spatial 
%%% dimensions and data will be rearranged before interpolation; 
%%%
%%% The default interpolation method is spline.
%%% Using interp_method input, a different method can be set: 
%%%                             linear, nearest, spline, cubic and makima.
%%%
%%% Be aware, makima, cubic and spline are not fast interpolation methods!
%%%
%%% Input
%%% data:           [#samples x (Kx) x (Ky) x (Kz)] data matrix.
%%% new_dim:        vector [1 x nD] describing the interpolated 
%%%                 data-dimensions.
%%% spat_dim:       indexes of the spatial dimensions of [1 x nD];
%%% interp_method:  the interpolation method as a string: linear, nearest, 
%%%                 spline (default), cubic and makima.
%%%
%%% Algorithm NFO:
%%% The interpn function in Matlab interpolates as one big volume whilst
%%% the spectroscopy data usually is not one correlated volume. It can
%%% include averages, different flip angles, multiple acquisitions, 
%%% channels and more. Therefor, the data is interpolated over the 
%%% spatial dimensions of the data, iterating over all extra dimensions 
%%% seperatly.
%%%
%%% Dr. Quincy van Houtum, ELH-Institute Essen, 2023/04.
%%% quincyvanhoutum@gmail.com

% Process input % --------------------------------------------- %
if nargin < 4, interp_method = 'spline'; end
if nargin < 3, spat_dim = []; end


% Rearrange data % -------------------------------------------- %
if ~isempty(spat_dim)
    dim = size(data);

    % Permute vector
    perm_vec = [1 spat_dim];

    % Remaining indexes for permute vector
    remainder_ind = find(ismember(1:numel(dim), perm_vec) == 0);
    perm_vec = [perm_vec remainder_ind];

    % Permute
    data = permute(data, perm_vec);
end


% Data dimensions % ------------------------------------------- %
dim = size(data);
dim_new = [dim(1) new_dim];

% Dimensionality
dim_type = numel(new_dim); % 1D/2D/3D


% Data to cell-array % ---------------------------------------- %
cell_layout_spat = num2cell(dim(1:dim_type+1));
cell_layout_other = arrayfun(@ones,...
    ones(1,size( dim(dim_type+2:end),2)), dim(dim_type+2:end), 'Uniform', 0);
data = mat2cell(data, cell_layout_spat{:}, cell_layout_other{:});


% Data grid points % ------------------------------------------ %
dim_spat = dim(1:dim_type+1); % #S x Kx x Ky x Kz

% Vectors of the index range in all dimensions - original
dim_vec = arrayfun(@linspace, ...
    ones(size(dim_spat)), dim_spat, dim_spat, 'uniform', 0);

% Gridpoints for all points in space - original dimensions
dim_grid = cell(1,size(dim_spat,2)); % Container for output ndgrid
[dim_grid{:}] = ndgrid(dim_vec{:});


% Vectors of the new index ranges in all dimensions
dim_vec_new = arrayfun(@linspace, ones(size(dim_new)), ...
                            dim_spat, dim_new, 'uniform',0);

% Gridpoints for alll points in space - new dimensions
dim_grid_new = cell(1,size(dim_new,2)); % Container for output ndgrid
[dim_grid_new{:}] = ndgrid(dim_vec_new{:}); 


% Interpolate % ---------------------------------------------- %

% List of voxels.
dim_pre_reshape = size(data);
data = reshape(data,[],1); 

nvox = size(data,1); 
for vi = 1:nvox    
    % Interpolate over grid
    data{vi} = ...
        interpn(dim_grid{:}, data{vi}, dim_grid_new{:}, interp_method);          
end

% Revert reshape
data = reshape(data, dim_pre_reshape);

% Revert mat2cell
data = cell2mat(data);

% Rearrange data % -------------------------------------------- %
if ~isempty(spat_dim), data = permute(data, perm_vec); end

