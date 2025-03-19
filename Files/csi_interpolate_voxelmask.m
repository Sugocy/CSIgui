function data = csi_interpolate_voxelmask(data, new_dim, spat_dim, interp_method)
%%% Interpolate 1D, 2D or 3D voxel mask to a new resolution specified with 
%%% new_dim. Expected is the spatial dimensions on the first 3 indexes. 
%%% Dimensionality is found by the #elements in new_dim. If this condition 
%%% is not met, the spat_dim input must be used to set the spatial 
%%% dimensions and data will be rearranged before interpolation; 
%%%
%%% The default interpolation method is nearest - for best [0; 1] interp.
%%% Using interp_method input, a different method can be set: 
%%%                             linear, nearest, spline, cubic and makima.
%%%
%%% Be aware, makima, cubic and spline are not fast interpolation methods!
%%%
%%% Input
%%% data:           [(Kx) x (Ky) x (Kz)] data matrix.
%%% new_dim:        vector [1 x nD] describing the interpolated 
%%%                 data-dimensions.
%%% spat_dim:       indexes of the spatial dimensions of [1 x nD];
%%% interp_method:  the interpolation method as a string: linear, spline, 
%%%                 nearest (default), cubic and makima.
%%%
%%% Algorithm NFO:
%%% The interpn function in Matlab interpolates as one big volume whilst
%%% the spectroscopy data usually is not one correlated volume. It can
%%% include averages, different flip angles, multiple acquisitions, 
%%% channels and more. Therefor, the data is interpolated over the 
%%% spatial dimensions of the data, iterating over all extra dimensions 
%%% seperatly.
%%%
%%% This function is a slightly adapted version of csi_interpolate,
%%% correcting for sample-index size of 1.
%%%
%%% Dr. Quincy van Houtum, ELH-Institute Essen, 2024/12.
%%% quincyvanhoutum@gmail.com


% Process input % --------------------------------------------- %
if nargin < 4, interp_method = 'linear'; end
if nargin < 3 || isempty(spat_dim), spat_dim = [1 2 3]; end

% Data dimensions % ------------------------------------------- %
dim = size(data); dim_new = new_dim;

% Dimensionality
dim_type = numel(new_dim); % 1D/2D/3D


% Rearrange data % -------------------------------------------- %
if ~isempty(spat_dim)
    dim = size(data);

    % Permute vector
    perm_vec = spat_dim; % CHANGE HERE NOT [1 spat_dim]

    % Remaining indexes for permute vector
    remainder_ind = find(ismember(1:numel(dim), perm_vec) == 0);
    perm_vec = [perm_vec remainder_ind];

    % Permute
    data = permute(data, perm_vec);
end

% Data to cell-array % ---------------------------------------- %
cell_layout_spat = num2cell(dim(spat_dim));
cell_layout_other = arrayfun(@ones,...
    ones(1,size( dim(dim_type+1:end),2)), dim(dim_type+1:end), 'Uniform', 0);
data = mat2cell(data, cell_layout_spat{:}, cell_layout_other{:});


% Data grid points % ------------------------------------------ %
dim_spat = dim(1:dim_type); % #S x Kx x Ky x Kz
% HERE SOME CHANGE - NOT + 1

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
        interpn(dim_grid{:}, single(data{vi}), dim_grid_new{:}, interp_method);          
end

% Revert reshape
data = reshape(data, dim_pre_reshape);

% Revert mat2cell
data = cell2mat(data);

% Rearrange data % -------------------------------------------- %
if ~isempty(spat_dim), data = permute(data, perm_vec); end

% Create [0/1] boolean
data(data <= 0.5) = 0; data(data > 0.5) = 1; data = logical(data);
