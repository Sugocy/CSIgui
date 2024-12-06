function [data, permv, szr] = csi_combine_reshape(data, channel_index)
% This function reshapes the data for use with WSVD and other coil 
% combinations functions. It permutes the channel index to the second 
% dimension and creates a cell of the FID/SPEC for every channel per voxel.
% Returns all required variables to undo this operation  after any 
% calculations.
%
% Input:
%   data array 
%   index of the channels
%
% Output:
%   data  = {nDim x nChan} x nVoxels
%   permv = permute vector, used to revert this operation.
%   szr   = size of reshaped matrix
%
% Quincy van Houtum, PhD; 07/2023
% quincyvanhoutum@gmail.com
    
% Permute vector - channel index on index 2.
sz = size(data);                                      

% Create permute vector and permute data
rm_ind = 2:numel(sz); % Remainder indices vector
rm_ind(rm_ind == channel_index) = []; % Remove channel-index         
permv = [1 channel_index]; permv = cat(2, permv, rm_ind); % Permute vector.
data = permute(data, permv); % Permute

% Reshape to dimensional cell array: dt x channels x nVoxels
% Array size per dimension
szr = size(data); 
% Add a dimension if only one voxel.
if numel(szr) < 2, szr(3) = 1; end 
% Cell layout: {dt x nchan} x nCells
cell_layout = arrayfun(@ones,...
    ones(1,size(szr(3:end),2)),szr(3:end),'UniformOutput',0);
% Convert to a cell matrix with {dt x nChan} x [other indexes];
data = squeeze(mat2cell(data, szr(1), szr(2), cell_layout{:})); 

% Create a cell-list {dt x nChan} x nVoxels
data = reshape(data, [], 1); % Create a list of all voxels.