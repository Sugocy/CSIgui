function [data, restore_permv] = csi_combine_reshape_revert(data, permv, szr)
% This function reverst the data back to its original dimensional layout as
% before csi_combine_reshape w.r.t. the change of the coil/channel index.
%
% Input
% data  = array or cell-matrix with layout {nS x nCh} x nVox
% permv = permute vector created in csi_combine_reshape to move the channel
%         index to the second position
% szr   = array size before mat2cell-fnc in csi_combine_reshape
%
% NB. The channel/coil dimension is expected to be reduced to one by
% combining channels.
%
% Quincy van Houtum, PhD; 08/2023
% quincyvanhoutum@gmail.com 

% Revert cell to matrix % ----------------------------------------- %
if iscell(data)
    % Calculate permute vector to set exlcuded-cell indexes at correct
    % dimension.
    perm_index = [numel(size(data)) + [1 2] 1:numel(size(data))];
    data = permute(data,perm_index);
    
    % Apply cell to mat fcn
    data = cell2mat(data);
end

% Reshape to pre-cell-matrix size, changing channel-index
szr(2) = 1; szr(1) = size(data,1); 
data = reshape(data, szr);

% Create undo-permute-permute-vector
nindex = numel(permv); restore_permv = NaN(1,nindex);
for kk = 1:nindex, restore_permv(kk) = find(kk == permv); end

% Permute the data to original index order
data = permute(data, restore_permv);