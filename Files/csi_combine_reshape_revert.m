function [data, restore_permv] = csi_combine_reshape_revert(data, permv, szr)
% This function reverst the data back to its original dimensional layout as
% before csi_combine_reshape w.r.t. the change to the coil/channel index.
%
% NB. The channel/coil dimension is expected to be reduced to one by
% combining channels.

% Reshape the data
szr(2) = 1; szr(1) = size(data,1); 
data = reshape(data, szr); 

% Create undo-permute-permute-vector
nindex = numel(permv); restore_permv = NaN(1,nindex);
for kk = 1:nindex, restore_permv(kk) = find(kk == permv); end

% Permute the data
data = permute(data, restore_permv);