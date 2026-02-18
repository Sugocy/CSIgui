function output = allCombinations(input, varargin)
% Create a combination of given input lists.
% Example: {[0,1],[0,1],[0,5]} or 
% input = allCombinations([1,2,3], [3,2], [1:12]);
%
% Usefull to find combinations of indexes.

% List view
% input = flipud(input(:));

if nargin > 1, input = cat(2, {input}, varargin); end

% Convert to grid matrix (cell)
c = cell(1, numel(input)); [c{:}] = ndgrid(input{:});
% Replace
output = cell2mat(cellfun(@(v) v(:), fliplr(c), 'UniformOutput',false));

% This creates equal output as combvec-fcn
output = fliplr(output)';
