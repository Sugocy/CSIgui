function output = allCombinations(input)
% Create a combination of given input lists.
% Example: {[0,1],[0,1],[0,5]}
%
% Usefull to find combinations of indexes.

% List view
% input = flipud(input(:));

% Convert to grid matrix (cell)
c = cell(1, numel(input)); [c{:}] = ndgrid(input{:});
% Replace
output = cell2mat(cellfun(@(v)v(:), fliplr(c), 'UniformOutput',false));

% This creates equal output as combvec-fcn
output = fliplr(output)';
