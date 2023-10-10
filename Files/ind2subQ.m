function subs = ind2subQ(matrix_size, index)
% Returns the converted subindex for each dimensions of a matrix from a
% linear index.
%
% Q. van Houtum

sz = repmat(matrix_size, [1 numel(matrix_size)]);
subs = cell([1 numel(matrix_size)]);  
[subs{:}] = ind2sub(sz, index);
