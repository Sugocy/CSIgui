function [zc, mxmn] = findZeroCross(data)
% Find zero crossing by searching the data-value closest to zero and
% returning this index as zc, zero-crossing.
%
% mxmn = logical array describing zero-crossing over index, with zero
% before and one at and after zerocrossing.

% Spectrum real data and abs-values
pval = abs(real(data)); 

% Find primary zero-crossing: closest to zero.
[~, zc] = min(pval);

% Find maximum value
[~, max_val_ind] = max(real(data)); 

% Correct if the zero-crossing is found in a region before the maximum.
breakWhile = 0;
while zc < max_val_ind
    zc_prev = zc; [~, zc] = min(pval(zc_prev+1:end));     
    zc = zc + zc_prev - 1;        
    if zc == zc_prev, breakWhile = 1; break; end
end
if zc == 0, zc = 1; breakWhile = 1; end
if zc == size(data,2)-1, breakWhile = 1; end


% Create max-min logic array to split spectra before and after
% zero-crossing. Can be used to phase spectra such that after
% zero-crossing the spectra go through zero.
sz = size(data,2); mxmn = zeros(1,sz);
if ~breakWhile, mxmn(zc+1:end) = 1; end
mxmn = logical(mxmn);