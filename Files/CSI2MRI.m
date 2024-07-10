function mz_range = CSI2MRI(cz, mz, sz_cz, sz_mz)
% Returns min and max index of MRI-slices located in the CSI slice
% coordinate given.


% Number of MZ slices in a CZ slice
nMz_per_Cz = sz_cz/sz_mz;
% Step for ...from center of CZ, add and subtract half number in 1 CZ
nMz_per_Cz_step = floor(nMz_per_Cz/2);

% Index of slice MZ closest to CZ
[~, close_ind] = min(abs(mz-cz));

nMz_per_Cz_step_min = nMz_per_Cz_step - 1;
if nMz_per_Cz_step_min < 0, nMz_per_Cz_step_min = 0; end
min_of_slice = (close_ind - nMz_per_Cz_step_min); % bc if even...
if min_of_slice <= 0, min_of_slice = 1; end

max_of_slice = (close_ind + nMz_per_Cz_step);
if max_of_slice > size(mz,1), max_of_slice = size(mz,1); end

% MRI slices in range.
if (max_of_slice < min_of_slice) || (max_of_slice < 0)
    max_of_slice = min_of_slice + 1;
end

mz_range = [min_of_slice max_of_slice];
    