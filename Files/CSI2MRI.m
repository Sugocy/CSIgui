function mz_range = CSI2MRI(cz, mz, sz_cz, sz_mz)
% Returns min and max index of MRI-slices which lay around the CSI slice
% coordinate given.


% Number of MZ slices in a CZ slice
nMz_per_Cz = sz_cz/sz_mz;
nMz_per_Cz_step = floor(nMz_per_Cz/2);

% Index of slice MZ closest to CZ
[~, close_ind] = min(abs(mz-cz));

if close_ind == 1 % First slice

    % MRI slices in range.
    mz_range = [close_ind close_ind+nMz_per_Cz_step-1];
    
elseif close_ind > size(mz,1)-nMz_per_Cz_step
    
    % Close_ind + step = larger than Mz;
    mz_range = [size(mz,1)-nMz_per_Cz_step+1 size(mz,1)];
    
else
    
   % MRI slices in range.
   mz_range = [(close_ind-nMz_per_Cz_step+1) close_ind+nMz_per_Cz_step]; 
    
end


