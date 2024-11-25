function spec = csi_pca_denoising(spec, chan_ind, patch_size, svd_method)
% For a more detailed explanation of this method see sources below.
% 
% Input:        volume of spectra (3D!), channel index, odd patch-size.
% Output:       volume of spectra PCA denoised.
% 
% Data will be reshaped to nSamples x nVoxels x nChan for calculations,
% therefor it requires the index of the channels/coils. Other
% data-dimensions should be excluded and calculated separately! 
%
% This function uses parallel computing to speed up calculations.
%
% This function can use a faster SVD method with 10e-10 lower accuracy:
% svd_method [0], default Matlab or svd_method [1], fast custom.
%
% Explanation of PCA denoising method:
% Generate a local-patch with patch-size and convolve over the data set.
% For each voxel in the volume, all voxels within this patch are 
% concatenated into matrix M: rows represent the voxels in the local-patch,
% the columns represent the real and imaginary values (concatenated). 
% Singular Value Decomposition is applied to matrix M to find eigenvalues 
% using SVD-output matrix S: its elements are the square root of the
% eigenvalues. All eigenvalues within the Marchenko-Pastur distribution
% (MPD) are set to zero. The threshold for what its within the MPD can be 
% estimated according Veraart et al: when the variance after summing 
% N eigenvector-components is larger than the variance of the difference 
% between the principle components and the summated eigenvectors, then N-1
% is the threshold. Matrix M is reconstructed and the denoised local voxels
% are added to a new matrix, summating the denoised patched for every voxel
% patch-size cubed. After one run through all voxels, the data is
% averaged by patch-size cubed. This is performed for every channel 
% separately using a parallel computing method to speed up calculations.
%
% Sources:
% PCA denoising, Froeling et. al. doi.org/10.1002/mrm.28654
% Denoising, Veraart et. al. 10.1016/j.neuroimage.2016.08.016
%
% Q. van Houtum, 2023;
% quincyvanhoutum@gmail.com

 
% Default patch-size
if nargin < 3, patch_size = 5; svd_method = 0; end
if nargin < 4, svd_method = 1; end

% Data dimensions
dim = NaN(1,10); % This makes it compatible with Matlab <R2021
for kk = 1:10, dim(kk) = size(spec,kk); end

dim_spat = dim; dim_spat(chan_ind) = []; dim_spat(1) = [];

% Reshape data nS x nVox x nCh % --------------------------------------- %

% CH to outer dimensions
if chan_ind < numel(dim)                                       
    % Create permute vector and permute data
    permv = 1:numel(dim);               % Remainder indices vector
    permv(permv == chan_ind) = [];      % Remove channel-index         
    permv(end+1) = chan_ind;            % Put channel-index at end
    spec = permute(spec, permv);        % Permute
end

% Convert data to nS x nVox x nCh,
spec = reshape(spec, dim(1), [], dim(chan_ind));
                
% Calculate linear indices for a neighbourhood at a voxel:
lin_ind = csi_pca_denoising_neighbourhoods(dim_spat, patch_size);

psz = patch_size; 

% Parallel Pool % ------------------------------------------------------ %

mlyear = version('-release'); mlyear = str2double(mlyear(1:end-1));
parpoolName = 'Processes'; if mlyear <2023, parpoolName = 'local'; end

p = gcp('nocreate'); 
if isempty(p)
    nCores = floor(feature('numcores') .* 0.75);
    parpool(parpoolName, nCores);
end

tic
% Loop over each channel (chi) % --------------------------------------- %
nS = dim(1);

parfor chi = 1:dim(chan_ind)                % PARALLEL LOOP %
% for chi = 1:dim(chan_ind)
    data_single_chan = spec(:,:,chi);
    data_single_chan_denoised = zeros(size(data_single_chan));
    
    vox_patch = NaN(psz.^3, nS*2);
    % Loop over each voxel % ------------------------------------------- %
    for vi = 1:size(data_single_chan,2)                
        % Get all neighbourhood voxels at voxel vi of size patch_size in
        % directions (xyz).
        
        % Linear indices of the neighbourhood for this voxel.
        lin_nbh = lin_ind{vi};
        
        % -------------------------------- %

        % Get the local-voxels neighbourhood around voxel vi
        nbh = data_single_chan(:,lin_nbh);
            
        % Make a large matrix with real and imaginary on the col-index and
        % voxels on row-index
        vox_patch(:,1:nS) = real(nbh)';
        vox_patch(:,nS+1:end) = imag(nbh)';
        
        % ---
        % Apply Denois Matrix function (where the actual magic happens)
        % As a seperate function, needs the voxels in the local-volume.
        denoised_patch = ...
            csi_pca_denoising_marchenko_pastur(vox_patch, svd_method);
        % ---

        % Back to complex
        denoised_patch = complex(denoised_patch(:,1:nS)',...
                                 denoised_patch(:,nS+1:end)');
        
        % This is  wrong and flips the FID over the x-axis.
        % denoised_patch = complex(denoised_patch(:,1:nS),...
        %                          denoised_patch(:,nS+1:end))';

        % Summate the denoised data
        data_single_chan_denoised(:,lin_nbh) = ...
            data_single_chan_denoised(:,lin_nbh) + denoised_patch;
    end
    
    
    % Average and store denoised-data for this channel.
    spec(:,:,chi) = (data_single_chan_denoised ./ (psz.^3));

    

end
dt = toc; fprintf('PCA duration without parallel pool startup: %fs\n', dt);

% Undo Reshape and Permute % ------------------------------------------- %

% Reshape to nSamples x n[Spatial Dimensions] x nChan
spec = reshape(spec, [dim(1), dim_spat, dim(chan_ind)]);
% Permute back to original matrix size
if chan_ind < numel(dim) 
    % Create undo-permute-permute-vector
    nindex = numel(permv); restore_permv = NaN(1,nindex);
    for kk = 1:nindex, restore_permv(kk) = find(kk == permv); end
    % Permute
    spec = permute(spec, restore_permv);
end

% --- Used by PCA_Denoising
function lin_ind = csi_pca_denoising_neighbourhoods(dim_spat, patch_size)
% Calculate the linear indices of a neigbourhood of size "patch_size" for
% every spatial location in a (CSI) volume.
%
% PCA_Denoising convolves over full volume using a local neighbourhood or
% local patch for every channel/coil. This functions returns the linear 
% indices of a neighbourhood for every voxel allowing for a lookup approach 
% instead of calculating this every voxel/channel-iteration.
%
% uses sub2neighbourhood(sub, psz, dim);

% Patch size variable short naming
psz = patch_size;

% Number of voxels
N = prod(dim_spat);

% Prep var-containers
lin_ind = cell(N,1);

% Loop over every voxel (linear index)
for vi = 1:N
    % Get local-voxels neighbourhood indexes.

    % Step 1 - subindex of this linear index voxel vi
    sub = cell2mat(ind2subQ(dim_spat,vi));

    % Step 2 - subindexes of the neighbourhood (local indices/patch)
    sub_nbh = sub2neighbourhood(sub, psz, dim_spat)';

    % Step 3 - Convert to linear index for data-handling;
    sub_nbh = num2cell(sub_nbh);  lin_nbh = NaN(psz^3,1);
    for ni = 1:(psz^3)
        lin_nbh(ni) = sub2ind(dim_spat,sub_nbh{ni,:});
    end

    % Step 5 - Store in output var
    lin_ind{vi} = lin_nbh;
end

% --- Used by CSI_PCA_Denoising_Neighbourhoods
function ind_nbh = sub2neighbourhood(sub, psz, dim)
% Given a sub-index of specific volume of size dim, returns all
% neighbouring linear indexes with a neighbourhood size of psz such that 
% sub-index is the center; assuming psz is of odd-size.
%
% Applies a circular shift to coordinates outside of the dim-limits.

% Step 1 - calculate "surrounding" voxel-vector using patch-size
ind_nbh = NaN(numel(dim), psz);
for kk = 1:numel(dim)
    % Index
    ind_nbh(kk,:) = ...
        linspace(sub(kk)-(psz-1)/2,sub(kk)+(psz-1)/2, psz);
    
    % Correction for outside grid < (circular)
    ind_loc_vox_cor = (ind_nbh(kk,:) <= 0);
    ind_nbh(kk,ind_loc_vox_cor) = ...
        ind_nbh(kk,ind_loc_vox_cor) + dim(kk);
    
    % Correction for outside grid > (circular)
    ind_loc_vox_cor = (ind_nbh(kk,:) > dim(kk));
    ind_nbh(kk,ind_loc_vox_cor) = ...
        ind_nbh(kk,ind_loc_vox_cor) - dim(kk);
end

% Get all combinations
ind_nbh = allCombinations({ind_nbh(1,:),ind_nbh(2,:),ind_nbh(3,:)});