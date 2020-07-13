function [data, phase_shift] = csi_voxelshift(data, vox_shift, doi)
%%%% Description:                    Shift CSI space N voxels.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2018-08
%%% --------------------------------------------------------------------
%%% Shift the CSI space in the frequency domain by adding a phase shift to
%%% each voxel seperatly.
%%%
%%% Input:      data        - Array with each FID on the first dimension.
%%%                          (freq x M x N x P x ...)
%%%
%%%             vox_shift   - Number of voxels to shift
%%%
%%%             doi         - Dimensions of interest e.g. k-space
%%%
%%%                 
%%% Phase to voxel shift calculated from (45deg == 0.125 voxels).
%%%
%%% Contact: qhoutum2@umcutrecht.nl    

% Prepare general factors ---- %

% Get dimensions 
sz = size(data); dim = sz(doi); 
% General phase change for one voxel displacement
gen_phase = (1./dim) * (45./0.125);

% Phase change in radians
gen_phase_rad = gen_phase*(pi/180);


% Spatial Array ---- %

% Create a matrix with all indexes
range = arrayfun(@(x) 1:x, dim,'UniformOutput', 0);
[mx, my, mz] = meshgrid(range{:});

% Calculate sum( vx_sh(1) * (x-1), vx_sh(2) * (y-1), vx_sh(3) * (z-1))
% Sum of spatial index and voxel shift
tmp1 = [(mx(:)-1) (my(:)-1) (mz(:)-1)].*vox_shift;

% Multiply by phase change
phase_shift = reshape(sum(tmp1,2),dim).*gen_phase_rad(1);


% Repeat for data size ---- %

% Permute vector:
% Move the spatial index to the original dim index dimensions
perm_vec = 1:numel(sz); dim_pos = zeros(size(perm_vec)); dim_pos(doi) = 1;
perm_vec(doi) = 1:numel(doi);
perm_vec(~dim_pos) = numel(doi)+1:numel(sz);

% Permute
phase_shift = permute(phase_shift, perm_vec);

% Repeat the matrix at other indexes
rep_vec = sz; rep_vec(doi) = 1;
phase_shift = repmat(phase_shift,rep_vec);


% Complex --------- %

% Create complex array
phase_shift_complex = complex(cos(phase_shift), sin(phase_shift));

% Multiply -------- %
data = data.*phase_shift_complex;


