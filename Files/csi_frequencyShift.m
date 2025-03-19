function fid = csi_frequencyShift(fid, dt)
% Shift the frequency of each FID to correct linear phasing
%
% Input
% fid - Free induction decay of MRS data set. e.g. time domain.
% dt  - Time correction for delay in FID acquisition (-8 for pi?)

% Data details
sz = size(fid);

% Subtract linear phase dependent on dt
lin_phase = ((0:sz(1)-1).*(2.*pi.*dt) ./1000)';
% Phased
lin_phase = complex( cos(lin_phase), sin(lin_phase) );

% Convert to cell
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
fid = mat2cell(fid, sz(1), cell_layout{:}); 

% Division cell array
% div_array = repmat(lin_phase,size(fid));
% fid = fid./div_array;
div_array = repmat({lin_phase}, [1 sz(2:end)]);

% Apply operation to each cell e.g. FID
fid = cellfun(@(x,y) x./y, fid, div_array,  'UniformOutput', 0);                     

% Create array
fid = cell2mat(fid);
% Correct if 1D|2D
if numel(sz) <= 2, fid = squeeze(fid); end
