function nfo = twix_extractNFO(twix_hdr)
% Extract all required NFO from the twix-header for different Siemens
% platform versions: nucleus, dwelltime, fieldstrength, frequency, OS
% factor, dimensions, FOV, 
%
% Uses getFieldValues.m, findFieldnames.m

nfo = struct;

% Nucleus
vals = getFieldValues(twix_hdr, 'nucleus'); 
nfo.nucleus = getCommonValue(vals, 0);

% Dwelltime
vals = getFieldValues(twix_hdr, 'dwelltime'); 
nfo.dwelltime = getCommonValue(vals) * 1e-9;

% Fieldstrength
nfo.tesla = NaN;
stroi = {'fieldstrength', 'flNominalB0'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}); 
    nfo.tesla = getCommonValue(vals);
    if ~isnan(nfo.tesla), break; end
end

% Imaging Frequency
vals = getFieldValues(twix_hdr, 'frequency'); 
nfo.trans = getCommonValue(vals);

% OS Factors
nfo.OS = NaN;
stroi = {'ReadOSFactor', 'ReadoutOSFactor', 'OS'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}); 
    nfo.OS = getCommonValue(vals);
    if ~isnan(nfo.OS), break; end
end

% --- Geometric parameters --- %

% Matrix resolutions)
stroi = {'ReadResolution', 'MatrixSizeRead', 'Read'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}, 0); 
    nfo.dim(1) = getCommonValue(vals);
    if ~isnan(nfo.dim(1)), break; end
end

stroi = {'PhaseEncodingLines', 'MatrixSizePhase', 'Phase'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}, 0); 
    nfo.dim(2) = getCommonValue(vals);
    if ~isnan(nfo.dim(2)), break; end
end

% Dimensions: [AP LR FH]
stroi = {'SliceResolution', 'MatrixSizeSlice', 'Slice'};
stroi_exact = [1 0 0];
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}, 0, stroi_exact(kk)); 
    nfo.dim(3) = getCommonValue(vals);
    if ~isnan(nfo.dim(3)), break; end
end

% FOV
vals = getFieldValues(twix_hdr, 'ReadoutFOV'); 
nfo.fov(1) = getCommonValue(vals);
vals = getFieldValues(twix_hdr, 'PhaseFOV'); 
nfo.fov(2) = getCommonValue(vals);
vals = getFieldValues(twix_hdr, 'SliceThickness'); 
nfo.fov(3) = getCommonValue(vals);

% Orientation
stroi  = {'VoI_Normal_Tra', 'VoI_Normal_Sag', 'VoI_Normal_Cor'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}); 
    nfo.image_orientation = getCommonValue(vals);
    if ~isnan(nfo.image_orientation)
        nfo.image_orientation = stroi{kk}(end-2:end);
        break;
    end
end

% Also located in Config.VoI_Position_Cor/Sag/Tra
% Sagital means slices in LR dir
stroi = {'SliceResolution', 'MatrixSizeSlice', 'Slice'};
stroi_exact = [1 0 0];
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}, 0, stroi_exact(kk)); 
    nfo.dim(3) = getCommonValue(vals);
    if ~isnan(nfo.dim(3)), break; end
end

% Order of stroi matters.
% Other fields of interest: VoI_Position_Cor/Sag/Tra
stroi = {'VoiPositionCor', 'VoiPositionSag', 'VoiPositionTra'};
for kk = 1:numel(stroi)
    vals = getFieldValues(twix_hdr, stroi{kk}); 
    nfo.offcenter(kk) = getCommonValue(vals);
    if isnan(nfo.offcenter(kk)), nfo.offcenter(kk) = 0; end
end



end % End of main function

function vals = getCommonValue(vals, type)
% Given a cell-array of values/content, find the most frequent value and
% return this. Does NOT work with arrays - only single values or strings.
%
% type describes if query is a value (1, default) or a string (0)

if nargin == 1, type = 1; end

% Safety variable-type checks: empty, character, cell and array-size
ind_empty = cellfun(@isempty, vals); % Removing empty entries
ind_cell = cellfun(@iscell, vals);   % Removing cell-entries
    
if type % If a value-query
    
    ind_char = cellfun(@ischar, vals); 
    ind_big = ~cell2mat(cellfun(@(x) sum(size(x)) == 2, vals, 'uniform', 0));  
    ind_str = cellfun(@isstruct, vals);
    
    % Final index of interest
    ind = ones(size(vals)) - ind_char - ind_empty - ind_cell - ind_big - ind_str;        
    ind(ind < 0) = 0;  ind = logical(ind); 

    % Extract values within set properties (value vs. string)
    vals = cellfun(@(x) double(x), vals(ind));
    
    
    % Get most common value in array.
    if ~isempty(vals), vals = mode(vals); else, vals = NaN; end
    
else % If value is a string-query
    
    % Final index of interest
    ind = ones(size(vals)) - ind_empty - ind_cell;        
    ind(ind < 0) = 0;  ind = logical(ind); 

    % Extract values within set properties (value vs. string)
    vals = vals(ind);

    % Get most common value in array if not a single unique value.        
    if size(unique(vals),2) > 1
        vals = string(mode(categorical( vals)));
    else
        vals = cell2mat(unique(vals));
    end
        
end

end