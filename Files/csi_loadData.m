%% Description:               Load CSI .data-file into memory and index
% Creator: Ir. Q. van Houtum       Version: 1.5          Date: 2017-06
% --------------------------------------------------------------------
%
% Reads the raw CSI-data data-file and store all data into an array. Data
% from .list file is used to index each voxel.
% 
% Input:          none      (-)      - Opens up .data-file selection GUI
%                 path      (string) - Path to .list/.data-files
%                 list-data (struct) - Data structure of list-file
%                                      as returned from csi_loadList.
%
% NB: Untouched data. No operations are performed on the indexed
%     data itself.
%     
% Data output:
% 1. data.filepath        - path of data-file
%        .filename        - name of data-file
%        .noise           - noise entries from data-file
%        .raw             - raw complex values (excl. noise). Split/FID:
%                           reshape(raw,[size(raw,2)/#FIDsz FIDsz]);
%        .labels          - label for each dimensions after indexing.
%
%      ### ! INDEXING SUCCESFULL: data.raw == indexed data. !  ###
%      Reduces data size of output-variable. 
%
% Known issue: list-struct from csi_loadData not useable to load data
%              again....Due splitting of noise in list.data.
%
% Contact: qhoutum2@umcutrecht.nl

function [data, list] = csi_loadData(listpath, varargin)


tic;

switch nargin
    case 0 
        [fn, fp] = uigetfile('*.data','Select .data-file.');
        if fn == 0, return; end; fn = fn(1:end-5); fp = [fp '\'];
        save_on = 0; 
    case 1
        if ischar(listpath), [fp, fn] = fileparts(listpath); fp = [fp '\'];
        elseif isstruct(listpath), list = listpath;
        end
        save_on = 0; 
    case 2
        if ischar(listpath), [fp, fn] = fileparts(listpath); fp = [fp '\'];
        elseif isstruct(listpath), list = listpath;
        end
        save_on = varargin{1};
end

%% READ LIST-file
if ~exist('list','var')
    list = csi_loadList2([fp fn]); 
    if isempty(list)
        warning('Loading of list file failed.'); data = []; list = []; 
        return;
    end
end
% Add filepath.
if ~isfield(list, 'filepath'), list.filename = fn; list.filepath = fp; end


%% Read DATA-file

% Open data-file identifier, load data, split real/imag part and 
% close file ID. 
try fid_r = fopen([list.filepath list.filename '.data'], 'r');
    data_tmp = reshape(fread(fid_r, 'float32', 'ieee-le'),2,[]); 
    fclose(fid_r);
catch 
    data = []; list = [];
    warning('Unable to load data-file. File not found.');
    return;
end


%% EXTRACT noise from LIST-file and DATA-file
% Find noise entries in list-data and seperate from actual list.data.
% Do the same for the data-file data. Data.raw and data.noise are created.

% LIST.DATA.ROW-TYPE: filter noise component
nind = find(strcmp(list.data.row_type,'NOI') == 1);

% LIST.DATA.CELL-char backup and delete noise from list.data
list.data.noise.char = list.data.char(nind); list.data.char(nind) = [];
% LIST.DATA.ARR-double backup and delete noise from list.data
list.data.noise.double = list.data.array(nind,:); 
list.data.array(nind,:)= [];

% Data structure;
data = struct; data.filepath = list.filepath; data.filename = list.filename;
% Data start point in raw-data: time-resolution plus number of channels 
% required. Excludes noise.
ndimf      = list.F_resolution; nchan = size(nind,1);
data_start = nchan*ndimf+1; 

% -- NOISE: Get noise FIDs: stored as first nchan-lines in data-file.
noise_tmp  = complex(data_tmp(1,1:data_start-1),data_tmp(2,1:data_start-1));    
data.noise = reshape(noise_tmp, [ndimf nchan]);

% -- DATA: Get values of FID by excluding noise.
data.raw       = complex(data_tmp(1,data_start:size(data_tmp,2)),...
                         data_tmp(2,data_start:size(data_tmp,2)));    
clear('data_tmp');

    %%% NOISE and DATA are seperated in the list-and data-structs %%%
    

%% ANALYSE indexing-values in LIST-file  
% Get correct indexing values from list.data table to index data.raw.
% Search correct columns for indexing.

try % SEE END OF FILE!

% COL-NAME: Get max. index per dimensions using column names. .
% Noise excluded. Labels of interest for indexing:
loi = {'chan', 'kx', 'ky', 'kz', 'aver', 'mix', 'dyn', 'card',...
       'echo', 'loca', 'sign', 'extr1', 'extr2','rf','grad','enc','rtop',...
       'rr'};
% If more are required, just append the list here! Keep "chan" at 1! :)   

% LOI_IND:      Index of loi-columns in list-data.double;
% LOI_IND_NEG:  Set to one if this loi-column has negative index values
% LOI_MAX:      Maximum # of unique values in the loi-column in list-data;
loi_ind = zeros(size(loi,2),1); loi_max = loi_ind; loi_ind_neg = loi_ind; 

% Loop each label of interest (loi) to find the column-index of this label in
% list.data.array, its maximum nr of values and possible negative indexing.
% QH - 201906: edit for non-consecutive channel ID e.g. [54:69 16] is the
% ID for the 31P RX 16ch array and body coil for receive.
for loii = 1:size(loi,2)
    % Get data-index (column) of the label in list.data.array, index 
    % minus one due exclusion of type-column in list.data.array!
    loi_ind(loii) = find(strcmp(loi{loii}, list.data.col_names) == 1)-1;    
    
    % Get max # of unique indices in this loi-column in list.data.array
    % Matlab does not allow zero-indexing, therefor numel(unique).
    loi_max(loii) = numel(unique(list.data.array(:,loi_ind(loii))));
   
    % Matlab does not allow negative indexing: boolean of columns which
    % require special calculations.
    if sum(list.data.array(:,loi_ind(loii)) < 0 ) > 0 
        loi_ind_neg(loii)  = 1; 
    end
    
    % Correct for non-consecutive numbers.
    if (max(list.data.array(:,loi_ind(loii))) > loi_max(loii))
        fprintf('Corrected for high-value indexing: %s\n', loi{loii})
               
        % Create temp backup of the column of interest
        tmp = list.data.array(:,loi_ind(loii));
        % Loop each unique value and increment the new indexing ID 
        % automatically.
        uni = unique(list.data.array(:,loi_ind(loii))); % All unique values
        for ui = 1:size(uni,1)
            tmp(list.data.array(:,loi_ind(loii)) == uni(ui)) = (ui-1);
        end
        list.data.array(:,loi_ind(loii)) = tmp; clear('tmp');
    end
    
end

% Include only index-dimensions larger than one.
loi_max_one = loi_max == 1; % update all to this boolean indexing.

% Update loi_max and loi_ind_neg to loi_max_one:
loi_ind(loi_max_one~=0) = []; loi_ind_neg(loi_max_one~=0) = [];
% Update LOI, loi_max to loi_max_one:
loi_max(loi_max_one ~= 0) = []; 
% Update loi_ind_neg from boolean to actual index numbers in LOI
loi_ind_neg = find(loi_ind_neg == 1); % IN LOI

% Add labels used for indexing to output data-struct
data.labels = loi(loi_max_one==0);
% Add first dimension to labels
data.labels = {'sec', data.labels{:}};

% Set normal and negative indexing-column of loi
loi_ind_norm = 1:size(loi_ind,1); loi_ind_norm(loi_ind_neg) = [];
% loi_ind_norm need correction for zero indexing values.
% loi_ind_neg  need correction for both zero and negative indexing values.


%% EXTRACT indexing-values in LIST-data

% List-data: The required indexing values data per FID 
ldat = list.data.array(:,loi_ind);

% Special, correct for negative indexing
if ~isempty(loi_ind_neg)
    ldat(:,loi_ind_neg) = ldat(:,loi_ind_neg) + ...
       repmat((floor(loi_max(loi_ind_neg)/2) + 1)',size(ldat,1),1);
end

% Normal, correct for non-zero indexing Matlab only
if ~isempty(loi_ind_norm), ldat(:,loi_ind_norm) = ldat(:,loi_ind_norm) + 1;   
end


%% CONVERT subscript to linear indexing

% Update max value per index dimension with time-index
% size(data) equals ldat_max_dim
ldat_max_dim_db  = [ndimf loi_max'];                %double
ldat_max_dim_cl  = num2cell(ldat_max_dim_db,1);     %cell

% Create container for the total array 
data.indexed = complex(zeros(ldat_max_dim_cl{:}));

% Repeate every line in ldat ndimf time below itself!
t = repmat(1:size(ldat,1),ndimf,1); t = t(:)'; tmp = ldat(t,:);
% Vertically add time-index 1:ndimf.
data.sub_index_raw = [repmat((1:ndimf)' ,size(ldat,1),1) tmp]; %dbl
clear('tmp'); % Clear memory.
% Convert to cell.
data.sub_index_raw = num2cell(data.sub_index_raw',2);   	   %cell

% Linear index for every data-point
data.lin_index_raw = sub2ind(size(data.indexed),data.sub_index_raw{:});

%% INDEX DATA.RAW

% Store raw dat using the linear indices
data.indexed(data.lin_index_raw) = data.raw(:);
% </Done>.


%% Display info

% Display information of the loaded data and list file.
fprintf('#FIDS: %i | #Vector: %i  \n', size(ldat,1), ndimf);
fprintf('Dimensions loaded: \n'); fprintf('| %4s ', loi{loi_max_one==0}); 
fprintf('\nMax size: \n'); fprintf('| %4i ', loi_max); fprintf('\n');

% Display run-time
t1 = toc; 
fprintf('Processed data-file using list-file. Elapsed: %4.3fs \n', t1);



%% Overwrite raw to indexed
% And remove linear and sub index arrays.
rem_fields_str = {'lin_index_raw','sub_index_raw', 'indexed'};
data.raw = data.indexed; data = rmfield(data, rem_fields_str);



%% Save indexed data
if save_on ==1
    try
        save([fp 'indexed_' fn  '.mat'], 'data', '-v7.3');
        fprintf('Indexed data saved to: \n    path: %s\n    name: %s\n',...
                fp, ['indexed_' fn  '.mat']);
    catch
    end
end


% Catch if anything goes wrong during indexing. Store error. Output already
% set in code outside of try-block: data-struct and list-struct. 
catch err
    save(['mrs_loadData_errorlog_' datestr(rem(now,1),'HHMMSS') '.mat'],'err');
    warning('Raw data loaded. Indexing failed. See log.');
end

% </end>
    
