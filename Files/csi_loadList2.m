%% Description:       Read .list-file to return header content and data
% Creator: Ir. Q. van Houtum       Version: 1.4          Date: 2017-06
% --------------------------------------------------------------------
%
% Reads the raw list-file and store header parameters in structure
% "list". Voxel data from the list-file is stored in field "list.data"
% as a structure. See below for list.data field-content.
%
% BE AWARE: Reading and analysing list-file is a slow process (minute) 
% for data-sets larger than 8MB.
%
% Compatible with raw CSI data (Philips). Exported raw format 
% *.data & *.list.
%
% List structure output:
% I. Data: contains the data array of the list-file
%          i. Char         - cell arrray;  
%         ii. Array        - data as double excl. typ-column
%         iv. Column names - name for each column 
%          v. Row type     - type for each row (STD, NOI, etc.)
% II+*.* : all parameter-lines found in list-file.
%
% Update notes:
% v1.4 - Clean up and output ordering.
%
% Contact: qhoutum2@umcutrecht.nl

function [list] = csi_loadList2(fpn_list)
tic;

% Process input.
if nargin == 0 
    [fn,fp, fi] = uigetfile('*.list','Select .list-file.');
    if fi == 0, return; end; fn = fn(1:end-5);
else
    [fp, fn] = fileparts(fpn_list); fp = [fp '\'];
end


% Open file ID 
fid_r = fopen([fp fn '.list'], 'r'); 
if fid_r == -1
    warning('csi_loadList2:Wrong_Path',...
            'Unable to load list-file. File not found.'); 
    list = []; return; 
end


%% Loop each line in list-file
% Number of lines are unknown - therefor rough memory allocation.
% Create storage variables for list-data (ldata) and count iterations 
% for fdata indexing.
ldata = cell(1,1); rtype = cell(1,1); list = struct; n = 1; 
while ~feof(fid_r)                                                           
    % Get temporary line of file.
    tline = fgetl(fid_r); tline_split = strsplit(tline,' ');
    
    % DATA    
    if ~strcmp(tline_split(1), '#') && ~strcmp(tline_split(1), '.')         % If line starts with hash (35) or dot (46) else data line.
        % Store data-line and typ column.
        if length(tline_split)<2
            warning(['MATLAB:csi_loadList2:corrupt list-file. ' ,...
                     'Could not properly read lines.']);
            list = NaN;
            return;
        end
        ldata{n,1} = tline_split(2:end); 
        rtype{n,1} = tline_split{2};
        n = n+1;

    % HEADER
    elseif strcmp(tline_split(1), '.')
        % Every header-line has parameters stored in format as below.
        hdr    = textscan(tline,'. %f %f %f %s : %f %f', 1);
        % For older Matlab use genvarname(hdr{4}{:});
        try
            list.(matlab.lang.makeValidName(hdr{4}{:})) = [hdr{5} hdr{6}]; 
        catch
            list.(strrep(genvarname(hdr{4}{:}),'0x2D','_')) ...
                = [hdr{5} hdr{6}];
        end
        
    % COLUMN-NAMES
    else
        if size(tline_split,2)>6 && strcmp(tline_split{2},'typ') 
            data_col_names = tline_split(2:end); 
        end
    end % end of if-not-#/.
end % end of file-while loop


%% Convert data-type

% Container and convert each column excluding the first.
try
    ldata_dbl = NaN(size(ldata,1), size(data_col_names,2)-1);
    for li = 1:size(ldata,1)
        ldata_dbl(li,:) = str2double(ldata{li}(2:end)); 
    end
catch, warning('MATLAB:csi_loadList:char_2_double_error', ...
        'Could not convert list-data to double');
end


%% Set output fieldorder

% Get fieldnames
fieldn = fieldnames(list); 

% Add data: char, double, column names, row type and file-info.
list.data.array = ldata_dbl; list.data.char = ldata; 
list.data.col_names = data_col_names; list.data.row_type = rtype; 
list.filename = fn; list.filepath = fp; 

% Order
list = orderfields(list,{'data','filename','filepath',fieldn{:}});

% Display info.
t1 = toc; fprintf('Processed .list-file. Elapsed: %4.3fs \n', t1);

% </end>.


