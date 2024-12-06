function csi_writeSPAR(fpn, data_size, data_dim_labels)
%%%% Description:                                   Create new SPAR-file.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-07
%%% --------------------------------------------------------------------
%%% Called by csi_writeSDAT
%%%
%%% Writes SPAR-file for new SDAT-file. Not intended for stand-alone-use.
%%% Manually edit other parameters such as bandwidth/sample_frequency and 
%%% others in a text-editor or use the function csi_editSPAR.
%%%
%%% fpn: path
%%% data_size: size of data-array written to SDAT-file
%%% data_dim_labels: label for each index-dimension in data_size.
%%%
%%% Contact: qhoutum2@umcutrecht.nl

% JMRUI     - Required parameters in SPAR: Row, spec_num_row, spec_num_col.
% MATLAB    - Required parameters in SPAR: dimN_pnts.

%% 1. Load template SPAR
if exist('csi_writeSPAR_template.SPAR','file') 
    % Open file ID
    fid = fopen('csi_writeSPAR_template.SPAR','r'); 
    if fid == -1, warning('SPAR template missing!'); return; end
    % Loop each line in SPAR-file
    spar = cell(208,1); n = 1; 
    while ~feof(fid), spar{n,1} = fgetl(fid); n=n+1; end
    % Close file ID
    fclose(fid);
else
   fprintf(['csi_writeSPAR_template.SPAR not found.'...
            ' Unable to write SPARfile. Returning.\n']); 
end



%% 2. All required parameters required to match to be able to properly 
%     write an SPAR and SDAT file.
num_dimensions = numel(data_size); dims = data_size;



%% 3. Create labels for each index-dimension (if not given).
if nargin == 2
    data_dim_labels    = cell(1,num_dimensions);
    data_dim_labels{1} = 'sec'; data_dim_labels(2:end) = {'num'};
end



%% 4. Update SPAR parameters
% I.   Add the number of dimensions of data-array in SDAT-file. 
% II.  Add parameters for propper SDAT-file use in Matlab 
% III. Add parameters for propper SDAT-file use in JMRUI.

% I. Num_dimensions: global.
num_dim_line = find(cellfun(@isempty,strfind(spar,'num_dimensions')) == 0);
spar(num_dim_line) = {['num_dimensions : ' num2str(num_dimensions)]};
                    
                    
% II. MATLAB: dimN_pnts
for kk = 1:num_dimensions
    dimstr = sprintf('dim%i_', kk); % Get dimN string.
    
    % Create important lines - contain dimension size and label (dir/ext)
    tmp_paragraph = { [dimstr 'ext :['       data_dim_labels{kk} ']' ];...
                      [dimstr 'pnts : '      num2str(dims(kk))       ];...
                      [dimstr 'low_val : '   num2str(1.0)            ];...
                      [dimstr 'step : '      num2str(1.0)            ];...
                      [dimstr 'direction : ' data_dim_labels{kk}     ];...
                      [dimstr 't0_point : -'                         ];...
                     };
    
    % Indices of lines dimN_.
    line_index = find(cellfun(@isempty,strfind(spar,dimstr)) == 0);
    
    if ~isempty(line_index)
        % DimN exists in template - replace lines in SPAR.
        spar(line_index) = tmp_paragraph(:);
    else
        % Add new dimN paragraph to SPAR:
        % Layout: break line, tmp_paragraph with empty str/line, empty str;
        
        % A. Get previous dimN line-index
        line_index_prev = ...
        find(cellfun(@isempty,strfind(spar,sprintf('dim%i_', kk-1))) == 0);
        % B. Split spar after dimN, get new line index for tmp_paragraph in
        %    part A.
        line_index_AB = line_index_prev(end) + 1;  % split part A/B here.
        line_index_new= line_index_prev      + 14; % adds 14 lines.
        
        % C. Part A and B of SPAR
        partA = spar(1:line_index_AB); partB = spar(line_index_AB+1:end);
        
        % D. Add paragraph to new-line-indices
        partA(line_index_new) = tmp_paragraph(:);
        % E. Add break-line
        partA{line_index_new(1)-2} = ['!' repmat('-',1,53)];
        % F. Insert empty strings in empty cells
        partA([line_index_new-1; line_index_new(end)+1]) = ...
            repmat({''},1,size(line_index_new,1)+1); 
        
        % G. Replace SPAR - combine parts A and B
        spar = cat(1,partA,partB);
    end
end

% III. JMRUI: spec_num_row/col.

% Row-parameter SPAR
line_spec_num_row = cellfun(@isempty,strfind(spar,'spec_num_row')) == 0;
% Col-parameter SPAR
line_spec_num_col = cellfun(@isempty,strfind(spar,'spec_num_col')) == 0;

% ROW: Nr of FIDS (second dim).
spar(line_spec_num_row) = {['spec_num_row : ' num2str(dims(2))]};
% COL: Nr of Samples (first dim).
spar(line_spec_num_col) = {['spec_num_col : ' num2str(dims(1))]}; 


%% 5. Write SPAR
fid = fopen([fpn '.SPAR'],'w');
for kk = 1:size(spar,1), fprintf(fid, '%s\r\n', spar{kk}); end; fclose(fid); 
fprintf('SPAR-file written: %s\n', [fpn '.SPAR']);

