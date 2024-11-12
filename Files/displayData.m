function displayData(data, csi, tag, labels, opts)
% Display data as table, map or histogram. Data includes a data-matrix in
% combination with an image array - aimed at data-maps.
%
% This function is created for speed of plotting as the original function
% uses axes per voxel to plot the maps, creating very large
% graphics-objects. It causes a bottleneck in Matlabs processing regarding
% plot updates and thus some calculations.
%
% data   = [1 x X x Y x Z x ...]
% csi    = struct from CSIgui and includes the field csi.conv, with 
%          converted images to plot in the maps.
% labels = either custom labels of the data-dimensions or csi.data.labels.
% opts   = opts.doMask and opts.dataDisp, skipping user input to set output
%          of interest (maps, tables, histogram).


if nargin < 5, opts = struct; end
if nargin < 4, labels = csi.data.labels; end

% Add color-scheme
clr.main = [0 0 0]; clr.text_title = [0.8 0.8 0.8];
clr.hilight1 = [0.8 0 0]; clr.lines1 = [0 0 1]; 
clr.grid = [0.4 0.4 0.4]; csi.clr = clr;

% Permute data to correct format if necessary [ 1 x X x Y x Z ];
if size(data,1) ~= 1 
    sz = size(data); data = permute(data, [numel(sz)+1 1:numel(sz)]); 
end

% Initiate data-plot
% data: a data array of size [1 x X x Y x Z].
% data-tag: tag used to display in figures and input-queries.
% csi:  csi-structure from CSIgui including conv-field for image-plot.
CSI_dataAs_Initiate(data, csi, tag, labels, opts);

% --- Executes by map-scripts to start any visualiziation of data
function CSI_dataAs_Initiate(data, csi, data_tag, labels, opts)
% After calculating some maps or anything 3D, and one wants to display it.
% Call this function. It will ask the user the proper info; display type,
% filter by SNR and more.

if nargin < 4, labels = csi.data.labels; end
if nargin < 5, opts = struct; end

if ~isfield(opts, 'dataDisp')
    ename = repmat({'popup'}, 1, 2);
    qry = {'Display type: ','Apply excluding voxel mask:'};
    inp = {{'Map','Table','Histogram'}, {'Yes', 'No'}};

    % Display type from user
    uans = getInput(ename, qry, inp, data_tag);
    if isempty(uans)
        CSI_Log({sprintf('%s mapping skipped.', data_tag)},{''}); 
        return; 
    end
    dataDisp = lower(uans{1}); % Display method
    switch uans{2}, case 'Yes', doMask = 1; otherwise, doMask = 0; end
else
    dataDisp = opts.dataDisp;
    doMask = opts.doMask;
end

% If voxel mask is requested.
maskMsg = 'none';
if doMask
    if ~isfield(csi.data, 'voxelmask')
        csi = VoxelMask_Initiate(csi);        
    end
    
    % User can quit/cancel voxel-mask creation - this catches that error
    % and continues without creating masked-data
    if isfield(csi.data, 'voxelmask')
        mask = csi.data.voxelmask;  
    
        % Data-size - is expected to have spatial dimensions on ind(2:4);
        dsz = size(data); msz = size(mask);     
    
        % Retrieve mask for current size: if non-spatial dimensions are 
        % used to calculate a parameter, that index is not present in the 
        % data-volume.
        if numel(msz) ~= numel(dsz)
            dim = size(data);
            cell_ind = arrayfun(@(x) 1:x, dim, 'UniformOutput', false);  
            mask = mask(cell_ind{:});
        end
    
        % Apply mask to data
        data(mask) = NaN;
        maskMsg = 'applied';
    end
end

% Show statistics nfo
stats = csi_statistics_of_volume(data); stats.mask = maskMsg; 
stats.source = 'full-data-volume'; stats.nfo = data_tag; 
Statistics_Viewer(stats);


% Switch to data-display type.
switch dataDisp
    case 'table'      % Table % ----- %    
        % Send to tableData function, to show each slice as a table.
        % If higher dimensions are available, the data will be 
        % concentonated per slice. 
        CSI_dataAsTable(data, data_tag)
    
    case 'map'
        % Send data to dataAsTabs and create maps in a tabbed figure for 
        % all slices (or only the current slice plotted).
        
        % Prepare dataAsTabs input
        [data, color_scale, sloi] = CSI_dataAsTabs_Prepare(data);
        if isnan(data), return; end
        % Plot as tab
        CSI_dataAsTabs(csi, data, data_tag, labels, color_scale, sloi);
    
        % Show statistics nfo
        stats = csi_statistics_of_volume(data); stats.mask = maskMsg;
        stats.source = 'displayed-data'; stats.nfo = lower(data_tag);
        Statistics_Viewer(stats);
        
    case 'histogram'
        % Plot data as histogram - simple display.    
        fig = figure(); ax = axes(fig);
        mx = max(data(:), [], 'omitnan'); md = median(data(:),'omitnan');  
        tmp = data(:);
        histogram(ax, tmp(~isnan(tmp)), round(100+(mx./md)) );
        title([data_tag ' Histogram']); xlabel([data_tag ' Bins']);
end

% --- Executes by functions to display data as a table
function CSI_dataAsTable(data, datatag)
% Given data, data will be shown as a table per slice. Each higher index
% dimension in data will be concatonated below its origin slice, and
% seperated by an empty line including its original index of slice and
% higher dimensions.
%
% data,        (1,X,Y,Z,etc); Z  = slice. etc are higher dimensions.
% datatag,     Data origin string;

% Process input.
if nargin == 1, datatag = 'data'; end


% FIGURE: prepare % ----------------------------------------------- %
% Table figure
tfh = figure('Tag', ['CSIgui_table' datatag],'Name',...
            ['CSIgui - ' datatag],...
            'Color', 'Black','Toolbar', 'None', 'MenuBar', 'None',...
            'NumberTitle', 'Off', 'Resize','On'); 
pos = get(tfh,'Position');
tgui = guidata(tfh); tgui.fig = tfh;

% Tabgroup for all slices
tgui.tabgp = uitabgroup(tgui.fig,'Unit','normalized','Position',[0 0 1 1]);


[nZeros, nZeros_acc] = numzeros(median(data(:), 'omitnan'));

% Loop each slice - expected at dimension 4
for sl = 1:size(data,4)

    % PREP TAB % ------------------------------------------------------- %
    
    % Add tab-handle according nr of thumbnails!
    tgui.tab{sl}.tabh = uitab(tgui.tabgp ,...
          'Title',['Slice ' num2str(sl)],...
          'BackGroundColor','Black','ForegroundColor', 'Black');    

    % Table object
    tgui.tab{sl}.table = uitable(tgui.tab{sl}.tabh,...
        'Position', [30 20 pos(3:4)-60],...
        'Units','Normalized'); 
    
    % ADD BUTTONS HERE: To save all data
    tgui.tab{sl}.savebutton = uicontrol(tgui.tab{sl}.tabh,...
        'Position', [30 2.5 100 15],'String','Save all',...
        'BackgroundColor','Black','ForegroundColor',[0.94 0.94 0.94],...
        'Callback', @CSI_dataAsTable_SaveButton);
    
    % To save selected data
    tgui.tab{sl}.savebutton_sel = uicontrol(tgui.tab{sl}.tabh,...
        'Position', [132.5 2.5 100 15],'String','Save selected',...
        'BackgroundColor', 'Black', 'ForegroundColor', [0.94 0.94 0.94],...
        'Callback', @CSI_dataAsTable_SaveButton);
    
    % To copy to clipboard
    tgui.tab{sl}.savebutton_clip = uicontrol(tgui.tab{sl}.tabh,...
        'Position', [235 2.5 100 15],'String','Copy to Clipboard',...
        'BackgroundColor', 'Black', 'ForegroundColor', [0.94 0.94 0.94],...
        'Callback', @CSI_dataAsTable_SaveButton);
    
    % DATA: prepare % -------------------------------------------------- %

    % Slice-selection and higher dimension selection.
    % Size of data to a size vector cell
    sz = size(data); 
    if numel(sz) >=5, nDimC = num2cell(sz(5:end));
    else,             nDimC = {1};
    end
    % To linear vector per cell.
    nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

    % Slice data + other dimensions
    sliceData = data(:,:,:,sl, nDimC{:});

    % Convert to a 2D sliced array: index 2 and 3 are considered the 2D 
    % plane of the original ND data array.
    [tableData_num, tableData_Index] = slices2rowArray(sliceData,1);
    
    

    % Calculate NaN row index
    NaNstep   = size(data,3)+1; dSz = size(sliceData); dSz = dSz(5:end);
    nanrow = 1: NaNstep : (NaNstep * prod(dSz));

    % DATA: String % -------------------------------------------------- %
    
    % Create proper format to vieuw the value
    if nZeros <= 0 || nZeros == nZeros_acc        
        % Convert double-array to cell-array and then to string
        tableData_str = cellfun(@sprintf, repmat({'%3.3f'},...
            [size(tableData_num,1) size(tableData_num,2)]), ...
            num2cell(tableData_num), 'UniformOutput',0);
    else 
        prefix = '0'; if nZeros >= 9, prefix = ''; end
        formatstr = repmat({'%3.3fe-%s%i'},...
                    [size(tableData_num,1) size(tableData_num,2)]);
        prefix = repmat({prefix}, ...
                    [size(tableData_num,1) size(tableData_num,2)]);
        tableData_num = num2cell(tableData_num .* (10.^(nZeros+1)));
        nZeros_cell = repmat({nZeros},...
                    [size(tableData_num,1) size(tableData_num,2)]);

        tableData_str = cellfun(@sprintf, formatstr, tableData_num,...
            prefix, nZeros_cell, 'UniformOutput',0);
        isnan_ind = cellfun(@(x) isnan(x), tableData_num,'uniform', 0);
        tableData_str(cell2mat(isnan_ind)) = {'NaN'};
    end



    % Replace NaN with slice and higher index dimension information
    for ri = 1:length(nanrow)
        tableData_str{nanrow(ri),1} = ...
            ['Ind: ' sprintf('%i %i ', sl, tableData_Index(ri,2:end))];
        tableData_str(nanrow(ri),2:end) = {''}; % Remove other NaNs
    end

    % DATA: to Table % ------------------------------------------------ %
    tgui.tab{sl}.table.Data = tableData_str;
    
    % Add selection callback
    set(tgui.tab{sl}.table,...
        'CellSelectionCallback',@CSI_dataTable_CellSelection);

    % DEV: Is each first slice equal to first 2D array part of the t data?
    % squeeze(data(:,:,:,sl,1,1,1,1)) == ...
    %                          tableData_num(2:sz(3)+1,:) 
    % DEV: Is each scnd slice equal to scnd 2D array part of the t data?
    % squeeze(data(:,:,:,sl,2,1,1,1)) == ...
    %                           tableData_num(sz(3)+3:((sz(3)+1)*2),:) 

end % End table/slice loop.

% Save GUI-data to figure.
guidata(tgui.fig, tgui);

%tgui.fig.ToolBar = 'figure';
toolbar_create(tgui.fig)

% --- Executes by CSI_tableData.
function [sliceArray, sliceIndex] = slices2rowArray(data, split)
% Convert every slice and its higher index dimensions to a 2D array with
% every slice concatenated in the row direction.
%
% Data = 
% (1 x X x Y x 1 x Higher Index) data.
% Thus higher dimensions index are 5 and above.
%
% Split = 
% 1 or 0; if 1, every 2D array will be split with a NaN row at the
% top of the slice e.g. row number (1+X)*(1:prod(higher dims)) will be NaN; 
% Default set to off (0). 
%
% E.g. a 2x2 slice with higher index dimensions 4, 5, 3; will result in a
% 2D array of 2 columns with 4*5*3 rows. Every 2 rows, shows one slice.
% Including NaN-split, the resulting array willbe (4*5*3) +  (4+5+3)
%
% Index per slice is calculated using dim2index.

if nargin == 1, split = 0; end

% 1. Size of data to a size vector cell
sz    = size(data);
if numel(sz) >=4, nDimC = num2cell(sz(4:end));
else,             nDimC = {1};
end
% Convert to linear vector cell per higher index dimension.
nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

% 2. Get the ordered indexing e.g. all possible combinations in increasing
% order.
sliceIndex = slice2rowIndex(nDimC);

% 3. Concatenate 
sliceArray = []; % Container.
for sli = 1:size(sliceIndex,1) % Loop each index representing a slice.
    % Create cell index and get this higher index dim 2D array
    ind = num2cell(sliceIndex(sli,:)); 
    % Why not squeeze?
    % Well we need the transpose:
    % To move all to 2D array, Matlab will transpose the other higher
    % dimensions. To correct for this, apply a transpose again.
    % This is because the x/y dimension in 2/3 index must remain ofc.
    slice = permute(data(:,:,:,ind{:}),[2:numel(size(data)) 1])';
    
    % Concatenate
    if split == 1
        sliceArray = cat(1,sliceArray,NaN(1,size(slice,2)),slice);
    else
        sliceArray = cat(1,sliceArray,slice);
    end
end

% --- Executes by slice2rowArray.
function indArray = slice2rowIndex(ndimCell)
% Input:
% A cell array (1,N) with a 1 to "index size" vector ...
% 
% Outout:
% ...will return an array with the corresponding indexing for all elements
% in the initial array.
%
% E.g. Array of size (2,3):
% ndumCell = {{1:2}, {1:3}} --> returns: [1 1; 1 2; 1 3; 2 1; 2 2; 2 3]

% Get nr of linear vector dimension index cells to process
dims_sz   = size(ndimCell);

% Get the max per cell.
dimMax = NaN(1,dims_sz(2));
for kk = 1:size(ndimCell,2), dimMax(kk) = max(ndimCell{kk}); end
nRows = prod(dimMax);

% Loop each column. Loop the product of previous ndimCell maxima. Loop the
% product of the remaining ndimCell maxima. 
indArray = NaN(nRows, dims_sz(2));
for col = 1:size(ndimCell,2)                        % Every column
    for di = 1:prod(dimMax(1:col-1))                % Repeat prod prev max
        colrep   = prod(dimMax(col+1:end));         % Repeat prod rem  max
        colstart = prod(dimMax(col:end))*(di-1);    % Start at rowindex
        for kk = 1:dimMax(col)                      % 1 to max this col
            rstrt = colstart + (((kk-1)*colrep)+1); % Start at row..
            rstop = colstart + (colrep*kk);         % End of row..
            indArray( rstrt:rstop, col) = kk;       % Fill these rows @ col
        end   
    end
end

% --- Executes by selecting data in tabels
function CSI_dataTable_CellSelection(hObj, evt)

% Table GUI handle
tgui = guidata(hObj.Parent);
% Cell of interest
coi = evt.Indices;

tgui.selectedCells = coi; guidata(hObj.Parent, tgui);

% Abort if more than 1 cell is selected
if size(coi,1) > 1, return; end

% Get current selected tab
tab_title = tgui.tabgp.SelectedTab.Title;
tab_title_space = strfind(tab_title,' ');
tab_ind = str2double(tab_title(tab_title_space+1:end));

% Find distance between each row-title (Ind:)
bool_ind = contains(tgui.tab{tab_ind}.table.Data,'Ind');
[r, ~] = find(bool_ind == 1);
if size(r,1) > 1, deltaNaN = r(2)-r(1)-1; else,  deltaNaN = 1; end

% Location of cell in actual data array before table conversion
row_slice = sum(r<coi(1));
coi(1) = coi(1) - ((row_slice-1) .* deltaNaN) - row_slice;
coi = fliplr(coi); % Revert for row - column / x - y inversion!

% Replace titlebar name with selected voxel index
str_vertLine_ind = strfind(tgui.fig.Name, '|');
if isempty(str_vertLine_ind)
    tgui.fig.Name = [tgui.fig.Name...
                               ' | Selected voxel [x, y]: ' int2str(coi)];
else
    tgui.fig.Name = [tgui.fig.Name(1:str_vertLine_ind-2)...
                               ' | Selected voxel [x, y]: ' int2str(coi)];
end

% --- Executed by save button press in table
function CSI_dataAsTable_SaveButton(hObj, evt)
% Save data shown in table.

% Get app-structs and parent of this function
tab_obj = hObj.Parent; tgui = guidata(tab_obj); 
fparent = evt.Source.String;



% Save according to function-parent e.g. button.
saveToFile = 1;
switch fparent
    case 'Save all'
        dsize = [size(tgui.tab{1}.table.Data) size(tgui.tab,2)];
        data = NaN(dsize); datastr = cell(dsize);
        % Get data as string


        % Get data as double and string
        for sli = 1:size(tgui.tab,2)
            data(:,:,sli) = ...
                cellfun(@str2double, tgui.tab{sli}.table.Data);
            datastr(:,:,sli) = tgui.tab{sli}.table.Data;
        end

        % Find all header-rows: starts with "Ind"
        ind = contains(datastr, 'Ind');
        ind = find(ind == 1);
        ind = ind2subQ(dsize,ind);

        % row/col/slice indexes
        row = unique(ind{1}); 

        % Number of row/col/slices
        extr = size(row,1); % Number of additional slices
        nrow = (size(data,1) - extr)/extr; % nRows in slice
        ncol = size(data,2); % #Col in slice
        nsli = size(data,3); % #Number of actual slices
        
        outp = NaN(nrow,ncol,nsli, extr);
        for si = 1:nsli
            % Take slice of interest
            sloi = data(:,:,si);
            % Remove header-rows
            ind_header = ind{1}(ind{3} == si);
            sloi(ind_header,:) = [];  
            sloi = reshape(sloi,[nrow*extr, ncol]);
            % Loop each extra-dim
            for xi = 1:extr
                outp(:,:,si,xi) = sloi( ((xi-1)*nrow)+1:xi*nrow,:);        
            end
        end
        
        data2exp = data; data = outp; ind = [];
        
    case 'Save selected'
         % Tab index
        tab_title = tgui.tabgp.SelectedTab.Title;
        tab_title_space = strfind(tab_title,' ');
        tab_ind = str2double(tab_title(tab_title_space+1:end));
        
        % Get data
        ind = tgui.selectedCells;

        data = cellfun(@str2double, tgui.tab{tab_ind}.table.Data);
        data2exp = NaN(size(ind,1),1);
        for kk = 1:size(ind,1)
            data2exp(kk,1) = data(ind(kk,1),ind(kk,2));
        end
    case 'Copy to Clipboard' 
        saveToFile = 0;

        % Tab index
        tab_title = tgui.tabgp.SelectedTab.Title;
        tab_title_space = strfind(tab_title,' ');
        tab_ind = str2double(tab_title(tab_title_space+1:end));

        % Get selected indexes
        ind = tgui.selectedCells;

        % Get data
        data = cellfun(@str2double, tgui.tab{tab_ind}.table.Data);

        % Get array-config of selected data
        sz_sel_mx = max(ind); sz_sel_mn = min(ind);
        sz_sel_mx = sz_sel_mx - (sz_sel_mn - 1);
        ind_aim = ind - (sz_sel_mn-1);

        % This could be a checksum: but user may select different shapes
        % size(ind,1) == prod(sz_sel_mx);

        data2exp = NaN(sz_sel_mx);
        for kk = 1:size(ind,1)
            data2exp(ind_aim(kk,1), ind_aim(kk,2)) = ...
                data(ind(kk,1), ind(kk,2));
        end

        % Copy it as text - under construction
        % data2exp = string(data2exp) % Creates string array
        % Copy clipboard requires string-vector!

        clipboard("copy", data2exp)
end

if saveToFile

% Get file destination from user.
[fn, fp, fi] = ...
    uiputfile({'*.txt', 'Text Files (*.txt)';...
               '*.mat', 'MATLAB File (*.mat)'},...
               'Save table data...');
if fi == 0, return; end
ext = fn(end-3:end);

% Save to file.
switch lower(ext)
    case '.mat'
        selected = data2exp; index = ind;
        save([fp fn],'selected','data','index'); 
    case '.txt'
        csi_writeText(data2exp,[fp fn]);
end

end

% --- Executes by CSI_dataAs_Initiate to prepare data for tabbed-maps
function [data, color_scale, sloi] = CSI_dataAsTabs_Prepare(data)
% This function will get userinput, prep the data and return required
% variables to plot data-array in a tabbed-figure.
%
% data        = the data selection of interest for mapping
% color_scale = the color scales of the map either [min max] or histogram
%               optimized (max set from 98% of all data-points, minimizing  
%               outlier influence).
% sloi        = the slices of interest including higher dimensions. If all
%               spatial slices are chosen (4th dimension of the data), the
%               sloi is set to NaN.

% USER INPUT: Data of Interest
elm = repmat({'popup'}, 1, 2);
qry = {'Data range to show: ','Color scale range: '};
def = {{'All', 'Specific Slice'}, {'Min to Max', 'Histogram optimized'}};
uans = getInput(elm, qry, def, 'Data Display - Maps');                                                 
if isempty(uans), data = nan; color_scale = nan; sloi = nan; return; end

% \\ Get part of data-array to visualize
if strcmpi(uans{1},'Specific Slice')
    % --- Get slice of interest
    
    % // Prep queries and default answers

    % Default answers: list of 1 to index-size for every dimension
    % including 'all' string, as a cell.
    sl = size(data); sl = sl(4:end); 
    def = arrayfun(@(x) 1:x, sl, 'UniformOutput',0);
    def = cellfun(@(x) arrayfun(@num2str, x, 'UniformOutput', 0),...
        def, 'UniformOutput',0);
    def = cellfun(@(x) ['All' x], def, 'UniformOutput', 0);
    % Default query:
    nqry = numel(def);
    elm = repmat({'popup'}, 1, nqry); 
    qry = repmat({'Index of interest:'}, 1, nqry); 
    for kk = 1:nqry, qry{kk} = sprintf('Dimension %i - %s',kk, qry{kk});
    end
    
    % // User input
    uans_sloi = getInput(elm, qry, def, 'Data Display - Maps');                                                     
    if isempty(uans_sloi)
        data = nan; color_scale = nan; sloi = nan; return; 
    end   
    sloi = str2double(uans_sloi); all_ind = find(isnan(sloi));
    sloi = num2cell(sloi); 
    sloi(all_ind) = arrayfun(@(x) 1:x, sl(all_ind), 'uniform', 0);
    
    % // Process user input

    % Get data of interest.
    data_dim = size(data);    
    data_dim_range = arrayfun(@(x) 1:x, data_dim, 'uniform',0);
    data_dim_range(4:end) = sloi;

    % Get data of interest 
    data_cut = data(data_dim_range{:});

    % Replace SNR-all with SNR-cut
    data = data_cut;

    % If all slices are selected - slice of interest is NaN;.
    if ~isempty(all_ind) && all_ind(1) == 1, sloi = NaN; end
else
    sloi = NaN;
end

% \\ Calculate color-range of maps.
switch uans{2}
    case 'Min to Max'
        color_scale = [min(data(:)) max(data(:))];
    case 'Histogram optimized'
        [N, edges, bin] = histcounts(data(:));
        bool = cumsum(N) <= 0.98*size(bin,1);
        maxval = max(edges(bool));
        color_scale = [min(data(:)) maxval];                         
end

% --- Plot data as tabs and create maps/tab    
function fh_all = CSI_dataAsTabs(csi, data, tag, labels, color_range, sloi)
% Data is shown as a color map and also a numeric representation in the
% figure window. 
%
% clr           = Color-scheme.
% data          = data in of size N-D with Val x X x Y x Z for the 
%                 first 4 indexes.
% tag           = tag/name for the figure window.
% labels        = labels of the data for naming >4-indexes in the window.
% color_range   = range for color-map colors; 'auto' or [low-lim up-lim];
% sloi          = Slice of interest for image-plotting correct image.
%
% Uses the following function: (in order of use)
%
% CSI_dataAsTabs_create_figure: creates the figure with a background color
% and set specific options to the figure.
%   
% CSI_dataAsTabs_create_griddedTabs: Create a tab per slice including a 
% voxel-grid overlay and create plot_par variable. 
%
% CSI_dataAsTabs_addVoxelAxis: Add a transparent axis to each voxel/slice
%
% Then plots the data as colormaps by changing the background color of the
% voxel-axis. Uses a transparency. If a value is equal to NaN, will set
% color to black [0 0 0];

if nargin < 5, color_range = 'auto'; sloi = NaN; end
if nargin == 5, sloi = NaN; end

% Data-Dimensions safety
dim = size(data); 
if numel(dim) <= 2
    fprintf('Aborted data-as-tabs plot, data is not 2D or 3D.\n'); return; 
end

% Number of tab-windows
nwindows = 1; index_range = {1};
if numel(dim) > 4
    nwindows = dim(5); 
    if numel(dim) > 5
        % Data ranges
        index_range = cellfun(@(x) 1:x, num2cell(dim(6:end)), 'uniform', 0); 
    end
end

% Restore data
data_main = data; tag_main = tag;
if numel(labels) > 4
    loi = labels{5};
else
    loi = int2str(sloi);
end

fh_all = cell(1,nwindows);
for kk = 1:nwindows
        
% Get data for nwindow-kk (one dimension above slices) and all
% subsequent dimension.
data = data_main(:,:,:,:,kk, index_range{:});
tag = sprintf('%s - %s - %i', tag_main, loi, kk);

            % -------- % Figure: Create window % ----------------------- %

% Create a figure with specific settings for tabs
fh = CSI_dataAsTabs_create_figure(tag, csi.clr.main);


            % -------- % Figure: Create tabs % ------------------------- %
           
% Create a tab per slice including a voxel-grid overlay and
% create plot_par            
CSI_dataAsTabs_create_griddedTabs(fh, data, csi.clr) ;       

            % -------- % Figure: Create image axis % ------------------- %                

% Get images
plot_img = 0;
if isfield(csi,'conv')
    conv = csi.conv;
    img = MRI_matchSlices(csi, conv); 
    plot_img = 1;
end

if plot_img
    fh = CSI_dataAsTabs_create_ImageAxis(fh);
end

            % -------- % Figure: Create axis/voxel % ------------------- %                

% Add a transparent axis to each voxel per slice
tgui = CSI_dataAsTabs_addVoxelAxis(fh);


            % --------- % Plot Data As Map % --------------------------- %

% Color map data
clr_map = jet(128);
if ischar(color_range)   % Automatic color-range
    clr_val = linspace(min(data(:)),max(data(:)),size(clr_map,1));            
else
    clr_val = linspace(color_range(1),color_range(2),size(clr_map,1));      
end
tgui.plot_par.clr_val = clr_val; tgui.plot_par.clr_rng = color_range;            
tgui.plot_par.clr_map = clr_map;
% transparency value
tgui.plot_par.alpha = 0.33;

plot_par = tgui.plot_par;
for tabi = 1:plot_par.tabs_total                % Sli/tab loop.
    % loadBar(tabi./plot_par.tabs_total , 'Plotting data...');

    % Current tab - its index in tab format
    % This index lacks the slice+1 index (or the dim 1 above the slice
    % dimension). Therefor it requires correction regards pointing to data.
    tab_index = plot_par.tabs_index_table(tabi,:);
    tab_index_cell = num2cell(tab_index);
    
    % For data indexing - convert tab_index_cell to correct index-dims.
    sli = tab_index_cell{1}; % Slice from data.   
    if size(tab_index_cell,2) >= 2
        tab_index_for_data = tab_index_cell(2:end);
    else
        tab_index_for_data = {1};
    end
    
    % Plot images if available
    if plot_img
        if ~iscell(sloi) && isnan(sloi)
            sli_img = tab_index_cell{1}; % Slice from data.
        else
            % Slice of interest - applicable when specific slice is chosen
            % for displaying data. Correct image-slice is taken.
            sli_img = sloi{1}(1); 
        end
        
       
        img2plot = img(:,:,sli_img); 

        if size(img2plot,1) > 1 % Safety against NaN-values
        imagesc(img2plot, 'parent', tgui.himg{tab_index_cell{:}}); 


        % Image Contrast.
        if isfield(conv, 'contrast')
            v = version('-release'); v(end) = []; v = str2double(v);
            if v < 2023
                caxis(tgui.himg{tab_index_cell{:}}, conv.contrast);
            else
                clim(tgui.himg{tab_index_cell{:}}, conv.contrast);
            end
            tgui.plot_par.contrast_img =  conv.contrast;
        else
            contrast_min = min(img2plot(:));
            contrast_max = max(img2plot(:))*0.75;
            if contrast_max <= contrast_min
                contrast_max = contrast_min +1; 
            end 
           
            % Store image contrast
            tgui.plot_par.contrast_img(:,sli_img) = ... 
                [contrast_min contrast_max];
            
            v = version('-release'); v(end) = []; v = str2double(v);
            if v < 2023
                caxis(tgui.himg{tab_index_cell{:}},...
                      tgui.plot_par.contrast_img(:,sli_img));
            else
                clim(tgui.himg{tab_index_cell{:}},...
                     tgui.plot_par.contrast_img(:,sli_img));
            end
                       
        end
        colormap(tgui.himg{tab_index_cell{:}}, gray(255));
        set(tgui.himg{tab_index_cell{:}}, 'box', 'off');

        end
    end
     
    
    % Plot data
    tmp_data = squeeze(data(:,:,:,sli,1,tab_index_for_data{:}))';
    tgui.plot_par.plotobj{tab_index_cell{:}} = ...
        imagesc(tgui.plot_par.ax{tab_index_cell{:}}, tmp_data);
    
    % Colormap and transparency
    set(tgui.plot_par.ax{tab_index_cell{:}},...
        'colormap', clr_map, 'clim', color_range);
    transmap = double(~isnan(tmp_data)); 
    transmap(transmap == 1) = tgui.plot_par.alpha;
    set(tgui.plot_par.plotobj{tab_index_cell{:}}, 'AlphaData', transmap)
    tgui.plot_par.transmap{tabi} = transmap;

    % Axes cosmetics
    set(tgui.plot_par.ax{tab_index_cell{:}},...
        'LineWidth', 0.1, 'Xtick',[], 'Ytick', [],...
        'Box', 'off'); 
    
    % Turn of outer-axis i.e. data-imagesc
    tgui.plot_par.ax{tab_index_cell{:}}.Visible = 'off';    

end
% loadBar(NaN);

% Save tgui-data to gui itself 
guidata(fh, tgui);

% Save current object to tab-gui
fh_all{kk} = tgui;

% Create toolbar in figure
toolbar_create(tgui.fig);


end % end of nsub loop. %


% Plot colorbar
colorbarQ(clr_map, clr_val);

% --- Create axis for image plot in tab-figure
function fh = CSI_dataAsTabs_create_ImageAxis(fh)
% Add an axis meant for displaying overlay images. This should be done
% before adding voxel-axis to the figure or tab.

% Get GUI data of figure
tgui = guidata(fh);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

% Loop each tab of figure
for sli = 1:plot_par.tabs_total                % Sli loop/tabs loop
    
    tab_index = plot_par.tabs_index_table(sli,:);
    tab_index_cell = num2cell(tab_index);
    
    % Create an axis at the size of the figure.
    tgui.himg{tab_index_cell{:}} = ...
        axes('parent',tgui.tabh{tab_index_cell{:}},...
                'Position',[0 0 1 1], 'Color', 'None');
    
    set(tgui.himg{tab_index_cell{:}},...
                   'Color','None',...
                   'XColor', plot_par.colors.main,...
                   'YColor', plot_par.colors.main,...
                   'LineWidth', 1.7, 'Xtick',[], 'Ytick', [],...
                   'TickLength',[0 0.00001], 'Box', 'off'); 
     
end
guidata(fh, tgui);

% --- Create figure for CSI_dataAsTabs
function fh = CSI_dataAsTabs_create_figure(tag, clr)
% Create a window for the dataAsTab functions
% tag: name of the figure
% clr: background color

% Input handling
if nargin < 2, clr = [0,0,0]; end

% -------- % Figure: Create window % -------- %

% Create figure
fh = figure('Tag', tag ,'Name', tag ,...
            'Color', clr, 'Toolbar', 'None',...
            'MenuBar', 'None', 'NumberTitle', 'Off');                   

% 1. Default figure size and screen pixel size
def_sz = 720; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
% 2. Ratio to define figure height to def_size
fig_sz = [def_sz def_sz.*(scr_sz(2)/scr_sz(1))];
% 4. Position of figure.
fig_ps = [40 scr_sz(2)-(1.15*fig_sz(2))];
% 5. Apply
set(fh, 'Position', [fig_ps fig_sz]); 

% --- Create grid overlay in every tab (CSI_dataAsTabs)
function tgui = CSI_dataAsTabs_create_griddedTabs(fh, data, clrs)
% Add tabs and a grid to the figure
% 
% fg is the figure from CSI_dataAsTabs_create_figure()

% Create tab group
tgui = struct;
tgui.tabg = uitabgroup(fh); 


% Plotting parameters
plot_par = struct;
plot_par.colors   = clrs;              
plot_par.dim      = size(data);           % Data dimensions
plot_par.dim(1)   = [];                   % Remove data index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 


% Tab total parameters.
tabs_to_plot = prod(plot_par.dim(5:end)) * size(data,4);
% index_table = cell(1,tabs_to_plot);
index_range = cellfun(@(x) 1:x, ...
    num2cell([size(data,4) plot_par.dim(5:end)]), 'uniform', 0); 

% This array contains all indexes for each slice (first column) and all
% subsequent indexes except the slice+1 index - which is depicted as a
% seperate window.
indArray = slice2rowIndex(index_range);

plot_par.tabs_index_table = indArray;
plot_par.tabs_total = tabs_to_plot;


            % -------- % Figure: Calculate Axis Grid % ----------------- %
% Calculate the mesh-grid for the slice representing the voxels

% Resolution without any correction of linewidth
plot_par.res = 1./plot_par.dim(1:2); 

% Loop X, Y and Z
% X, Y and Z == data/csi-space. Position in figure starts at
% point zero (0). Resolution equals the nr of CSI-voxels to plot in the row 
% (y in image)and column (x in image) dimension of the figure. 
for kk = 1:numel(plot_par.res)
    plot_par.range{kk} = ...
    0 : plot_par.res(kk) : ...
       (plot_par.res(kk) * plot_par.dim(kk)-plot_par.res(kk));
end

% Caluclate a grid for each axis in the figure representing the voxels in
% CSI data.
[x,y] = meshgrid(plot_par.range{1},plot_par.range{2});
plot_par.grid.x = (x); plot_par.grid.y = flipud(y);

% 1D Correction 
% Transpose x/col and y/row. Creats Nx1 vs 1xN lists of the
% axis grid coordinates.
if plot_par.data_dim == 1, plot_par.grid.y = y'; plot_par.grid.x = x'; end



            % --------- % Figure: Loop slices/tabs % ------------------- %
% Adds a tab to the figure for all slices and other dimension (excluding 
% the index slice+1; which is seperated into different windows.
for tabi = 1:plot_par.tabs_total % Loop each tab 
    
    % Plotting index
    tmp_plotindex_tabs = plot_par.tabs_index_table(tabi,:);
    tmp_plotindex_tabs_cell = num2cell(tmp_plotindex_tabs);
    
    % Create a tab in the tabgroup for this slices
    tgui.tabh{tmp_plotindex_tabs_cell{:}} = ...
        uitab(tgui.tabg, 'Title', int2str(tmp_plotindex_tabs),...
        'BackgroundColor', plot_par.colors.main,...
        'ForegroundColor', plot_par.colors.grid);

    % Plot voxel grid
    % Input: target figure, target figure size, data dimensions, range and color.
    CSI_2D_grid(tgui.tabh{tmp_plotindex_tabs_cell{:}},...
        fh.Position(3:4), plot_par.dim, ...
        plot_par.range, plot_par.colors.grid);

end % End of slice/tabs loop

tgui.plot_par = plot_par; tgui.fig = fh;

% Add 
guidata(fh, tgui);

% --- Add axis to voxels in CSI_dataAsTabs
function tgui = CSI_dataAsTabs_addVoxelAxis(fh)
% Add to a transparent axis to each voxel per slice

% Get GUI data of figure
tgui = guidata(fh);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

% Storage for each axis.
ax = cell(1, plot_par.tabs_total);

% Loop each tab of figure
for sli = 1:plot_par.tabs_total                % Sli loop/tabs loop
    % loadBar(sli./plot_par.tabs_total , 'Adding voxel-axis...');

    tab_index = plot_par.tabs_index_table(sli,:);
    tab_index_cell = num2cell(tab_index);
    
    % Create axis
    ax{tab_index_cell{:}} = axes('parent',tgui.tabh{tab_index_cell{:}});     
    
    tbsz = tgui.tabh{tab_index_cell{:}}.Position;
    pos = [0 0 tbsz(3) tbsz(4)];
    set(ax{tab_index_cell{:}}, 'Unit', 'Normalized',...
           'Position', pos,...
           'LineWidth', 0.1, 'Xtick',[], 'Ytick', [],...
           'Box', 'off');   

end
tgui.plot_par.ax = ax; guidata(fh, tgui);

% --- Initiate loading/creating voxel-mask 
function csi = VoxelMask_Initiate(csi, clr)

if nargin < 2
    clr.main = [0 0 0]; clr.text_title = [0.941 0.941 0.941];
    clr.hilight1 = [0.8 0 0];
end


% Get image-data
plot_img = 0;
if isfield(csi, 'conv')
    conv = csi.conv; plot_img = 1;
    img = MRI_matchSlices(csi, conv); 
    if isfield(conv,'contrast')
        contrast = conv.contrast;    
    else
        contrast = [min(img(:)) max(img(:))];
        if contrast(2) <= contrast(1), contrast(2) = contrast(1)+1; end
    end
end

% Ask user to load or use editor
uans = getInput({'popup'}, {'Create voxel-mask:'}, {{'Editor', 'Import'}} , 'Voxelmask');
if isempty(uans), return; end

switch uans{1}
    case 'Editor'

    % Set plot parameters % ------------------- %
    % Plot_par requires data-dimensions and colors. Other parameters are
    % calculated by VoxelMask_Editor.
    plot_par = struct; 
    
    plot_par.colors = clr; % clr... main, text_title, hilight1
    plot_par.plot_img = plot_img; % Image-nfo (boolean, image-data, contrast)
    if plot_par.plot_img, plot_par.img = img; plot_par.contrast = contrast;
    else,                 plot_par.img = NaN; plot_par.contrast = NaN;
    end
    
    % Volume and voxel dimensions
    plot_par.dim      = size(csi.data.raw);   % Data dimensions
    plot_par.dim(1)   = [];                   % Remove time index e.g. 1
    plot_par.dim = plot_par.dim(1:3); plot_par.data_dim = 3; % 3D Only.
    
    % Add previous-voxel mask if available.
    if isfield(csi, 'voxelmask')
        msz = size(csi.voxelmask);
        if numel(msz) > 3
            ind = num2cell(ones(1,numel(msz)-4)); 
            plot_par.mask = squeeze(csi.voxelmask(:,:,:,:,ind{:})); 
        else
            plot_par.mask = squeeze(csi.voxelmask); 
        end
    end
    
    % Start VoxelMask_Editor - mask is returned on exit of editor.
    mask = VoxelMask_Editor(plot_par);

    case 'Import'
                
        % Get file selection
        if isfield(csi,'filepath'), fp = csi.filepath; else, fp = []; end
        [fi, fp, idx] = uigetfile({'*.mat'}, 'Import mask-data', fp);
        if idx == 0, return; end

        % Load mask
        load([fp fi], 'mask');
end


% Convert mask to data-size % ----------------------------------------  %
csz = size(csi.data.raw); msz = size(mask);

% Permute vector to move mask to spatial dimension
permv = [numel(msz)+1 1:numel(msz)];
mask = permute(mask, permv); msz = size(mask);    

if numel(msz) ~= numel(csz) % If there are higher order indexes    
    % Repeat mask over non-spatial dimensions    
    dimrem = csz(numel(msz)+1:numel(csz));
    dimrem = [ones(1,numel(msz)) dimrem];

    mask = repmat(mask, dimrem); msz = size(mask);
    if sum(msz(2:end) ~= csz(2:end))
        mask = zeros(size(csi.data.raw));
    end
end

% Clean up
csi.data.voxelmask = mask;

% --- Match MRI-slices to CSI-volume
function [img, img_all, img_all_slice_range] = MRI_matchSlices(csi, conv)
% Using the conv struct, it calculates images corresponding to the MRSI
% data indexing. Number of slices in img is equal to number of slices in
% MRSI data.

        % ----- % Create MRS matching images % ----- %

nSlices = size(csi.data.raw, 4);
szData = size(conv.data);
img = NaN(szData(1),szData(2),nSlices);

cz = unique(csi.ori.mesh.z); mz = unique(conv.mesh.z);
[v, i] = min(abs(cz - repmat(mz, size(cz,1))')');
nMz_per_Cz = csi.ori.res(3)./conv.res(3);
nMz_per_Cz_step = floor(nMz_per_Cz ./2);

min_of_slice = i - (nMz_per_Cz_step-1);
min_of_slice(min_of_slice <= 0) = 1;

max_of_slice = i + nMz_per_Cz_step;
max_of_slice(max_of_slice > size(conv.data,3)) = size(conv.data,3);

img_all_slice_range = [min_of_slice' max_of_slice'];

for sli = 1:nSlices              % For every CSI slice

    % Slice coordinates for CSI and CONV
    cz = unique(csi.ori.mesh.z); mz = unique(conv.mesh.z);

    % Find CONV slices matching to CSI slices.
    % img_range = CSI2MRI(cz(sli), mz, csi.ori.res(3), conv.res(3));
    % if sli == 1, img_all_slice_range = NaN(nSlices, size(img_range,2)); end
    % img_all_slice_range(sli,:) = img_range;
    
    img_range = img_all_slice_range(sli,:);

    % Minimum and maximum index
    indMn = img_range(1); indMx = img_range(2);
    
    % 08/01/24 QH
    if indMx > size(conv.data,3), indMx = size(conv.data,3); end

    % Type of converted CSI-slice image
    pstr = 'Projection';

    % Calculate matching image
    switch pstr        
        case 'Summation'
            img(:,:,sli) = sum(conv.data(:,:,indMn:indMx),3); 
        case 'Projection'
            img(:,:,sli) = mean(conv.data(:,:,indMn:indMx),3);  
        case 'Minimum'
            img(:,:,sli) = conv.data(:,:,indMn); 
        case 'Maximum'
            img(:,:,sli) = conv.data(:,:,indMx); 
        case 'Middle'
            imRange = indMn:indMx; 
            midsl = imRange(round(size(imRange,2)/2));
            img(:,:,sli) = (conv.data(:,:,midsl)); 
    end    


end % End of slice loop

img_all = conv.data;

% --- Executes to plot a 2D grid in map-plots
function CSI_2D_grid(target, target_sz, dim, range, grid_clr)

% Get figure size for normalization of grid thickness
w = target_sz(1); h = target_sz(2);
% Define using vertical line the thickness of horizontal ones 
wv = 1.08; wh = wv; wv = wv./w; wh = wh./h; 
% Line height or width for both vertical and horizontal grid lines;
h = 1; 

% Vertical lines (e.g. the x direction -> nr columnes - 1)
for kk = 1:dim(1)-1
    % Position of grid
    xpos = range{1}(kk+1);
    % Set text
    uicontrol(target,'Style','Text','Unit','Normalized',...
             'Position',[xpos 0 wv h],'BackgroundColor',grid_clr);    
end
    
% Horizontal (e.g. the y direction -> nr rows - 1)
for kk = 1:dim(2)-1
    % Position of grid
    ypos = range{2}(kk+1);
    % Set text
    uicontrol(target,'Style','Text','Unit','Normalized',...
             'Position',[0 ypos h wh],'BackgroundColor',grid_clr);
end


%% TOOLBAR FUNCTIONS % -------------------------------------------------- %


% --- Executed by data visualization with seperate figures functions
function toolbar_create(figure_object)
% Add a toolbar to a figure with datacursor mode and save-figure as button.
% Uses functions: toolbar_dataPointer, toolbar_saveFigure.


% Create toolbar in figure
tb = uitoolbar(figure_object);

% Data pointer toggle button for toolbar
tt_data = uitoggletool(tb); % Toggle button - pointer
tt_data_icon = imread('Images\mouse_pointer.png');
tt_data_icon = imresize(tt_data_icon,[20,20]);
tt_data_icon(tt_data_icon == max(tt_data_icon(:))) = ...
    round(max(tt_data_icon(:))*0.94);
tt_data.CData = tt_data_icon;
tt_data.ClickedCallback = @toolbar_dataPointer;
tt_data.TooltipString = 'Toggle data cursor.';

% Save figure button for toolbar
tt_save = uipushtool(tb); % Push button - save data
tt_save_icon = imread('Images\floppy_disk.png');
tt_save_icon = imresize(tt_save_icon,[20,20]);
tt_save_icon(tt_save_icon ==  max(tt_save_icon(:))) = ...
    round(max(tt_save_icon(:))*0.94);
tt_save.CData = tt_save_icon;
tt_save.ClickedCallback = @toolbar_saveFigure;
tt_save.TooltipString = ...
    'Save figure. NB. Transparency is lost when saving as matlab-figures.';

% Save all tabs as image button for toolbar
tt_save_all = uipushtool(tb); % Push button - save data
tt_save_all_icon = imread('Images\floppy_disk_save_all.png');
tt_save_all_icon = imresize(tt_save_all_icon,[20,20]);
tt_save_all_icon(tt_save_all_icon ==  max(tt_save_all_icon(:))) = ...
    round(max(tt_save_all_icon(:))*0.94);
tt_save_all.CData = tt_save_all_icon;
tt_save_all.ClickedCallback = @toolbar_saveFigure_allTabs;
tt_save_all.TooltipString = ...
    'Save all tabs to image-file automagically.';

% Set the colorscale range of the plot
tt_color_scale = uipushtool(tb); % Push button - color range scale
tt_color_scale_icon = imread('Images\colorbar.png');
tt_color_scale_icon = imresize(tt_color_scale_icon,[20,20]);
tt_color_scale_icon(tt_color_scale_icon ==  max(tt_color_scale_icon(:))) = ...
    round(max(tt_color_scale_icon(:))*0.94);
tt_color_scale.CData = tt_color_scale_icon;
tt_color_scale.ClickedCallback = @toolbar_SetColorScale;
tt_color_scale.TooltipString = ...
    'Adjust color-scale range of plotted data.';


% Change transparency button for toolbar
tt_alphaD = uipushtool(tb); % Push button - change alpha
tt_alphaD_icon = imread('Images\transparency4.png');
tt_alphaD_icon = imresize(tt_alphaD_icon,[20,20]);
tt_alphaD_icon(tt_alphaD_icon ==  max(tt_alphaD_icon(:))) = ...
    round(max(tt_alphaD_icon(:))*0.94);
tt_alphaD.CData = tt_alphaD_icon;
tt_alphaD.ClickedCallback = @toolbar_SetAlphaValue;
tt_alphaD.TooltipString = ...
    'Change the transparency of the color-map i.e. the voxels.';


% Restore transparency button for toolbar
tt_alpha = uipushtool(tb); % Push button - change alpha
tt_alpha_icon = imread('Images\fix_transparency6.png');
tt_alpha_icon = imresize(tt_alpha_icon,[20,20]);
tt_alpha_icon(tt_alpha_icon ==  max(tt_alpha_icon(:))) = ...
    round(max(tt_alpha_icon(:))*0.94);
tt_alpha.CData = tt_alpha_icon;
tt_alpha.ClickedCallback = @toolbar_fixTransparency;
tt_alpha.TooltipString = ...
    'Fix transparency of color-map. Alpha-maps arent stored by Matlab.';


% Grid on/off
tt_grid = uipushtool(tb); % Push button - change alpha
tt_grid_icon = imread('Images\grid2.png');
tt_grid_icon = imresize(tt_grid_icon,[20,20]);
tt_grid_icon(tt_grid_icon ==  max(tt_grid_icon(:))) = ...
    round(max(tt_grid_icon(:))*0.94);
tt_grid.CData = tt_grid_icon;
tt_grid.ClickedCallback = @toolbar_toggleGrid;
tt_grid.TooltipString = 'Toggle visibility of the voxel-grid.';

% Data2Table
tt_data2table = uipushtool(tb); % Push button - data 2 table
tt_data2table_icon = imread('Images\table3.png');
tt_data2table_icon = imresize(tt_data2table_icon,[20,20]);
tt_data2table_icon(tt_data2table_icon ==  max(tt_data2table_icon(:))) = ...
    round(max(tt_data2table_icon(:))*0.94);
tt_data2table.CData = tt_data2table_icon;
tt_data2table.ClickedCallback = @toolbar_MapToTable;
tt_data2table.TooltipString = ...
    'Show data in a table. Also allows saving data to file.';

% Current-plot Statistics
tt_statistics = uipushtool(tb);
tt_statistics_icon = imread('Images\statistics_icon.png');
tt_statistics_icon = imresize(tt_statistics_icon,[20,20]);
tt_statistics_icon(tt_statistics_icon ==  max(tt_statistics_icon(:))) = ...
    round(max(tt_statistics_icon(:))*0.94);
tt_statistics.CData = tt_statistics_icon;
tt_statistics.ClickedCallback = @toolbar_Statistics;
tt_statistics.TooltipString = ...
    'Calculate statistics of plotted data.';

% --- Executes on press of toolbar's save Figure button
function toolbar_saveFigure(src, ~)
% Save figure when clicking save button in toolbar.
%
% ISSUE: Figure's of matlab are not saved with their transparency value.
% This means any transparent background or field will disappear when
% opened afterward.

% Get figure object
fig = src.Parent.Parent;

% Get save-location
[fn,fp,index] = uiputfile(...
    {'*.fig','Matlab'; '*.png','PNG'; '*.jpeg','JPEG'},'Save figure.');
if  index == 0, return; end

% Save figure
[~,~,ext] = fileparts(fn); ext = ext(2:end);
export_fig([fp '\' fn], '-transparent',['-' ext],'-nocrop','-m1', fig);

% --- Executes on press of toolbar's data Pointer toggle.
function toolbar_dataPointer(src, ~)
% Toggles the data cursor in a figure using the src-button.
datacursormode(src.Parent.Parent, 'toggle');

% --- Executes on press of toolbar's fix transparency button.
function toolbar_fixTransparency(src,~)
% Loop all axes in the figure and set transparency element.

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

if ~isfield(tgui,'plot_par')
    fprintf('Error: no plot-data stored in figure'); return; 
end
plot_par = tgui.plot_par;

% \\ ADJUST: axis alpha value
ax = plot_par.plotobj(:);
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        
        tab_index = plot_par.tabs_index_table(axi,:);
        tab_index_cell = num2cell(tab_index);
    
        % Colormap and transparency     
        plot_par.transmap(plot_par.transmap == plot_par.alpha) =...
            plot_par.alpha;
        set(tgui.plot_par.plotobj{tab_index_cell{:}}, ...
            'AlphaData', plot_par.transmap);
    end    
end
guidata(figobj, tgui);

% --- Executes on press of toolbar's save all tabs button.
function toolbar_saveFigure_allTabs(src,~)

% Get tab-group object
figobj = src.Parent.Parent;
for kk = 1:size(figobj.Children,1)
    % Find tab group child of figure object.
    if strcmp(figobj.Children(kk).Type,'uitabgroup')
        tabg = figobj.Children(kk); ind = kk;
    end
end
if ~exist('tabg','var'), return; end

% Get save-location
[fn, fp, index] = uiputfile(...
    {'*.png','PNG'; '*.jpeg','JPEG'},'Save figure.');
if  index == 0, return; end
[~,fn,ext] = fileparts(fn); ext = ext(2:end);
    
% Set active-tab and print to file.
tabs = tabg.Children; scrsz = get(0,'screensize'); scrsz = scrsz(3:4);
for kk = 1:size(tabs,1)
    % Set tab-kk to active
    tabg.SelectedTab = tabs(kk);
     
    pause(0.5); % Make sure display is updated.       

    check = 0; n = 0;
    while check ~= 1
        tabg = figobj.Children(ind);
        if isequal(tabg.SelectedTab, tabs(kk)), check = 1; end
        n = n + 1; if n == 100, break; end
    end

    pause(0.5);

    tab_title_str = str2double(strsplit(tabs(kk).Title));
    tab_title_str = sprintf('%03i', tab_title_str);

    % Save figure
    fntmp = [fp fn '_' tab_title_str '.' ext];
    
    pause(0.75);

    % Solution is screenshot!
    figpos = figobj.Position; % [l b w h]
    figpos(2) = abs(figpos(2) + figpos(4) - scrsz(2));
    img = screenshot(figpos); % [l t w h]
    
    % Save to file
    imwrite(img, fntmp);    
end

% --- Executes on press of toolbar's set color scale button
function toolbar_SetColorScale(src,~)
% Adjust color scale range of current tab-plot.
% Works for maps i.e. CSI_dataAsTabs

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

if ~isfield(tgui,'plot_par')
    fprintf('Error: no plot-data stored in figure'); return; 
end
plot_par = tgui.plot_par;
if ~isfield(plot_par,'clr_val')
    fprintf('Error: no color-values stored in figure'); return; 
end

% Current tab
toi = cell2mat(cellfun(@isequal,tgui.tabh, ...
        repmat({tgui.tabg.SelectedTab}, size(tgui.tabh)), 'uniform', 0));
toi = tgui.tabh{toi == 1};
tab_dims = strsplit(toi.Title, ' ');
tab_dims = num2cell(str2double(tab_dims));

% \\ USER INPUT: new color range

% Current color-range
cur_clr_rng = plot_par.clr_rng;

% Check if images are plotted
contrast = NaN; 
if isfield(tgui.plot_par, 'contrast_img'), contrast = 1; end

% Get image-contrast
if ~isnan(contrast) && numel(tgui.plot_par.contrast_img) == 2
    contrast = tgui.plot_par.contrast_img;
elseif ~isnan(contrast)
    contrast = tgui.plot_par.contrast_img(:,tab_dims{1})';
end

uans = getInput({'edit', 'edit'},...
                {'Color-scale range for data-map:','Image contrast:'}, ...
                {num2str(cur_clr_rng), num2str(contrast)}, ...
                'Data Display - Maps');            
if isempty(uans), return; end
new_clr_rng = str2double(strsplit(uans{1}));

% Calculate new color-range table
clr_map = plot_par.clr_map;
clr_val = linspace(new_clr_rng(1),new_clr_rng(2),size(clr_map,1));      

% Store into gui-data-obj
tgui.plot_par.clr_val = clr_val; tgui.plot_par.clr_rng = new_clr_rng;            
tgui.plot_par.clr_map = clr_map;

% New Image contrast
tgui.plot_par.contrast_img = str2double(strsplit(uans{2}));
if ~isnan(tgui.plot_par.contrast_img)
    himg = [tgui.himg{:}]; set(himg, 'CLim', tgui.plot_par.contrast_img);
end

% \\ ADJUST: axis colors
ax = [plot_par.ax{:}];
set(ax,'colormap', clr_map, 'clim', new_clr_rng);

% Update GUI
guidata(figobj, tgui)

% Plot colorbar with new values.
colorbarQ(clr_map, clr_val);

% --- Executes on press of toolbar's set alpha value button
function toolbar_SetAlphaValue(src,~)
% Change the transparency of the voxels, i.e. the alpha value.

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

if ~isfield(tgui,'plot_par')
    fprintf('Error: no plot-data stored in figure'); return; 
end
plot_par = tgui.plot_par;
if isfield(plot_par,'alpha')
    cur_alpha_val = plot_par.alpha;
else
    cur_alpha_val = 1;
end

% \\ USER INPUT: new color range
uans = getInput({'edit'}, {'New transparency value:'},...
    {num2str(cur_alpha_val)}, 'Data Display - Maps');
if isempty(uans), return; end
new_alpha = str2double((uans{1}));
tgui.plot_par.alpha = new_alpha;

% \\ ADJUST: axis alpha value
ax = plot_par.plotobj(:);
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        
        tab_index = plot_par.tabs_index_table(axi,:);
        tab_index_cell = num2cell(tab_index);
    
        % Colormap and transparency     
        plot_par.transmap{axi}(plot_par.transmap{axi} == cur_alpha_val) = ...
            new_alpha;
        set(tgui.plot_par.plotobj{tab_index_cell{:}}, ...
            'AlphaData', plot_par.transmap{axi});
    end    
end
guidata(figobj, tgui);

% --- Executes on press of toolbar's convert to table button
function toolbar_MapToTable(src,~)
% Get data from maps and display in table

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

% \\ GET: Data
ph = [tgui.plot_par.plotobj{:}];
data = cat(3,ph.CData);
 
% data = reshape(data, size(tgui.plot_par.ax));
sz = size(data);
data = permute(data, [numel(sz)+1 1:numel(sz)]);

CSI_dataAsTable(data,tgui.fig.Name);

% --- Executes on press of toolbar's magnifying-glass-stats-button
function toolbar_Statistics(src,~)
% do something
% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

% Which data to visualize
uans = getInput({'popup'}, {'Calculate stastics for:'},...
                         {{'Current', 'All', 'All/Slice'}}, ...
                           'Data Display - Maps');
if isempty(uans), return; end


% \\ GET: Data
ph = [tgui.plot_par.plotobj{:}];
data = cat(3,ph.CData);

% Correct X/Y
permvec = 1:numel(size(data)); permvec([2 1]) = permvec([1 2]);
data = permute(data, permvec);
 
% Reshape to [1 x X x Y x Z];
data = reshape(data, tgui.plot_par.dim); sz = size(data);
data = permute(data, [numel(sz)+1 1:numel(sz)]);


% Tab of interest
if strcmp(uans{1}, 'Current')
    toi = cell2mat(cellfun(@isequal,tgui.tabh, ...
        repmat({tgui.tabg.SelectedTab}, size(tgui.tabh)), 'uniform', 0));
    toi = tgui.tabh{toi == 1};
    tab_dims = strsplit(toi.Title, ' ');
    tab_dims = num2cell(str2double(tab_dims));
    data = squeeze(data(:,:,:,tab_dims{:}));

    % Calculate and show results.
    stats = csi_statistics_of_volume(data);
    stats.index = toi.Title;
    Statistics_Viewer(stats);

elseif strcmp(uans{1}, 'All')
    % Calculate and show results.
    stats = csi_statistics_of_volume(data);
    Statistics_Viewer(stats);

elseif strcmp(uans{1}, 'All/Slice')
    for kk = 1:numel(tgui.tabh)
        % Get data of tab
        data = tgui.plot_par.plotobj{kk}; data = data.CData;
        toi = tgui.tabh{kk};
        % Statistics
        stats = csi_statistics_of_volume(squeeze(data));
        % NFO
        stats.tab = kk; stats.index = toi.Title;
        % Display to user
        Statistics_Viewer(stats);
    end
end

% --- Toggle the voxel-grid in the tabs
function toolbar_toggleGrid(src,~)

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

if ~isfield(tgui, 'tabh'), return; end

% Toggle visibility of axes on and off.
for tabi = 1:numel(tgui.tabh)
    kids = tgui.tabh{tabi}.Children;
    for kk = 1:numel(kids)
        if strcmp(kids(kk).Type, 'uicontrol')
            if strcmp(kids(kk).Style, 'text')
                if strcmp(kids(kk).Visible, 'on')
                    kids(kk).Visible = 'off';
                else
                     kids(kk).Visible = 'on';
                end
            end
        end
    end
end