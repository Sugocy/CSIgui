function [varargout] = VoxelMask_Editor(plot_par)
% Voxel Editor - V2 - Dr. Quincy van Houtum, 2024
%
% Input ---------------------------------------------------------------- %
% plot_par struct with fields...
% dim:          the dimensions of the voxel-grid [x, y, z];
% plot_img:     boolean for plotting image-data too (1), or not (0)
% img:          image data to plot.
% contrast:     contrast for the image-plot, or NaN to auto-calculate via
%               min and max of image data
% colors:       struct with colors for GUI-color coding expecting fields: 
%               main/text_title/hilight1
% mask:         an initial mask to highlight. Can be used to reload a 
%               previous mask.
%
%
% Output --------------------------------------------------------------- %
%
% Mask of size dim with value 1 where voxels were selected.
%
%
% Version 2: updated plot-method to use single-axis per tab to speed up 
%            GUI startup (by a factor of 15!).

if ~isfield(plot_par, 'colors')
    plot_par.colors.main         = [0.1 0.1 0.1];
    plot_par.colors.text_title   = [0.502 0.502 0.502];
    plot_par.colors.hilight1     = [0.8 0 0];
end

% --- % Create Figure for plotting 2D
plot_par.fh = VoxelMask_createFigure(...
              'Voxelmask-Editor', plot_par.colors.main);

% --- % Calculate plot-parameters
% Voxel-resolution, voxel axes meshgrid.
plot_par = VoxelMask_PlotParameters(plot_par);

% --- % Add tabs with voxel-grid
plot_par = VoxelMask_setTabs(plot_par);

% ---------------------------------------------------------------------- %
% ------------- Guidata(fh) contains the plot_par structure ------------ %
% ---------------------------------------------------------------------- %
% Includes tab-group and handles. All other GUI-elements will be added to
% the guidata.

% --- % Add image-axis.
VoxelMask_addImageAxes(plot_par.fh);

% --- % Add voxels per tab
VoxelMask_addVoxelAxes(plot_par.fh);

% --- % Add images/slice
VoxelMask_addImages(plot_par.fh);

% --- % Add menubar with options
VoxelMask_addMenu(plot_par.fh)

% --- % Add visual mask
VoxelMask_PlotMask(plot_par.fh);

tgui = guidata(plot_par.fh);
plot_par = tgui.plot_par;

set(tgui.plot_par.fh,'WindowButtonMotionFcn',@mouseHover);
set(tgui.plot_par.fh,'WindowbuttonDownFcn',@mouseClick);

% Get GUI data of figure
guidata(tgui.plot_par.fh, tgui);

% Set output
fh = plot_par.fh;


% ---------------------------------------------------------------------- %

% Wait until user closes.
uiwait(fh); % UIresume by closeRqstFnc.

if nargout
    % Update gui-data
    tgui = guidata(fh);
    % create mask-output
    varargout{1} = tgui.plot_par.mask;

    delete(fh);
else
    delete(fh);
end
delete(fh);

function fh = VoxelMask_createFigure(tag, clr)
% Create a figure-window with specific settings.
% Create a window for the dataAsTab functions
% tag:  Name of the figure
% clr:  Background color - if empty, default set to black.

% Input handling
if nargin < 2, clr = [0,0,0]; end

% -------- % Figure: Create window % -------- %

% Create figure
fh = figure('Tag', tag ,'Name', tag ,...
            'Color', clr, 'Toolbar', 'None',...
            'MenuBar', 'None', 'NumberTitle', 'Off',...
            'CloseRequestFcn', @closeReqFcn);                   

% 1. Default figure size and screen pixel size
def_sz = 720; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
% 2. Ratio to define figure height to def_size
fig_sz = [def_sz def_sz.*(scr_sz(2)/scr_sz(1))];
% 4. Position of figure.
fig_ps = [40 scr_sz(2)-(1.15*fig_sz(2))];
% 5. Apply
set(fh, 'Position', [fig_ps fig_sz]); 

function plot_par = VoxelMask_PlotParameters(plot_par)
% The voxel-specific plot_parameters are calculated for all slices.
% Calculate the mesh-grid for the slice representing the voxels, including
% voxel resolution in the figure.
%
% Applies corrections for 1D data and corrections for X/Y - COL/ROW
% matrix convention in MATLAB.

% Resolution without any correction of voxel-grid linewidth
plot_par.res = 1./plot_par.dim(1:2); 

% Loop X, Y and Z
% X, Y and Z == data/csi-space. Position in figure starts at
% point zero (0). Resolution equals the nr of CSI-voxels to plot in the row 
% (y in image) and column (x in image) dimension of the figure. 
for kk = 1:numel(plot_par.res)
    plot_par.range{kk} = ...
    0 : plot_par.res(kk) : ...
       (plot_par.res(kk) * plot_par.dim(kk)-plot_par.res(kk));
end

% Calculate a grid for each axis in the figure representing the voxels in
% CSI data.
[x,y] = meshgrid(plot_par.range{1},plot_par.range{2});
plot_par.grid.x = (x); plot_par.grid.y = flipud(y);

% 1D Correction 
% Transpose x/col and y/row. Creats Nx1 vs 1xN lists of the
% axis grid coordinates.
if numel(plot_par.dim) == 1, plot_par.grid.y = y'; plot_par.grid.x = x'; end

function plot_par = VoxelMask_setTabs(plot_par)
% Add tabs and voxel-line-grid to figure for voxel-mask selection.
% 
% plot_par.fh is the figure from CSI_VoxelMask_createFigure() but can be
% any figure-handle.

% Create tab group
tgui = struct; tgui.tabg = uitabgroup(plot_par.fh); 

% Tab total parameters: only 3rd dimension data
tabs_to_plot = plot_par.dim(3);
index_range = {1:plot_par.dim(3)};

% This array contains all indexes for each slice (first column) 
indArray = allCombinations(index_range)';
plot_par.tabs_index_table = indArray;
plot_par.tabs_total = tabs_to_plot;


% --------- % Figure: Loop slices/tabs % ------------------- %
% Adds a tab to the figure for all slices
for tabi = 1:plot_par.tabs_total % Loop each tab 
    
    % Plotting index
    tmp_plotindex_tabs = plot_par.tabs_index_table(tabi,:);
    tmp_plotindex_tabs_cell = num2cell(tmp_plotindex_tabs);
    
    % Create a tab in the tabgroup for this slices
    tgui.tabh{tmp_plotindex_tabs_cell{:}} = ...
        uitab(tgui.tabg, 'Title', int2str(tmp_plotindex_tabs),...
        'BackgroundColor', plot_par.colors.main,...
        'ForegroundColor', plot_par.colors.text_title);

    % Plot voxel grid
    % Input: target figure, target figure size, data dimensions, range and color.
    VoxelMask_addGrid(tgui.tabh{tmp_plotindex_tabs_cell{:}},...
        plot_par.fh.Position(3:4), plot_par.dim, ...
        plot_par.range, plot_par.colors.text_title);

end % End of slice/tabs loop
tgui.plot_par = plot_par; tgui.fig = plot_par.fh;

% Add 
guidata(plot_par.fh, tgui);

function VoxelMask_addImageAxes(fh)
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

function VoxelMask_addVoxelAxes(fh)
% Add to a transparent axis to each voxel per slice

% Get GUI data of figure
tgui = guidata(fh);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

plot_par.alpha = 0.35; 

% Storage for each axis.
ax = cell(1,plot_par.tabs_total);

% Loop each tab of figure
for sli = 1:plot_par.tabs_total                % Sli loop/tabs loop
    %     loadBar(sli./plot_par.tabs_total , 'Adding voxel-axis...');

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
plot_par.ax = ax;  tgui.plot_par = plot_par; guidata(fh, tgui); loadBar(NaN);

function VoxelMask_addImages(fh)
% Add image to each tab/slice of the figure.

% Get GUI data of figure
tgui = guidata(fh);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

% If no images need to be plotted, return.
if plot_par.plot_img == 0, return; end

tgui = guidata(plot_par.fh);


for tabi = 1:plot_par.tabs_total
    loadBar(tabi./plot_par.tabs_total , 'Plotting images...');

    % Current tab - its index in tab format
    tab_index = plot_par.tabs_index_table(tabi,:);
    tab_index_cell = num2cell(tab_index);
    
    % Image-index
    sli_img = tab_index_cell{1}; % Slice from data.

    % Image-data
    img2plot = plot_par.img(:,:,sli_img); 
    
    % Plot image
    if size(img2plot,1) > 1 % Safety against NaN-values
    imagesc(img2plot, 'parent', tgui.himg{tab_index_cell{:}}); 

    % Get image Contrast.
    if isfield(plot_par, 'contrast')
        contrast = plot_par.contrast;        
    else
        contr_min = min(img2plot(:)); contr_max = max(img2plot(:))*0.75;
        if contr_max <= contr_min, contr_max = contr_min+1; end 
        contrast = [contr_min contr_max];  plot_par.contrast = contrast;
    end
    % Set image Contrast (Matlab version specific)
    v = version('-release'); v(end) = []; v = str2double(v);
    if v < 2023, caxis(tgui.himg{tab_index_cell{:}}, contrast);
    else,        clim(tgui.himg{tab_index_cell{:}}, contrast);
    end
    % Colormap to gray
    colormap(tgui.himg{tab_index_cell{:}}, gray(255));
    
    end
end % End tab-loop
tgui.plot_par = plot_par; guidata(plot_par.fh, tgui); loadBar(NaN);

function VoxelMask_addGrid(target, target_sz, dim, range, grid_clr)

% Get figure size for normalization of grid thickness
w = target_sz(1); h = target_sz(2);
% Define using vertical line the thickness of horizontal ones 
wv = 2; wh = wv; wv = wv./w; wh = wh./h; 
% Line height or width for both vertical and horizontal grid lines;
h = 1; 

% Vertical lines (e.g. the x direction -> nr columns - 1)
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

function VoxelMask_addMenu(fh)
% Add menu bar with options to the figure

% Get GUI data of figure
tgui = guidata(fh);

% --- % File
tgui.mb.file.main = uimenu(fh, 'Label', 'File');
tgui.mb.file.export = uimenu(tgui.mb.file.main,'Label', 'Export', ...
                     'Callback', @VoxelMask_MB_export);
tgui.mb.file.import = uimenu(tgui.mb.file.main,'Label', 'Import', ...
                     'Callback', @VoxelMask_MB_import);
tgui.mb.file.close = uimenu(tgui.mb.file.main,'Label', 'Close', ...
                     'Callback', @VoxelMask_MB_close);
                 
% --- % Selection
tgui.mb.selection.main = uimenu(fh, 'Label', 'Selection');


% Selection > Auto
tgui.mb.selection.auto.main = ...
    uimenu(tgui.mb.selection.main, 'Label', 'Auto');
tgui.mb.selection.auto.all =  uimenu(tgui.mb.selection.auto.main,...
    'Label', 'All', 'Callback', @VoxelMask_MB_auto);                 
tgui.mb.selection.auto.slice =  uimenu(tgui.mb.selection.auto.main,...
    'Label', 'Slice', 'Callback', @VoxelMask_MB_auto);                 


% Selection > Copy
tgui.mb.selection.copy.main = ...
    uimenu(tgui.mb.selection.main, 'Label', 'Copy');
tgui.mb.selection.copy.all =  uimenu(tgui.mb.selection.copy.main,...
    'Label', 'Replace','Callback', @VoxelMask_MB_copy);            
tgui.mb.selection.copy.slice =  uimenu(tgui.mb.selection.copy.main,...
    'Label', 'Add', 'Callback', @VoxelMask_MB_copy);                    
tgui.mb.selection.copy.memory =  uimenu(tgui.mb.selection.copy.main,...
    'Label', 'to Memory', 'Callback', @VoxelMask_MB_copy);   
tgui.mb.selection.copy.paste =  uimenu(tgui.mb.selection.copy.main,...
    'Label', 'Paste', 'Callback', @VoxelMask_MB_copy); 

% Selection > Clear ...                 
tgui.mb.selection.clear.main = uimenu(tgui.mb.selection.main, ...
    'Label', 'Clear');         
% Selection > Clear > All  
tgui.mb.selection.clear.all = uimenu(tgui.mb.selection.clear.main, ...
    'Label', 'All', 'Callback', @VoxelMask_MB_clear);         
% Selection > Clear > Slice
tgui.mb.selection.clear.slice = uimenu(tgui.mb.selection.clear.main, ...
    'Label', 'Slice', 'Callback', @VoxelMask_MB_clear);                          
            
% Selection > Inverse ...
tgui.mb.selection.inverse.main = uimenu(tgui.mb.selection.main, ...
    'Label', 'Inverse');         
% Selection > Inverse > All  
tgui.mb.selection.inverse.all = uimenu(tgui.mb.selection.inverse.main, ...
    'Label', 'All', 'Callback', @VoxelMask_MB_inverse);         
% Selection > Inverse > Slice
tgui.mb.selection.inverse.slice = uimenu(tgui.mb.selection.inverse.main, ...
    'Label', 'Slice', 'Callback', @VoxelMask_MB_inverse);     

% --- % Image
tgui.mb.image.main = uimenu(fh, 'Label', 'Image');

% Image > Contrast ...
tgui.mb.image.contrast.main = uimenu(tgui.mb.image.main, ...
    'Label', 'Contrast', 'Callback', @VoxelMask_MB_contrast); 



% Update gui-data
guidata(fh, tgui);                 

function tgui = VoxelMask_PlotMask(fh)
% Plot a surface to use as a click-mask

% Get GUI data of figure
tgui = guidata(fh);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

if ~isfield(plot_par, 'mask')
    tgui.plot_par.mask = zeros(tgui.plot_par.dim);
end

for tabi = 1:plot_par.tabs_total 
    
    % Current tab - its index in tab format
    % This index lacks the slice+1 index (or the dim 1 above the slice
    % dimension). Therefor it requires correction regards pointing to data.
    tab_index = tgui.plot_par.tabs_index_table(tabi,:);
    tab_index_cell = num2cell(tab_index);
    
    % For data indexing - convert tab_index_cell to correct index-dims.
    sli = tab_index_cell{1}; % Slice from data.   
    
    % Plot data
    tmp_data = squeeze(tgui.plot_par.mask(:,:,sli))';
    tgui.plot_par.plotobj{tab_index_cell{:}} = ...
        imagesc(tgui.plot_par.ax{tab_index_cell{:}}, tmp_data);

    % Colormap and transparency
    set(tgui.plot_par.ax{tab_index_cell{:}},...
        'colormap', jet(128), 'clim', [0 1]);
    transmap = double((tmp_data) == 1); 
    transmap(transmap == 1) = tgui.plot_par.alpha;
    set(tgui.plot_par.plotobj{tab_index_cell{:}}, 'AlphaData', transmap)
    tgui.plot_par.transmap = transmap;

    % Axes cosmetics
    set(tgui.plot_par.ax{tab_index_cell{:}},...
        'LineWidth', 0.1, 'Xtick',[], 'Ytick', [],...
        'Box', 'off'); 

    % Turn of outer-axis i.e. data-imagesc
    tgui.plot_par.ax{tab_index_cell{:}}.Visible = 'off';    

end
guidata(fh, tgui);

function mouseHover(hobj,~)
% Handles figure-window mouse-position and value string plus tracking.

% Get data handle
tgui = guidata(hobj); 

% Find current tab
inp_selectedTab = repmat({tgui.tabg.SelectedTab}, size(tgui.tabh));
% Tab of interest
tab_nr = find(cellfun(@isequal, tgui.tabh, inp_selectedTab)==1);

% Current mouse-position
C = tgui.plot_par.ax{tab_nr}.CurrentPoint;

% Get value of mouse position
pos = round(C(1,1:2)); % X and Y (COL / ROW)

tag = ''; 
if isfield(tgui, 'tracking') && tgui.tracking
    tag = '| Tracking!'; 
end

hobj.Name = (sprintf('Voxelmask-Editor: %02i | %02i | %s', pos, tag));

% Manage tracking in different ways
if isfield(tgui, 'tracking') && tgui.tracking
    if isfield(tgui, 'list')
        
        if sum(~ismember(tgui.list, pos, 'row')) == size(tgui.list,1)...
                && ~tgui.tracking_type
            tgui = MaskTile_enable(tgui,pos);
            tgui.list(end+1,:) = pos;         
        else
            if max(ismember(tgui.list, pos, 'row') == 1)
                if MaskTile_query(tgui,pos)
                    MaskTile_enable(tgui,pos);
                else
                    MaskTile_disable(tgui,pos);
                end
            else
                tgui = MaskTile_toggle(tgui, pos);
                tgui.list(end+1,:) = pos; 
            end
            
            
        end
        guidata(hobj, tgui);
        tgui = VoxelMask_PlotMask(tgui.fig);                
                        
    else
        tgui.list = [];
        tgui.list(end+1,:) = pos;  
        
        tgui = MaskTile_toggle(tgui, pos); guidata(tgui.fig, tgui);
        tgui = VoxelMask_PlotMask(tgui.fig);        
    end
    
end

guidata(hobj, tgui);

function mouseClick(hobj, ~)
clicktype = get(hobj, 'selectiontype');
% 'normal' for left moue button
% 'alt' for right mouse button
% 'extend' for middle mouse button
% 'open' on double click

gui = guidata(hobj); 

% Find current tab
inp_selectedTab = repmat({gui.tabg.SelectedTab}, size(gui.tabh));
% Tab of interest
tab_nr = find(cellfun(@isequal, gui.tabh, inp_selectedTab)==1);

% Current mouse-position: X and Y (COL / ROW)
click_pos = round(gui.plot_par.ax{tab_nr}.CurrentPoint(1,1:2));       
if sum(click_pos <= 0) >=1 , return; end

switch clicktype
    case 'normal'
        if MaskTile_query(gui, click_pos) % isEnabled
            gui = MaskTile_disable(gui,click_pos);
        else % isDisabled
            gui = MaskTile_enable(gui,click_pos);
        end
        guidata(hobj, gui);
        VoxelMask_PlotMask(hobj);      
    case 'alt'       
        gui.tracking_type = 1;
        if ~isfield(gui, 'tracking') || gui.tracking == 0 
            gui.tracking = 1; % Enable tracking            
        else          
            gui.tracking = 0; % Disable tracking  
            if isfield(gui,'list'), gui = rmfield(gui,'list'); end
        end                   
        guidata(hobj, gui); 
        mouseHover(hobj);
   
    case 'extend'
        gui.tracking_type = 0;
        if ~isfield(gui, 'tracking') || gui.tracking == 0 
            gui.tracking = 1; % Enable tracking   
        else          
            gui.tracking = 0; % Disable tracking  
            if isfield(gui,'list'), gui = rmfield(gui,'list'); end
        end                   
        guidata(hobj, gui); 
        mouseHover(hobj);         
end

function gui = MaskTile_disable(gui, pos)
% Find current tab
inp_selectedTab = repmat({gui.tabg.SelectedTab}, size(gui.tabh));
% Tab of interest
tab_nr = find(cellfun(@isequal, gui.tabh, inp_selectedTab)==1);
% Change mask
if ~isempty(tab_nr) && (sum(pos > 0) == numel(pos))
    sz = size(gui.plot_par.mask);
    if ~(sum(pos(1:2) > sz(1:2)) > 0)
        gui.plot_par.mask(pos(1),pos(2), tab_nr) = 0;
    end
end

function gui = MaskTile_enable(gui, pos)
% Find current tab
inp_selectedTab = repmat({gui.tabg.SelectedTab}, size(gui.tabh));
% Tab of interest
tab_nr = find(cellfun(@isequal, gui.tabh, inp_selectedTab)==1);
% Change mask
if ~isempty(tab_nr) && (sum(pos > 0) == numel(pos))  
    sz = size(gui.plot_par.mask);
    if ~(sum(pos(1:2) > sz(1:2)) > 0)
        gui.plot_par.mask(pos(1),pos(2), tab_nr) = 1;
    end
end

function isEnabled = MaskTile_query(gui, pos)
% Find current tab
inp_selectedTab = repmat({gui.tabg.SelectedTab}, size(gui.tabh));
% Tab of interest
tab_nr = find(cellfun(@isequal, gui.tabh, inp_selectedTab)==1);
% Change mask
isEnabled = 0;
if ~isempty(tab_nr) && (sum(pos > 0) == numel(pos))
    isEnabled = gui.plot_par.mask(pos(1),pos(2), tab_nr);
end

function gui = MaskTile_toggle(gui, pos)
if MaskTile_query(gui, pos) % isEnabled
    gui = MaskTile_disable(gui,pos);
else % isDisabled
    gui = MaskTile_enable(gui,pos);
end
guidata(gui.fig, gui);

function closeReqFcn(hobj, ~, ~)
% Custom close-gui fcn.

try % FAST
    % Resume UI, see initialisation, to create output and close figure.
    tgui = guidata(hobj); uiresume(tgui.plot_par.fh);
catch % SLOW
   fobj = findobj('Type','Figure','Tag','Voxel Mask');
   if ~isempty(fobj), uiresume(fobj);  end
end


% ---------------------------------------------------------------------- %

% --- Executes by menubar: file > close
function VoxelMask_MB_close(hobj, ~, ~)
closeReqFcn(hobj);

% --- Executes by menubar: file > export
function VoxelMask_MB_export(hobj, ~, ~)
VoxelMask_Selection_Export(hobj);

% --- Executes by menubar: file > import
function VoxelMask_MB_import(hobj, ~, ~)
VoxelMask_Selection_Import(hobj);

% --- Executes by menubar: selection > copy
function VoxelMask_MB_copy(hobj, ~, ~)
VoxelMask_Selection_Copy(hobj);

% --- Executes by menubar: selection > clear
function VoxelMask_MB_clear(hobj,~,~)
VoxelMask_Selection_Clear(hobj);

% --- Executes by menubar: selection > auto
function VoxelMask_MB_auto(hobj, ~, ~)
VoxelMask_Selection_Auto(hobj);

% --- Executes by menubar: selection > Inverse
function VoxelMask_MB_inverse(hobj, ~, ~)
VoxelMask_Selection_Inverse(hobj);

% --- Executes by menubar: image > Contrast
function VoxelMask_MB_contrast(hobj, ~, ~)
VoxelMask_Contrast(hobj);

% ---------------------------------------------------------------------- %

function VoxelMask_Selection_Copy(hobj)
% Copy selection in current tab to all other tabs.

% Get gui-data
tgui = guidata(hobj); plot_par = tgui.plot_par;

% Use object text to replace or add "copy"-selection.
setMemory = 0; useMemory = 0; doReplace = 0;
switch hobj.Text
    case 'Replace',    doReplace = 1; 
    case 'Add',        doReplace = 0; 
    case 'to Memory',  setMemory = 1; 
    case 'Paste',      useMemory = 1;
end

% Find current tab-obj
inp_selectedTab = repmat({tgui.tabg.SelectedTab}, size(tgui.tabh));
% Tab of interest
toi = cellfun(@isequal, tgui.tabh, inp_selectedTab)==1;     

if ~useMemory
    % Get selected of current slice
    sel2copy = plot_par.mask(:,:,toi);    
end

% Store selection if copy-to-memory
if setMemory, tgui.copied = sel2copy; guidata(hobj, tgui); return; end

% Copy to all slices 
if ~setMemory && ~useMemory

    if doReplace
        % Copy selection to other tabs
        plot_par.mask = repmat(sel2copy, [1 1 size(plot_par.mask,3)]);
    else
        tmp = repmat(sel2copy, [1 1 size(plot_par.mask,3)]);
        plot_par.mask = arrayfun(@(x,y) x + y,tmp, plot_par.mask);
        plot_par.mask(plot_par.mask == 2) = 1;
    end

% Copy from memory
else
    % Use copy from memory
    sel2copy = tgui.copied;
    plot_par.mask(:,:,toi) = sel2copy;
end

% Clean up and plot
tgui.plot_par = plot_par; guidata(hobj, tgui);
VoxelMask_PlotMask(hobj);

function VoxelMask_Selection_Clear(hobj)
% Clear selection on all or current tab(s).

% Get gui-data
tgui = guidata(hobj); plot_par = tgui.plot_par;

% Use object text to clear.
switch hobj.Text
    case 'All'               
        % Set axes-background-color to none for all axes.
        plot_par.mask = zeros(size(plot_par.mask));
    case 'Slice'
        % Find current tab
        inp_selectedTab = repmat({tgui.tabg.SelectedTab}, size(tgui.tabh));
        % Tab of interest
        toi = find(cellfun(@isequal, tgui.tabh, inp_selectedTab)==1);
        % Clear data
        sz =  size(plot_par.mask);
        plot_par.mask(:,:,toi) = zeros(sz(1:2));
end

% Save changes.
tgui.plot_par = plot_par; guidata(hobj, tgui);
% Replot
VoxelMask_PlotMask(hobj);

function VoxelMask_Selection_Auto(hobj)
% Automatically select voxels on background of images.


% Get gui-data
tgui = guidata(hobj); plot_par = tgui.plot_par;


% All or current slice
switch hobj.Text
    case 'All'
        doCurrent = 0;
    case 'Slice'
        doCurrent = 1;    
        % Find current tab
        inp_selectedTab = repmat({tgui.tabg.SelectedTab}, size(tgui.tabh));
        % Tab of interest
        toi = find(cellfun(@isequal, tgui.tabh, inp_selectedTab)==1);
end


% If no image-data, return;
if ~plot_par.plot_img, return; end

% IMAGE PROCESSING % -------------------------------------------------- %

% Binary image - separating BG/FG
img = plot_par.img; cutoff_perc = 0.7;
% Cutoff value
[N, edges] = histcounts(plot_par.img(:), 256);
csum = cumsum(N)./sum(N);
[~, ind] = min(abs(csum - cutoff_perc)); cutoff = edges(ind);
% Binary image
img(img < cutoff) = 0; img(img >= cutoff) = max(img(:));
% display3D(img, 'tag', 'Binary IMG')

tmp = img; img(isnan(img)) = 0;

% Dilate and erode binary image
SE = strel('disk', 10); img = imdilate(img, SE); %img = imdilate(img, SE); 
SE = strel('disk', 10); img = imerode(img, SE); %img = imerode(img, SE); 
% display3D(img, 'tag', 'Dilated/Eroded')

% Mask the image - cut at its boundaries.
sz = size(img);
for si = 1:sz(3)
    ii = img(:,:,si); ii = imclose(ii, true(64)); ii = bwconvhull(ii);
    bw = bwboundaries(ii); 
    if ~isempty(bw)
    [~, mxind] = ...
        max(cell2mat(cellfun(@(x) numel(x), bw, 'UniformOutput', 0)));
    img(:,:,si) = poly2mask(bw{mxind}(:,2), bw{mxind}(:,1), sz(1), sz(2));
    end
end
% display3D(img, 'tag', 'Boundary mask');


% CUT-UP IMAGE %  ----------------------------------------------------- %

% Cut image into voxel-sized chunks: vectors for indexes to cut
% !! Here; row and col of plot_par.dim are swapped as with Matlab
% convention. !!
im_sz = size(plot_par.img); im_vec = cell(1,2);
im_sz_new = plot_par.dim([2 1]); 
for kk = 1:2
    im_vec{kk} = round(linspace(1, im_sz(kk), im_sz_new(kk)+1));
end

% Cut image up: swap of r/c wrt grid
im_avg = NaN([im_sz_new  plot_par.dim(3)]); 
im_med = im_avg; im_max= im_avg; im_std = im_avg;
for sli = 1:plot_par.dim(3)
    for ri = 1:im_sz_new(1)
        for ci = 1:im_sz_new(2)
            r = im_vec{1}(ri):(im_vec{1}(ri+1));
            c = im_vec{2}(ci):(im_vec{2}(ci+1));  
                  
            tmp = NaN(size(r,2), size(c,2));
            for rii = 1:size(r,2), tmp(rii,:) = img(r(rii), c, sli); end
                        
            im_avg(ri, ci, sli) = mean(tmp(:)); 
            im_med(ri, ci, sli) = median(tmp(:));
            im_std(ri, ci, sli) = std(tmp(:)); 
            im_max(ri, ci, sli) = max(tmp(:));

            clear tmp;
        end
    end
end
% display3D(im_avg, 'tag', 'Cut mean'); 
% display3D(im_med, 'tag', 'Cut median');
% display3D(im_std, 'tag', 'Cut std');


% Threshold new binary image
[N, edges] = histcounts(im_avg(:), 128); 
csum = cumsum(N)./sum(N); cutoff_perc = 0.2; 
[~, ind] = min(abs(csum - (cutoff_perc))); cutoff = edges(ind);
im_avg(im_avg <= cutoff) = 0; im_avg(im_avg > cutoff) = 1;

% Voxels to highlight
vox2light = (im_avg == 0);
vox2light = permute(vox2light,[2 1 3]);
% display3D(vox2light);

% If current slice only...
if doCurrent
    tmp = vox2light(:,:,toi); vox2light = plot_par.mask;
    vox2light(:,:,toi) = logical(tmp);
end

% Magic!
N = numel(size(plot_par.mask));
plot_par.mask = vox2light; % permute(vox2light,[N+1 1:N]);
tgui.plot_par = plot_par;
guidata(hobj, tgui);

% Update
VoxelMask_PlotMask(hobj); 

function VoxelMask_Selection_Inverse(hobj)
% Inverse selection of voxels

% Get gui-data
tgui = guidata(hobj); 

% Use object text to replace or add "copy"-selection.
switch hobj.Text
    case 'All', doSlice = 0;
    case 'Slice', doSlice = 1;
end

% Get current mask
selection = tgui.plot_par.mask;

% Create inverse
if doSlice
    % Find current tab-obj
    inp_selectedTab = repmat({tgui.tabg.SelectedTab}, size(tgui.tabh));
    % Tab of interest
    toi = cellfun(@isequal, tgui.tabh, inp_selectedTab)==1;       
    inversed = selection; inversed(:,:,toi) = ~inversed(:,:,toi);
else
    inversed = ~selection;
end

% Apply and update
tgui.plot_par.mask = inversed; guidata(hobj, tgui);

% Plot
VoxelMask_PlotMask(hobj);

function VoxelMask_Selection_Export(hobj)

gui = guidata(hobj);

% Get selection
mask = gui.plot_par.mask;

% Save selection
[fi, fp, idx] = uiputfile({'*.mat'}, 'Export mask-data', 'VoxelMask.mat');
if idx == 0, return; end

save([fp fi], 'mask');

function VoxelMask_Selection_Import(hobj)

% Get file selection
[fi, fp, idx] = uigetfile({'*.mat'}, 'Import mask-data');
if idx == 0, return; end

% Load mask
load([fp fi], 'mask');

% Get GUI data of figure
tgui = guidata(hobj); 

% Store mask
tgui.plot_par.mask = mask;

% Update gui-data
guidata(hobj, tgui);

% Plot mask
VoxelMask_PlotMask(hobj);


function VoxelMask_Contrast(hobj)

tgui = guidata(hobj); plot_par = tgui.plot_par;

% If no images need to be plotted, return.
if plot_par.plot_img == 0, return; end

% Get contrast from user
qry = {'Contrast:'}; def = {num2str(plot_par.contrast)}; 
type = {'edit'}; tstr = 'VoxelMask'; 
uans = getInput(type, qry, def, tstr);
if isempty(uans), return; end
plot_par.contrast = str2double(strsplit(uans{1}));
contrast = plot_par.contrast;

% Apply contrast to all tabs
for tabi = 1:plot_par.tabs_total
    % Current tab - its index in tab format
    tab_index = plot_par.tabs_index_table(tabi,:);
    tab_index_cell = num2cell(tab_index);
        
    % Set image Contrast (Matlab version specific)
    v = version('-release'); v(end) = []; v = str2double(v);
    if v < 2023, caxis(tgui.himg{tab_index_cell{:}}, contrast);
    else,        clim(tgui.himg{tab_index_cell{:}}, contrast);
    end
    % Colormap to gray
    colormap(tgui.himg{tab_index_cell{:}}, gray(255));
    
end % End tab-loop
tgui.plot_par = plot_par; guidata(plot_par.fh, tgui);


% ---------------------------------------------------------------------- %

function output = allCombinations(input)
% Create a combination of given input lists.
% Example: {[0,1],[0,1],[0,5]}
%
% Usefull to find combinations of indexes.

% List view
% input = flipud(input(:));

% Convert to grid matrix (cell)
c = cell(1, numel(input)); [c{:}] = ndgrid(input{:});
% Replace
output = cell2mat(cellfun(@(v)v(:), fliplr(c), 'UniformOutput',false));

% This creates equal output as combvec-fcn
output = fliplr(output)';
