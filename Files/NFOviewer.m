function NFOviewer(str)
% GUI to view a structure in a tree-format with ability to search for
% fields with a specific text.
%
% Dr. Q. van Houtum, quincyvanhoutum@gmail.com
% Version 1.0 - 04/2024


% Close if opened
fh = findall(groot, 'Type', 'figure', 'Name', 'Struct Viewer'); close(fh);

% Create figure
fh = uifigure('ToolBar','none','MenuBar','none','Name','Struct Viewer');
w = 480; h = 480; scrsz = get(0, 'ScreenSize');
pos_xy = round((scrsz(3:4)./2) - [w./2 h./2]);
fh.Position = [pos_xy w h];

% GUI-data
gui = guidata(fh);

% Create uitree
gui.tree = uitree(fh, 'Position',[0 2 round(w./2) h-2]);
matlab_version = version('-date'); 
matlab_version = str2double(matlab_version(end-3:end));
if matlab_version >= 2022
    gui.tree.ClickedFcn = @clickTree;
else
    gui.tree.SelectionChangedFcn = @clickTree;
end

% Create edit-box: search bar
editw = round(w./2); edith = 24; editx = round(w - editw); edity = h-edith;
gui.edit = uieditfield(fh, 'Position', [editx edity editw edith]);
gui.edit.ValueChangedFcn = @(src, evt) searchTree(src, evt);

% Create list-box: search results
lbw = editw; lbh = (h/2) - (2 * edith) - edith; 
lbx = editx; lby = h - edith - lbh;
gui.lboxresult = uilistbox(fh, 'Position',[lbx lby lbw lbh]);
gui.lboxresult.Items = {}; gui.lboxresult.Value = {};
gui.lboxresult.ValueChangedFcn = @(src, evt) searchResult(src, evt);

% Create label TITLE
lblx = editx; lbly = (edity - edith - lbh); lblw = editw; lblh = edith;
gui.lblval = uilabel(fh,'Position', [lblx lbly lblw lblh]);
gui.lblval.Text = ' Size | Bytes | Class'; 
gui.lblval.HorizontalAlignment = 'Left';

% Create label NFO
nfolblx = editx; nfolbly = lbly - lblh; 
gui.lblnfo = uilabel(fh,'Position', [nfolblx nfolbly lblw lblh]);
gui.lblnfo.Text = ' - | - | -'; 
gui.lblnfo.HorizontalAlignment = 'Left';

% Create list-box: tree-value
lbw = editw; lbh = h - edith - (2 * lblh) - lbh; 
lbx = editx; lby = gui.tree.Position(2);
gui.lboxvalue = uilistbox(fh, 'Position',[lbx lby lbw lbh]);
gui.lboxvalue.Items = {}; gui.lboxvalue.Value = {};

% Fill the tree
fillTree(str, gui.tree);

% Set structure as app-data
setappdata(fh, 'struct', str); guidata(fh, gui);

end

function tree = fillTree(str, tree)
% Fill up a uitree with all fieldnames in a structure, including any
% substructures. Lazy me - list of nodes grows within recursive loop.

fn = fieldnames(str); 
for fi = 1:numel(fn)
    % Add node to tree
    tmp = uitreenode(tree, 'Text', fn{fi});
    % Recursive over substructures
    if isstruct(str.(fn{fi})), fillTree(str.(fn{fi}), tmp); end
end

end

function clickTree(hobj, ~)
% Activate when user clicks node in tree

% Selected node
selNode = hobj.SelectedNodes;
% Node address in structure
fieldsToNode = getNodeFieldAddress(selNode);

% Get main structure and guidata
fh = getFigureHandle(hobj); str = getappdata(fh, 'struct');
gui = guidata(fh);

subfields = strsplit(fieldsToNode,'.'); val = str;
for si = 1:numel(subfields)
    val = val.(subfields{si});
end

% Set value
setValue(gui, val)


end

function setValue(gui, val)

% Value info string 
nfo = whos('val');
size_input_str = [ '[' repmat('%i ', 1, numel(nfo.size))];
size_input_str(end) = ']'; 
gui.lblnfo.Text = ...
    sprintf([size_input_str ' | %4i | %s'], nfo.size, nfo.bytes, nfo.class);

if ischar(val)
    gui.lboxvalue.Items = {val};
elseif isfloat(val)
    % Is a non-2D array
    isND = numel(size(val)) > 2;
    % Is it a list in row or col direction
    if ~isND
        if size(val,1) == 1 && size(val,2) > 1
            isList = 1; isList_row = 1;
        elseif size(val,2) == 1 && size(val,1) > 1
            isList = 1; isList_row = 0;
        else
            isList = 0;
        end
    else
        isList = 0;
    end
    
    % If column list, transpose
    if isList && ~isList_row, val = val'; end
    

    if isList
        val = num2cell(val);
        val = cellfun(@num2str, val, 'UniformOutput', false);
        gui.lboxvalue.Items = val;
    elseif isND || ~isList
        tmp = whos('val'); megabytes = tmp.bytes./1024.^2;
        if megabytes <= 50
            val = num2cell(reshape(val,[],1));
            val = cellfun(@num2str, val, 'UniformOutput', false);
            gui.lboxvalue.Items = val;
        else
            % Calculate statistics and convert to string-array for display.
            stat_str = val2statistics(val);

            % Display
            tmp = sprintf('Large data-array: %.2fMB', megabytes);
            stat_str = cat(1, tmp, {''}, stat_str);
            gui.lboxvalue.FontName = 'courier';
            gui.lboxvalue.Items = stat_str;
        end
    else        
        gui.lboxvalue.Items = {num2str(val(:))};
    end
elseif isstruct(val)
    gui.lboxvalue.Items = fieldnames(val);
elseif iscell(val)
    if numel(val{:}) <= 1
        setValue(gui, val{:});
    else
        bool_char = cellfun(@ischar, val);
        if sum(bool_char) == numel(val)
            gui.lboxvalue.Items = val;
        else
            disp tbd-non-char-cell
        end        
    end
else
    disp tbd-other-data-format
end


end

function stat_str = val2statistics(val)
% Calculate statistics of val (array) and create a string from the 
% resulting structure.   

% Get statistics
stats = csi_statistics_of_volume(val);

% Fieldnames
stat_names = fieldnames(stats);

% Values
stat_vals = cellfun(@getfield, ...
    repmat({stats},size(stat_names)), stat_names,...
    'UniformOutput', 0);

% Convert any cell in cell-values
isAcell = cellfun(@iscell, stat_vals);
tmp = cellfun(@cell2mat, ...
    stat_vals(isAcell),'UniformOutput',0);
stat_vals(isAcell) = tmp;

% Convert to string
stat_vals_str = cellfun(@num2str, stat_vals,...
    'UniformOutput', 0);

% Combine with alignment
stat_names_sz = cellfun(@numel, stat_names);
align_val = max(stat_names_sz)+4;
align_dff = align_val - stat_names_sz; 
stat_names_wspace = arrayfun(@(x) repmat(' ', 1, x), ...
    align_dff, 'UniformOutput', false);
stat_str = cellfun(@(x,y,z) [x ':' y z], ...
    stat_names, stat_names_wspace, stat_vals_str,...
    'UniformOutput', false);

end

function fh = getFigureHandle(hobj)
fh = hobj;
while ~strcmp(fh.Type, 'figure'), fh = fh.Parent; end
end

function [field_address] = getNodeFieldAddress(node)

field_address = node.Text; upNode = node.Parent;
while strcmpi(upNode.Type, 'uitreenode')
    field_address = strjoin({upNode.Text, field_address}, '.');
    upNode = upNode.Parent;
end

end

function searchTree(hobj, ~)

% Figure handle
fh = getFigureHandle(hobj);
% Guidata
gui = guidata(fh);
% Tree-struct
str = getappdata(fh, 'struct');
% Search query
qry = hobj.Value;
% Search tree-fields
flds = searchFields(str, qry, 0);
% Display in result-listbox
gui.lboxresult.Items = flds;

end

function searchResult(hobj, ~)

% Figure handle
fh = getFigureHandle(hobj);
% Guidata
gui = guidata(fh);
% Tree-struct
str = getappdata(fh, 'struct');
% Selected item
sel = gui.lboxresult.Value;

sel = strsplit(sel, '.'); tmp = str; noi = gui.tree;
for kk = 1:numel(sel)
    ind = strcmp({noi.Children.Text}, sel{kk});
    noi = noi.Children(ind);    
    tmp = tmp.(sel{kk}); 
end
% Select sel in uitree
gui.tree.SelectedNodes = noi;


% Set value in value-listbox
setValue(gui, tmp)


end

function flds = searchFields(struct_inp, query, case_sens)
% Find all fieldnames in a structure with multiple substructs and returning
% every fieldname containing the query.
%
% Input:
% struct_inp =  structure
% query =       string with (fieldname) string of interest.
% case_sens=    case sensitive search on [0] or off [1], default is off.
%
% Output:
% flds =        fields with "query" in it including any nested
%               substructures.
%
% quincyvanhoutum@gmail.com

if nargin < 3, case_sens = 1; end

% Fieldnames of structure
fn = fieldnames(struct_inp); 

% If any fieldnames contain query, add to output.
flds = fn(contains(fn, query,'IgnoreCase',case_sens)); 

% Loop each field to find nested structures and find query in those
% fieldnames.
for kk = 1:size(fn,1)    
    if isstruct( struct_inp.(fn{kk}) )
        % Field in struct_inp is a structure - read its fields
        tmp = searchFields(struct_inp.(fn{kk}), query, case_sens);
        
        if ~isempty(tmp)
            for ff = 1:size(tmp,1)
                flds = [flds; strjoin([fn{kk} tmp(ff)],'.')];
            end
        end
    end
end

end