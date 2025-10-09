function CSI_Log(varargin)
%%% Input: info-label + info-value as input.
%%%
%%% Updates listbox in main CSIgui window displaying csi-data information
%%% and processing steps.

obj = findobj('Tag','CSIgui_main', 'Type', 'Figure');
if isempty(obj), return; end

% Get GUIdata from CSIgui_main;
gui = guidata(obj);

% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi')
    csi_available = 0;
else
    csi_available = 1;
end

% Old string in CSI info listbox
str_old = get(gui.listbox_CSIinfo,'String');

% Check for lines to save
if isempty(str_old),                  str_save = {};
else
    if size(str_old,1) == 2 % Two lines: is there DIM-lines?
        if strfind(str_old{1},'Dim'), str_save = {};
        else,                         str_save = str_old;
        end
    elseif size(str_old,1) > 2 % More lines: is there DIM-lines?
        if strfind(str_old{1},'Dim'), str_save = str_old(3:end);
        else,                         str_save = str_old;
        end
    else
        str_save = str_old;
    end
end

% Process input info arguments
if nargin >= 1
    infolabel = varargin{1}; 
    if nargin < 2
        infoval = repmat({''}, 1,numel(infolabel)); 
    else
        infoval  = varargin{2}; 
    end

    % 1. PROCESS USER INPUT
    new_info = cell(1,size(infolabel,2));
    for kk = 1:size(infolabel,2)
       val = infoval{kk}; 
       if ~ischar(val), val = num2str(val); end % Convert to a double.
       new_info{kk} = ...
           sprintf('%s - %s %s', datetime('now','Format','HH:mm:ss'), ...
                                 infolabel{kk}, val);
    end
else
    new_info = {};
end

% 2. PROCESS STANDARD INFO
% Data size and labels.
if csi_available
    csi = getappdata(gui.CSIgui_main,'csi');
    
    % Get data dimensions
    sz = size(csi.data.raw); 

    % Get dimension labels
    if isfield(csi.data, 'labels'), lb = csi.data.labels;
    else,                           lb = 'Unknown';
    end

    % Standard first two lines: size and labels
    stand_info = { 
    ['Dim size  : '  sprintf('%s   ',num2str(sz))];
    ['Dim label : '  sprintf(repmat('%s | ', 1, size(lb,2)), lb{:})]};
else
    stand_info = {};
end


% 3. CREATE NEW
str_new = [stand_info(:)' new_info(:)' str_save(:)'];

try 
    % Print to info-box.
    set(gui.listbox_CSIinfo, 'String', str_new, 'Value', 1);
catch err
    % Warn if something goes wrong.
    warning('Error updating the CSIgui log listbox. Execution halted.');
    fprintf('ID: \n%s',err.identifier);
    fprintf('Message:\n%s',err.message);   
end

% Flush java-visuals
drawnow;