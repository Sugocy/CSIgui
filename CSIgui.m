function varargout = CSIgui(varargin)
% Spectroscopy GUI; initial purpose to merge MRI and MRS data in Matlab.
%
% Possible labels arguments:
%            'data','list','csi','spec','image', 'mrs', 'labels'
%            {filepath}, {filepathi}
% Input:
% CSIgui(datafield, label);
%
% UNDER DEVELOPMENT - 20181001

% Last Modified by GUIDE v2.5 09-Apr-2019 21:52:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CSIgui_OpeningFcn, ...
                   'gui_OutputFcn',  @CSIgui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    % [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{1}); Ori
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before CSIgui is made visible.
function CSIgui_OpeningFcn(hObject, eventdata, gui, varargin)
% This function has no output args, see OutputFcn.
% Time ID is created here. Defines GUI its own ID.

 
%% ROOT CSIGUI

% Add CSIgui root and Files-folder to search path Matlab.
fpn_csigui = mfilename('fullpath'); froot = fileparts([fpn_csigui '.m']);

% Add CSIgui path to Matlab path-root.
if exist([froot '\Files'],'dir') == 7
    addpath(froot); addpath([froot '\Files']); savepath;
else
    % Disable this when compiling to executable. 
    warndlg('Missing directory "File" in CSIgui root.', 'CSIgui - error');
end


%% CSIGUI VERSION

% Choose default command line output for CSIgui
gui.output = hObject;
% This field create ID for this CSIgui-instance.
gui.ID = sprintf('%.0f', now.*10^10);
% Define CSIgui version here.
gui.version = '2.0';



%% Add menu-bar

loadBar(0.1, 'Creating menubar...');
gui = CSIgui_setMenuBar(gui);

loadBar(0.3, 'Executing JavaScript...');
gui = CSIgui_setMenuBar_JavaScript(gui);



%% Process user-input
loadBar(0.6, 'Processing input...');

% Following setup applies:  CSIgui(value, label);
% 
% Example: CSIgui(data, 'csi');
% Supported: 'csi', 'mri', 'spec', 'data', 'list', 
%            'filepath', 'filepathi','labels', 'fp', 'fpi'
% The filepathi is meant for images. Other for CSI data.

% Struct for input.
inp = struct; label_list = {'data','list','csi','spec','image','labels',...
                            'mrs', 'filepath', 'filepathi', 'fp','fpi'};
for inpi = 2:2:size(varargin,2)
    % Find corresponding label for accepted input labels
    lab_ind = strcmpi(label_list, varargin{inpi}); 
    if sum(lab_ind) == 0
        fprintf('Input label %s is not a supported input argument!\n',...
                varargin{inpi});
    else
        % If abbreviation is used, set to full input label.
        if strcmp(varargin{inpi}, 'fp'),  varargin{inpi} = 'filepath'; end
        if strcmp(varargin{inpi}, 'fpi'), varargin{inpi} = 'filepathi'; end
        
        % Save input argument as field in inp-struct;
        inp.(lower(varargin{inpi})) = varargin{inpi-1};       
    end
end

% Add input structure to GUI-data if some input is indeed added!
if ~isempty(fieldnames(inp))
    gui.inp = inp;                       % Add to GUI-data structure
    guidata(hObject, gui);               % Save updated GUI-data to figure
    drawnow;
    openFile(hObject, eventdata, gui);   % Instantly process in openFile.
    % Delete processed user input from memory
    gui = guidata(hObject);
    delete_UserInput(hObject, gui);
else
    % Update gui structure
    guidata(hObject, gui);
end

% Run GUI coloring theme function
setGUIcolor(hObject); % This stores the right color palet in gui-handle

% Close menubar
loadBar(1, 'Launching GUI...'); pause(0.0005); loadBar(NaN);

% --- Executes at launch
function gui = CSIgui_setMenuBar(gui)
% Set CSIgui menubar items and corresponding callbacks
% See gui.menubar for all available items.

if ~isfield(gui,'menubar')
    
% 1.File % ------------------------------------------------------------- %
gui.menubar.file.main = uimenu(gui.CSIgui_main,'Label','File');

% 1.File > Open..
gui.menubar.file.open = ...
    uimenu(gui.menubar.file.main,'Label', 'Open', 'Enable', 'on',...
    'Callback', @openFile);

% 1.File > Close.
gui.menubar.file.close = ...
    uimenu(gui.menubar.file.main, 'Separator','on','Label', 'Exit',...
    'Enable', 'on','Callback', @CSIgui_main_CloseFromMenu);            


% 2.MRSI % ------------------------------------------------------------- %
gui.menubar.MRSI.main = uimenu(gui.CSIgui_main, 'Label', 'MRSI');

% 2.MRSI > Undo
gui.menubar.MRSI.undo = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Undo', 'Enable', 'on',...
    'Accelerator', 'z', 'Callback', @CSIgui_Undo);

% 2.MRSI > Plot
gui.menubar.MRSI.plot = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Plot', 'Enable', 'on',...
    'Callback', @CSI_2D_initiate2D);

% 2.MRSI > Operators
gui.menubar.MRSI.operations.main = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Operators', 'Enable', 'on');

% 2.MRSI > Operations >  Multiply
gui.menubar.MRSI.operations.mul = ...
    uimenu(gui.menubar.MRSI.operations.main,'Label', 'Multiply',...
    'Enable', 'on','Callback', @button_CSI_Multiply_Callback);

% 2.MRSI > Operations >  Divide
gui.menubar.MRSI.operations.div = ...
    uimenu(gui.menubar.MRSI.operations.main,'Label', 'Divide',...
    'Enable', 'on','Callback', @button_CSI_Divide_Callback);

% 2.MRSI > Operations >  Sum
gui.menubar.MRSI.operations.sum = ...
    uimenu(gui.menubar.MRSI.operations.main,'Label', 'Sum',...
    'Enable', 'on','Callback', @button_CSI_Sum_Callback);

% 2.MRSI > Operations >  Mean
gui.menubar.MRSI.operations.avg = ...
    uimenu(gui.menubar.MRSI.operations.main, 'Label', 'Average',...
    'Enable', 'on','Callback', @button_CSI_Average_Callback);


% 2.MRSI > Color Scaling (2D)
gui.menubar.MRSI.ColorScale.main = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Color Scaling', 'Enable', 'On');

% 2.MRSI > Color Scaling (2D) > Slice
gui.menubar.MRSI.ColorScale.Slice = ...
    uimenu(gui.menubar.MRSI.ColorScale.main,'Label', 'Slice',...
    'Check', 'Off','Callback', @CSI_2D_Scaling_Color_Set, 'Enable', 'On');

% 2.MRSI > Color Scaling (2D) > Volume
gui.menubar.MRSI.ColorScale.Volume = ...
    uimenu(gui.menubar.MRSI.ColorScale.main,'Label', 'Volume',...
    'Check', 'On','Callback', @CSI_2D_Scaling_Color_Set, 'Enable', 'On');

% 2.MRSI > Color Scaling (2D) > Static
gui.menubar.MRSI.ColorScale.Static = ...
    uimenu(gui.menubar.MRSI.ColorScale.main, 'Label', 'Static',...
    'Check', 'Off','Callback', @CSI_2D_Scaling_Color_Set, 'Enable', 'On');

% 2. MRSI > Color Scaling (2D) > Colorbar
gui.menubar.MRSI.ColorScale.Colorbar = ...
    uimenu(gui.menubar.MRSI.ColorScale.main, 'Label', 'Colorbar',...
    'Check', 'Off','Separator','on',...
    'Callback',@CSI_2D_Scaling_plotColorbar, 'Enable', 'On');

% 2.MRSI > Axis Scaling (2D)
gui.menubar.MRSI.AxisScale.main = ...
    uimenu(gui.menubar.MRSI.main, 'Label', 'Axis Scaling', 'Enable', 'On');

% 2.MRSI > Axis Scaling (2D) > Slice
gui.menubar.MRSI.AxisScale.Voxel = ...
    uimenu(gui.menubar.MRSI.AxisScale.main,'Label', 'Voxel',...
    'Check', 'Off','Callback', @CSI_2D_Scaling_Axis_Set, 'Enable', 'On');

% 2.MRSI > Axis Scaling (2D) > Volume
gui.menubar.MRSI.AxisScale.Slice = ...
    uimenu(gui.menubar.MRSI.AxisScale.main,'Label', 'Slice',...
    'Check', 'On', 'Callback', @CSI_2D_Scaling_Axis_Set, 'Enable', 'On');

% 2.MRSI > Axis Scaling (2D) > Volume
gui.menubar.MRSI.AxisScale.Volume = ...
    uimenu(gui.menubar.MRSI.AxisScale.main,'Label', 'Volume',...
    'Check', 'Off', 'Callback', @CSI_2D_Scaling_Axis_Set, 'Enable', 'On');
        
% 2.MRSI > Set Domain
gui.menubar.MRSI.domain.main = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Set Domain', 'Enable', 'on');
% 2.MRSI > Set Domain > Time
gui.menubar.MRSI.domain.time = ...
    uimenu(gui.menubar.MRSI.domain.main,'Label', 'Time', ...
    'Enable', 'on','Callback', @CSI_setDomain);
% 2.MRSI > Set Domain > Frequency
gui.menubar.MRSI.domain.frequency = ...
    uimenu(gui.menubar.MRSI.domain.main,'Label', 'Frequency',...
    'Enable', 'on','Callback', @CSI_setDomain);



% 2.MRSI > Noise
gui.menubar.MRSI.noise = ...
    uimenu(gui.menubar.MRSI.main, 'Label', 'Noise','Enable', 'on',...
    'Callback',@CSI_viewNoise);        
        
% 2.MRSI > Save > ...
gui.menubar.MRSI.export.main = ...
    uimenu(gui.menubar.MRSI.main,'Separator','on','Label', 'Export',...
    'Enable', 'on');        

% 2.MRSI > Save > SDAT

gui.menubar.MRSI.export.Data = ...
    uimenu(gui.menubar.MRSI.export.main,'Separator','off',...
    'Label', 'Data', 'Enable', 'on','Callback', @CSI_saveData);  

% 2.MRSI > Save > Figure 2D Plot
gui.menubar.MRSI.export.Fig = ...
    uimenu(gui.menubar.MRSI.export.main,'Separator','off',...
    'Label', 'Figure', 'Enable', 'on', 'Callback', @CSI_saveFig);  
        
% 2.MRSI > Save > LOG (CSI_info)
gui.menubar.MRSI.export.Log = ...
    uimenu(gui.menubar.MRSI.export.main,'Separator','on','Label','Log',...
    'Enable', 'on', 'Callback', @CSI_saveLog);          

% 3.Image % ------------------------------------------------------------ %
gui.menubar.Image.main = uimenu(gui.CSIgui_main, 'Label', 'Image');

% 3. Image > Plot
gui.menubar.Image.plot = ...
    uimenu(gui.menubar.Image.main, 'Label', 'Plot','Enable', 'on',...
    'Callback', @button_plotIMG_Callback);

gui.menubar.Image.export.main = ...
    uimenu(gui.menubar.Image.main, 'Label', 'Export','Enable', 'on',...
    'Separator','On');

gui.menubar.Image.export.Images = ...
    uimenu(gui.menubar.Image.export.main, 'Label', 'Images',...
    'Enable', 'on', 'Callback', @MRI_saveIMG);

% 4.View % ------------------------------------------------------------- %
gui.menubar.View.main = uimenu(gui.CSIgui_main, 'Label', 'View');

% 4.View > Theme
gui.menubar.View.theme.main = ...
    uimenu(gui.menubar.View.main,'Label', 'Theme', 'Enable', 'On');

% 4.View > Theme > Day
gui.menubar.View.theme.day = ...
    uimenu(gui.menubar.View.theme.main,'Label', 'Day','Check', 'Off',...
        'Callback', @setGUIcolor_bymenu, 'Enable', 'On');
    
% 4.View > Theme > Night
gui.menubar.View.theme.night = ...
    uimenu(gui.menubar.View.theme.main,'Label', 'Night','Check', 'On',...
        'Callback', @setGUIcolor_bymenu, 'Enable', 'On');
    
% 4.View > Theme > Custom
gui.menubar.View.theme.custom = ...
    uimenu(gui.menubar.View.theme.main,'Label', 'Custom','Check', 'Off',...
        'Callback', @setGUIcolor_bymenu, 'Enable', 'On');    
    
% 4.View > Set Custom
gui.menubar.View.set_custom = ...
    uimenu(gui.menubar.View.main,'Label', 'Set Custom',...
        'Callback', @setGUIcolor_custom, 'Enable', 'On');      

% 5.Help % ------------------------------------------------------------- %
gui.menubar.About.main = uimenu(gui.CSIgui_main, 'Label', 'Help');

% 5.About > Help
gui.menubar.About.help = ...
    uimenu(gui.menubar.About.main,'Label', 'Help',...
        'Callback', @openHelp, 'Enable', 'On');
% 5.About > About
gui.menubar.About.about = ...
    uimenu(gui.menubar.About.main,'Label', 'About', 'Separator','on',...
        'Callback', @openAbout, 'Enable', 'On');

end

% --- Executes at launch
function gui = CSIgui_setMenuBar_JavaScript(gui)
% Sets tooltips to menu-bar entries.

% ----------- %%% Some Jave Menu ToolTip String Magic %%% -------------- %
try
    % Flush java.
    drawnow; pause(0.01); 
    % Surpress obsolete function message
    warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    
    % Get figure java frame
    jFrame = get(handle(gui.CSIgui_main),'JavaFrame');
    % Get menubar component
    jMenuBar = jFrame.fHG2Client.getMenuBar; 
    
    % Get CSI menubar; java starts at zero, 3-1 = 2;
    jMRSI = jMenuBar.getComponent(1); % CSI menu 
    
    % Check if right component
    if strcmp(jMRSI.getText,'MRSI')
        % Click once to open the menu 
        jMRSI.doClick;  pause(0.01); 
        
        menu_name = {'Color Scaling', 'Axis Scaling','Noise','Set Domain'};
        menu_tip = {...
         ['Set spectra color scaling in 2D MRSI plots relative to'... 
          'maximum per slice, volume or as static e.g. single color.'],...
         ['Set spectra y-axis scaling in 2D MRSI plots to maximum '...
          'per voxel, slice or volume.'],...
         ['View noise data from list/data file. Use the show '...
          'CSI-button or plot in menu > CSI to display. ' ...
          'To revert back to the original data, click here again.'], ...
          'Set data domain of MRSI data to frequency or time.'};
        
        nMenus = size(fieldnames(gui.menubar.MRSI),1)-1; % Remove main CSI
        for kk = 0:nMenus-1 % Start indexin @ zero - Matlab vs Java.
            jtmp = jMRSI.getMenuComponent(kk);
            jtmp_str = jtmp.toString;
            dot_loc = strfind(jtmp_str, '.');
            % 6    12    53    68  ... % seperator
            % 4    14    17    22  ... % MenuItem
            
            % Safety to skip specific menu items without getText such
            % as seperators.
            if dot_loc(1) ~= 6
                jtmp_name = jtmp.getText;
                loc = contains(menu_name, char(jtmp_name));
                ind = find(loc == 1);
                if ~isempty(ind)
                    jtmp.setToolTipText(menu_tip{ind});
                end
            end
          
        end    
    end
    
catch err
    fprintf('Jave Menubar setToolTipText error: %s', err.message);
    save([datestr(now,'HHMMSS') '_JavaError.mat'],  'err');
end

% --- Outputs from this function are returned to the command line.
function varargout = CSIgui_OutputFcn(~, ~, gui)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% gui    structure with gui and user data (see GUIDATA)

% Get default command line output from gui structure
varargout{1} = gui.output;

% --- Executes on button press in button_DeleteInput.
function button_DeleteInput_Callback(hObj, ~, gui)
% Delete initial user input.
delete_UserInput(hObj, gui);

% --- Delete user input
function delete_UserInput(hObj, gui)
% Delete the userinput field from GUI structure - Only required during
% launch of CSIgui.
if isfield(gui,'inp')
    gui = rmfield(gui, 'inp'); guidata(hObj,gui); 
end

% --- Open file or process input.
function openFile(hObj, ~, ~)
% 1. Open up selection UI or analyze input arguments
% 2. Open file according file-extension or user input
% 3. CREATE APPDATA:    CSI: CSI and SPECTRO data (3D/2D/1D)
%                       MRI: MR Image data        (3D/2D/1D)
%
% If user sends 2 filepath inputs; both will be loaded from here.
% Userinput is deleted at the end of this script.

% Get fresh gui-data.
gui = guidata(hObj); clear_log = 1;

% FILE INPUT + USER INPUT % -------------------------------------------- %
if ~isfield(gui, 'inp')
    
    % Get default path if available
    if isappdata(gui.CSIgui_main,'csi')
        csi = getappdata(gui.CSIgui_main,'csi');
        if isfield(csi,'filepath'), fp = csi.filepath; else, fp = []; end
    else, fp = []; 
    end
    
    % Select-file UI.
    [fn, fp, fi] = uigetfile({'*.*','All Files (*.*)';...
     '*.list;*.data;*.LIST;*.DATA',...
                'RAW MRS files (*.list, *.data)';...
     '*.dcm;*.DCM;*.par;*.rec;*.PAR;*.REC',...
                'Image files (*.dcm, *.par, *.rec)';...
     '*.spar;*.SPAR;*.sdat;*.SDAT',...
                'MRS files (*.sdat, *.spar)';...
     '*.txt;*.TXT;',...
                'Text files (*.txt)'},...
     'Select a file', fp); 
    % Canceled file selection
    if fi == 0, return; end 
    
    % Get file-extension for further processing, see below.
    [fp, fn, ext] = fileparts([fp, fn]); ext = lower(ext);
elseif isfield(gui,'openThisFile')
    % Get FP/FN/EXT from openThisFile field.
    fp = gui.openThisFile.fp;
    fn = gui.openThisFile.fn; 
    ext = lower(gui.openThisFile.ext);
    clear_log = gui.openThisFile.clear_log;
else
    % Input from the user is present.
    % If given a filepath or filepathi, processed here:
    %       Will be set to gui.openThisFile to process.
    % Other labels are processed at the ext = userinput code below.
    
    % 1. Filepath MRS + filepath MRI
    if isfield(gui.inp, 'filepath')  && isfield(gui.inp, 'filepathi') 
        
        % A. Process MRI and MRS filepaths
        [fp{1}, fn{1}, ext{1}] = fileparts(gui.inp.filepath{1}); 
        [fp{2}, fn{2}, ext{2}] = fileparts(gui.inp.filepathi{1}); 
        % Loop the MRS and MRI filepaths
        for fi = 1:2
            
            % Add to gui structure.
            gui.openThisFile.fp = fp{fi};  gui.openThisFile.fn = fn{fi};
            gui.openThisFile.ext= ext{fi}; 
            % Clear log once.
            if     fi == 1, gui.openThisFile.clear_log = 1;
            elseif fi == 2, gui.openThisFile.clear_log = 0;
            end
            guidata(hObj,gui);
            
            % Open this file.
            openFile(hObj,[], []);   
        end
        
        % Clear the openThisFile field in GUI and save.
        gui = rmfield(gui, 'openThisFile'); guidata(hObj,gui);
        
        % Stop the iteration of openFile();
        return; 
        
    % 2. Filepath MRS only.
    elseif isfield(gui.inp, 'filepath') 
        % Further processing as if user opened a file using the select file 
        % UI 
        [fp, fn, ext] = fileparts(gui.inp.filepath{1});
        ext = lower(ext);
        
    % 3. Filepath MRI only.   
    elseif isfield(gui.inp, 'filepathi') 
        % Further processing as if user opened a file using the select file 
        % UI 
        [fp, fn, ext] = fileparts(gui.inp.filepathi{1});
        ext = lower(ext);
        
     % 4. Data array, data struct or image array - etc.
    else   
        ext = 'userinput'; 
        % See last elseif statement below processing different file types
    end
end


                       % --- % Parse the input % --- %

if clear_log
    set(gui.listbox_CSIinfo, 'String', {}); % Clear info in listbox.    
end

% Update LOG
CSI_Log({['Loading ' ext(2:end) ' file.']}, {'Please wait.'});

if strcmp(ext,'.dcm')                                          % DICOM

    % Parse dicom file
    success = parse_dicom(fp, fn, gui);

elseif strcmp(ext,'.par') || strcmp(ext,'.rec')                % PAR/REC

    % Parse par/rec file
    success = parse_parrec(fp, fn, gui);

elseif strcmp(ext,'.list') || strcmp(ext,'.data')              % LIST/DATA
    
    % Parse the list/data file
    success = parse_listdata(fp, fn, gui);
    % Reorder dimensions of data
    CSI_ReorderDim_Auto(gui);
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj, [], 1); % Automatic
     
elseif strcmp(ext,'.mat')                                      % MAT
    
    % Parse mat-file
    success = parse_mat(fp, fn, gui);
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);
    
elseif strcmpi(ext,'.sdat') || strcmpi(ext,'.spar')            % SDAT/SPAR

    % Parse the sdat/spar file
    success = parse_sdatspar(fp, fn, gui);
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);

elseif strcmpi(ext, '.txt')                                    % TEXT
    
    % Parse the text file
    success = parse_text(fp, fn, gui);
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);
   
elseif strcmpi(ext,'userinput')                                % USER
    
    % Parse the userinput
    success = parse_userinput(gui);

else                                                           % ERROR 
    % Update LOG and return
    CSI_Log({'Warning. File format not supported.'},{ext});
    warning('The selected file format is not supported.'); return;
end

% Update LOG
if success
    CSI_Log({['Loading ' ext(2:end) ' file succeeded.']},{''});
else
    CSI_Log({['Loading ' ext(2:end) ' file failed.']},{''});
end

% Delete any user input from app-memory
delete_UserInput(hObj, gui);

% Fresh gui data
gui = guidata(hObj);    

                           % --- GUI updating --- %
              
maxL = 40; % Max name-length
if isappdata(gui.CSIgui_main,'csi') 
    csi = getappdata(gui.CSIgui_main,'csi');
    if isfield(csi, 'filename')
        name = csi.filename;                             % Get filename
        if size(name,2) > maxL, name = name(1:maxL); end % Max length name
        set(gui.txt_fnCSI, 'String', name);              % Set in GUI
    end
end
if isappdata(gui.CSIgui_main,'mri')
    mri = getappdata(gui.CSIgui_main,'mri');
    if isfield(mri, 'filename')
        name = mri.filename;                             % Get filename
        if size(name,2) > maxL, name = name(1:maxL); end % Max length name
        set(gui.txt_fnIMG, 'String', name);              % Set in GUI
    end
end

% --- Parse LIST/DATA file
function success = parse_listdata(fp, fn, gui)
% Load list & data file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

% Load list/data file into memory.
[csi.data, csi.list] = csi_loadData([fp '\' fn], 0);

% If error while loading.
if isempty(csi.data) 
    CSI_Log({'Loading LIST/DATA file failed.'},{''});
    success = 0; return; 
else
    success = 1;
end

% Add dimension info.
csi.data.dim = size(csi.data.raw);
% Add filename and extension info.
csi.ext = '.data'; csi.filename = fn; csi.filepath = fp;

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi);   

% --- Parse SDAT/SPAR file
function success = parse_sdatspar(fp, fn, gui)
% Load spar & sdat file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

% Read SDAT (Own function!)
[csi.data.raw, csi.list] = csi_readSDAT([fp '\' fn]);

% Check if loading succeeded
if isnan(csi.data.raw)
    success = 0; 
    CSI_Log({'Loading SDAT/SPAR file failed.'},{''}); 
    return;
else
    success = 1;
end

% Add dimension info.
csi.data.dim = size(csi.data.raw);
% Add filename extension info.
csi.ext = 'sdat'; csi.filename = fn; % Save filename
% Add dimension labels.
% Create labels 
csi.data.labels = csi.list.dim_labels;

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi);  

% --- Parse TEXT file
function success = parse_text(fp, fn, gui)
% Load text file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

% Load textfile
[csi.data.raw] = csi_readText([fp '\' fn '.txt']);

% Check if loading succeeded.
if isnan(csi.data.raw)
    CSI_Log({'Loading text-file failed.'},{''}); 
    success = 0; return; 
else
    success = 1;
end

% Add dimension info.
csi.data.dim = size(csi.data.raw);
% Add filename extension info.
csi.ext = 'txt'; csi.filename = fn; % Save filename

% Add dimension labels: Create labels 
for di = 1:numel(csi.data.dim)
    csi.data.labels{di} = [char(96+di) char(96+di)]; 
end

 % Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi);   

% --- Parse MAT file
function success = parse_mat(fp, fn, gui)
% Load mat file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

% Expected fields of interest
% foi = {'data','ext','filename','filepath',...
%        'xaxis','csi','raw','dim','labels', 'conv','mri'};

% Check mat-file integrety   
mat_cont = whos('-file',[fp '\' fn '.mat']);
if ~strcmp(mat_cont(1).name,'csigui')
    success = 0;
    CSI_Log(...
    {'Incorrect mat-file. Use by CSIgui generated mat-file.',...
     'Expected fields:'},...
    {'Required structure: csigui.',...
     'raw, dim, filepath and name, noise, split, conv, mri,'}); return;
else, success = 1; 
end

% Integrety verified, read matfile
inp = load([fp '\' fn '.mat'], 'csigui');
csigui = inp.csigui; clear inp;

              % ------- % Process MRI struct % ------- %

% Find MRI struct and store             
if isfield(csigui,'mri')
    mri = csigui.mri; csigui = rmfield(csigui,'mri');
    setappdata(gui.CSIgui_main,'mri',mri); 
end


              % ------- % Process CONV struct % ------- %

% Find conv struct and store
if isfield(csigui,'conv')
    converted = csigui.conv; csigui = rmfield(csigui,'conv');
    setappdata(gui.CSIgui_main,'conv',converted); 
end

              % ------- % Process CSI struct % ------- %

% Set remaining CSI input to structure;
csi.data = csigui;

% Set extensions
csi.ext = 'mat';

% Set mat-filename
csi.data.filename = fn; csi.data.filepath = fp;
csi.filename = fn; csi.filepath = fp;

% Find xaxis stuct and store
if isfield(csigui, 'xaxis')
    csi.xaxis = csigui.xaxis; csi.data = rmfield(csi.data,'xaxis');
end

if isfield(csigui, 'ori')
    csi.ori = csigui.ori; csi.data = rmfield(csi.data,'ori');
end
                 % ------- % Process LOG % ------- %

% Resubmit log if available
if isfield(csigui,'log')
    CSI_Log({'Copying log from mat-file.'},{'See below.'})
    log_data = cellstr(csigui.log);
    % Exclude first two lines
    log_data = log_data(3:end);
    
    % Add to log listbox
    gui.listbox_CSIinfo.String = ...
    cat(1, gui.listbox_CSIinfo.String,{''} ,log_data);

    % Remove log-field
    csi.data = rmfield(csi.data,'log');
end

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi); 
      
% --- Parse USERINPUT
function succes = parse_userinput(gui)
% Analyse userinput file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct and/or mri-struct

succes = 0;

% Check if userinput is available.
if ~isfield(gui,'inp'), CSI_Log({'No userinput available.'},{''}); return; end

                           % SOME INFORMATION %
% Expected userinput labels: 
% 'data','list','csi','spec','image','mrs','labels'
% Filepath(i) labels are already excluded from processing here; 

% Label DATA % ----------------------------------------------------- %
% If data-struct is given (from csi_loadData)
if isfield(gui.inp,'data')
    userInfo = 'Data-struct loaded.';

    try        
        csi.data = gui.inp.data;             % Set data struct       
        csi.data.dim = size(csi.data.raw);   % Add dimension info.
        csi.ext = '.data'; csi.filename = data.filename; % File nfo

        % Save CSI data in app-data
        setappdata(gui.CSIgui_main,'csi',csi);   

    catch err
        CSI_Log({'Failed processing user input!'},{err.message});   
        succes = 0; return;
    end

    % Update process info for user.
    CSI_Log({'User input processed: '},{userInfo});

    % Calculate xaxis from available data
    CSI_2D_Scaling_calc_xaxis(hObj,[],1); 
    
    % Parse succes
    succes = 1;
end

% Label LIST % ----------------------------------------------------- %
% If list-struct is given (from csi_loadList/csi_loadData)
if isfield(gui.inp,'list')
    userInfo = 'List-struct loaded.';

    % Set input.
    csi.list = gui.inp.list; csi.filename = csi.list.filename;
    % Save as appdata
    setappdata(gui.CSIgui_main,'csi', csi);   

    % Update process info for user.
    CSI_Log({'User input processed: '},{userInfo});
    
    % Parse succes
    succes = 1;
end

% Label CSI/SPEC % ------------------------------------------------- %
% If csi, spec or mrs is given (spectrum array)
if isfield(gui.inp,'csi') || ...
        isfield(gui.inp,'spec') || ...
                isfield(gui.inp, 'mrs')
    
    % Info for LOG
    userInfo = 'MRSI data loaded.';

    % Save the appropriate input data
    if     isfield(gui.inp,'mrs'),  csi.data.raw = gui.inp.mrs;
    elseif isfield(gui.inp,'csi'),  csi.data.raw = gui.inp.csi;
    elseif isfield(gui.inp,'spec'), csi.data.raw = gui.inp.spec;
    end

    % Add dimension info.
    csi.data.dim = size(csi.data.raw);
    % Create labels: if labels are given as input this will be removed
    % automatically.
    for di = 1:numel(csi.data.dim)
        csi.data.labels{di} = [char(96+di) char(96+di)]; 
    end

    % Add extension info.
    csi.ext = '.spec_only'; csi.filename = 'User input';

    % Save CSI data in app-data
    setappdata(gui.CSIgui_main,'csi',csi);   

    % Clear user info.
    set(gui.listbox_CSIinfo, 'String', {});

    % Update process info for user.
    CSI_Log({'User input processed: '},{userInfo});

    % Calculate xaxis from available data
    CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);
    
    % Parse succes
    succes = 1;
end

% Label LABELS % --------------------------------------------------- %
% Add dimension labels to the data structure.
if isfield(gui.inp,'labels')
    userInfo = 'Userinput: Saved labels for each index.';

    % Add the labels to the csi.data.labels field. Other random labels
    % will be discarded: see spec/csi as input labels.
    csi.data.labels = gui.inp.labels;

    % Save CSI data in app-data
    setappdata(gui.CSIgui_main,'csi',csi);   

    % Update processed info to user.
    CSI_Log({'User input processed: '},{userInfo});
    
    % Pase succes
    succes = 1;
end

% --- Parse DICOM file
function success = parse_dicom(fp, fn, gui)
% Analyse dicom file, parse the data and store it into CSIgui appdata.
%
% Created: mri-struct

% Read and order dicom data
[m,i,g] = dicomread7T({[fp '\' fn '.dcm']}); mri = struct;
if isnan(m)
    CSI_Log({'Loading DICOM file failed.'},{''}); 
    success = 0; return; 
else
    success = 1;
end

% Order dicom data
[mri.data, mri.par] = dicomorder3(m, i); 

% If scout is detected
% cell_ind = find(cellfun(@isempty, b.M) == 0);
mri.examinfo = g; mri.ext = 'dcm'; 
mri.filename = fn; mri.filepath = fp;


% Save MRI data in app-data
setappdata(gui.CSIgui_main,'mri',mri);   

% --- Parse PARREC file
function success = parse_parrec(fp, fn, gui)
% Analyse parrec file, parse the data and store it into CSIgui appdata.
%
% Created: mri-struct

% Load parrec.
[m,p,o,s,e] = readparrec4([fp '\' fn]); 

if isempty(m) && isempty(p) && isempty(o)
    success = 0; CSI_Log({'Loading PAR/REC file failed.'},{''});
    return; 
else
    success = 1;
end

% Create mri-struct and store image data
mri = struct; mri.examinfo = e; mri.par = s; mri.data.M = m; 
mri.data.P = p; mri.data.P = o;  mri.filename = fn; % Save filename

% Due par/rec - rotate 90deg
fina = fieldnames(mri.data);
for permi = size(fina,2)
    mri.data.(fina{permi}) = permute(mri.data.(fina{permi}),[2 1 3 4 5 6 ]); 
end

mri.ext = '.par';

% Save MRI data in app-data
setappdata(gui.CSIgui_main,'mri',mri);   

% --- Parse TEXT Protocol file
function success = parse_protocolText(fp, fn, gui)


% Read the text protocol file
nfo = csi_readTextProtocol([fp '\' fn]);
if ~isstruct(nfo) && isnan(nfo), success = 0; return; end

% Get CSI data
if ~isappdata(gui.CSIgui_main,'csi'), success = 0; return; end
csi = getappdata(gui.CSIgui_main,'csi');

              % --------- %  Process NFO % --------- % 

% NFO field                                 CSIgui variabel 
% Nucleus                                   Nucleus
% SpectralBW_Hz                             BW
% VoxelSizeRL_mm VoxAP_mm VoxFH_mm          Resolution
% SliceOffc_AP_P__mm RL_L__mm FH_H__mm'};   Offcenter
       
% Copy some frequency parameters

% Covert nucleus 
nucleus = cat(2,...
regexp(nfo.Nucleus,'\d*','Match'), regexp(nfo.Nucleus,'[A-Z]','Match'));
% Save nucleus and BW to axis info
csi.xaxis.nucleus = cat(2,nucleus{:});
csi.xaxis.BW = nfo.SpectralBW_Hz;

% Res AP
if isfield(nfo,'VoxAP_mm')
    csi.ori.res(1) = nfo.VoxAP_mm;
    csi.ori.offcenter(1) = [nfo.SliceOffc_AP_P__mm]; 
end

% Res RL
if isfield(nfo,'VoxRL_mm') 
   csi.ori.res(2) = nfo.VoxRL_mm;     
   csi.ori.offcenter(2) = [nfo.RL_L__mm]; 
end
if isfield(nfo, 'VoxelSizeRL_mm')
   csi.ori.res(2) = nfo.VoxelSizeRL_mm;
end

% Res FH
if isfield(nfo,'VoxFH_mm')
   csi.ori.res(3) = nfo.VoxFH_mm;
   csi.ori.offcenter(3) = [nfo.FH_H__mm]; 
end       
    


                % --------- % Clean Up % --------- %

% Save data
setappdata(gui.CSIgui_main,'csi',csi); 

% Calculate xaxis
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main);

% Calculate coordinates
button_CSI_setCoordinates_Callback([], [], gui);

success = 1;

% --- Parse SPAR Protocol file
function success = parse_protocolSpar(fp, fn, gui)


% Read the text protocol file
nfo = csi_readSPAR([fp '\' fn]);
if ~isstruct(nfo) && isnan(nfo), success = 0; return; end

% Get CSI data
if ~isappdata(gui.CSIgui_main,'csi'), success = 0; return; end
csi = getappdata(gui.CSIgui_main,'csi');

              % --------- %  Process NFO % --------- % 

% NFO field                                 CSIgui variabel 
% Nucleus                                   Nucleus
% SpectralBW_Hz                             BW
% VoxelSizeRL_mm VoxAP_mm VoxFH_mm          Resolution
% SliceOffc_AP_P__mm RL_L__mm FH_H__mm'};   Offcenter
       
% Copy some frequency parameters

% Covert nucleus 
% nucleus = 

% Save nucleus and BW to axis info
csi.xaxis.nucleus = nfo.nucleus;
csi.xaxis.BW = nfo.sample_frequency;

% Copy some coordnate parameters
% csi.ori.res = [nfo.slice_distance nfo.RL_L__mm nfo.slice_distance];
% csi.ori.offcenter = ...
%     [nfo.si_ap_off_center nfo.si_lr_off_center nfo.si_cc_off_center];



                % --------- % Clean Up % --------- %

% Save data
setappdata(gui.CSIgui_main,'csi',csi); 

% Calculate xaxis
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main);

% Calculate coordinates
button_CSI_setCoordinates_Callback([], [], gui);

success = 1;

% --- Executes on button press (+) in button_CSI_ReadInfo.
function button_CSI_ReadInfo_Callback(~, ~, gui)
% Read the MRS protocol text file and set specific MRS parameters


% Get app-data
csi = getappdata(gui.CSIgui_main,'csi'); 
if isempty(csi), CSI_Log({'No MRSI data loaded.'},{''}); return; end


                % --------- % Userinput % --------- %
                

% Get default file path if available
if isappdata(gui.CSIgui_main,'csi')
    csi = getappdata(gui.CSIgui_main,'csi');
    if isfield(csi.data,'filepath'), fp = csi.data.filepath; else, fp = []; end
else, fp = []; 
end

% Select-file UI.
[fn, fp, fi] = uigetfile(...
    {'*.txt','Text files (*.txt)';'*.spar','Spar files (*.spar)'},...
    'Select a file', fp); 
% Canceled file selection
if fi == 0, return; end 

% Get file-extension for further processing, see below.
[fp, fn, ext] = fileparts([fp, fn]); ext = lower(ext);

                % --------- %  Read file % --------- % 
          
switch ext
    case '.txt'
        success = parse_protocolText(fp, [fn ext], gui);
    case '.spar'
        success = parse_protocolSpar(fp, [fn ext], gui);
end


                 % --------- % Clean up % --------- %        
                  
if success, CSI_Log({'Read the protocol nfo file.'},{''});
else,       CSI_Log({'Failed reading protocol nfo file.'},{''});
end
                    

% Help and About % ---------------------------------------------------- %
% --------------------------------------------------------------------- %

% --- Executes by user click in menubar
function openHelp(~,~, ~)
% Open the CSIgui.pdf which includes function information and a help
% section.
pdf_fp = 'CSIgui - Help.pdf';

if exist(pdf_fp, 'file')
    winopen(pdf_fp); % Open using default PDF reader.
else 
    CSI_Log({'Help file not found.'},...
            {'Please check root directory of CSIgui'})
end

% --- Executes by user click in menubar
function openAbout(hObject,~ ,~)
% Open a small about window, showing miscellaneous GUI information.
gui = guidata(hObject);

% GUI Colors.
clr_bg = gui.colors.main; clr_tx = gui.colors.text_main;
clr_hl = gui.colors.hilight1;

% Load about.txt via loadAbout function.
about_text  = loadAbout(); about_text = sprintf(about_text);
about_title = 'About CSIgui';

% Figure: CSIgui About
aboutfig = figure('Tag', 'CSIgui_About',...
               'Name', ['CSIgui v' gui.version ],...
               'Color', clr_bg, 'MenuBar','none', 'Toolbar', 'none', ...
               'resize', 'off','numbertitle', 'off'); 

% Get PC screensize for panel size
scrsz = get(0,'Screensize'); 
% Width, height and x and y position on screen
w = 480; h = 280; 
x = ceil( ( scrsz(3)/2 ) - (0.5*w) ); % Center position of screen.
y = ceil( ( scrsz(4)/2 ) - (0.5*h) );
% Set figure position
set(aboutfig, 'position', [x y w h],'Unit','Normalized');

 uicontrol(aboutfig, 'Style','Text','Unit', 'Normalized',...
          'Position', [0.1 0.9 0.5 0.05], 'String', about_title,...
          'ForegroundColor', clr_hl, 'BackgroundColor', clr_bg,...
          'HorizontalAlignment','Left','FontWeight','Bold');

 uicontrol(aboutfig, 'Style', 'Listbox', 'Unit', 'Normalized',...
        'position', [0.1 0.1 0.8 0.8], 'String', about_text,...
        'ForegroundColor', clr_tx, 'BackgroundColor', clr_bg,...
        'HorizontalAlignment','Left');
            

    
    
% getUserInput EDITS GUI % --------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI to get user input
function userInput = getUserInput(question, defans, clrs)
% Opens up a menu for user to input data.
% question: cell with each required user input
% defans  : default answers.
%
% If user canceled, userInput is returned empty.


% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name','CSIgui - UserInput','Tag', 'CSIgui_UI_edits'); 
subdat      = guidata(fig_userinp); axis off;

% Number of questions to ask user.
nrofquest = size(question,2);
% Total size of figure.
dw = 10;  dy = 10; w = 280+2*dw; 
% Define values for EDIT and TEXT sizes.
sz_edit = [w-2*dw 25]; sz_text = [w-2*dw 20]; sz_button = [(w-2*dw)/2 20];
% Define height
h = nrofquest * (sz_edit(2)*2) + (sz_button(2)*2) + 2*dy;

% Size tune parameters and set figure size and window-position
scrsz = get(0,'Screensize'); figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

for kk = 1:nrofquest
    % Question.
    subdat.h_query{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Text','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_edit(2)*kk + sz_text(2)*(kk-1)) ...
                    sz_text(1) sz_text(2)],...
        'Tag', 'h_edit','String',question{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor',  clrs(1,:),...
        'HorizontalAlignment', 'Left');
    
    % Default answer.
    if ~ischar(defans{kk}), def_ans_tmp = num2str(defans{kk});
    else, def_ans_tmp = defans{kk};
    end
    subdat.h_edit{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Edit', 'Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_edit(2)*kk + sz_text(2)*(kk))...
                    sz_edit(1) sz_edit(2)],...
        'Tag', 'h_edit','Fontsize', 8,...
        'String',def_ans_tmp,...
        'Foregroundcolor',  clrs(2,:),'BackgroundColor', clrs(1,:),...
        'KeyPressFcn', @getUserInput_editKeyPress);
end

% Set button CONTINUE
subdat.h_buttonContinue = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw        dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_CONTINUE','String','Continue','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_setOutput);

% Set button SKIP
subdat.h_buttonSkip = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw+sz_button(1)  dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_SKIP','String','Skip','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_setOutput);

% --------------------- % Wait for user input

% Update gui-data of this function' figure.
guidata(fig_userinp, subdat); 
% Wait for user t2bdone
uiwait(fig_userinp); pause(0.1);

% --------------------- % Process user input

% Get updated gui-data
if ishandle(fig_userinp)
    subdat = guidata(fig_userinp);
    if isfield(subdat, 'ans')
        userInput = subdat.ans;
    else
        userInput = [];
    end
    % Close figure
    close(fig_userinp);
else
    userInput = [];
end

% -. Executed by button press in edits
function getUserInput_editKeyPress(H, E) 
drawnow;
if strcmp(E.Key, 'return')
    getUserInput_setOutput(H, E);    
else
    return;
end

% -. Executed by button in getUserInput to save answer
function getUserInput_setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

%userInputSource = get(hObject,'Style');
userInputString = get(hObject,'String'); % What did get us here? 

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','CSIgui_UI_edits');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

if strcmp(userInputString, 'Skip')
    val = {};
else
    % Get strings in edit-boxes!
    val = cell(1,size(subdat.h_edit,2));
    for kk = 1:size(subdat.h_edit,2)
        val{kk} = get(subdat.h_edit{kk},'String');
    end
end
% Set userinput value and Update storage
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);

% getUserInput: POPUP GUI % -------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI with popup menu's to get user input
function userInput = getUserInput_Popup(popup_title, popup_input, clrs)
% Create a simple gui with dropdown e.g. popup menu to choose certain 
% options.


% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name','CSIgui - UserInput','Tag', 'CSIgui_UI_popup'); 
subdat = guidata(fig_userinp); axis off;

% Nr of dropdown menu's
Ninput = size(popup_title,2);

% Total size of figure.
dw = 10;  dy = 10; w = 280+2*dw; 
% Define values for EDIT and TEXT sizes.
sz_popup  = [w-2*dw 25]; sz_text = [w-2*dw 20]; 
sz_button = [(w-2*dw)/2 20];

% Define height
h = Ninput * (sz_popup(2)*2) + (sz_button(2)*2) + 2*dy;

% Size tune parameters and set figure size and window-position
scrsz = get(0,'Screensize'); 
figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

% Set each title text and popup menu
for kk = 1:Ninput
    %Popup menu title
    subdat.h_title{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Text','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_popup(2)*kk + sz_text(2)*(kk-1)) ...
                    sz_text(1) sz_text(2)],...
        'Tag', 'h_popup','String',popup_title{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left');
    
    % Dropdown menu title    
    subdat.h_popup{kk} = ...
        uicontrol(fig_userinp, 'Style', 'popupmenu','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_popup(2)*kk + sz_text(2)*(kk))...
                    sz_popup(1) sz_popup(2)],...
        'Tag', 'h_popup','String',popup_input{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left');
   
end

% Set button CONTINUE
subdat.h_buttonContinue = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw        dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_CONTINUE','String','Continue','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Popup_setOutput);

% Set button SKIP
subdat.h_buttonSkip = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw+sz_button(1)  dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_SKIP','String','Skip','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Popup_setOutput);

% --------------------- % Wait for user input

% Update gui-data of this function' figure.
guidata(fig_userinp, subdat); 
% Wait for user t2bdone
uiwait(fig_userinp); pause(0.1);

% --------------------- % Process user input

% Get updated gui-data
if ishandle(fig_userinp)
    subdat = guidata(fig_userinp);
    if isfield(subdat, 'ans')
        userInput = subdat.ans;
    else  
        userInput = [];
    end
    % Close figure
    close(fig_userinp);
else
    userInput = [];
end

% -. Executed by button in getUserInput_Popup to save answer
function getUserInput_Popup_setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

%userInputSource = get(hObject,'Style');
userInputString = get(hObject,'String'); % What did get us here? 

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','CSIgui_UI_popup');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

if strcmp(userInputString, 'Skip')
    val = {};
else
    % Get strings input from popup menu
    val = cell(1,size(subdat.h_popup,2));
    for kk = 1:size(subdat.h_popup,2)
        pval = get(subdat.h_popup{kk},'Value');
        pstr = get(subdat.h_popup{kk},'String');
        val{kk} = pstr{pval};
    end
end
% Set userinput value and Update storage
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);


% getUserinput: BUTTON GUI % ------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI with popup menu's to get user input
function userInput = getUserInput_Buttons(qst, button_text, clrs)
% Create a simple gui with two buttons to choose between two options.
% qst:          One strin/char array for the question.
% Button text:  Two string/char arrays for both buttons.

% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name','CSIgui - UserInput','Tag', 'CSIgui_UI_button'); 
subdat = guidata(fig_userinp); axis off;

% Total size of figure.
dw = 10;  dy = 10; w = 280+2*dw; 
% Define sizes for BUTTON and TEXT.
sz_button = [(w-2*dw)/2 20]; sz_text = [w-2*dw 20]; 
% Define figure height.
h = (sz_button(2)*2) + 2*dy;

% Set figure size and window-position
scrsz = get(0,'Screensize'); figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

% Set question
 subdat.h_title = ...
   uicontrol(fig_userinp, 'Style', 'Text','Unit', 'Pixels',...
   'Position',[dw h/2 sz_text(1) sz_text(2)],...
   'Tag', 'h_title','String', qst, 'Fontsize', 8,...
   'Foregroundcolor',clrs(2,:),'BackgroundColor',clrs(1,:),...
   'HorizontalAlignment', 'Left');

% Set button A
subdat.h_buttonA = ...
  uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
  'Position',[dw dy  sz_button(1) sz_button(2) ],...
  'Tag', 'buttonA','String',button_text{1},'Fontsize', 12,...
  'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
  'Callback', @getUserInput_Buttons_setOutput);

% Set button B
subdat.h_buttonB = ...
  uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
  'Position',[dw+sz_button(1)  dy  sz_button(1) sz_button(2) ],...
  'Tag', 'buttonB','String',button_text{2},'Fontsize', 12,...
  'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
  'Callback', @getUserInput_Buttons_setOutput);

% --------------------- % Wait for user input

% Update gui-data of this function its figure.
guidata(fig_userinp, subdat); 
% Wait for user t2bdone
uiwait(fig_userinp); pause(0.1);

% --------------------- % Process user input

% Get updated gui-data and extract the answer
%if exist(fig_userinp,
if ishandle(fig_userinp)
    subdat = guidata(fig_userinp);
    if isfield(subdat, 'ans')
        userInput = subdat.ans;
    else
        userInput = [];
    end
    
    % Close figure
    close(fig_userinp);
    
else
    userInput = [];
end

% -. Executed by button in getUserInput_Popup to save answer
function getUserInput_Buttons_setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

%userInputSource = get(hObject,'Style');
userInputString = get(hObject,'String'); % What did get us here? 

% Set output value e.g. string of button
val = userInputString;

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','CSIgui_UI_button');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

% Set userinput value and update gui-data
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);


% getUserinput: RADIO GUI % -------------------------------------------- %
% ---------------------------------------------------------------------- %

% --- % Executed to get userinput from radiobuttons
function userInput = getUserInput_Radio(question, defans, clrs)
% Opens up a menu for user to select radio buttons
% question: cell with each required user input
% defans  : default answers: on and off (1/0)
%
% If user canceled, userInput is returned empty.


% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name','CSIgui - UserInput','Tag', 'CSIgui_UI_radio'); 
subdat      = guidata(fig_userinp); axis off;

% Number of questions to ask user.
nrofquest = size(question,2);
% Total size of figure.
dw = 15;  dy = 10; w = 280+2*dw; 
% Define values for EDIT and TEXT sizes.
sz_radio = [w-2*dw 25]; sz_button = [(w-2*dw)/2 20];
% Define height
h = nrofquest * (sz_radio(2)) + (sz_button(2)*2) + 2*dy;

% Size tune parameters and set figure size and window-position
scrsz = get(0,'Screensize'); figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

for kk = 1:nrofquest
    % Radio button with default answer
    subdat.h_radio{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Radiobutton','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_radio(2)*kk) ...
                    sz_radio],...
        'Tag', 'h_edit','String',question{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor',  clrs(1,:),...
        'HorizontalAlignment', 'Center','Value', defans(kk));
end

% Set button CONTINUE
subdat.h_buttonContinue = ...
    uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw        dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_CONTINUE','String','Continue','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Radio_setOutput);

% Set button SKIP
subdat.h_buttonSkip = ...
    uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw+sz_button(1)  dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_SKIP','String','Skip','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Radio_setOutput);

% --------------------- % Wait for user input

% Update gui-data of this function' figure.
guidata(fig_userinp, subdat); 
% Wait for user t2bdone
uiwait(fig_userinp); pause(0.1);

% --------------------- % Process user input

% Get updated gui-data
if ishandle(fig_userinp)
    subdat = guidata(fig_userinp);
    if isfield(subdat, 'ans')
        userInput = subdat.ans;
    else
        userInput = [];
    end
    % Close figure
    close(fig_userinp);
else
    userInput = [];
end

% -. Executed by button in getUserInput_radio to save answer
function getUserInput_Radio_setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

%userInputSource = get(hObject,'Style');
userInputString = get(hObject,'String'); % What did get us here? 

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','CSIgui_UI_radio');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

if strcmp(userInputString, 'Skip'), val = [];
else
    % Get strings in edit-boxes!
    val = NaN(1,size(subdat.h_radio,2));
    for kk = 1:size(subdat.h_radio,2)
        val(kk) = get(subdat.h_radio{kk},'Value');
    end
end
% Set userinput value and Update storage
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);




% CSI BUTTONS % -------------------------------------------------------- %
% ---------------------------------------------------------------------- %

% --- Executes on button press in button_CSI_Average: 
function button_CSI_Average_Callback(~, ~, gui, backup)
% Average over a specific dimension in MRS data. If no aver or nsa index is
% present, user input is required.
% This functions aims to average over dimension "aver" or "nsa".
%
% Uses CSI_average();
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % ------------------------------- %
% Create backup
if backup, CSI_backupSet(gui, 'Before averaging.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% AVERAGE % ---------------------------------------- %

% Get dimension to average from labels..
if isfield(csi.data, 'labels')
    nsa_ind = csi_findDimLabel(csi.data.labels,{'aver','nsa'});
    nsa_ind(isnan(nsa_ind)) = [];
else
    nsa_ind = [];
end
    
% Average
[data_averaged, nsa_ind] = CSI_average(csi.data.raw,nsa_ind);
if isempty(data_averaged)
    CSI_Log({'No index to average.'},{'Cancelled averaging.'}); return; 
end

% Replace data with averaged and set averaged index size to one.
csi.data.raw = data_averaged; csi.data.dim(nsa_ind) = 1;

% CLEAN UP % ---------------------------------------- %

% Store appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% Update CSI info panel
if nsa_ind <= size(csi.data.labels,2)
    CSI_Log({'Averaged over dimensions'},...
                   {strjoin(csi.data.labels(nsa_ind), ' ')});
else
    CSI_Log({'Averaged over dimensions'},{nsa_ind});
end

% --- Executes on button press in button_CSI_Mean: 
function button_CSI_Mean_Callback(~, ~, gui)
% Calculcate the mean over a specific dimension given by user.
% This functions allows manual averaging of MRSI data.
%
% Uses CSI_average();

% BACKUP + APPDATA % ------------------- %

% Create backup
CSI_backupSet(gui, 'Before apodization.');

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Mean % ------------------------------- %

% Ask user which dimension
index = getUserInput({'Enter the index to average over:'},{'4'});
if isempty(index), CSI_Log({'Skipped averaging.'},{''}); return; end
index = str2double(index{1});

% Calculate mean.
[data_averaged, index] = CSI_average(csi.data.raw,index);

% Average data
csi.data.raw = data_averaged;
% Change size of averaged index
csi.data.dim(index) = 1;

% CLEAN UP % --------------------------- % 
% Store appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% Update CSI info panel
if index <= size(csi.data.labels,2)
    CSI_Log({'Averaged MRSI data over dimension'},...
        {strjoin(csi.data.labels(index), ' ')});
else
    CSI_Log({'Averaged MRSI data over dimension'},{index});
end

% --- Executed by button CSI_average_Callback.
function [data_averaged, index] = CSI_average(data,index)
% Average data over specific index. If the index not given or empty 
% empty, the user is asked to enter an index e.g. dimension to average 
% over.

if nargin == 1 || isempty(index) % If not or an empty index is given
    index = getUserInput({'Enter index to average over:'},{'4'});
    if isempty(index), data_averaged = []; return; end
    index = str2double(index{1});
end

% Average
data_averaged = mean(data,index);



% --- Executes on button press in button_CSI_FFT_Kspace.
function button_CSI_FFT_Kspace_Callback(~, ~, gui, backup)
% Apply spatial FFT over MRSI data. If no spatial labeled dimensions e.g.
% indexes are available user input is required.
%
% Uses csi_rawfft;
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % ------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before spatial FFT.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% GET OPTION INPUT % ------------------- %

% quest = {'Circular shift (0) or Fourier shift (1)'};
% defan = {{'0','1'}};
% uans = getUserInput_Popup(quest,defan);
% if isempty(uans), return; end
% 
% % Set option
% shift_opt = str2double(uans{1});
shift_opt = 2; % SET TO AUTOMATIC

% K-space index: find the spatial dimensions/indexes: ask if not found.
spat_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
spat_dim(isnan(spat_dim)) = [];
if isempty(spat_dim)
    spat_dim = getUserInput({'Spatial index in MRSI data? (kx/ky/kz): '},...
                            {'2 3 4'});
    if isempty(spat_dim), CSI_Log({'Skipped FFT.'},{''}); return; end
    spat_dim = str2double(strsplit(spat_dim{1},' '));
end

% SPATIAL FFT % ------------------------ %

% Raw data fourier - CSI specific. 1D K-space data does NOT require this. 
% If FID to SPectra is required, use csi_fft
csi.data.raw = csi_rawfft(csi.data.raw, spat_dim, shift_opt);
% shift_opt

% CLEAN UP % -------------------------- %

% Store appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% After update info

% Spatial info to user.
if spat_dim <= size(csi.data.labels,2) 
    CSI_Log({'Applied FFT over spatial dimension(s)'},...
                    {strjoin(csi.data.labels(spat_dim), ' | ')});
else
    CSI_Log({'Applied FFT over spatial dimension(s) '},{spat_dim});
end

% Shift method to user.
if shift_opt == 1
    tmp = 'FFT shift'; 
elseif shift_opt == 2
    tmp = 'Automatic (fft- circular- for odd/even dimensions)';
else
    tmp = 'Circular shift';
end
CSI_Log({'Used FFT shift method: '}, {tmp});

% --- Executes on button press in button_CSI_FFT.
function button_CSI_FFT_Callback(~, ~, gui, backup)
% Apply forward fourier to transform the MRSI data from the spatial domain
% to the frequency domain.
% 
% Uses csi_fft();
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % ------------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before forward FFT.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% FFT % -------------------------------------------- %

% Get user input
uans = getUserInput_Popup(...
           {'Correct for N term?','Shift before and after FFT? (Echo)'},...
           {{'No','Yes'},{'No','Yes'}});
if isempty(uans), CSI_Log({'Skipped FFT.'},{''}); return; end

if strcmpi(uans{1},'yes'),correct_N = 1; else, correct_N = 0; end
if strcmpi(uans{2},'yes'),dbl_shift = 1; else, dbl_shift = 0; end

% Apply FFT
csi.data.raw = csi_fft(csi.data.raw,correct_N,dbl_shift);

% CLEANUP % ---------------------------------------- %

% Save appdata
setappdata(gui.CSIgui_main, 'csi', csi);
% After update info
CSI_Log({'Forward FFT applied.'},{''});

% --- Executes on button press in button_CSI_iFFT.
function button_CSI_iFFT_Callback(~, ~, gui, backup)
% Apply inverse fourier to transform MRS data from the frequency domain to
% the spatial domain.
%
% Uses csi_ifft();
if nargin < 4 , backup= 1; end

% BACKUP + APPDATA % ------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before inverse FFT.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% iFFT % ------------------------------- %

% Apply iFFT
csi.data.raw = csi_ifft(csi.data.raw);

% CLEAN UP % --------------------------- %

% Store data
setappdata(gui.CSIgui_main, 'csi', csi);

% After update info
CSI_Log({'Inverse FFT applied.'},{''});



% --- Executes on button press in button_CSI_Apodization_Kspace.
function button_CSI_Apodization_Kspace_Callback(~, ~, gui, backup)
% Apodization: Apply a 3D filter over spatial dimensions (indexes) 
%              of the data. 
%
% Uses csi_apodization3D();
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % --------------------------------- %
% Create backup
if backup, CSI_backupSet(gui, 'Before apodization.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USERINPUT % ---------------------------------------- %

% If required e.g. no spatial dimensions are found.
kspc_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
kspc_dim(isnan(kspc_dim)) = [];
if isempty(kspc_dim)
    kspc_dim = getUserInput({'Apply apodization over index dimensions: '},...
                            {'2 3 4'});
    if isempty(kspc_dim), return; end
    kspc_dim = str2double(strsplit(kspc_dim{1},' '));
end   

% Ask window type
uans = getUserInput_Popup({'Apodization filter-window:',...
                           'Show filter?'},...
                          {{'Hamming', 'Hann', 'Blackman','Flattop'},...
                           {'No','Yes'}});
if isempty(uans), CSI_Log({'Skipped Apodization.'},{''}); return; end

% Selected type
type = uans{1};   

% Show filter yes or no?
switch uans{2}
    case 'Yes', show_filter = 1;
    case 'No',  show_filter = 0;
end

% APODIZATION 3D % ---------------------------------- %

% Apply apodization
[csi.data.raw, win] = csi_apodization3D(csi.data.raw, kspc_dim, type);


% Show filter
if show_filter
    % Get data
    win_4plot = permute(win,[3 2 4 1 5 6 7 8]);
    % Display the 3D array.
    display3D(win_4plot,'colormap','jet','limit',[0 1]);    
end
    

% CLEAN UP % ---------------------------------------- % 
% Store appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% Update CSI info panel
if kspc_dim <= size(csi.data.labels,2)
    CSI_Log({['Apodization (' type ') applied over dimensions']},...
        {strjoin(csi.data.labels(kspc_dim), ' | ')});
else
    CSI_Log({['Apodization (' type ') applied over dimensions']},...
        {kspc_dim});
end
 
% --- Executes on button press in button_CSI_Apodization_FID.
function button_CSI_Apodization_FID_Callback(~, ~, gui, backup)
% Apply apodization to all spectra. Multiple methods are available:
% Hamming, Hanning, Blackman, Exponential, Gaussian.
% gaussmf(x, [std center]);
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % --------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before 1D apodization.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USER INPUT % ---------------------------------- %

% Get filter to apply.
quest = {'Filter algorithm:', 'FID or Echo:'};
defan = {{'Gaussian', 'Hamming','Hann', 'Exponential','Blackman',...
          'Flattop', 'Sinebell'}, {'FID','Echo'}};
uanstype = getUserInput_Popup(quest,defan);
if isempty(uanstype), CSI_Log({'Skipped apodization.'},{''}); return; end

% Get options for specific algorithms
switch uanstype{1}
    case 'Gaussian'
        % ((1/(1/BW x nSamples) * (ans/nSamples)) / pi) == (bw/ans)/pi 
        % Hz = bw / (N*pi) && N = (Hz * Pi)/BW
        if isfield(csi.xaxis,'BW')
            uopts = getUserInput({'Apodization factor: (Hz)'},{20});            
            hz = str2double(uopts{1});
            uopts{1} = num2str(csi.xaxis.BW ./ (hz * pi));
        else
            uopts = getUserInput(...
            {'Standard deviation e.g. half length at FWHM: (Samples)'},...
            {(csi.xaxis.N/4)});

        end
        
    case 'Exponential'
        uopts = getUserInput(...
            {'Exponential decay target value: '},{0.2});
    case 'Sinebell'
        uopts = getUserInput({'Sinebell shift: (Samples): '},{0});
end

% Set options if applicable
if exist('uopts', 'var')
    if isempty(uopts), CSI_Log({'Skipped apodization.'},{''}); return; end
    opts = str2double(strsplit(uopts{1},' '));
else, uopts = {''}; opts(1) = 0; 
end

% Set FID or Echo options
if strcmp(uanstype{2},'FID'), opts(2) = 1; else, opts(2) = 2; end

% APPLY FILTER % ---------------------------------- %

% Check data domain type availability
domain = CSI_getDomain(gui);
if strcmp(domain,'freq'), csi.data.raw = csi_ifft(csi.data.raw); end

% Create window and apply filter
[csi.data.raw, win] = CSI_filterSpectra(csi.data.raw, uanstype{1}, opts);

% Revert to frequency domain if required.
if strcmp(domain,'freq'), csi.data.raw = csi_fft(csi.data.raw); end

% SAVE % ------------------------------------------ %

% Save data
setappdata(gui.CSIgui_main, 'csi', csi)

% Set opts to show in Hz if available
if exist('hz','var'), uopts{1} = hz; end 
% Update LOG
if  (size(win,1) > 1) % No apodization if window size equals 1.
    CSI_Log({['Filter ' uanstype{1} ' applied. Possible opts:']},...
               uopts);
end


% --- Normalize CSI data.
function CSI_Normalize(gui)
% Normalize data to maximum in a spectrum or the volume.
% The real-part of the spectrum is used.
%
% Additional: peak specific?

% Get userinput: normalize how?
uans = getUserInput_Popup({'Normalize by: '},...
    {{'Maximum per voxel','Maximum in volume','Specific peak per voxel'}});
if isempty(uans), CSI_Log({'Skipped normalize data.'},{''}); return; 
end  


% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get data-array
data = csi.data.raw;

% Split units: use real for normalization
dataR = CSI_getUnit(data,'Real'); dataI = CSI_getUnit(data,'Imaginary');


% Data as cell index layout
sz = size(dataR); cell_layout = ...
arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

% Create cell of data.
datac = mat2cell(dataR, sz(1), cell_layout{:});



switch uans{1}
    case 'Specific peak per voxel'
        % Get peak of interest
        doi = CSI_getDataAtPeak(data, csi.xaxis);
        if isnan(doi), return; end
        % Maximum to normalize to.
        max_val = max(real(doi),[],1); 
        datac = cellfun(@(x,y) x./y, ...
            datac,repmat({max_val},size(datac)), 'Uniform', 0);
        
    case 'Maximum per voxel'
        
        % Normalize to maximum of each voxel
        datac = cellfun(@(x) x./max(x(:)), datac, 'Uniform', 0);
        
    case 'Maximum in volume'
        
        % Normalize to maximum in volume
        datac = cellfun(@(x,y) x./max(y), ...
            datac, repmat({max(dataR(:))},size(datac)), 'Uniform', 0);
        
end

% Replace with original data
dataR = cell2mat(datac); data = complex(dataR, dataI); % Complex data set
csi.data.raw = data;

% Save appdata.
setappdata(gui.CSIgui_main, 'csi', csi);


% --- Executed by Apodization FID button
function [data, window] = CSI_filterSpectra(data, ftype, opts)
% Calculates the correct ftype window and applies it to each spectrum in 
% data.

% Get filter size e.g. nSamples.
fsz = size(data,1);

% Get filter
window = CSI_filter(fsz, ftype, opts);
if isnan(window), CSI_Log({'Cancelled apodization.'},{''}); return; end

% Rearrange data to cell 
sz = size(data); 
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
data = mat2cell(data, sz(1), cell_layout{:}); 


% Apply window 
data = cellfun(@times, data, repmat({window}, size(data)),...
            'UniformOutput', 0);
        
       
% Rearrange data to matrix. 
% Convert.
data = cell2mat(data);
% Undo if 1D data.
if numel(sz) <= 2, data = squeeze(data); end

% --- Executed by CSI_filterSpectra to get filter window
function filter_window = CSI_filter(fsz, ftype, opts)
% Get anapodization filter: hamming, hann, blackman, gaussian, flattop or
% exponential.

switch lower(ftype)
    case 'gaussian'
        % Size vector and options center and standard deviation
        filter_window = gaussmf(0:fsz-1, [opts(1) 0]);
    case 'hamming'
        % Size, type and orienation
        filter_window = window1D(fsz.*2, 'hamming', 'Symmetric');
        filter_window = filter_window(fsz+1:end,1);
    case 'hann'
        % Size, type and orienation
        filter_window = window1D(fsz.*2, 'hann', 'Symmetric');
        filter_window = filter_window(fsz+1:end,1);
    case 'blackman'
        % Size, type and orienation
        filter_window = window1D(fsz.*2, 'blackman', 'Symmetric');
        filter_window = filter_window(fsz+1:end,1);
    case 'flattop'
        % Size, type and orienation
        filter_window = window1D(fsz.*2, 'flattop', 'Symmetric');
        filter_window = filter_window(fsz+1:end,1);
    case 'exponential'
        % Size, exponent length and stretch
        filter_window = window1D_exp(fsz, opts(1))';
    case 'sinebell'
        % Size and shift
        filter_window = sinebell(fsz, opts(1));
    otherwise
        CSI_Log({'Filter type unknown.'},{'Returning'});
end

% Transpose.
if size(filter_window,2) > size(filter_window,1)
    filter_window = filter_window';
end

% FID or Echo: for echoes mirror and append to front
if size(opts,2) > 1 && opts(2) == 2
    sz = size(filter_window);
    % Flip and append
    filter_window = cat(1,flip(filter_window,1), filter_window);
    % Correct window length
    filter_window = filter_window(...
        round(linspace(1,size(filter_window,1),sz(1))));
end

% Display the filter window.
fh = figure(); plot(filter_window); title(['Applied ' ftype ' filter']);

uans = getUserInput_Buttons('Apply shown filter?',{'Yes', 'No'});
if strcmp(uans,'No'), filter_window = NaN; end

% Close the figure
close(fh);

% --- Executes on button press in button_CSI_Linewidth.
function button_CSI_Linewidth_Callback(~, ~, gui)
% Calculate linewidth at FWHM of a specific peak.
%
% Uses: csi_LineWidth();

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Userinput % --------------------------------- %
uans = getUserInput_Popup({'Display type: ',...
                           'Save line width data (.txt): '},...
                         {{'Table', 'Graph'},{'No','Yes'}});
if isempty(uans)
    CSI_Log({'Skipped line width calculations.'},{''}); return; 
end

% Display type
dataDisp = uans{1};
% Save data boolean
switch uans{2}, case 'Yes', dataSave = 1; case 'No', dataSave = 0; end

% Get data at peak of interest
[doi, ~, range] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);

% Prepare Data % --------------------------------- %
tmp = csi.data.raw;

% Maximum in this range
[~,mi] = max(real(doi),[],1);

% Set maximum value as peak centre.
peak_pos = mi+range(1);

% Set peak width 
peak_width = ceil(csi.xaxis.N/100);
if peak_width < 3, peak_width = 3; end
    
% Peak range used to find FWHM
peak_pLow = num2cell(peak_pos - peak_width);
peak_pHig = num2cell(peak_pos + peak_width);
poi_range_perVox = cellfun(@(x,y) x:y, peak_pLow, peak_pHig,'uniform',0); 

% Input cell data
sz = size(tmp); 
cell_layout = arrayfun(...
    @ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
tmp = mat2cell(tmp, sz(1), cell_layout{:}); 

% Input axis
ax_inp = repmat({csi.xaxis.ppm}, size(tmp));
% Set plotting off - Overload matlab otherwise!
plot_off = repmat({0},size(tmp)); % plot_off{1,2,4,1} = 1;


% Calculate line width % --------------------------------- %
linewidth = cellfun(@csi_LineWidth, tmp, ax_inp, poi_range_perVox,...
    plot_off, 'Uniform',0);


% Data display % --------------------------------- %
switch dataDisp
    case 'Table'        
        CSI_dataAsTable(cell2mat(linewidth), 'LineWidth')
    case 'Graph'
        CSI_dataAsGraph(cell2mat(linewidth), gui, 'LineWidth')
end

% Save data % --------------------------------- %
if dataSave
   % Default filepath to start UI in.
    if isfield(csi.data,'filepath'), fp_def = csi.data.filepath; 
    else, fp_def = []; 
    end

    % Get file path and extension from user
    [fn, fp, fi] = ...
        uiputfile({'*.txt', 'Text file'},'Save data...', fp_def);
    if fi == 0, return; end
    csi_writeText(cell2mat(linewidth), [fp fn]); 
end


% --- Executes on button press in button_CSI_setFrequency.
function button_CSI_setFrequency_Callback(hObject, eventdata, gui)
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
CSI_2D_Scaling_calc_xaxis(hObject, eventdata);


% --- Executes on button press in button_CSI_setLabels.
function button_CSI_setLabels_Callback(~, ~, gui)
% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get labels
old_labels = csi.data.labels; 

% Allow editing of labels by user
uans = getUserInput({'Change labels of each index:'},...
                    {strjoin(csi.data.labels,' ')});
if isempty(uans),CSI_Log({'Skipped label change.'},{''}); return; end

% Split the new labels.
csi.data.labels = strsplit(uans{1},' ');

% Save data
setappdata(gui.CSIgui_main, 'csi',csi);

% Update info
CSI_Log({'New labels:','Previous labels:',},...
         {strjoin(csi.data.labels,' | '),strjoin(old_labels,' | ')});

                              
% --- Executes on button press in button_CSI_ReorderDim.
function button_CSI_ReorderDim_Callback(~, ~, gui)
% Permute e.g. reorder csi.data.raw.

% BACKUP + APPDATA % ------------------- %
% Create backup
CSI_backupSet(gui, 'Before reordering.');

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% PERMUTE % ---------------------------- %
% Get order from user
uans = getUserInput(...
    {'New order of MRS data dimensions?:'},{1:size(csi.data.labels,2)});
% If user pressed skip.
if isempty(uans),CSI_Log({'Skipped data re-ordering.'},{''}); return; end

% Str-ans to dbl-ans: new data order.
new_order = str2double(strsplit(uans{1}));

% I. Reorder data 
csi.data.raw = permute(csi.data.raw ,new_order);

% II. Reorder labels
% To allow increase in dimensionality, as labels(new_order) will not work,
% iterate each label and assign to new_order-index accordingly.
% If empty labels due higher dimensionality are present, replace by empty
% string.
new_label = cell(1,max(new_order));
for li = 1:size(csi.data.labels,2)
    if sum(new_order == li)
        new_label{new_order == li} = csi.data.labels{li};
    else
        new_label{li} = []; % Label not found, empty labeling.
    end
end
empty_lab = cellfun(@isempty, new_label) == 1;
if sum(empty_lab) ~= 0, new_label(empty_lab == 1) = {'-'}; end
csi.data.labels = new_label;

% III. Reorder dimensions
csi.data.dim = size(csi.data.raw);
% Might be required to change this way and use new_label approach as
% dim-field of csi.data may contain data size equal to 1 for example after
% averaging or combining channels.

% CLEAN UP % ----------------------------- %

% Update appdata
setappdata(gui.CSIgui_main,'csi',csi);

% Update CSI info panel
if size(new_order) <= size(csi.data.labels,2)
    CSI_Log({'Reordered dimensions of MRS data to '},...
        {strjoin(csi.data.labels, ' | ')});
else
    CSI_Log({'Reordered dimensions of MRS data to '},{new_order});
end

% --- Executes on button press in button_CSI_ReorderDim_Auto.
function button_CSI_ReorderDim_Auto_Callback(hObject, ~, ~)
% Get gui data and run CSI_ReorderDim_Auto script
gui = guidata(hObject); CSI_ReorderDim_Auto(gui);

% --- Executes after loading CSI data and by several other scripts
function CSI_ReorderDim_Auto(gui)
% Automatically permute e.g. reorder CSI data for easy 2D display.
% The index will be ordered highest to lowest and looks for spatial labels
% to order these to the 2nd, 3rd and 4th index e.g. x, y and z in MRS data
% display.

% Return if no csi gui available.
if ~isappdata(gui.CSIgui_main,'csi'), return; end

% Get CSI data
csi = getappdata(gui.CSIgui_main, 'csi');

% 1. Find any spatial dimensions
kspc_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
kspc_dim(isnan(kspc_dim)) = [];

% 2.If kspace dims are not found...
if isempty(kspc_dim)
    % Create a new dim order from max to min dim-sizes
    sz = csi.data.dim;
    [~, sort_ind] = sort(sz(2:end), 'descend');
    new_order = [1 sort_ind+1];
else
   dim_order = 1:numel(csi.data.dim);
   % Set the new order with kspace dimensions first
   new_order = [1 kspc_dim];
   % Find excluded dimensions
   excl_dim_order = ismember(dim_order, new_order);
   % Append new_order.
   new_order = [new_order dim_order(excl_dim_order == 0)];
end

% Reorder using new_order.
csi.data.raw    = permute(csi.data.raw, new_order);
csi.data.labels = csi.data.labels(new_order);
csi.data.dim    = csi.data.dim(new_order);

% Update appdata; 
setappdata(gui.CSIgui_main,'csi',csi);

% After update info
if new_order <= size(csi.data.labels,2)
    CSI_Log({'Automatically reordered CSI data. New order:'},...
        {strjoin(csi.data.labels, ' | ')});
else
    CSI_Log({'Automatically reordered CSI data. New order:'},...
        {new_order});
end

% --- Executes on button press in button_CSI_Bin.
function button_CSI_Bin_Callback(~, ~, gui, backup)
% Bin a specific dimension, splitting up the csi data.

if nargin < 4, backup = 1; end


% BACKUP + APPDATA % ------------------------------- %
% Create backup
if backup, CSI_backupSet(gui, 'Before applying bins.'); end

% Check if csi app data exists
if ~isappdata(gui.CSIgui_main, 'csi'),return; end


% User Input % -------- %
% Dimension to bin over and number of bins?
uans = getUserInput({'Dimension to bin: ','Number of bins: '},{1,2});
if isempty(uans),CSI_Log({'Skipped binning.'},{''}); return; end
dim = str2double(uans{1}); nbin = str2double(uans{2});


% Create backup
CSI_backupSet(gui, 'Before binning.');
% Get CSI data
csi = getappdata(gui.CSIgui_main, 'csi');


% Process data
tmp = csi.data.raw;  sz = size(tmp);

stepsz = sz(dim)/ nbin;
if ~(mod(sz(dim),nbin)==0) % BIN size does not fit the dimension
    CSI_Log({'Number of bins does not fit the data. Bin size:'},{stepsz});
    return;
end
% Calculate vector with index in tmp for each bin.
bin_vec = 1:stepsz:sz(dim); bin_vec(1,end+1) = sz(dim)+1;
% Index cell for all other dimensions
dim_vec = cellfun(@(x) 1:x, num2cell(sz),'uniform',0);
dim_vec{dim} = NaN;

sz_temp = sz; sz_temp(dim) = stepsz;
bin = NaN([sz_temp(:)',nbin]);
for bi = 1:size(bin_vec,2)-1
    
     % Index to correct arrays
     dim_vec_temp = dim_vec; 
     dim_vec_temp{dim} = bin_vec(bi):bin_vec(bi+1)-1;
     % Data
     tmp_bin = tmp(dim_vec_temp{:});
     % Store
     dim_vec_temp = dim_vec; dim_vec_temp{dim} = 1:stepsz;
     dim_vec_temp{end+1} = bi;
     bin( dim_vec_temp{:} ) = tmp_bin; 
    
end

% Store data
csi.data.raw = bin; csi.data.labels = cat(2,csi.data.labels,'bin');
csi.data.dim = size(csi.data.raw);

% Save appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% Log
msg = ...
sprintf('Binned dimension "%s" over %i bins', csi.data.labels{dim}, nbin);
CSI_Log({msg},{''});


% --- Executes on button press in button_CSI_Sum.
function button_CSI_Sum_Callback(~, ~, gui, backup)
% Summate MRSI data over a specific dimensions

% BACKUP + APPDATA % ------------------- %

if nargin < 4 , backup= 1; end
if backup, CSI_backupSet(gui, 'Before summation.');  end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');




% FUNCTION % --------------------------- %

% Apply summation
[data_summed, index] = CSI_summate(csi.data.raw);

% Replace with summed data
csi.data.raw = data_summed;
% Change data dimension size
csi.data.dim(index) = 1;

% CLEAN UP % --------------------------- %

% Store appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% After update info
if index <= size(csi.data.labels,2)
    CSI_Log({'Summated MRS data over dimension'},...
        {strjoin(csi.data.labels(index), ' ')});
else
    CSI_Log({'Summated MRS data over dimension'},{index});
end

% --- Executed by button CSI_Sum_Callback
function [data_summed, index] = CSI_summate(data, index)
% Sum over a specific index in data; if no index to summate over is given,
% the user is asked to give one.

if nargin == 1 % Ask user which dimension if no index given.
    index = getUserInput({'Enter data index to summate over:'},{'4'});
    if isempty(index),CSI_Log({'Skipped summation.'},{''}); return; end
    index = str2double(index{1});
end

% Summate over index
data_summed = sum(data,index);


% --- Executes on button press in button_T1_MRS.
function button_T1_MRS_Callback(hObject, eventdata, gui)
% Calculate T1 for each voxel over a specific dimension in the data.
% Requires TR input from user.
%
% 1. Get TR from user
% 2. Prepare data
% 3. Prepare plot-figure for each slice: axis/voxel.
% 4. Calculate T1 for each voxel plot.
% 5. Display data, fitted-data, T1, R2 and CI for each voxel.
% 6. Save T1 data to mat-file.
%
% Fit using S = S0 (1-exp (-TR/T1)

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get TR/TI of the data. % --------------------------------------------- %
uans = getUserInput(...
    {'Repetition time (TR):','Inversion time start (TI):',...
     'Inversion time step (dTI):','Fitted over the fifth index.'},...     
    {2000,5,10,'Skip to manually reorder the data first if necessary.'}); 
% User pressed skip or something else.                
if isempty(uans)
    CSI_Log({'Skipped T1 calculations.'}, {''}); return; 
end             
% Define numerical TE!
TRinit = str2double(uans{1}); 
TIstep = str2double(uans{3}); TIstart = str2double(uans{2});

% Prepare the data % --------------------------------------------------- %

% Switch to correct data unit.
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};
switch data_unit
    case 'Real',      data = real(csi.data.raw);
    case 'Imaginary', data = imag(csi.data.raw);
    case 'Magnitude', data = abs(csi.data.raw);
    case 'Phase',     data = angle(csi.data.raw);
    otherwise,        data = abs(csi.data.raw);
end

% Calculate the maximum values.
data_max = max(data,[],1); 

% CALCULATE T1/VOXEL/SLICE % ------------------------------------------- %
% NB: other dimensions besides the TR-index (5) are ignored.
T1 = struct;
CSI_Log({'Calculating T1 for each voxel and slice.'},...
               {'Please be patient.'});
for sli = 1:size(data_max,3)
    
% FIGURE PER SLICE % --------------------------------------------------- %
% Figure for all axis objects in this slice.
fh = figure('Tag', 'CSIgui_plotT1','Name',...
            ['CSIgui - T1 plot for slice: ' num2str(sli)],...
            'Color', 'Black','Toolbar', 'figure', 'MenuBar', 'None',...
            'NumberTitle', 'Off');  
        
% 1. Default figure size and screen pixel size
def_sz = 720; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
% 2. Ratio to define figure height to def_size
fig_sz = [def_sz def_sz.*(scr_sz(2)/scr_sz(1))];
% 3. Position of figure.
fig_ps = [40 scr_sz(2)-(1.15*fig_sz(2))];
% 4. Apply
set(fh, 'Position', [fig_ps fig_sz]);

% PLOT DIMENSIONS % ---------------------------------------------------- %
% Explanation:
% Get dimension index sizes and process
% Use dim(2) and dim(3) as XY, Use dim(3) as Z/SLICE, Use dim(...) 
% as dimensionality (Line). Plot_para.dim = [x y slice others...];
% 1. Create plot grid for each voxel e.g. spectrum
plot_par.dim      = size(data_max);     % Data dimensions
plot_par.dim(1)   = [];                 % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);% 3D/2D/1D volume. 

% Set in plot_par the 1D dimension to the second index e.g. column.
if plot_par.data_dim == 1
    % Add a Y dimension.
    plot_par.dim(2) = 1; plot_par.data_dim = 2;  
end

% AXIS LAYOUT % -------------------------------------------------------- %
% Layout of the axis in the figure for each voxel.

% Axis Resolution: Relative to figure (Normalized) e.g. [1 1];
plot_par.res = 1./plot_par.dim; 
% Grid with locations for each axis to plot to in the figure.
for axi = 1:numel(plot_par.dim)
    plot_par.range{axi} =   ...
    0 : plot_par.res(axi) : ...
       (plot_par.res(axi) * plot_par.dim(axi) - plot_par.res(axi));
end
[x,y] = meshgrid(plot_par.range{1},plot_par.range{2});
plot_par.grid.x = (x); plot_par.grid.y = flipud(y);

% 1D CORRECTION -------------------------------------------------------- %
% Flip UD of the y/row-grid values and transpose x/col 
if plot_par.data_dim == 1
    plot_par.grid.y = flipud(y)'; plot_par.grid.x = x'; 
end


% Create a cell array for other higher dimension index. 
% ---> Line/Echo dimension, after slice. All taken accumulated!
if numel(plot_par.dim) >3, nDimC = num2cell(plot_par.dim(4:end));
else, nDimC = {1};
end
% To linear vector per cell.
nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

% PLOT LIMIT COLORING % ------------------------------------------------ %
vol_scaling = get(gui.menubar.Opts_ColorScale_Volume, 'Checked');
switch vol_scaling
    case 'on'
        % Take full volume of the data and calc min and max
        data_ylimits = [min(data_max(:)), max(data_max(:))];  
    case 'off'
        % Plot data of the whole slice for ylimit coloring.
        plot_data_slice = squeeze(data_max(:,:,:,sli,nDimC{:}));
        % Calc limit
        data_ylimits = [min(plot_data_slice(:)), max(plot_data_slice(:))];           
end
% Check limits agrees with rules: lim(1) < lim(2)
if data_ylimits(2) <= data_ylimits(1),data_ylimits(2) =  data_ylimits(1)+1; 
end
% Get plot color
[clrs, clrs_data_range] = CSI_2D_Scaling_calc_ColorOfPlots(data_ylimits);  

% CALC & PLOT VOXEL % --------------------------------------------------- %
% Be aware of the row/column switch to get proper x/y axis display in the
% figure. E.g. dim(1) = column = x, dim(2) = row = y; Otherwise this would
% be Matlab style reversed.
for ci = 1:plot_par.dim(1) % Column loop.
    for ri = 1:plot_par.dim(2) % Row loop.
        
        % X and Y position in the figure of the axis object for this voxel
        x = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
        % Position of axis to plot
        pos = [x y plot_par.res(1) plot_par.res(2)];
        % Create axis with pos(3,4) size at pos(1,2) position
        plot_par.ax{ri,ci} = axes('parent',fh,'position',pos);
        
        % Get plot data.
        plot_data = squeeze(data_max(:,ci,ri,sli,nDimC{:}));
       
        % CALCULATE T1 %%%%%%%%%%%%%%%
        % T1 returns: T1, R2, Confidence interval and amplitude.
        TI_vec = NaN(1,size(plot_data',2)); TI_vec(1) = TIstart+TRinit;
        for kk = 2:size(plot_data',2)
            TI_vec(kk) = TI_vec(kk-1)+TIstep;
        end
        
%         [T1.T2(ri,ci,sli,1), T1.R2(ri,ci,sli,1:2),T1.CI(ri,ci,sli,1:2),...
%             T1.A(ri,ci,sli,1)] = T1_expIR(plot_data', TR_vec, []);
%         
%         % Calculate y-values from fit parameters.
%         % Exponential function of S = A*exp(-TE/T2);
%         fitfunc = @(param,x) param(1).*(1-exp(-param(2)./param(3))); 
%         yfitval = fitfunc(...
%             [T1.T2(ri,ci,sli,1) T1.A(ri,ci,sli,1)], TR_vec);

        % Y LIMIT calculated coloring
        % Actual display ylimit changed below at Axis Cosmetics.
        ylimit = [min(plot_data(:)) max(plot_data(:))];
        if ylimit(2) == ylimit(1)
            ylimit = ...
            [-0.98*abs(max(plot_data(:))) 1.02*abs(max(plot_data(:)))];
        end
        
        % PLOT COLOR 
        % Set relative to maximum Y-data. See PLOT LIMIT COLORING 
        [~, clr_ind] = min(abs(clrs_data_range-ylimit(2)));
        plot_color = clrs(clr_ind,:); % See before ri/ci for-loops.
  
        % PLOT DATA MAX
        plot(plot_par.ax{ri,ci}, TI_vec, plot_data, 'color', plot_color,...
            'Linewidth', 1.5); hold on;
%         % PLOT DATA FIT
%         plot(plot_par.ax{ri,ci}, TR_vec, yfitval, 'Color', 'r',...
%             'LineWidth', 1.5,'LineStyle',':');
%         % Plot T2 and CI and R2 as text.
%         text(TR_vec(2), plot_data(1).*0.95, ...
%             sprintf('T2: %3.1f - R^2: %1.2f\n95%%CI: %3.2f | %3.2f',...
%              T1.T2(ri,ci,sli,1),T1.R2(ri,ci,sli,1),... 
%              T1.CI(ri,ci,sli,1),T1.CI(ri,ci,sli,2)),...
%             'Color', [0.6 0.6 0.6],'FontSize',8,'FontWeight','Bold');            
        
        % Axis Cosmetics
        ylimit = ylimit * [0.975 0; 0 1.025]; ylim(ylimit);
        set(plot_par.ax{ri,ci},'Color', 'None','YLim', ylimit,...
                               'XColor', [0.4 0 0],'YColor', [0.4 0 0],...
                               'LineWidth', 1.75, 'Xtick',[], 'Ytick',[]);
        
    end % for row loop
end % for column loop
end % for each slice loop

% Display info to user.
CSI_Log({'T1 calculated. Saving T1 results to mat-file.'}, {''});
% Save the T2 data struct as a MAT-file.
% save([datestr(now,'yyyymmdd_HHMMSS') '_T2Calc_' data_unit '.mat'], 'T2');
disp('')

% --- Executes on button press in button_CSI_T2.
function button_CSI_T2_Callback(hObject, eventdata, gui)
% Calculate T2 for each voxel over a specific dimension in the data.
% Requires TE/echotime input from user.
%
% 1. Get TE from user
% 2. Prepare data
% 3. Prepare plot-figure for each slice: axis/voxel.
% 4. Calculate T2 for each voxel plot.
% 5. Display data, fitted-data, T2, R2 and CI for each voxel.
% 6. Save T2 data to mat-file.
%
% Uses T2_exp();

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');


                % ----------- % Userinput % ----------- %

% ECHO TIME % --- %
uans = getUserInput({'Echotime spacing:', 'Echotime of the first echo:'},...
                    {50,0});
% User pressed skip or something else.                
if isempty(uans)
    CSI_Log({'Skipped T2 calculations,'},...
            {'Echotime information required.'}); 
    return; 
end 
% Define numerical TE
TEinit = str2double(uans{2}); TEstep = str2double(uans{1});

% IMAGE PLOT % --- %
uans = getUserInput_Popup({'Plot Images: '},{{'Yes','No'}});
if isempty(uans), CSI_Log({'Skipped T2 calculations.'},{''}); return; end
switch uans{1}, case 'Yes', imgPlot = 1; case 'No', imgPlot = 0; end

% Peak of interest: active only
% The doi_axis is used to get peak of interest in echo data if available 
% and requested by the user: active type in csi.data.split.
[doi, doi_xaxis] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
if isnan(doi), return; end % Skipped

% FID & Echo detection
combine_fe = 0;
if isfield(csi.data,'split')
    % From user: Calculate using FID and echoes?
    uans = getUserInput_Popup({'Use both FID and Echo data?'},...
        {{'Yes', 'No'}});
    if isempty(uans)
        CSI_Log({'Skipped splitting FID & Echo.'},{''}) ; return; 
    end
    
    % Combine fid and echo yes or no?
    switch uans{1}, case 'Yes', combine_fe = 1; end
end

               % ----------- % Prepare Data % ----------- %

% Get data unit setting from gui
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};

if combine_fe
    % Get calculated maxima from FID and Echo
    switch csi.data.split.type
        case 'FID'

            % FID: Convert doi e.g. data at peak of interest
            data_fid = CSI_getUnit(doi, data_unit);
            % FID: Maximum of FIDs
            max_fid = max(data_fid,[],1); 
            
            % ECHO: get peak in echo data
            doi_echo = CSI_getDataAtPeak_Stored(doi_xaxis, gui);
            
            % ECHO: Convert echo data and calculate maximum
            data_echo = CSI_getUnit(doi_echo, data_unit);
            max_echo = max(data_echo,[],1);
            
        case 'Echo'
            
            % ECHO: Convert doi e.g. data at peak of interest
            data_echo = CSI_getUnit(doi, data_unit);
            max_echo = max(data_echo,[],1); 
            
            % FID: get peak in FID data
            doi_fid = CSI_getDataAtPeak_Stored(doi_xaxis, gui);
            
            % FID: Convert FID data and calculate maximum
            data_fid = CSI_getUnit(doi_fid, data_unit);
            max_fid = max(data_fid,[],1);
            
    end
    
    new_data = squeeze(data_fid(:,4,4,4,1));
    for kk = 1:size(data_echo,5)
        new_data = cat(1,new_data, data_echo(:,4,4,4,kk));
    end
    figure(); plot(new_data); title('FID + All echoes');
    
    % As T2 fit echoes are expected the 5th dimension...
    % Concatenate at index 5.
    data_max = cat(5, max_fid, max_echo);
    
else
    
    % Convert data to unit of interest
    data = CSI_getUnit(doi, data_unit);

    % Calculate the maximum values.
    data_max = max(data,[],1); 
end


                % --------- % CALCULATE T2 % --------- %

% Set echoes (index 5) to first dimension
sz = size(data_max); szd = numel(sz); perm_vec = 1:szd;
perm_vec(1) = 5; perm_vec(5) = 1; data_main = permute(data_max,perm_vec);

% Cell layout
sz = size(data_main); 
cell_layout = ...
    arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

% Creat cell of data.
data_cell = mat2cell(data_main, sz(1), cell_layout{:}); 
data_cell = cellfun(@transpose, data_cell,'UniformOutput',0);

% Echotime vector
TE_vec = NaN(1,size(data_main,1)); TE_vec(1) = TEinit;
for kk = 2:size(data_main,1),TE_vec(kk) = TE_vec(kk-1)+TEstep; end

% Calc T2 for each cell
[data_T2, data_R2, data_CI, data_A] = ...
    cellfun(@T2_exp, data_cell, ...
                     repmat({TE_vec},size(data_cell)), ...
                     'UniformOutput',0);

% Y-Fit values from fit-parameters
% Exponential function: S = A*exp(-TE/T2);
data_yfit = cellfun(@T2_formula, ...
    data_T2, data_A, repmat({TE_vec},size(data_cell)),'UniFormOutput',0);
data_yfit = cellfun(@transpose, data_yfit,'UniformOutput',0);    


            % --------- % CREATE DISPLAY DATA % ---------- %
% For correct sprintf split the CI and RI in two seperate cell arrays.

% SPLIT R2
sz = size(data_R2);
perm_vec = 1:numel(size(data_R2)); perm_vec(1) = 2; perm_vec(2) = 1;
a = permute(data_R2,perm_vec); b = cell2mat(a); c = permute(b,perm_vec);
extr_vec = arrayfun(@(x) 1:x, sz(2:end),'UniformOutput',0);
data_R2a = c(1,extr_vec{:});
data_R2a = mat2cell(data_R2a,1,cell_layout{:}); 

% SPLIT CI
a = permute(data_CI,perm_vec); b = cell2mat(a); c = permute(b,perm_vec);
data_CIa = c(1,extr_vec{:}); data_CIb = c(2,extr_vec{:});
data_CIa = mat2cell(data_CIa,1,cell_layout{:});
data_CIb = mat2cell(data_CIb,1,cell_layout{:});

% Text to display.
txt = repmat({'T_2: %3.1f | R^2: %1.2f\n95%%CI: %3.0f | %3.0f'},...
              size(data_cell));         
data_text = cellfun(@sprintf, txt, data_T2,...
                data_R2a, data_CIa, data_CIb,'UniFormOutput',0);

% Create output data for displaying graphs

% NORMALIZE
sz = size(data_main); 
cell_layout = ...
arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
data_main = mat2cell(data_main, sz(1), cell_layout{:});

data_main = cellfun(@(x) x./max(x), data_main,'Uniform',0);


data_yfit = cellfun(@(x) x./max(x), data_yfit,'Uniform',0);


data2plot{1} = cell2mat(data_main); data2plot{2} = cell2mat(data_yfit);

% Send to special dataAsGraph function
CSI_dataAsGraph2(data2plot,data_text,'T2',gui, imgPlot);


                % -------- % SAVE T2 DATA % ---------- %

% Display info to user.
str_time = datestr(now,'yyyymmdd_HHMMSS');
CSI_Log({'T2 calculated. Saving T2 results to mat-file.'}, ...
               {['Filename prefix: ' str_time]});

% Create arrays with equal "spatial" relation/indexing to T2 plot.
% e.g. permute the row and column index, convert cell to array.
T2 = permute(cell2mat(squeeze((data_T2))),[2 1 3:numel(size(data_T2))]);
R2 = permute(cell2mat(squeeze((data_R2a))),[2 1 3:numel(size(data_R2))]);
CIl = permute(cell2mat(squeeze((data_CIa))),[2 1 3:numel(size(data_CIa))]); 
CIh = permute(cell2mat(squeeze((data_CIb))),[2 1 3:numel(size(data_CIb))]);  

% Create output struct
T2data.T2 = T2; T2data.R2 = R2; T2data.CI_low = CIl; T2data.CI_high = CIh;

% Save the T2 data struct as a MAT-file.
save([ str_time '_T2Calc_' data_unit '.mat'], 'T2data');


% --- Executes on button press in button_CSI_ZeroFill.
function button_CSI_ZeroFill_Callback(~, ~, gui, backup)
% Zero fill the data in the time domain. Requires userinput to get lentgh
% of the FID after zero filling.
if nargin < 4, backup = 1; end;
% BACKUP + APPDATA % ------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before zero filling.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USER INPUT % ------------------------- %

uans = getUserInput({'Requested #samples after zero filling: '},...
                    {size(csi.data.raw,1).*2});
if isempty(uans), CSI_Log({'Skipped zero filling.'},{''}) ; return; end

% Convert to length after zero filling.
N = str2double(uans{1});

% Get zerofill direction from FID/Echo e.g post/both
uans = getUserInput_Popup({'FID or Echo: (Post/Both)'},{{'FID','Echo'}});
if isempty(uans), CSI_Log({'Skipped zero filling.'},{''}) ; return; end

switch uans{1} % Direction for zerofilling of fids or echoes
    case 'FID', dir = 'Post'; 
    case 'Echo',dir = 'Both';
end

% ZERO FILL % --------------------------- %

% Correct data domain (I/II)
domain = CSI_getDomain(gui);
if strcmp(domain,'freq'), csi.data.raw = csi_ifft(csi.data.raw); end

% Use the csi_zeroFill function.
% First dimension represents the data to zerofill. The integer N
% will be the resulting sample size.
csi.data.raw = csi_zeroFill(csi.data.raw, N, dir);

% Change other parameters
csi.data.dim(1) = N;

% Correct data domain (I/II)
if strcmp(domain,'freq'), csi.data.raw = csi_fft(csi.data.raw); end


% CLEAN UP % ---------------------------- %

% Save appdata
setappdata(gui.CSIgui_main,'csi',csi);

% Recalculate xaxis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% adjwandkawdnawnd
% if function above is set to auto, it sets the wrong limits 
% NEEDS UPDATING.

% Update info
CSI_Log({'Applied zero filling. Sample size increased to'},{N});

% --- Executes on button press in button_CSI_AutoPhase.
function button_CSI_AutoPhase_Callback(hObject, eventdata, gui,backup)
% Apply zero order phase correction to all voxels in the data set.
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % ---------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before auto phase correction.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');


% USER INPUT % ---------------------------------- %

% Get peak of interest.
[poi] = CSI_getPeakOfInterest(csi.xaxis);
if isempty(poi), return; end

% Get method of auto-phasing
uans = getUserInput_Popup({'Auto phasing method:'},...
    {{'Match real part to maximum absolute signal.',...
      'Maximize real part of signal.'}});
if isempty(uans), CSI_Log({'Skipped zero-phasing.'},{''}) ; return; end

switch uans{1}
    case 'Match real part to maximum absolute signal.', phase_method = 2;
    case 'Maximize real part of signal.',phase_method = 1;
end


% APPLY CORRECTION % ----------------------------- %

% POI from user.
if length(poi) > 1, poi = poi(1):poi(2); end

% Correct data domain (I/II)
domain = CSI_getDomain(gui);
if strcmp(domain,'time'), csi.data.raw = csi_fft(csi.data.raw); end

% MRSI data to cell format
sz = size(csi.data.raw); 
cell_layout = arrayfun(@ones,...
    ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
cell_mrsi = mat2cell(csi.data.raw, sz(1), cell_layout{:});

% Apply auto zerophase to each cell
cell_mrsi_phased = ...
cellfun(@csi_autoZeroPhase, ...
            cell_mrsi, ...                              % data
            repmat({poi},    size(cell_mrsi)),...       % range
            repmat({phase_method},size(cell_mrsi)),...  % method
            repmat({0},      size(cell_mrsi)),...       % plot
            'UniformOutput', 0);            

% MRSI data to array        
array_mrsi = cell2mat(cell_mrsi_phased);

% Write to csi-struct
csi.data.raw = array_mrsi;

% Correct data domain (II/II)
if strcmp(domain,'time'), csi.data.raw = csi_ifft(csi.data.raw); end

% CLEAN UP % ------------------------------------ %

% Save
setappdata(gui.CSIgui_main,'csi',csi);

% Update info to user.
CSI_Log({'Applied automatic zero order phase correction:'},uans);

% --- Executes on button press in button_CSI_Flip.
function button_CSI_Flip_Callback(~, ~, gui)
% Flip order of a specific dimension in the spectral data set. 
if nargin < 4, backup = 1; end

% Check for CSI data
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

if backup, CSI_backupSet(gui, 'Before flipping.'); end

% Get dimension to flip over
useans = getUserInput({'Dimension:'}, {2});
if isempty(useans), CSI_Log({'Skipped flipping data'},{''}) ; return; end
dim = str2double(useans{1});

% Flip over this dimension
csi.data.raw = flip(csi.data.raw,dim);

% Update appdata
setappdata(gui.CSIgui_main,'csi',csi); 

% Display information to user.
CSI_Log({'CSI data flipped. Dimension:'},{dim});

% --- Executes on button press in button_CSI_Rotate.
function button_CSI_Rotate_Callback(~, ~ , gui)
% Rotate over a specific plane.
% Also required rotating labels and dimensions of this plane!

% Get csi data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USERINPUT % ------------------ %

% Ask user which axis to rotate over
uans = getUserInput({'Dimensional plane to rotate?'},{'2 3'});
if isempty(uans), CSI_Log({'Skipped rotating data.'},{''}) ; return; end
% Plane to rotate.
rotpl = str2double(strsplit(uans{1},' ')); 

% ROTATE % --------------------- %
% Rot90 works only if rotating plane is on index 1 and 2.

% Rotate
csi.data.raw = CSI_rotate(csi.data.raw, rotpl, 1);

% Apply rotation to labels
csi.data.labels(rotpl) =  csi.data.labels(fliplr(rotpl));
% And dimensions
csi.data.dim(rotpl) = csi.data.dim(fliplr(rotpl));

% 4. Save appdata
setappdata(gui.CSIgui_main, 'csi', csi);

% ------------------------------------------------- %

% Display information to user.
CSI_Log({'Rotated plane'},...
    {[csi.data.labels{rotpl(2)} ' | ' csi.data.labels{rotpl(1)}]});


function [spec, permv] = CSI_rotate(spec, plane, nturns)
% Rotates the MRS plane by 90 degrees counter clockwise for every turn.
% permv = the permute vector to make rot90 to work.

if nargin < 3, nturns = 1; end


% Plane dimensions which rotate
rot_plane = plane;

% Create permute vector: plane to index 1 and2.
tot_vec = 1:numel(size(spec)); lia = ismember(tot_vec, rot_plane);
permv = [rot_plane tot_vec(lia == 0)];

% Permute
spec = permute(spec, permv);

% Rotate 
spec = rot90(spec, nturns);

% Permute back
spec = ipermute(spec, permv);




% --- Executes on button press in button_CSI_Multiply.
function button_CSI_Multiply_Callback(~, ~, gui, backup)
% Multiply the MRS data

if nargin < 4 , backup= 1; end

% Create backup
if backup, CSI_backupSet(gui, 'Before multiplication.'); end

% Get CSI data
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');
    

% Userinput % --------------------- %

% Get factor
uans = getUserInput({'Multiply by '},{2});
if isempty(uans), CSI_Log({'Skipped multiplication.'},{''}) ; return; end
fac = str2double(uans{1});

% Multiply
csi.data.raw = csi.data.raw.*fac;

% Clean Up % ---------------------- %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% Update LOG.
CSI_Log({'MRS data multiplied by '},{fac});

% --- Executes on button press in button_CSI_Divide.
function button_CSI_Divide_Callback(hobj, ~, gui)
% Divide the MRS data


% If no GUI as input
if nargin ~= 3, gui = guidata(hobj); end

% Create backup
CSI_backupSet(gui, 'Before dividing');   


% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');
  

% Userinput % --------------------- %

% Get factor
uans = getUserInput({'Divide by '},{2});
if isempty(uans), CSI_Log({'Skipped division.'},{''}) ; return; end
fac = str2double(uans{1});

% Multiply
csi.data.raw = csi.data.raw./fac;

% Clean Up % ---------------------- %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% Update LOG.
CSI_Log({'MRS data divided by '},{fac});


% --- Executes on button press in button_Normalize.
function button_Normalize_Callback(~ , ~, gui, backup)
% Normalize data to specific peak maximum: multiple methods available. See
% CSI_Normalize();
if nargin < 4 , backup = 1; end

% Create backup
if backup, CSI_backupSet(gui, 'Before normalization.'); end

% Normalize
CSI_Normalize(gui);




% --- Executes on button press in button_CSI_Combine.
function button_CSI_Combine_Callback(~, ~, gui, backup)
% Combine the channels of all coils.
% 1. Get CSI data
% 2. Get requirements for channel combination from user.
% 3. Process data and apply chosen algorithm or settings.
if nargin < 4, backup = 1; end

% BACKUP + APPDATA % --------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before combining channels.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% WSVD APPLIED % ----------------------- %
% ----- Get user input and channel index.
% Get options from user.
% qst    = {'WSVD Combination? (y/n): (Other options are ignored) ',...
%           'Exclude channels: (Leave empty if none)',...
%           'Summate only: (y/n)'}; defans = {'y','', 'n'};
% uans = getUserInput(qst, defans);
% if isempty(uans), return; end % User pressed skip!


% This should be a gateway only
% 1. WSVD:   noise scans or noise mask, exclude channels.
% 2. Manual: exclude channels, summate only.

% Figure: CSI Combine channels options
fig_cO = figure('NumberTitle', 'Off', 'resize', 'off',...
                'Color', 'Black', 'MenuBar','none', 'Toolbar', 'none', ...
                'Tag', 'CSI_CombOpt', 'Name', 'Combine');
gui_combCoil = guidata(fig_cO); gui_combCoil.fig = fig_cO;

% Figure size and position
CSImain_pos = get(gui.CSIgui_main,'Position'); w = 240; h = 60; 
szdiff  = (CSImain_pos(3:4)-[w h])./2;      % Get figure in middle of main
fig_pos = [CSImain_pos(1:2)+szdiff w h];    % Use dSize of w/h to set.
set(fig_cO,'Position',fig_pos);


% Add buttons
% Buttons, handles and their info description.
bName = {'Manual','WSVD'};
bCall = {@CSI_Combine_Manual, @CSI_Combine_WSVD};
bInfo = {'Manual combination of channels.',...
         'Whitened singular voxel decomposition for channel combinations.'};
         
% Add buttons     
nButtons = size(bName,2); bw = 60; bh = 20;
bpos = [(w-(bw*nButtons))/3, (h-bh)/2]; % Buttons have gap;
for bi = 1:nButtons
   gui_combCoil.button{bi} = ...
   uicontrol(fig_cO, 'Style','pushbutton','String', bName{1,bi},...
         'Tag', ['button_combOpts_' bName{1,bi}], 'Callback',bCall{bi},...
         'Position',[(bpos(1,1)*bi)+((bi-1)*bw) bpos(1,2) bw bh],...
         'ForegroundColor', [0.94 0.94 0.94],'BackgroundColor', 'Black',...
         'ToolTipString',bInfo{bi});
end
% Update gui information
guidata(gui_combCoil.fig , gui_combCoil);

% --- Executes if No WSVD for coil combination is selected.
function CSI_Combine_Manual(hobj, ~)
% Comine channels using a manual method. Include specific channels and
% summate the combined channels or calculate the mean.
%
% Requires prior user input. (String).
% userInp{1} = 'Channels to exclude'
% userInp{2} = 'Summation only'
%
% The above will be in this function --> See to do list Combine update.


% GUI prepwork %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% If exists, close the combine options menu
if exist('hobj', 'var')
    gui_combCoil = guidata(hobj);close(gui_combCoil.fig);
end

% Get GUI and object: CSIgui_main
CSIgui_obj = findobj('Tag', 'CSIgui_main');gui = guidata(CSIgui_obj);

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
% Get CSI data
csi = getappdata(gui.CSIgui_main, 'csi');



% User input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get UserInput
qst = {'Exclude channels:','Summate only: (y/n)'}; defans = {5,'n'};
uans = getUserInput(qst, defans);
if isempty(uans), CSI_Log({'Skipped combing channels.'},{''}) ; return; end


% Set UserInput
% 1. Exclude channels
ch_excl = str2double(strsplit(uans{1},' '));
if isempty(ch_excl), ch_excl = NaN; end
% 2. Summate only
if strcmp(uans{2}, 'y'), sum_only = 1; else, sum_only = 0; end

% Data prepwork %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find channel dimension index
if isfield(csi.data, 'labels')
    chind = csi_findDimLabel(csi.data.labels,{'chan'});
end

% If no channel dimension/index found - ask user.
if ~exist('chind', 'var') || isempty(chind)
    chind = getUserInput({'Channel index: '},{'2'});
    if isempty(chind), CSI_Log({'Skipped combining channels.'},{''}) ; return; end
    chind = str2double(chind{1});
end

% ------- Prepare the exclusion of channels

% Create included channels array
ch_incl = 1:size(csi.data.raw, chind); 
if ~isnan(ch_excl), ch_incl(ch_excl) = []; end
% Create cell with all data per dimension of this spectral data
dims = csi.data.dim; dim_cell = cell(1,size(dims,2)); 
for kk = 1:size(dims,2), dim_cell{kk} = 1:dims(kk); end
% Replace the channel dimension with the included channels.
dim_cell{chind} = ch_incl;

% ------- Summate the channel dimensions over all included channels.
csi.data.raw = sum(csi.data.raw(dim_cell{:}),chind);

% ------- Create mean from summated if requested
if ~sum_only , csi.data.raw = csi.data.raw./size(ch_incl,2);
end

% ------- Update APP and APP data
setappdata(gui.CSIgui_main, 'csi', csi);

% After update, show info to user.
CSI_Log({'Channels combined over dimension',...
                'Excluded channels', 'Included channels'},...
               {chind, ch_excl, ch_incl});
if ~sum_only 
    CSI_Log({'Divided the summed data by'},{size(ch_incl,2)});
else
    CSI_Log({'No division applied to data.'},{'Summated only.'});
end

% --- Executes if WSVD for coil combination is selected.
function CSI_Combine_WSVD(hobj,~)
% Combine channels using WSVD algorithm:
%            Whitened singular voxel decomposition steps.
% 1. Prepare noise component: use mask or noise pre-scans.  See ref.
%           - User input required.
% 2. Average data: best results in WSVD after averaging.    See ref.
%           - If average dimension not
% 3. Reshape data: create {dt x #channels} x #voxels.
% 4. WSVD: see ref or help of wsvd script.
% 5. Reshape data: revert back to original array dimensions.
%
% REFERENCE + Script creator:
%   1. Rodgers and Robson, Magnetic Resonance in Medicine, 2010.
%      (http://dx.doi.org/10.1002/mrm.22230).
%   2. Rodgers and Robson, Magnetic Resonance in Medicine, YYYY.
%      (http://dx.doi.org/INSERT FINAL DOI HERE).
%   
%   Copyright Chris Rodgers, University of Oxford, 2014.


% GUI prepwork % ------------------------------------------- %
 
% If exists, close the combine options menu
if exist('hobj', 'var')
    gui_combCoil = guidata(hobj);close(gui_combCoil.fig);
end

% Get GUI and object: CSIgui_main
CSIgui_obj = findobj('Tag', 'CSIgui_main');gui = guidata(CSIgui_obj);

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
% Get CSI data
csi = getappdata(gui.CSIgui_main, 'csi');


% Average and channel index
ind_avg = csi_findDimLabel(csi.data.labels, {'aver'});
ind_cha = csi_findDimLabel(csi.data.labels, {'chan'});
if isnan(ind_cha)
    CSI_Log({'Aborted coil channel combination.'},{'No "chan" dimension.'})
    return; 
end

% Create backup % ------------------------------------------- %

% Create a backup of the current data set.
button_setBackup_Callback([], [], gui)
CSI_Log({'WSVD; Created a backup before combining channels.'},...
               {'Use the backup get button to revert back.'});

% User input % ----------------------------------------------- %

% Get user input for noise component
uans = getUserInput({'Use the noise prescans? (y/n):','Noise mask size:',...
                     'Exclude channels:'},{'y',50,''}); 
if isempty(uans), CSI_Log({'Skipped combining channels.'},{''}) ; return; end
                
% Get num mask size and channels to exclude.
mask_size = str2double(uans{2}); 
% Use prescans or create a noise mask per voxel. 
% If mask == 1, a noise mask is used and the prescans are ignored.
if strcmp(uans{1}, 'y')
    if isfield(csi.data, 'noise')
        mask = 0; 
        CSI_Log({'WSVD; Using the full noise pre-scans'},...
                       {'to calculate noise covariance matrix'});
    else
        mask = 1;
        CSI_Log({'WSVD; Prescans unavailable.'},...
                       {'Using mask to calculate noise covariance matrix.'});
    end
else
    mask = 1;
end
CSI_Log({'WSVD; Size of noise mask:'},{mask_size});



% Number of samples and channels.
nchan = size(csi.data.raw,ind_cha); ndimf = size(csi.data.raw,1);

% Exclusion of channels % ------------------------------------------- %

% Convert user input.
ch_excl = str2double(strsplit(uans{3},' ')); 
if isnan(ch_excl),ch_excl = []; end
% Create included channels array
ch_incl = 1:size(csi.data.raw, ind_cha); ch_incl(ch_excl) = [];

% Averaging % ------------------------------------------------------ %

% If data not averaged - ask user if to do so.
if ~isnan(ind_avg)
if size(csi.data.dim, ind_avg) ~= 1
    uans = getUserInput(...
        {'Average channel data before WSVD? (y/n)'},{'y'});
    if isempty(uans), uans{1} = 'n'; end
    
    % Average data
    if strcmp(uans{1},'y')
        csi.data.raw = mean(csi.data.raw,ind_avg);   % Average
        CSI_Log({'WSVD: Averaged data over dimension: '},{ind_avg});
    end
end
end

% Reshaping (1/2) % -------------------------------------------------- %

% Reorder data setting the channel data to the second dimension index.
% Data format: [dt x channels x rest ... ]

% Permute vector - channel index on index 2.
sz = csi.data.dim;                                    % Size of data array
rm_ind = 2:numel(sz); rm_ind(rm_ind == ind_cha) = []; % Vector of #dims
permv = [1 ind_cha]; permv = cat(2,permv,rm_ind);     % Permute vector.
% Permute data
tmp = permute(csi.data.raw, permv); tmp_label = csi.data.labels(permv);

% Reshape to dimensional cell array: dt x channels x nVoxels
% Array size per dimension
szr = size(tmp); 
% Add a dimension if only one voxel.
if numel(szr) < 2, szr(3) = 1; end 
% Cell layout: {dt x nchan} x nCells
cell_layout = arrayfun(@ones,...
    ones(1,size(szr(3:end),2)),szr(3:end),'UniformOutput',0);
% Convert to a cell matrix with {dt x nChan};
tmp_cell = squeeze(mat2cell(tmp, szr(1), szr(2), cell_layout{:})); 

% WSVD % ------------------------------------------------------------ %

% Number of voxels to combine
tmp_cell = reshape(tmp_cell, [], 1); % Create a list of all voxels.
nvox     = size(tmp_cell,1);

% Containers for WSVD
comb = struct; comb.data = zeros(szr(1),nvox); 
comb.qual = zeros(nvox,1); comb.ampl = zeros(nvox, size(ch_incl,2)); 
comb.W = zeros(size(ch_incl,2), nvox); 

% phaseRefChannel
% Debug
% elseif isfield(options,'noiseCov') && strcmp(options.noiseCov,'diag')
% Use only the noise variances (not off diagonal elements)...
% noiseCov=diag(diag(cov(rawSpectra(noiseMask,:))));

% WSVD loop. 
% Apply for every indices excluding the channel index: e.g. every voxel.
for vi = 1:nvox
    % Get spectrum and exclude channels given by user.
    tmp_spec = tmp_cell{vi,1}; tmp_spec = tmp_spec(:,ch_incl);
    
    % Create noise coVariance matrix.
    if mask == 1 % Use mask
        noiseMask = ndimf-mask_size:ndimf-1;
        noiseCov = cov(tmp_spec(noiseMask,:)); 
    else         % Use noise pre-scans
        % FFT of noise prescans. (all channels)
        noise_spec = csi_fft(csi.data.noise);
        noiseMask = 1:ndimf; noiseCov = cov(noise_spec(:,ch_incl));              
    end

    % WSVD algorithm
    % Data, quality, coil amplitude and weights.
    [comb.data(:,vi), comb.qual(vi),comb.ampl(vi,:),comb.W(:,vi)] = ...
        wsvd(tmp_spec, noiseMask', 'noiseCov', noiseCov);
end


% Reshaping (2/2) % ------------------------------------------------- %
% Reshape the cell array back to the an array: dt x chan x residual.
% Reorder to original order as before reshape (1/2).

% Reshape to before cell array size with channel index size = 1.
% tmp = array of spectra with channel dim at 2nd index
sz_cell = size(tmp); sz_cell(2) = 1; sz_cell(1) = ndimf; 
tmp = reshape(comb.data, sz_cell); % prod(sz_cell) == numel(comb.data)!

% Reorder back to starting order of index dimensions.
% Using the label-change during the first reshaping permute, find the
% permute vector to go back to initial state.

% Check if no repeated labels are present
if size(unique(csi.data.labels)) ~= size(csi.data.labels)
    unilab = unique(csi.data.labels);
    for kk = 1:unilab
        
    end
    
end

cur2prev_permv = csi_findDimLabel(tmp_label, csi.data.labels);

csi.data.raw   = permute(tmp,cur2prev_permv);

% Update csi.data.dim data.
new_dim = NaN(1,size(csi.data.dim,2));
for kk = 1:size(csi.data.dim,2), new_dim(kk) = size(csi.data.raw,kk); end
csi.data.dim = new_dim; 

% REMAINING ISSUE: Quality; Amplitude; Weight.
% Save? Export? Use for what?

% Save data % ------------------------------------------------------- %

% Set CSI app data
setappdata(gui.CSIgui_main, 'csi', csi);

% Display info to user
CSI_Log({'WSVD; Channels are combined.',...
                'WSVD; Average quality:','WSVD; Included channels:'},...
                {'',mean(comb.qual,1),ch_incl});            
            

            
% --- Executes on button press in button_CSI_MaxValue.
function button_CSI_MaxValue_Callback(~, ~, gui)
% Calculate the max for each voxel per slice (Dim 2,3, and 4). Data is
% shown as graph per slice with higher index dimensions represented as a
% line for each voxel or as a table per slice including every higher 
% index dimension per slice-table, seperated by its index.

% Get option from user: per slice or in 3D.
uans = getUserInput_Popup({'Select display of maxima option: '},...
                         {{'Maximum per slice','Maximum in 3D'}});                     
if isempty(uans), CSI_Log({'Skipped calculating data maximum.'},{''}) ; return; end                
 
% Launch maximum function
switch uans{1}
    case 'Maximum per slice'
        CSI_Max_Per_Slice([], [], gui);
    case 'Maximum in 3D'
        CSI_Max_In_3D([], [], gui);
end

% --- Executed by MaxValue button callback
function CSI_Max_Per_Slice(~, ~, gui)
% Calculate the maximum for each voxel and display data as either a graph
% or a table.
%
% Uses CSI_dataAsGraph(); CSI_dataAsTable();

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% USERINPUT % --------------------------------- %
uans = getUserInput_Popup({'Display type: ',...
                           'Select peak: '},...
                          {{'Graph', 'Table'},{'Yes','No'}});
if isempty(uans)
    CSI_Log({'Skipped max/slice.'},{''}) ; return; 
end

% User requested data dsisplay
dataDisp = lower(uans{1});

% User requested peak selection
switch uans{2}
    case 'Yes', selectPeak = 1;
    case 'No',  selectPeak = 0;
end

% FID & Echo detection
combine_fe = 0;
if isfield(csi.data,'split')
    % From user: Calculate using FID and echoes?
    uans = getUserInput_Popup({'Use both FID and Echo data?'},...
        {{'Yes', 'No'}});
    if isempty(uans), CSI_Log({'Skipped max/slice.'},{''}) ; return; end
    
    % Combine fid and echo yes or no?
    switch uans{1}, case 'Yes', combine_fe = 1; end
end




% Data % -------------------------------- %

% Get the data and apply data unit.
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};

% SELECT PEAK% -------------------------------- %

if selectPeak
    [doi, ~, range] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
else
    doi   = csi.data.raw;
    range = [csi.xaxis.none(1) csi.xaxis.none(end)]; 
end

% MAX/VOXEL % -------------------------------- %


if combine_fe
    
    doi = csi.data.raw; 
    switch csi.data.split.type
        case 'FID', active = 'fid';  stored = 'echo';
        case 'Echo',active = 'echo'; stored = 'fid';
    end
    
    % Get maximum for stored and active type
    data_max_active = max(CSI_getUnit(doi,data_unit),[],1);
    data_max_stored = ...
        max(CSI_getUnit(csi.data.split.(stored),data_unit),[],1);

    
    % Concatenate data
    switch active
        case 'fid',  data_max = cat(5, data_max_active, data_max_stored);
        case 'echo', data_max = cat(5, data_max_stored, data_max_active);
    end
    
else
    % Get data in unit
    data = CSI_getUnit(doi, data_unit);
    % Calculate the maximum values/voxel.
    data_max = max(data,[],1);

end



% DISPLAY INFO % ----------------------------- %
CSI_Log({['Maxima calculated. Displaying data as ' dataDisp '.']},...
               {'Please be patient.'});
           
% DISPLAY DATA % ----------------------------- %
switch dataDisp
    case 'graph'
        % Display data % ------- % Graph.
        CSI_dataAsGraph(data_max, gui, 'Maxima');
    case 'table'
        % Display data % ------- % Table.
        CSI_dataAsTable(data_max,'Maxima');
end

% --- Executed by MaxValue button callback
function CSI_Max_In_3D(~, ~, gui)
% Caululate maximum for each voxel and launch the max3D plot function. 
% The maximum values will be normalized to allow propper display of the
% scatter 3D plot.
%
% Uses CSI_max3D();

% Get CSI app-data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Slice data: Get data-unit to plot
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};
data = CSI_getUnit(csi.data.raw, data_unit);

% Launch max 3D function
CSI_max3D(data,[gui.colors.main;gui.colors.text_main]);

% --- To visualise 3D scatter plot of a N-dimensional data set.
function CSI_max3D(data,clrs)
% Plot the maximum of each point in data as a 3D scatter. Data values are
% normalized to the maximum in the full data set.
% Each dimension higher than index(3) is plotted as a seperate tab.
% 
% Input expected as x,y and z spatial indices to be on index/dimension 
% 1,2 and 3. 
% Clrs is optional; (2x3) with row(1) the background and row(2) text color.

% Check input arguments 
if nargin == 1,clr_txt = [0.94 0.94 0.94]; clr_bg = [0 0 0];
else, clr_txt = clrs(2,:); clr_bg = clrs(1,:);
end

% Get maximum for each fid/spec entry: dim Dt x .. x .. x .. etc.
data = (max(data,[],1)); 
% Remove first index only by permuting to the end.
data = permute(data,[2:numel(size(data)) 1]);
% Normalize
data = data./max(data(:));

% Spatial dimensions are index(1:3);
dim = size(data);
% Find any higher order index dimensions
if numel(dim) > 3, dimrem = dim(4:end); else, dimrem = 1; end

% For each spatial index (X,Y,Z) plot a 3D scatter of the max values.
rem_index = slice2rowIndex(num2cell(dimrem)); % All additional indices

figh = figure(); 
tabg = uitabgroup(figh);

% Create plotgrid: identical for all
y = 1:size(data,1); x = 1:size(data,2); z = 1:size(data,3);
[x,y,z] = meshgrid(x,y,z);

tabh = cell(1,size(rem_index,1)); axh = cell(1,size(rem_index,1));
for ri = 1:size(rem_index,1)
    
    % Get plot data
    ind = num2cell(rem_index(ri,:));
    tmp = data(:,:,:,ind{:}).*1000;
    
    % Add tab and axes in tab
    indm = rem_index(ri,:); 
    title_str = strjoin(strsplit(num2str(indm),' '),'/');
    tabh{ri} = uitab(tabg, 'Title',title_str,...
            'BackgroundColor',clr_bg, 'ForegroundColor',clr_txt./3);
    axh{ri}  = axes('parent',tabh{ri},'Position',[0.1 0.1 0.8 0.8]);
    
    % Plot scatter
    tmp(tmp == 0) = eps; tmp(~isfinite(tmp)) = eps;tmp(isnan(tmp)) = eps;
   % tmp = tmp + abs(min(tmp(:))).*1.01;
    scatter3(axh{ri}, x(:), y(:),z(:),abs(tmp(:)),(tmp(:)), 'filled');
    
    % Axes cosmetics
    lim_sz = size(tmp); lim_sz = lim_sz+1;
    set(axh{ri},'XColor',clr_txt,'YColor',clr_txt,...
                'ZColor',clr_txt,'Color',clr_bg, 'GridColor', clr_txt,...
                'YLim',[0 lim_sz(1)],'XLim',[0 lim_sz(2)],...
                'ZLim',[0 lim_sz(3)]);
    xlabel('Columns','Color',clr_txt);ylabel('Rows','Color',clr_txt); 
    zlabel('Slices', 'Color',clr_txt);
    

end

% --- Executes on button press in button_CSI_SNR.
function button_CSI_SNR_Callback(~, ~, gui)
% Calculate SNR for each dimension in the dimensional MRS data.
%
% Uses csi_SNR();

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

           % --------------- % USERINPUT % --------------- %
           % POI, mask-size, snr-method & display-method

% POI: Peak of SNR
range_none = CSI_getPeakOfInterest(csi.xaxis);
if isempty(range_none), return; end

% SNR Noise mask
uans = getUserInput({'Size of SNR noise mask: '},{50});
if isempty(uans), CSI_Log({'Skipped SNR calculations.'},{''}) ; return; end

% Convert user input
mask_size = str2double(uans{1});

% SNR-method (real/magnitude) and display method (table or graph)
uans = getUserInput_Popup({'SNR Signal unit: ','Display type: '},...
                         {{'Real', 'Magnitude'},{'Graph','Table'}});
if isempty(uans), CSI_Log({'Skipped SNR calculations.'},{''}) ; return; end

dataDisp = lower(uans{2}); % Display method
switch uans{1}             % SNR method
    case 'Real', SNRmethod = 1; case 'Magnitude', SNRmethod = 0; end

              % --------------- % SNR % --------------- %

% Display Info %
CSI_Log({'Calculating SNR per voxel, please be patient.'},{''});

% Calculate using noise mask
SNR_all = csi_SNR(csi.data.raw, mask_size, SNRmethod, range_none);
% Convert NaNs to zero
SNR_all(isnan(SNR_all)) = 0; 

% Display Info %
CSI_Log({'Starting data display:',...
         'SNR calculated for each voxel. Overall average: '},...
        {dataDisp,nanmean(SNR_all(:))});

           % --------------- % DISPLAY DATA % --------------- %

% Table or graph
switch dataDisp
    case 'table'      % Table % ----- %
    
    % Send to tableData function, to show each slice as a table.
    % If higher dimensions are available, the data will be concentonated
    % per slice. Graph may be better suited.
    %CSI_tableData(SNR_all, 'SNR');
    CSI_dataAsTable(SNR_all, 'SNR')
     
    % Save full table to file.
    %     save([datestr(now,'YYmmdd_HHMMSS') '_SNR.mat'],...
    %         'SNR_all','ppm_range','mask_size');
    
    case 'graph'      % Graph % ----- %
    % Send to displayData function, to show each slice as a graph.
    % If higher dimensions only represent only 1 value/voxel, a dot is 
    % shown. Adviced is using table in those cases.
    CSI_dataAsGraph(SNR_all, gui, 'SNR');
end

% --- Executes on button press in button_ws: export appdata to workspace.
function button_ws_Callback(~, ~, gui)
% Export appdata from CSIgui to workspace; mainly used for debugging.

appdat = getappdata(gui.CSIgui_main);
assignin('base', 'appdat', appdat);

% --- Executes on button press in button_Info.
function button_Info_Callback(hObject, eventdata, gui)

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Update LOG
log = get(gui.listbox_CSIinfo,'String');
csi.data.log = char(log);

% Send to workspace
assignin('base', 'csi',csi.data);
csi.data
fprintf('Variable "csi" accessible from workspace.\n');


% --- Executes on button press in button_CSI_Peak_Map.
function button_CSI_Peak_Map_Callback(~, ~, gui)
% Plot a map of a specific peak maximum including images and voxel grid.


% INITIATE % --- %

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get peak of interest
[doi, doi_axis, range] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
doir = real(doi); doim = imag(doi);


% MAP % --- %

% Maximum positions and values: e.g. Map
map = max(doir, [], 1);
% Normalize map
nfac =  max(map(:)); nmap = map./nfac;


% FIGURE wIMAGES % --- %
tmp = MRI_plotImage_tabbed(gui,'CSI_Map');
obj = tmp.fig;

% ADD VOXELS % --- %
MRI_plotImage_tabbed_addVoxels(obj);
tgui = guidata(obj); plot_par = tgui.plot_par;

% Map 2 Colors: Calculate map range
clr_map = jet(128);
clr_val = linspace(0,1,size(clr_map,1));

% Loop each tab of figure
for sli = 1:size(tgui.tabh,2)                   % Sli loop.
    for ci = 1:plot_par.dim(1)                  % Col loop.
        for ri = 1:plot_par.dim(2)              % Row loop.

            % VOXEL VALUE
            vox_val = nmap(:,ci,ri,sli);
            [~, clr_ind]= min(abs(clr_val-vox_val));
            plot_par.ax{ri,ci,sli}.Color = [clr_map(clr_ind,:) 0.33];            

        end 
    end 
end

% Plot colorbar
colorbarQ(clr_map, clr_val);


function MRI_plotImage_tabbed_addVoxels(obj)
% Use MRI_plotImage_tabbed(gui,tag) to create the appriopriate
% figure and allow this function to add an axis for each voxel in all tabs.

% ADD VOXELS % --- %

% Get GUI data of figure
tgui = guidata(obj);
% Plot data for each tab: voxel grid and more plot settings.
plot_par = tgui.plot_par;

ax = cell([plot_par.dim, size(tgui.tabh,2)]);
% Loop each tab of figure
for sli = 1:size(tgui.tabh,2)                   % Sli loop.
    % Plot csi voxel per axis in plot_par.grid
    for ci = 1:plot_par.dim(1)                  % Col loop.
        for ri = 1:plot_par.dim(2)              % Row loop.

            % AXIS DETAILS 
            % X and Y position in the figure of the axis to plot.
            x   = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
            % Position of axis to plot
            pos = [x y plot_par.res(1) plot_par.res(2)];
            % Create axis with pos(3,4) size at pos(1,2) position
            if ~ishandle(tgui.tabh{sli}), return; end
            ax{ri,ci,sli} = axes('parent',tgui.tabh{sli},'position',pos);           

            % AXIS COSMETICS
            set(ax{ri,ci,sli},...
                   'Color','None',...
                   'XColor', plot_par.colors.main,...
                   'YColor', plot_par.colors.main,...
                   'LineWidth', 1.7, 'Xtick',[], 'Ytick', [],...
                   'TickLength',[0 0.00001], 'Box', 'off');             

        end 
    end 
end
plot_par.ax = ax; tgui.plot_par = plot_par; guidata(obj,tgui);

% --- Executes on button press in button_CSI_phaseShift.
function button_CSI_phaseShift_Callback(hObject, eventdata, gui, backup)
% Add function csi_frequencyShift
if nargin < 4, backup = 1; end

% Create backup
if backup, CSI_backupSet(gui, 'Before phase shift'); end

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get shift from user
uans = getUserInput({'FID Acquisition delay (ms)'},{0.05});
if isempty(uans), CSI_Log({'Skipped phase shift.'},{''}) ; return; end
dt = str2double(uans{1});

% Shift FIDS
csi.data.raw = csi_frequencyShift(csi.data.raw, dt);

% CLEAN UP % ---------- %

% Save data
setappdata(gui.CSIgui_main,'csi', csi);

% Update LOG
CSI_Log({['Applied a phase shift of ' num2str(dt) 'ms to each FID']},{''});

% --- Executes on button press in button_CSI_VoxelShift.
function button_CSI_VoxelShift_Callback(~, ~, gui, backup)
% Shift CSI K-space to spatially shift the volume a number of voxels

if nargin < 4 , backup= 1; end
% Create backup
if backup, CSI_backupSet(gui, 'Before voxel shift'); end

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% USER INPUT % ------------------- %

% Get shift from user
uans = getUserInput({'Number of voxels to shift: '},{[0.5 0.5 0.5]});
if isempty(uans), CSI_Log({'Skipped voxel shift.'},{''}) ; return; end
shift = str2double(strsplit(uans{1},' '));

% Find the spatial dimensions/indexes: ask if not found.
k_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
k_dim(isnan(k_dim)) = [];
if isempty(k_dim)
    k_dim = getUserInput({'Spatial index in MRSI data? (kx/ky/kz): '},...
                            {'2 3 4'});
    if isempty(k_dim), CSI_Log({'Skipped voxel shift.'},{''}) ; return; end
    k_dim = str2double(strsplit(k_dim{1},' '));
end

% SHIFT DATA % ------------------- %

% Check data domain type availability: Need time domain data
domain = CSI_getDomain(gui);
if strcmp(domain,'freq'), csi.data.raw = csi_ifft(csi.data.raw); end

% Shift the data
[csi.data.raw, ~] = csi_voxelshift(csi.data.raw, shift, k_dim);

% Check data domain type availability: Need time domain data
domain = CSI_getDomain(gui);
if strcmp(domain,'freq'), csi.data.raw = csi_ifft(csi.data.raw); end


% CLEAN UP % -------------------- %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% LOG
CSI_Log({'Applied phase change to shift voxels. Shifted by:'},...
        {strjoin(strsplit(uans{1},' '),' | ')});
                
% --- Executes on button press in button_CSI_FidOrEcho.
function button_CSI_FidOrEcho_Callback(hObject, eventdata, gui)
% Split data to FID and Echo data, specially design for use with AMESING
% data. Echoes and FIDs require different processing resulting in different
% sample sizes.

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if ~isfield(csi.data,'split')
    
                    % ------- % Cut FID/ECHO % ----------%
    
    % Create backup
    CSI_backupSet(gui, 'Before FID/Echo split');                    
                    
    % Userinput % ------------ %

    % Get FID/Echo dimension
    uansDim = getUserInput_Popup(...
    {'Which dimension represents the FID and echoes?'},{csi.data.labels});
    if isempty(uansDim), CSI_Log({'Skipped FID & Echo split.'},{''}); return; 
    end
    
    % FID/Echo dimensions and size
    fidecho_dim = strcmp(csi.data.labels, uansDim{1});
    dim_sz = csi.data.dim(fidecho_dim); 

    % Get FID index @ FID/Echo dimension
    uansInd = getUserInput_Popup({'Which index in this dimension is the FID?'},...
                          {{1:dim_sz}});
    if isempty(uansInd), CSI_Log({'Skipped FID & Echo split.'},{''}); return; 
    end
    
    % FID index in FID/Echo dimension
    fid_ind = str2double(uansInd{1});                 

    % Split and save % ---------- %
    
    % Data and size
    data = csi.data.raw; sz = size(csi.data.raw);

    % All dimensions indexing
    split_ind = arrayfun(@(x) 1:x, sz,'UniformOutput',0);

    % To cut the fid
    split_ind_fid = split_ind; split_ind_fid{fidecho_dim} = fid_ind;

    % To cut the echo
    split_ind_echo = split_ind;tmp = split_ind{fidecho_dim}; tmp(fid_ind) = [];
    split_ind_echo{fidecho_dim} = tmp;

    % Cut FID and ECHO
    fid = data(split_ind_fid{:}); echo = data(split_ind_echo{:});
    
    % Save data
    csi.data.split.echo = echo; csi.data.split.fid = fid;
    
    % Save xaxis
    csi.data.split.xaxis.echo = csi.xaxis;
    csi.data.split.xaxis.fid = csi.xaxis;
    
    % Set FID as active split-data set.
    csi.data.raw = fid; csi.data.split.type = 'FID';

    % Log message
    log_msg = 'FID and echoes are split up. In memory:';
else

               % --------- % Switch FID/ECHO  % ---------- %
               
    switch csi.data.split.type
        case 'FID' % FID is loaded
            
            % Backup FID data + axis
            csi.data.split.fid = csi.data.raw; 
            csi.data.split.xaxis.fid = csi.xaxis;
            
            % Set ECHO data
            csi.data.raw = csi.data.split.echo; 
            csi.data.split.type = 'Echo'; 
            
            % Set ECHO xaxis
            csi.xaxis = csi.data.split.xaxis.echo;
            
        case 'Echo' % Echo is loaded
            
            % Backup ECHO data + axis
            csi.data.split.echo = csi.data.raw; 
            csi.data.split.xaxis.echo = csi.xaxis;
            
            % Set FID data
            csi.data.raw = csi.data.split.fid;
            csi.data.split.type = 'FID'; 
            
            % Set FID xaxis
            csi.xaxis = csi.data.split.xaxis.fid;
            
    end

    % Update data dimensions
    csi.data.dim = size(csi.data.raw);
               
    % Log message
    log_msg = ['FID and echo data swapped. In memory: '];
end


              % ------------ % CLEAN UP % ------------ %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% Update LOG.
CSI_Log({log_msg}, {csi.data.split.type});

% --- Executes on button press in button_Log_DeleteLine.
function button_Log_DeleteLine_Callback(hObj, eventdata, gui)
CSI_Log_deleteLine(hObj);



% --- Executes on button press in button_CSI_AutoProcessing
function button_CSI_AutoProcessing_Callback(hObj, ~, gui)
% Automatically process MRS data.
CSI_AutoProcessing_initiate(gui);

% --- Executes on button press in button_CSI_AutoProcessing
function CSI_AutoProcessing_initiate(gui)
% Automatically process data using the following methods
%        'Average k-space.', ...                1
%        'Apodize k-space.', ...                2
%        'Spatial FFT.', ...                    3
%        'Load protocol text- or spar-file.'    4 % OFF by default
%        'Set parameters: frequency.', ...      5 
%        'Set parameters: geometry.', ...       6 % OFF by default
%        'Apodize FID.', ...                    7
%        'Zero fill FID', ...                   8
%        'FID to Spectrum (Forward FFT).', ...  9
%        'Automatic zero-phasing', ...          10
%        'Combine channels.', ...               11 
%        'Automatic zero-phasing'               12
%        'Convert MRI to CSI space.'            13 % OFF by default.
% 
% Before starting automatic processing, specific processing steps can be
% enabled and/or disabled. Order shown is order applied.

% Backup.
CSI_backupSet(gui,'Initiated auto-processing.');

% ---------------------------------------------- %
%                 Do Not Edit Order              %
% The following options are executed automatically.
str = {'Average k-space.', ...
       'Apodize k-space.', ...
       'Spatial FFT.', ...
       'Load protocol text- or spar-file.',...
       'Set parameters: frequency.', ...      
       'Set parameters: geometry.', ...  
       'Apodize FID.', ...
       'Zero fill FID', ...
       'FID to Spectrum (Forward FFT).', ...
       'Automatic zero-phasing', ...
       'Combine channels.', ...
       'Automatic zero-phasing.',...
       'Convert MRI to CSI space.'};
%                 Do Not Edit Order              %
% ---------------------------------------------- %
 
   
% Find any open CSI_auto windows and close.
obj = findobj('Type', 'Figure', 'Tag', 'CSI_Auto');
if ~isempty(obj), delete(obj); end

% Create CSI_Auto
fh = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                      'NumberTitle', 'Off', 'Color', gui.colors.main,...
                      'Name','CSIgui - Automatic Processing','Tag', 'CSI_Auto'); 
gdat = guidata(fh); axis off;


% FIGURE POSITION 
% Based on nr of elements in str; e.g. every automatic processing step.

% Each process stepp has 2x text heights which include a dh (offset)
% per text and a line for seperation.
dw = 5; dh = 4; szt = [230 18]; szl = round(dh/2); % DEFINES FIGURE SIZE

% pos x and y and w h for each object in a single element
% [txt; nfo; line];
step = szt(2)*2 + dh*2 + szl; button_space = 15;
w = szt(1) + dw*2; h = step*size(str,2) + button_space;
szc = gui.CSIgui_main.Position; 

% Set figure position
figpos = round((szc(3:4)./2) - ([w h]./2)) + szc(1:2);
figpos(1) = figpos(1) + szc(3)/2 + w/2;
fh.Position = [figpos w h];

% Add processing items % ---------------- %
nelem = size(str,2);
txt = cell(1,nelem);nfo = cell(1,nelem);
rad = cell(1,nelem);lyn = cell(1,nelem);
for kk = 1:nelem
 
    % Text label
    txt{kk} = uicontrol('Style', 'text', 'String', str{kk},...
    'ForegroundColor',gui.colors.text_title,...
    'BackgroundColor',gui.colors.main,...
    'HorizontalAlignment', 'Left');
    % step = ((dh + szt(2)) * (kk)) + ((dh + szt(2))*(kk-1));
    step = 2 * kk * ( dh + szt(2) ) - dh - szt(2);
    txt{kk}.Position = [ dw (h - step)  szt];
    txt{kk}.FontWeight = 'Bold';txt{kk}.FontSize = 10;
    
    % Radio
    rad{kk} = uicontrol('Style', 'radiobutton', 'String', '',...
    'BackgroundColor',gui.colors.hilight2,...
    'ForegroundColor',gui.colors.main,...
    'HorizontalAlignment', 'Left');
    step = 2 * kk * ( dh + szt(2) ) - dh - szt(2);
    rad{kk}.Position =  [ w-18 (h - step)  18 18];
    if kk == 4 || kk == 6 ||  kk == 13 % DEFAULT OFF OPTIONS
        rad{kk}.Value = 0; 
    else
        rad{kk}.Value = 1;
    end
    
    % Info label
    nfo{kk} = uicontrol('Style', 'text',...
    'String', '...',...
    'ForegroundColor',gui.colors.text_main,...
    'BackgroundColor',gui.colors.main,...
    'HorizontalAlignment', 'Right');
    step = (dh + szt(2))* 2 * (kk);
    nfo{kk}.Position = [ dw (h - step)  szt];
    
    % Line
    lyn{kk} = uicontrol('Style', 'text', 'String', '',...
    'BackgroundColor',gui.colors.hilight2,...
    'ForegroundColor',gui.colors.main,...
    'HorizontalAlignment', 'Left');
    step = ((dh + szt(2))* 2 * (kk))+szl*2;
    lyn{kk}.Position = [ dw (h - step)  szt(1) szl];

end

% Start button.
uicontrol('style', 'pushbutton', 'string', 'Start',...
    'Position', [dw dh 75 button_space],...
    'BackgroundColor',gui.colors.main,...
    'ForegroundColor',gui.colors.text_main,...
    'Callback', @CSI_AutoProcessing_execute);

% Store GUI stuff,.
gdat.lyn = lyn; gdat.nfo = nfo; gdat.txt = txt; gdat.rad = rad;
gdat.str = str; gdat.fh = fh;gdat.csigui = gui.CSIgui_main;
guidata(fh,gdat);

% --- Executes by Start-button in CSI_AutoProcessing
function CSI_AutoProcessing_execute(hobj, ~)
%        'Average k-space.', ...                1
%        'Apodize k-space.', ...                2
%        'Spatial FFT.', ...                    3
%        'Load protocol text- or spar-file.'    4 % OFF by default
%        'Set parameters: frequency.', ...      5 
%        'Set parameters: geometry.', ...       6 % OFF by default
%        'Apodize FID.', ...                    7
%        'Zero fill FID', ...                   8
%        'FID to Spectrum (Forward FFT).', ...  9
%        'Automatic zero-phasing', ...          10
%        'Combine channels.', ...               11 
%        'Automatic zero-phasing'               12
%        'Convert MRI to CSI space.'            13 % OFF by default.
gdat = guidata(hobj); nfo = gdat.nfo; str = gdat.str; rad = gdat.rad;
hObj = gdat.csigui; gui = guidata(hObj);

bu = 0; nelem = size(str,2);
for kk = 1:nelem
    if rad{kk}.Value == 1
        if     kk == 1
            % Average k-space
            button_CSI_Average_Callback([], [], gui, bu);
        elseif kk == 2
            % Apodization k-space
            button_CSI_Apodization_Kspace_Callback([], [], gui, bu);
        elseif kk == 3
            % FFT Spatial
            button_CSI_FFT_Kspace_Callback([], [], gui, bu);
        elseif kk == 4
            button_CSI_ReadInfo_Callback([], [], gui)
        elseif kk == 5
            % Set frequency parameters
            CSI_2D_Scaling_calc_xaxis(hObj,[]);
        elseif kk == 6
            % Set geometry parameters
            button_CSI_setCoordinates_Callback([], [], gui)
        elseif kk == 7
            % Apodization FID
            button_CSI_Apodization_FID_Callback([],[], gui, bu);      
        elseif kk == 8
            % Zerofill
            button_CSI_ZeroFill_Callback([], [], gui, bu);
        elseif kk == 9
            % FFT time2spec
            button_CSI_FFT_Callback([], [], gui, bu);
        elseif kk == 10
            % Autophase
            button_CSI_AutoPhase_Callback([], [], gui, bu);
        elseif kk == 11
            % Combine channels
            CSI_Combine_WSVD;
        elseif kk == 12
            % Autophase
            button_CSI_AutoPhase_Callback([], [], gui, bu);
        elseif kk == 13
            % Convert MRI to CSI space.
            MRI_to_CSIspace(gui);
        end
    
        % Show executed in window
        nfo{kk}.String = 'Succesfully executed.';
        nfo{kk}.ForegroundColor = [0 1 0];
    else
        nfo{kk}.String = 'Skipped.';
        nfo{kk}.ForegroundColor = [0.65 0.1 0.1];
    end
end


% Display data functions % --------------------------------------------- %
% Data can be visualized either as a graph or as a table.
% Specific function to display calculated data for all slices, in a
% single figure, with a tab for each slice in the array. 
% 
% E.g. Visualize data wrt spatial information.


% ---- Executed by functions to graphically display data/slice
% Uses cell input to overlay multiple lines
function CSI_dataAsGraph2(data, data_text, datatag, gui, img)
% Display data as a line graph according using the indexing as a spatial
% grid for plotting.
%
% data = a cell with a N-dimensional array. The 1st,2nd and 3rd dimension
% are interperted as the X,Y and Tab index, other higher indices are
% plotted as seperate lines. 

% If no data tag is present.
if nargin == 2, datatag = 'Data'; end

% Plot image data
if nargin == 5, plotImg = img; else, plotImg = 0; end

% Get colors for plotting
clr_bg  = gui.colors.main;       % Background
clr_tx1 = gui.colors.text_main;  % In plot
clr_tx2 = gui.colors.text_title; % In tab
clr_li1  = gui.colors.lines1;    % Line color
clr_li2  = gui.colors.lines2;    % Unused

            % ----------- % Prepare Figure % ----------- %

% Figure for all axis objects in this slice.
fh = figure('Tag', ['CSIgui_plot' datatag],'Name',...
            ['CSIgui - ' datatag ],'Color', clr_bg,...
            'Toolbar', 'figure', 'MenuBar', 'None','NumberTitle', 'Off');  
        
% Height and width of figure linked to screen ratio! 
fig_W = 720; scr_sz = get(0, 'screensize'); 
% Figure size
fig_sz = [fig_W fig_W.*(scr_sz(4)/scr_sz(3))];
% Figure position on screen
fig_ps = round(scr_sz(3:4)./2) - (fig_sz./2);
% Apply
set(fh, 'Position', [fig_ps fig_sz]);

% GUI data to figure
tgui = guidata(fh); tgui.fig = fh;
        
% Tabgroup for all slices
tgui.tabgp = uitabgroup(tgui.fig,'Unit','normalized','Position',[0 0 1 1]);

% Check line data size e.g. cell content
lines2plot = size(data,2);
line_types = repmat({'none',':','--'},[1 10]);
line_mark  = repmat({'s','.','o','+','x','d'},[1 10]);


                % -------- % Images to Plot % --------- %
if plotImg, img = MRI_matchSlices(gui.CSIgui_main); end



               % -------- % Axes Plot grid % -------- %
             
% Explanation: Get dimension/index sizes and process.
% Use dim(2) and dim(3) as X/W Y/H, Use dim(3) as Z/SLICE, Use dim(...) 
% as dimensionality (Line). Plot_para.dim = [x y slice others...];

% Create a plot/axis grid for each voxel
plot_par.dim      = size(data{1});      % Data dimensions
plot_par.dim(1)   = [];                 % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);% 3D/2D/1D volume. 

% 1D Correction: % ------------------------------------ %
% For 1D data, the correct axis locations arent calculated. To correct, a
% second dimension is inserted and the 1D data is relocated to the Y e.g.
% column dimension. 
% REMARK: 
%       The 2nd dimension shown, equals x thus the COLUMNS
%       The 3rd dimension shown, equals y thus the ROWS
% Set in plot_par.dim the 1D dimension to the second index e.g. column.
if plot_par.data_dim == 1
    % Add a Y dimension.
    plot_par.dim(2) = 1; plot_par.data_dim = 2;  
end

% Axis Resolution: Relative to figure (Normalized) e.g. [1 1];
plot_par.res = 1./plot_par.dim; 

% Grid with locations for each axis to plot to in the figure.
for axi = 1:numel(plot_par.dim)
    plot_par.range{axi} =   ...
    0 : plot_par.res(axi) : ...
       (plot_par.res(axi) * plot_par.dim(axi)-plot_par.res(axi));
end
[x,y] = meshgrid(plot_par.range{1},plot_par.range{2});
plot_par.grid.x = (x); plot_par.grid.y = flipud(y);

% 1D Correction: if 1D transpose x/col and y/row
if plot_par.data_dim == 1, plot_par.grid.y = y'; plot_par.grid.x = x'; end

% Create a cell array for other non-spatial dimension index. 
if numel(plot_par.dim) >3, nDimC = num2cell(plot_par.dim(4:end));
else, nDimC = {1};
end
% To linear vector per cell.
nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

% Get plot color method 
clrScale_vol = get(gui.menubar.MRSI.ColorScale.Volume, 'Checked'); 
clrScale_sli = get(gui.menubar.MRSI.ColorScale.Slice,  'Checked'); 
clrScale_sta = get(gui.menubar.MRSI.ColorScale.Static, 'Checked');


                  % -------- % Slice Loop% --------- %
                  
for sli = 1:size(data{1},4)
    
              % ------- % TAB: create tab/slice % -------- %
                
% Add tab-handle according nr of thumbnails!
tgui.tab{sli}.tabh = uitab(tgui.tabgp ,'Title',['Slice ' num2str(sli)],...
      'BackGroundColor',clr_bg ,'ForegroundColor', clr_tx2);   
   
  
                % ------- % Plot color method % -------- %
  
% Plot color method applies to first line in data-set only 
% E.g. first cell used to calculate the line color
if strcmp(clrScale_vol,'on')       
    
    % Take full volume/data 
    data_ylimits_color = [min(data{1}(:)) max(data{1}(:))];
    
elseif strcmp(clrScale_sli,'on')  
    
    % Get minimum and maximum of slice sli
    tmp = data{1}(:,:,:,sli,nDimC{:}); 
    data_ylimits_color = [min(tmp(:)) max(tmp(:))];
    
elseif strcmp(clrScale_sta, 'on')
    
    % Color scaling limited to 1 color: static line color
    data_ylimits_color = NaN;
    
end

% Get plot colors range for different max-limits.
if ~isnan(data_ylimits_color)
    
    % Check limits agrees with rules: lim(1) < lim(2)
    if data_ylimits_color(2) <= data_ylimits_color(1),...
            data_ylimits_color(2) = data_ylimits_color(1)+1; 
    end

    [clrs, clrs_data_range] = CSI_2D_Scaling_calc_ColorOfPlots(data_ylimits_color);
else
    clrs = clr_li1; 
    tmp = data{1}(:,:,:,sli,nDimC{:}); clrs_data_range = max(tmp(:));
end

       % ---------- % Plot y-axis limit method  % ---------- %

% Y-limit of axis scaling by VOXEL, SLICE OR VOLUME.
axScale_vox = get(gui.menubar.MRSI.AxisScale.Voxel, 'Checked');
axScale_sli = get(gui.menubar.MRSI.AxisScale.Slice, 'Checked');
axScale_vol = get(gui.menubar.MRSI.AxisScale.Volume,'Checked');

% If voxel, set scale to zero, slice to one and volume to two. For slice
% and volume, already calculate the limits. 
% For coloring - also voxel scaling needs ylimit of slice as this is
% relative to the slice minima and maxima.
if strcmp(axScale_vox,'on')
    scale = 0; 
elseif strcmp(axScale_sli,'on') 
    scale = 2;
    % Plot data of the whole slice for ylimit scaling and coloring
    % From each slice for each data set
    tmp_mi = cell(1,size(data,2)); tmp_ma = cell(1,size(data,2));
    for kk = 1:size(data,2)
        tmp_mi{kk} = min(data{kk}(:,:,:,sli,nDimC{:}));
        tmp_ma{kk} = max(data{kk}(:,:,:,sli,nDimC{:}));
    end
    % Get minimum and maximum
    tmp_mi = cell2mat(tmp_mi); tmp_ma = cell2mat(tmp_ma);
    data_ylimits_axis = [min(tmp_mi(:)) max(tmp_ma(:))];
elseif strcmp(axScale_vol,'on') 
    % Take full volume of the data and calc min and max
    scale = 2; tmp = cell2mat(data);
    data_ylimits_axis = [min(tmp(:)) max(tmp(:))];  
end

             % ---------- % Plot Image % ---------- %
if plotImg && ~isnan(sum(img(:)))
    
    conv_data = getappdata(gui.CSIgui_main,'conv'); % Conv data
    img2plot = img(:,:,sli);                        % Get image to plot.
    
    % Create axis for image
    hold on; 
    imax = axes('parent',tgui.tab{sli}.tabh,...
        'Position',[0 0 1 1], 'Color', 'None');
    
    % Plot Images
    if (sum(img2plot(:)) == 0) % Image is only zeroes.
        colormap(gray(255)); set(imax,'Color', 'Black'); alpha(1); 
    else
        % Image plotting:
        % Imscale as it plots over the entire figure and does not
        % imply any border issues as with imshow-function.
        imagesc(img2plot, 'parent', imax); 

        % Image Contrast.
        if isfield(conv_data, 'contrast')
            caxis(imax, conv_data.contrast);
        else
            caxis(imax,[min(img2plot(:)) max(img2plot(:))*0.5]);
        end
        colormap(gray(255));
    end 
    
    % Change text-color and line-color
    clr_tx1 = clr_tx2;

end
    
   
                  % ------- %  Plot Voxel/Axis % ------- %
                  
% Be aware of the row/column switch to get proper x/y axis display in the
% figure. E.g. dim(1) = column = x, dim(2) = row = y; Otherwise this would
% be Matlab style reversed.
for ci = 1:plot_par.dim(1)     % Column loop.
    for ri = 1:plot_par.dim(2) % Row loop.
        
        % Axis Location
        % X and Y position in the figure of the axis object for this voxel
        x   = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
        % Position of axis to plot
        pos = [x y plot_par.res(1) plot_par.res(2)];
        % Create axis with pos(3,4) size at pos(1,2) position
        plot_par.ax{ri,ci} = axes('parent',tgui.tab{sli}.tabh, ...
                                  'position',pos);
        
        % For plotting - get all lines at this voxel in a variable 
        plot_data_merged = cell(1,lines2plot);
        for kk = 1:lines2plot
            plot_data_merged{kk} = data{kk}(:,ci,ri,sli,nDimC{:});
        end
        plot_data_merged = cell2mat(plot_data_merged);
                              
        % Plot each line of the voxel in this axis.
        for kk = 1:lines2plot
            % DATA 2 PLOT % ----------------------- %
            % Get voxel 2 plot data.
            plot_data = data{kk}(:,ci,ri,sli,nDimC{:});
            plot_data = reshape(plot_data,size(plot_data,1),[]);
            
            % COLOR SCALING % --------------------- %
            % This is NOT the Ylimit eventually used for axis-display.
            % Color scaling by voxel/slice/volume already in the 
            % color range variable
            ylimit = [min(plot_data) max(plot_data)];
            % Relative to maximum Y-data.
            [~, clr_ind] = min(abs(clrs_data_range-ylimit(2)));
            plot_color = clrs(clr_ind,:); % See before ri/ci for-loops.
            
            
            % Only first line/points have scaled color.
            if kk > 1, plot_color = clr_li2; end
            
            
            %%%% PLOT CSI data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot(plot_data,'color', plot_color, 'LineWidth', 1.5,...
                       'Marker',line_mark{kk},'LineStyle',line_types{kk});        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;

            
            % TEXT AND AXIS COSMETICS % ---------------- %
            if kk == 1
                

                % Axis Y-limit method
                if  scale == 0                  % Voxel
                    ylimit = [min(plot_data_merged(:)) ...
                              max(plot_data_merged(:)) ];
                elseif scale == 1 || scale == 2 % Slice or Volume
                    % Use pre-calculated limit relative to slice or volume
                    ylimit = data_ylimits_axis;    
                end

                % Axis Y-limit value
                ylimit = [ylimit(1)-(0.15*abs(ylimit(2)))  ...
                          ylimit(2)+(0.5*abs(ylimit(2)))];
                if ylimit(2) <= ylimit(1), ylimit(2) = ylimit(1)+1; end
                
                % Axis Y-limit
                ylim(ylimit);

            
                % Axes X Limit
                xlimit = [0 size(plot_data,1)+1];
                plot_par.ax{ri,ci}.XLim = xlimit; 

                % Add Text
                xpos_text = plot_par.ax{ri,ci}.XLim(end)*0.95;
                nDimC2 = repmat({1},1,numel(size(data_text)));
                text(xpos_text, ylimit(2), ...
                     data_text{:,ci,ri,sli,nDimC2{:}},...
                    'Color', clr_tx1, 'FontSize',9, 'FontWeight','Bold',...
                    'HorizontalAlignment', 'Right',...
                    'VerticalAlignment', 'Top');    

                % Axis Cosmetics 
                set(plot_par.ax{ri,ci},...
                 'Color','None','XColor', [0.4 0 0],'YColor', [0.4 0 0],...
                 'LineWidth', 1.75, 'Xtick',[], 'Ytick',[]);
            end   
            
        end % End of data-set plot loop
       
    end % for row loop
end     % for column loop
end     % for slice loop

% --- Executes by functions to graphically display data/slice
% Display all non-spatial indexes as one line
function CSI_dataAsGraph(data, gui, datatag)
% Enable display of data calculated from CSI data set. The data is shown
% per slice and higher dimensions are used for line-plot generation per
% plot. The mean, min and max of the data is shown per voxel as a text in
% the upper right corner of the axis. Datatag enables setting correct data
% name in the plot-figure title bar.
%
%           Requires csi appdata retrieved from gui.CSIgui_main
% Data should have the same dimensionality/size as csi.data.raw for proper
% functioning: 
% data-index(1) = time or 1; (2) = X; (3) = Y; (4) = Z; (5...) = other;
% 
% Used for displaying SNR; Max Value


% If no data tag is present.
if nargin == 2, datatag = 'Data'; end

% IMPORTANT.
% If 1D data, the row dimension is shifted towards higher LINE index: >=5.
% E.g. only 1 slice. 1 Voxel. nRow line points for the other dimensions. 

% Get objects and data % --------- %

% REQUIRED APPDATA: data, samples, index to plot.
csi = getappdata(gui.CSIgui_main,'csi'); 
if isempty(csi), CSI_Log({'No CSI data loaded!'},{':O'}); return;end

% Copy data
plot_par.colors = gui.colors; plot_par.xaxis = csi.xaxis;

% Dimensions % ------------------- %

% Axis size to fit figure
plot_par.dim      = size(csi.data.raw);   % Data dimensions
plot_par.dim(1)   = [];                   % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 
        
% Correction 1D
% Set in plot_par the 1D dimension to the second index e.g. column.
if plot_par.data_dim == 1, plot_par.dim(2) = 1; plot_par.data_dim = 2; end

% Create Figure % ----------------- %
plot_par = CSI_2D_setFigure(plot_par, [], datatag);

            
% Data unit % --------------------- % 
% As data is given, no additional data unit is applied and real-data is
% expected.
plot_par.data_unit = 'Real';

% GUI data to figure
tgui = guidata(plot_par.fh); tgui.fig = plot_par.fh;

% Tabgroup for all slices
tgui.tabgp = uitabgroup(tgui.fig,'Unit','normalized','Position',[0 0 1 1]);

% Loop each slice
for sli = 1:size(data,4) % ----------------- %

    % TAB: creat slice tab % --------------- %

    % Add tab-handle according nr of thumbnails!
    tgui.tab{sli}.tabh = uitab(tgui.tabgp ,...
          'Title',['Slice ' num2str(sli)],...
          'BackGroundColor','Black','ForegroundColor', 'Black');   

    % Data: get slice-data % ---------------- %  
      
    % Create cell for indexing outside the spatial-dimensions 
    if numel(plot_par.dim) >3, nDimC = num2cell(plot_par.dim(4:end));
    else, nDimC = {1};
    end
    % To linear vector per cell.
    nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);
    
    % Get data.
    plot_par.data2D = data(:,:,:,sli,nDimC{:});
    
    % Plot Scaling Settings %  ------------ %
    
    % Add Scaling options: plot color, axis-y and x limit 
    % by volume/static/voxel.
    plot_par = CSI_2D_getPlotSettings(plot_par, gui, data);

    % PLOT PER VOXEL % ------------------------ %
    for ci = 1:plot_par.dim(1)     % Column loop.
        for ri = 1:plot_par.dim(2) % Row loop.
        
        % Safety
        if ~ishandle(plot_par.fh), return; end
        
        
        % Axis details % ----------------- %
        % X and Y position in the figure of the axis object for this voxel
        x = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
        % Position of axis to plot
        pos = [x y plot_par.res(1) plot_par.res(2)]; 
        % Create axis with pos(3,4) size at pos(1,2) position @ tab.
        plot_par.ax{ri,ci} = axes('parent',tgui.tab{sli}.tabh, ...
                                  'position',pos);
        
        % Get data % ------------------------------- %
        
        % Plot data
        plot_data = plot_par.data2D(:,ci,ri,:,:,:,:,:,:,:);
        % Reshape such other dimensions then spatial are represented as one
        % line.
        plot_data = reshape(plot_data, size(plot_data,5),[]);
        
        % Set color scale % ------------------------ %
        
        % Relative to maximum Y-data.
        ylimit = [min(plot_data(:)) max(plot_data(:))];
        [~, clr_ind] = min(abs(plot_par.clrs_data_range-ylimit(2)));
        plot_color = plot_par.clrs(clr_ind,:); 
  
        
        %%%% PLOT CSI data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(plot_data,'color', plot_color, 'LineWidth', 0.005,'Marker','s');        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Yaxis Scaling % ------------------ %
        
        % Scale Y axis by limits of voxel, slice or volume.
        if     plot_par.scale == 0 % to voxel, no change.             
        elseif plot_par.scale == 1 || plot_par.scale == 2 
                % Use pre-calculated limit relative to slice or volume
                ylimit = plot_par.axScale_ylimit;    
        end
            
        % Calculate y-axis limits and apply
        ylimit = [ylimit(1)-(0.05*abs(ylimit(2)))  ...
                  ylimit(2)+(0.05*abs(ylimit(2)))];
        % Safety for the Ylimits;
        if ylimit(2) <= ylimit(1), ylimit(2) = ylimit(1)+1; end
        % Set Ylimit
        ylim(plot_par.ax{ri,ci}, ylimit);
           
        
        % Add text % ---------------------- %
        if size(plot_data,2) > size(plot_data,1)
            xlimit = [0 size(plot_data,2)+1];
        else
            xlimit = [0 size(plot_data,1)+1];
        end
        
        plot_par.ax{ri,ci}.XLim = xlimit; 
        
        xpos_text = plot_par.ax{ri,ci}.XLim(end);
        text(plot_par.ax{ri,ci}, xpos_text, ylimit(2), ...
             sprintf('Mean: %2.2g\n Min|Max: %2.2g|%2.2g',...
             nanmean(plot_data(:)), max(plot_data(:)), min(plot_data(:))),...
             'Color', [0.6 0.6 0.6],'FontSize',8,'FontWeight','Bold',...
             'HorizontalAlignment', 'Right','VerticalAlignment', 'Top');    
        
        % Cosmetics % -------------------- %
        set(plot_par.ax{ri,ci},'Color', 'None',...
                       'XColor', [0.4 0 0],'YColor', [0.4 0 0],...
                       'LineWidth', 1.75, 'Xtick',[], 'Ytick',[]);
        
        end % row loop
    end % column loop
end % slice loop

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
            'NumberTitle', 'Off', 'Resize','Off'); 
pos = get(tfh,'Position');
tgui = guidata(tfh); tgui.fig = tfh;

% Tabgroup for all slices
tgui.tabgp = uitabgroup(tgui.fig,'Unit','normalized','Position',[0 0 1 1]);

% Loop each slice - expected at dimension 4
for sl = 1:size(data,4)

    % PREP TAB % ------------------------------------------------------- %
    
    % Add tab-handle according nr of thumbnails!
    tgui.tab{sl}.tabh = uitab(tgui.tabgp ,...
          'Title',['Slice ' num2str(sl)],...
          'BackGroundColor','Black','ForegroundColor', 'Black');    

    % Table object
    tgui.tab{sl}.table = uitable(tgui.tab{sl}.tabh,...
        'Position', [30 20 pos(3:4)-60]); 
    
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

    % Convert double-array to cell-array and then to string
    tableData_str = cellfun(@sprintf, repmat({'%4.4f'},...
        [size(tableData_num,1) size(tableData_num,2)]), ...
        num2cell(tableData_num), 'UniformOutput',0);

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

    % DEV: Is each first slice equal to first 2D array part if the t data?
    % squeeze(data(:,:,:,sl,1,1,1,1)) == ...
    %                          tableData_num(2:sz(3)+1,:) 
    % DEV: Is each scnd slice equal to scnd 2D array part if the t data?
    % squeeze(data(:,:,:,sl,2,1,1,1)) == ...
    %                           tableData_num(sz(3)+3:((sz(3)+1)*2),:) 

end % End table/slice loop.

% Save GUI data to figure.
guidata(tgui.fig, tgui);

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

% Get file destination from user.
[fn, fp, fi] = ...
    uiputfile({'*.txt', 'Text Files (*.txt)';...
               '*.mat', 'MATLAB File (*.mat)'},...
               'Save table data...');
if fi == 0, return; end
ext = fn(end-3:end);

% Save according to function-parent e.g. button.
switch fparent
    case 'Save all'
        for sli = 1:size(tgui.tab,2)
            data(:,:,sli) = ...
                cellfun(@str2double, tgui.tab{sli}.table.Data);
        end
        data2exp = data; ind = [];
        
    case 'Save selected'
         % Tab index
        tab_title = tgui.tabgp.SelectedTab.Title;
        tab_title_space = strfind(tab_title,' ');
        tab_ind = str2double(tab_title(tab_title_space+1:end));

        data = cellfun(@str2double, tgui.tab{tab_ind}.table.Data);
        
        % Get data
        ind = tgui.selectedCells;

        data = cellfun(@str2double, tgui.tab{tab_ind}.table.Data);
        data2exp = NaN(size(ind,1),1);
        for kk = 1:size(ind,1)
            data2exp(kk,1) = data(ind(kk,1),ind(kk,2));
        end
end

% Save to file.
switch lower(ext)
    case '.mat'
        selected = data2exp; index = ind;
        save([fp fn],'selected','data','index'); 
    case '.txt'
        csi_writeText(data2exp,[fp fn]);
end



% COORDINATES % -------------------------------------------------------- %
% MRI and CSI

% --- Executes on button press in button_CSI_setCoordinates.
function button_CSI_setCoordinates_Callback(~, ~, gui)
% Calculate coordinate of the CSI-volume
% Requires: resolution, #samples, offcenter, dimensions, fov.

% Return if no csi data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get possible ori struct
if isfield(csi,'ori'), ori = csi.ori; else, ori = struct; end


% Spatial Index % ----------- %

space_str = {'kx', 'ky', 'kz','x','y','z'};
space_dim = csi_findDimLabel(csi.data.labels,space_str);  
space_dim(isnan(space_dim)) = []; 
% Ask user if no kx/ky/kz are found
if isempty(space_dim)
    uans = getUserInput({'Spatial dimensions MRS data? (x/y/z)'},...
                        {[2 3 4]});
    if isempty(uans), CSI_Log({'Skipped FID & Echo split.'},{''}); return; 
    end
    % Convert userans to double
    space_dim = str2double(strsplit(uans{1},' '));
end
    

% Read Header Data % ----------- %
% For SPAR and SDAT data

if ~isfield(ori,'res') % Only if not already read the header file.
    if strcmp(csi.ext,'spar') || strcmp(csi.ext,'sdat')
        % Get from header file: Offcenter, slicethickness, slicegfov, fov 
        
        % Offcenter
        ori.offcenter(1) = csi.list.si_ap_off_center;
        ori.offcenter(2) = csi.list.si_lr_off_center;
        ori.offcenter(3) = csi.list.si_cc_off_center;
        % FOV
        ori.fov(3) = csi.list.slice_thickness;
        ori.fov(2) = csi.list.phase_encoding_fov;
        % RES
        ori.res(2) = ori.fov(2)./csi.data.dim(2);
        ori.res(3) = csi.list.slice_distance;
        ori.res(1) = ori.res(2);
    
    end
end

% UserInput % ----------- %

                    % RESOLUTION  +   OFFCENTER %

% Create default answers from previous input or header data
if isfield(ori,'res'), dans{1} = ori.res; 
else, dans{1} = '20 20 20';
end
if isfield(ori,'offcenter'), dans{2} = ori.offcenter;
else, dans{2} = '0 0 0';
end

% Offset and coordinates
qry = {'Resolution (AP LR FH) [mm]:', 'Offcenter (AP LR FH) [mm]: '};
uans = getUserInput(qry, dans);
if isempty(uans), CSI_Log({'Skipped calculating coordinates.'},{''}); return; 
end


                        % CORRECTION FACTORS %

% Default questions
qry2 = {'Apply odd/even size related half-voxel shift correction:'...
        'Apply FFT related half-voxel shift correction:'};

% Default answers
if isfield(ori,'vox_cor')
    if strcmp(ori.vox_cor, 'Yes'), dans2{1} = {'Yes','No'};
    else,                          dans2{1} = {'No', 'Yes'};       
    end
else, dans2{1} = {'Yes','No'}; 
end

if isfield(ori,'FFT_cor') 
    if strcmp(ori.FFT_cor, 'Yes'), dans2{2} = {'Yes','No'};
    else,                          dans2{2} = {'No', 'Yes'};       
    end
else, dans2{2} = {'No','Yes'}; 
end

% Ask user for corrections
uans2 = getUserInput_Popup(qry2, dans2);
if isempty(uans2)
    CSI_Log({'Skipped calculating CSI coorindates.'},{''}); return; 
end

switch uans2{1}, case 'Yes',vox_cor = 1; case 'No', vox_cor = 0; end
switch uans2{2}, case 'Yes',fft_cor = 1; case 'No', fft_cor = 0; end

% Calculate Coordinates % ----------- %

% Prep misc. parameters
ori.res       = str2double(strsplit(uans{1},' '));
ori.offcenter = str2double(strsplit(uans{2},' '));
ori.dim       = csi.data.dim(space_dim);

% Correction factors saved.
ori.vox_cor = uans2{1}; ori.fft_cor = uans2{2};       

% Coordinates of CSI data.
% Adds fields: coordinate vector (vector) & coordinate limits (limit) &
% volume limit (limit_vol).
csi.ori = csi_coordinates(ori,'center', vox_cor, fft_cor); 


% Calculate volume gird % ----------- %

% Meshgrid
[csi.ori.mesh.x, csi.ori.mesh.y, csi.ori.mesh.z] = ...
    meshgrid(csi.ori.vector{1} , csi.ori.vector{2}, csi.ori.vector{3});


% Clean up % ------------------------ %

% Update LOG
CSI_Log({   'CSI-parameters ---------------------',...
            'Direction:', 'Dimensions:', 'Resolution:',...
            'Offcenter:', 'FOV:', ...
            'Voxel shift due odd/even: ','Voxel shift due FFT-method:',...
            'Voxel limit (Min):', 'Voxel limit (Max):', ...
            'Volume limit (Min):', 'Volume limit (Max)', '',''},...
           { '',                      '[AP/LR/FH]',...
             csi.ori.dim,             csi.ori.res,...
             csi.ori.offcenter(1,:),  csi.ori.fov,...
             csi.ori.vox_cor,         csi.ori.fft_cor,...
             csi.ori.lim(:,1)',       csi.ori.lim(:,2)',...
             csi.ori.lim_vol(:,1)',   csi.ori.lim_vol(:,2)',...
             '----------------------------------------',''});

% Save to appdata
setappdata(gui.CSIgui_main,'csi',csi);   

% --- Executes on button press in button_MRI_setCoordinates.
function button_MRI_setCoordinates_Callback(~, ~, gui)
%
% Calculate MRI coordinates for PAR and DCM files.
% Using the appropriate functions.

% Get MRI data.
if ~isappdata(gui.CSIgui_main, 'mri'), return; end
mri = getappdata(gui.CSIgui_main, 'mri');

% Calculate for extensions DCM or REC
if     strcmp(mri.ext, 'dcm')
    MRI_coordinates_DCM(gui, mri); 
elseif strcmp(mri.ext, 'rec') || strcmp(mri.ext, '.par')
    MRI_coordinates_PAR(gui, mri);
end

% Get updated appdata.
mri = getappdata(gui.CSIgui_main,'mri');

% Show data to user.
if ~isfield(mri,'ori'), return; end
mid_slice = ceil(size(mri.ori.offcenter,1)/2);

% Show details to user: Res, dims, FOV, Limits.
CSI_Log({'MRI-parameters ---------------------', ...
    'Dimensions', 'Resolution: ', ...
    'Full FOV: ',...
    'Voxel limit (min): ','Voxel Limit (max): ',...
    'Offcentre (Middle): ', '', ''},...
    {'', mri.ori.dim, mri.ori.res, ...
    sprintf(' %3.2f',mri.ori.fov),...
    sprintf(' %3.2f',mri.ori.lim(:,1)),...
    sprintf(' %3.2f',mri.ori.lim(:,2)),...
    sprintf(' %3.2f',mri.ori.offcenter(mid_slice,:)),...
    '----------------------------------------',' '});

% --- Executes on button press in button_MRI_plotImage_Grid.
function button_MRI_plotImage_Grid_Callback(~, ~, gui)
% Plot converted images with the grid superposed on the images. No spectra
% are shown.

% Plot.
MRI_plotImage_tabbed(gui);

% --- Executes on button press in button_IMGinCSI.
function button_IMGinCSI_Callback(hobj, evt, gui)
% Show image in current displayed MRS array.
MRI_plotImage_current_CSI(hobj);

% --- Executes on button press in button_MRI_convert.
function button_MRI_convert_Callback(~, ~, gui)
%
% Convert MRI images to CSI data and create appdata CONV.
%
% Requires:
%       image data and the coordinates meshgrid, 
%       csi limit, resolution and dimensions.

MRI_to_CSIspace(gui);


% MRI Functions % ------------------------------------------------------ %

% --- Executes on button press in button_CSI_setCoordinates.
function MRI_coordinates_DCM(gui, mri)
% Calculate MRI coordinates for DICOM files.

imgori = struct;

% Get data dimensions, offcenter and resolution
try 
    
    % Get image types available.
    fn = fieldnames(mri.par); fnoi = fn{1}; 
    
    % Get stack details % ----------------------------------------------- %
    n_stacks = size(mri.data.(fnoi),3);
    if n_stacks > 1, stack_mode = 1; else, stack_mode = 0; end
    
    % Survey mode % ----------------------------------------------------- %
    % Enabled for multiple stack processing
    if stack_mode
        
        % Get each stack its information: 
        fn_stack = fieldnames(mri.examinfo.StacksInfo);
        
        % Loop each stack item_(1:Nstacks);
        stack_ori = cell(1,size(fn_stack,1));
        for fi = 1:size(fn_stack,1)
            % Orientation by viewaxis
            stack_ori{fi} = ...
            mri.examinfo.StacksInfo.(['Item_' int2str(fi)]).MRStackViewAxis;
        end
        
        % Convert stacks MR view axis to image orientation
        for kk = 1:size(stack_ori,2)
            switch stack_ori{kk}
                case 'FH', stack_ori{kk} = 'TRA';
                case 'RL', stack_ori{kk} = 'SAG';
                case 'AP', stack_ori{kk} = 'COR';
            end
        end
        
        % Select orientation of interest (by user)
        uans = getUserInput_Popup(...
            {'Select the stack with orientation of interest:'},{stack_ori});
        if isempty(uans)
            CSI_Log({'Skipped calculating image coordinates.'},{''}); return; 
        end
        
        % Get answer index in stack_ori.
        stack_ind = find(strcmp(stack_ori,uans)==1);
        stack_ori = stack_ori{stack_ind}; 
        
        % Get for the correct stack with prefered orientation, the 
        % orientation nr of slices and total nr slices.
        
        % Get orientation prefered stack images and par header.
        % Get stack details
        imginfo = mri.par.(fnoi)(stack_ind,:);
        imgdata = mri.data.(fnoi)(:,:,stack_ind,:);
        
        CSI_Log({'Detected survey image set.'},...
                       {['Using stack with orientation ' stack_ori]});
        
        % Replace stack.item_1 field with the stack of interest. Following
        % function will look at item_1 field only. Stack contains all stack
        % item_1:Nstacks for each slice.
        stack_item_oi = ...
            extractField(imginfo, ['Stack.Item_' num2str(stack_ind)]); 
        % Add it to stack.item_1 of each slice
        for sli = 1:size(stack_item_oi,1)
            imginfo{sli}.Stack.Item_1 = stack_item_oi{sli};
        end
        
        % Store the stack index! --> Used for converting images correctly.
        imgori.stack_of_interest = stack_ind;
        
    else
        % Get images and parameter header
        imginfo = mri.par.(fn{1});
        imgdata = mri.data.(fn{1});
    end
    
    % Get orientation options % ----------------------------------------- %
    opts = getOriOpts_dcm(imgdata, imginfo);
    
catch err
    err.message
    return;
end

% Process remaining ranges % -------------------------------------------- %

% RESOLUTION
imgori.offcenter = opts.ipp;
imgori.res       = opts.vox;
imgori.gap       = opts.gap;

% DIMENSIONS
imgori.dim       = opts.imSz;

if size(imgori.dim,2) == 2
    CSI_Log({'Cannot merge MR images to MRS using only 1 image.'},...
       {''}); return;
end
imgori.fov = imgori.dim.*imgori.res;

% RANGE
imgori.vec{1} = opts.range.col;
imgori.vec{2} = opts.range.row;
imgori.vec{3} = opts.range.slice;
imgori.lim    = cat(1,opts.lim.row,opts.lim.col, opts.lim.slice);

% MESH
imgori.mesh.x = opts.mesh.x; imgori.mesh.y = opts.mesh.y;
imgori.mesh.z = opts.mesh.z;

% Save data
mri.ori = imgori;
setappdata(gui.CSIgui_main, 'mri', mri);

% --- Executes on button press in button_IMGcoordinates.
function MRI_coordinates_PAR(gui, mri)
% Calculate MRI coordinates for Par/Rec-files.


imgori = struct; % Struct with orienation details.

% Check if orientation of all slices is equal
fns = fieldnames(mri.par);
stack_ori_all = extractField(mri.par.(fns{1}), 'acqori_tsc');
stack_ori_uni = unique(stack_ori_all);
if size(stack_ori_uni,1) > 1, stack_mode = 1; else, stack_mode = 0; end

if stack_mode
    % Stacks available, expect multiple orientations
    uans = getUserInput_Popup(...
        {'Select the stack with orientation of interest:'},{stack_ori_uni});
    if isempty(uans)
        CSI_Log({'Skipped calculating image coordinates.'},{''}); return; 
    end
    stack_ori = uans{1};
    
    % Get stack of interest indexes
    stack_ind = find(strcmp(stack_ori_all,stack_ori)==1);
    
    % Extract slices of stack of interest
    imginfo = mri.par.(fns{1})(stack_ind,:);
    imgdata = mri.data.(fns{1})(:,:,stack_ind,:);
    
    % Save stack of interest.
    imgori.stack_of_interest = stack_ind;
else
    imginfo = mri.par.(fns{1});
    imgdata = mri.data.(fns{1});
end

% Get data dimensions, offcenter and resolution
try 

    % OFFCENTER
    % Load from par as [ap fh lr]; 
    % Only view-axis index dimension differs. Flip columns to get LR/AP/FH.
    imgori.offcenter = cell2mat(extractField(imginfo, 'offcenter'));
    % Offcenter indexed as [ap lr fh];
    imgori.offcenter(:,[1 2 3]) = imgori.offcenter(:,[3 1 2]);               
    
    % RESOLUTION
    imgori.res = imginfo{1}.resolution;
    imgori.gap = imginfo{1}.gapsize;
    
    % DIMENSIONS
    imgori.dim = size(imgdata);
        
catch err
    err.message % Show the error message of development
    
    % Ask user input
    qry   = {'Offcentre [LR], [AP], [FH]','Resolution', 'Gap'};
    dfans = {[0 0 0], [20 20 20], 0};
    userinp = getUserInput(qry,dfans);
    if isempty(userinp)
        CSI_Log({'Skipped calculating image coordinates.'},{''}); return; 
    end
    
    imgori.offcenter = str2double(strsplit(userinp{1},' '));
    imgori.resolution= str2double(strsplit(userinp{2},' '));
    imgori.gap       = str2double(strsplit(userinp{3},' '));
end



imgori.fov = imgori.dim.*imgori.res;

% Create range vectors of the image-volume (mm);
imgori = PARcoordinates(imgori,'center');   

% Some notes about PARREC and coordinates.
%
% imgori        = PARcoordinates(imgori,'center',1);   %vx cor
% NOT FOR PARREC - TOPLEFT IS DEF WRONG.
% imgori        = PARcoordinates(imgori,'topleft');   
% imgori        = PARcoordinates(imgori,'topleft',1);  %vx cor 
% 1. Struct with offcenter/slice, dimensions and resolution
% 2. Topleft or Center - origin of offcenter values
% 3. Optional: voxel correction e.g half a voxel is corrected from the
%              offcenter points - if top left means left top of top left
%              voxel, you need vox correction.
% CREATE MESH GRID FOR IMAGE_DATA
% 1. XYZ-Coordinate for each voxel-center
% 2. Uses vectors of imgori created with PARcoordinates
% 3. FOV != abs(lim(min)) + lim(max); FOV == abs(lim(min max)) + voxel;


% Create grid of the image volume (mm);
[x,y,z] = meshgrid(imgori.vec{1} , imgori.vec{2}, imgori.vec{3});
imgori.mesh.x = x; imgori.mesh.y = y; imgori.mesh.z = z; 

% Save data
mri.ori = imgori;
setappdata(gui.CSIgui_main, 'mri', mri);

% DIsplay info
CSI_Log({'MRI-parameters order','Dimensions','Resolution',...
                'Gapsize','Offcenter','Full FOV',...
                'Center-Vox Limit (Min)','Center-Vox Limit (Max)'},...
               {'[AP/LR/FH]',imgori.dim,imgori.res,imgori.gap,...
                 imgori.offcenter(1,:), imgori.fov,...
                 imgori.lim(:,1)',imgori.lim(:,2)'});

             
% Display information to user.
CSI_Log({'Image coordinates calculated.'},...
               {'Conversion MR images to MRSI space enabled.'});

% --- Calculate images matching CSI slices
function [img, img_all, img_all_slice_range] = MRI_matchSlices(hObj)
% Using the conv struct, it calculates images corresponding to the MRSI
% data indexing. Number of slices in img is equal to number of slices in
% MRSI data.

% Gui-data struct.
gui = guidata(hObj);

% Get CONV struct
if isappdata(gui.CSIgui_main,'conv')
    conv = getappdata(gui.CSIgui_main,'conv');
else
    img = NaN; return;
end

if ~isappdata(gui.CSIgui_main,'csi') img = NaN; return; end

csi = getappdata(gui.CSIgui_main,'csi');
if ~isfield(csi,'ori'), img = NaN; return; end

        % ----- % Create MRS matching images % ----- %

nSlices = size(csi.data.raw,4);
img = NaN(size(conv.data,1),size(conv.data,2),nSlices);
for sli = 1:nSlices              % For every CSI slice

    % Slice coordinates for CSI and CONV
    cz = unique(csi.ori.mesh.z); mz = unique(conv.mesh.z);

    % Find CONV slices matching to CSI slices.
    img_range = CSI2MRI(cz(sli), mz, csi.ori.res(3), conv.res(3));
    img_all_slice_range(sli,:) = img_range;
    
    % Minimum and maximum index
    indMn = img_range(1); indMx = img_range(2);

    % Type of converted CSI-slice image
    pstr = get(gui.popup_plotIMG, 'String'); 
    pval = get(gui.popup_plotIMG, 'Value');

    % Calculate matching image
    switch pstr{pval}            
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

% --- Executed by MRIconvert2csi button callback
function MRI_to_CSIspace(gui)
% Convert MRI data to CSI space - requires mri and mrsi coordinates to be
% calculated.

            % ---------- % App-Data structs % ---------- %

% Get data: MRI
if ~isappdata(gui.CSIgui_main, 'mri')
    CSI_Log({'Requires MRI data!'},{';)'}); return; 
end
mri = getappdata(gui.CSIgui_main, 'mri');

% Get data: CSI
if ~isappdata(gui.CSIgui_main, 'csi')
    CSI_Log({'Requires MRSI data!'},{';)'}); return; 
end
csi = getappdata(gui.CSIgui_main, 'csi');


      % ---------- % Orientation Information MRSI/MRI % ---------- %

% Check for required orientation data field: CSI
if ~isfield(csi, 'ori') 
    button_CSI_setCoordinates_Callback([], [], gui);
	csi = getappdata(gui.CSIgui_main, 'csi');
    if ~isfield(csi, 'ori'), return; end
end

% Check for required orientation data field: MRI
if ~isfield(mri, 'ori')
    button_MRI_setCoordinates_Callback([], [], gui); 
    mri = getappdata(gui.CSIgui_main, 'mri');
    if ~isfield(mri, 'ori'), return; end
end    

% If existing conv data exists, get it.
if isappdata(gui.CSIgui_main, 'conv')
    conv = getappdata(gui.CSIgui_main, 'conv');
end


              % ---------- % Userinput % ---------- %

              
% Get available imagetypes
imtmp = fieldnames(mri.data);

% Get imagetype of interest; imtoi;
uans = getUserInput_Popup({'Image type to convert to MRSI space:'},...
                           {imtmp});
if isempty(uans), CSI_Log({'Skipped image convert.'},{''}); return; end
    
% Image type of interest
imtoi = uans{1}; 


             % ---------- % Convert IMG % ---------- %

% Converted resolution equals original MR images resolution. (3D)
% However, to fit correctly in csi grid, the resolution is changed (below).
conv.res = mri.ori.res; % Initial... May change!

% Calculate a resolution fitting the CSI space such that there is a integer
% amount of image pixels fitted in each CSI direction of space.
res_fit = csi.ori.res ./ conv.res;      % #MRpix / CSIpix
res_rem = res_fit - floor(res_fit);     % Pixel change
res_new = csi.ori.res ./ floor(res_fit);% New MRpix resolution 

% New resolution for each direction
conv.res = res_new;

% Volume limits of CSI but with half a voxel distance for voxel limits e.g.
% a total voxel (MRI res) vs the volume.
conv.fov       = csi.ori.fov;     % Does not change regards to CSI
conv.lim_vol   = csi.ori.lim_vol; % Volume of MRSI grid
conv.lim(:,1)  = conv.lim_vol(:,1) + (0.5.*conv.res)'; % Voxel limits
conv.lim(:,2)  = conv.lim_vol(:,2) - (0.5.*conv.res)'; % Used for coords

% Range of volume/grid of MRSI for MRI
for kk = 1:size(conv.lim,1)    
    % Number of pixels in kk-direction
    N = (conv.fov(kk)./conv.res(kk));
    
    % Grid vector defining volume
    conv.vec{kk} = linspace(conv.lim(kk,1), conv.lim(kk,2), N);
end

% Image grid in MRSI space 
% This grid lays in the CSI-space e.g. within limits of CSI FOV but is
% sampled in x, y and z as close to the resolution of the image as possible
[x,y,z] = meshgrid(conv.vec{1} ,conv.vec{2}, conv.vec{3});
conv.mesh.x = x; conv.mesh.y = y; conv.mesh.z = z; 

% Interp values @ CSI space % --------------------------------- %
% Interp3 is used to convert images to csi-grid in conv-struct.

% Check for stack availability - these images need seperate inclusion!
if isfield(mri.ori,'stack_of_interest')
    stack_ind = mri.ori.stack_of_interest;
else
    stack_ind = 1; 
end

if strcmp(mri.ext ,'.par')
    image_convert_data = mri.data.(imtoi)(:,:,:);
else
    image_convert_data = mri.data.(imtoi)(:,:,stack_ind,:);
end

% Interp
conv.data = interp3(mri.ori.mesh.x, mri.ori.mesh.y, mri.ori.mesh.z,...      % Original MRI coordinates
                    squeeze(image_convert_data),...                         % Original MRI values
                    conv.mesh.x, conv.mesh.y, conv.mesh.z,'Spline',0);      % Requested coordinates
conv.dim  = size(conv.data);

% Export any previous contrast info for image display.
if isfield(mri,'contrast'), conv.contrast = mri.contrast; end


% Post processing MRSI % ------------------------------------ %
% Apply required rotations or flips.

% % Apply additional rotation if it is list/data file
% if strcmp(csi.ext,'.list') || strcmp(csi.ext,'.data')
%     % Patient position: HFS and FFS
%     ppos = mri.par.(imtoi){1}.PrivatePerFrameSq.Item_1.PatientPosition;
%     
%     % Find spatial index: x,y, and z.
%     spat_ind = csi_findDimLabel(csi.data.labels,...
%                                {'kx','ky','kz','x','y','z'});
%     spat_ind(isnan(spat_ind)) = [];
%     if isempty(spat_ind)
%         % Ask user for k-space/spatial index.
%         uans = getUserInput(...
%                {'Spatial dimensions in MRS data? (x/y/z)'},{[2 3 4]});
%         if isempty(uans), return; end
%         spat_ind = str2double(strsplit(uans{1},' '));
%     end
%     
%     % Get, if available, unflipped original data.
%     if isfield(conv,'Original_MRSI')
%         mrsi_to_flip = conv.Original_MRSI;
%     else
%         mrsi_to_flip = csi.data.raw; conv.Original_MRSI = mrsi_to_flip;
%     end
% 
%     % Rotate according to HFS/FFS?
%     switch ppos
%         case 'HFS'
%         csi.data.raw = flip(flip(mrsi_to_flip,spat_ind(1)),spat_ind(3));
% 
%         case 'FFS'
%         csi.data.raw = flip(flip(mrsi_to_flip,spat_ind(1)),spat_ind(3));
%     end
% end

% Save data: CONV CREATED - IMAGE data in CSI-space.
setappdata(gui.CSIgui_main,'conv', conv);
setappdata(gui.CSIgui_main,'csi', csi);

% DIsplay info
CSI_Log({'Converted Images -------------------',...
         'Direction:','Dimensions:','Resolution',...
         'Voxel limit (Min)', 'Voxel limit (Max)',...
         'Volume limit (Min)', 'Volume limit (Max)',...
         '',''},...
        {'', '[AP/LR/FH]',conv.dim,conv.res,...
         conv.lim(:,1)',conv.lim(:,2)',...
         conv.lim_vol(:,1)',conv.lim_vol(:,2)',...
         '----------------------------------------',''});

% --- Save MRI and/or converted images
function MRI_saveIMG(hObj,~)
% Save the converted images to mat-file

                    % ------ % APP data % -------- %

% Gui-data struct.
gui = guidata(hObj);

% Get MRI struct
if ~isappdata(gui.CSIgui_main,'mri'), return; end
mri = getappdata(gui.CSIgui_main,'mri');

% Get CONV struct
if isappdata(gui.CSIgui_main,'conv')
    conv = getappdata(gui.CSIgui_main,'conv');
else
    conv = [];
end

                   % ------ % UserInput % -------- %

uans = getUserInput_Popup(...
    {'Which images to save: '},{{'All','Converted','MRS Matching'}});
if isempty(uans), CSI_Log({'Skipped exporting images.'},{''}); return; end


% UI for file destination from user
if isfield(mri,'filepath'), fp = mri.filepath; end
[fn, fp] = uiputfile({'*.mat','MATLAB file'},'Save images to file...',fp); 
if fn == 0, return; end
[~,fn,ext] = fileparts(fn); % Analyze file extension.

switch uans{1}
    case 'All'                          % Save MRI and converted %
        save([fp fn ext], 'mri','conv');
        
    case 'Converted'                    % Save converted %
        save([fp fn ext], 'conv');
        
    case 'MRS Matching'                 % Save MRS slice specific %


              % ----- % Create MRS matching images % ----- %
                
        % Create all the images in CSI space
        img = MRI_matchSlices(hObj);
        
        % Save the file.
        save([fp fn ext], 'img');
end

% --- Plot images in CSI slice
function MRI_plotImage_current_CSI(hObj)

gui = guidata(hObj);
% Get app-data
csi = getappdata(gui.CSIgui_main,'csi'); 
if isempty(csi), CSI_Log({'No CSI data in memory.'},{''}); return; end


[~, img_all, soi_range] = MRI_matchSlices(hObj);

panelobj = findobj('Tag','CSIpanel_2D_DataToDisplay'); 
if ~isempty(panelobj)
    % If it exists get guidata panel 2D.
    pan2D_gui = guidata(panelobj);
    % Get values from panel
    plotindex = cell(1,size(pan2D_gui.sliders,2));
    for kk = 1:size(pan2D_gui.sliders,2)
        plotindex{kk} = get(pan2D_gui.sliders{kk}, 'Value');
    end
else, plotindex = {1}; % Slice and other indices set to 1.
end

soi = soi_range(plotindex{1},:);
img_in_slice = img_all(:,:,soi(1):soi(2));

display3D(img_in_slice, 'limit',[0 max(img_in_slice(:))*0.75]);

% --- Plot images/tab and include voxel grid.
function tgui_data = MRI_plotImage_tabbed(gui, tag)
% Plot each image in csi-space and overlay a grid. The spectra are NOT
% shown. Each slice/array2D is plotted in a TAB. 
%
% Returns a figure: CSIgui - tabs with tag "tag" or CSIgui_tabs (default);
%
% OUTPUT: tgui_data. Includes figure object: tgui_data.fig;
%                    and tab-objects are stored in tgui.tabh;
% 
% Without output use:
% obj = figure('Tag', 'CSIgui_tabs'); if isempty(obj), return; end
% if size(obj,2)>1, obj = obj{1}; end

if nargin<2, tag = 'CSIgui - tabs'; end

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');


% Size of data
nSlices = size(csi.data.raw,4);

% GUI Colors
clr_bg = gui.colors.main;
clr_tx = gui.colors.text_main; clr_tb = gui.colors.text_title;

% PLOT_PAR: STRUCT WITH ALLL 2D PLOT OPTIONS
plot_par.colors = gui.colors;
                   % -------- % Figure: create window % -------- %

% Create figure
fh = figure('Tag', tag ,'Name', 'CSIgui - Tabs',...
            'Color', clr_bg, 'Toolbar', 'None', 'MenuBar', 'None',...
            'NumberTitle', 'Off');                   

% Create tab group
tabg = uitabgroup(fh); tabh = cell(1,nSlices);

% 1. Default figure size and screen pixel size
def_sz = 720; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
% 2. Ratio to define figure height to def_size
fig_sz = [def_sz def_sz.*(scr_sz(2)/scr_sz(1))];
% 4. Position of figure.
fig_ps = [40 scr_sz(2)-(1.15*fig_sz(2))];
% 5. Apply
set(fh, 'Position', [fig_ps fig_sz])        
 


              % -------- % Data: Dimensions % -------- %
              
plot_par.dim      = size(csi.data.raw);   % Data dimensions
plot_par.dim(1)   = [];                   % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 



% 1D Correction
% Set in plot_par the 1D dimension to the second index e.g. column. This 
% adds a y-dimension/height.
if plot_par.data_dim == 1, plot_par.dim(2) = 1; plot_par.data_dim = 2; end

              % -------- % Figure: Axis Grid % -------- %

% Axes linewidth: One point == 1/72 inch.
% axlw_pt = 1;

% Resolution without any correction of linewidth
plot_par.res = 1./plot_par.dim(1:2); 

% Loop X Y and Z
% X, Y and Z == csi-space. Position in figure starts at
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


              % --------- % Prepare images % --------- %
              
img = MRI_matchSlices(gui.CSIgui_main);

              % --------- % Loop all slices % --------- %

tgui_data = struct;
for sli = 1:nSlices 
    
plot_par.plotindex = {sli};
    
% Create a tab for this slice in in tabgroup   
tgui_data.tabh{sli} = uitab(tabg, 'Title',int2str(sli),...
    'BackgroundColor', clr_bg, 'ForegroundColor',clr_tb);

                 % ------ % Plot Images % ------- %
                
if isappdata(gui.CSIgui_main, 'conv')
    
    conv_data = getappdata(gui.CSIgui_main,'conv'); % Conv data
    img2plot = img(:,:,sli);                        % Get image to plot.
    
    % Create axis for image
    hold on; 
    imax = axes('parent',tgui_data.tabh{sli},'Position',[0 0 1 1], 'Color', 'None');
    
    % Plot Images
    if (sum(img2plot(:)) == 0) % Image is only zeroes.
        colormap(gray(255)); set(imax,'Color', 'Black'); alpha(1); 
    else
        % Image plotting:
        % Imscale as it plots over the entire figure and does not
        % imply any border issues as with imshow-function.
        imagesc(img2plot, 'parent', imax); 

        % Image Contrast.
        if isfield(conv_data, 'contrast')
            caxis(imax, conv_data.contrast);
        else
            caxis(imax,[min(img2plot(:)) max(img2plot(:))*0.5]);
        end
        colormap(gray(255));
    end 
    
    plot_par.colors.main = [0 0 0];
    
end

            % ------ %  PLOT GRID: Overlay of the axis. % ------ % 
CSI_2D_grid(tgui_data.tabh{sli},...
    fh.Position(3:4), plot_par.dim, plot_par.range, gui.colors.grid);

end
tgui_data.plot_par = plot_par; % add plot-par to output.
tgui_data.fig = fh;

% Save tgui-data to gui
guidata(fh, tgui_data);


% PLOT 2D CSI % -------------------------------------------------------- %


% --- Executes on button press in button_plotCSI.
function button_plotCSI_Callback(hObject, eventdata, gui)
% Show CSI button; to display CSI and possible converted MRI data.
% Launches:
%       CSIpanel_2D_DataToDisplay for navigation through all dimensions
%       CSI_plot2D to open up the CSIgui 2D figure and plot MRS data

if isappdata(gui.CSIgui_main, 'csi')
    panel_2D_DataSliders(hObject, eventdata, gui);
    CSI_2D_initiate2D(hObject,eventdata,gui);
end

% After loading update info
CSI_Log;

% --- % Executed by button: intiate plot 2D data. [show CSI]
function CSI_2D_initiate2D(~,~,~)
% ## Creates the CSIgui_plot figure ##
%
% 1. Get CSI data and analyse dimensions
% 2. Plot data in figure CSIgui_plot2D // CSIgui 2D-plot.
%
% The CSI data is plotted relative to the figure, thus the figure describes
% the FOV of the CSI data. Each axis shown in this 2D plot represents a
% pixel or voxel from the CSI data.
% 
% plot_par = struct with all plot information e.g. axis layouts, data and
%            more.

% GUIDATA main CSIgui app.
csgui_main_obj = findobj('Tag', 'CSIgui_main', 'Type', 'Figure');
gui            = guidata(csgui_main_obj);

% REQUIRED APPDATA: data, samples, index to plot.
csi = getappdata(csgui_main_obj,'csi'); 
if isempty(csi), CSI_Log({'No CSI data loaded!'},{':O'}); return;end

% Figure: close CSIgui_plot2D % ------------- % 

openFig.obj = findobj('Tag', 'CSIgui_plot2D');
if ~isempty(openFig.obj) && ~isempty(guidata(openFig.obj))
    openFig.gui = guidata(openFig.obj);     % Get figure gui-data.
    openFig.grid_sz = size(openFig.gui.ax); % Get grid size of axes;   
end
plot_par.colors = gui.colors;
                       
% Dimensions % ------------------- %

% Axis size to fit figure
plot_par.dim      = size(csi.data.raw);   % Data dimensions
plot_par.dim(1)   = [];                   % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 

% Correction 1D % ------------------- %
% For 1D data, the correct axis locations arent calculated. To correct, a
% second dimension is inserted and the 1D data is relocated to the Y e.g.
% column dimension. 
% REMARK: 
%       The 2nd dimension shown, equals x thus the COLUMNS
%       The 3rd dimension shown, equals y thus the ROWS

% Set in plot_par the 1D dimension to the second index e.g. column.
if plot_par.data_dim == 1
    % Add a Y dimension.
    plot_par.dim(2) = 1; plot_par.data_dim = 2;  
end


% Figure: create|reuse  % ------------------- %

% Create figure (2D object)
if isfield(openFig,'grid_sz') && ...            
        (sum(openFig.grid_sz == plot_par.dim([2 1]))==2)
    % Reuse figure %
    
    % Copy structs;
    plot_par = openFig.gui; plot_par.fh = openFig.obj; 
else 
    % Create figure %
    
    if isfield(openFig, 'gui') && ~isempty(openFig.gui)
        old_pos = openFig.gui.fh.Position; % Save position of figure.
    else
        old_pos = 0;
    end
    delete(openFig.obj); clear figOpen;
    
    % Create figure and figure-options
    plot_par = CSI_2D_setFigure(plot_par, old_pos);

end % End of "if create fig yes/no"

% ---- % Save xaxis to plot parameters: every time as can be updated.
plot_par.xaxis = csi.xaxis; 

% ------- % Get data for plotting
plot_par = CSI_2D_getData(plot_par, gui,  csi.data.raw);

% ------- % Get and set all other plot options
plot_par = CSI_2D_getPlotSettings(plot_par, gui, csi.data.raw);

% ------- % Plot images
plot_par = CSI_2D_plotImages(plot_par, csgui_main_obj);

% ------- % Plot datat using options from plot_par
CSI_2D_plotVoxels(plot_par,gui);

% ------ % Set the figure ratio to voxel size
% CSI_2D_setFigure_ratio(gui);

% --- % Executed to create a plot2D figure.
function plot_par = CSI_2D_setFigure(plot_par, init_pos, fig_tag)
% Create the 2D plot figure using dimension and color settings in plot_par.
%
% Input plot_par fields: .dim, .colors;
% Added plot_par fields: .res, .grid, .range, data_dim;
%
% Additional input:
% init_pos - initial position of the figure.
% fig_tag  - tag for the figure.

% Process input-arguments
if     nargin == 1, init_pos = 0; fig_tag = 'CSIgui_plot2D'; 
elseif nargin == 2,               fig_tag = 'CSIgui_plot2D';
end

% Create a new figure
fh = figure('Tag',fig_tag, 'Name', ['CSIgui 2D-plot: ' fig_tag],...
        'Color',plot_par.colors.main, 'Toolbar', 'None', 'MenuBar', 'None',...
        'NumberTitle', 'Off');                   
set(fh, 'CloseRequestFcn', @CSI_close2D);
        
% Axis: Resolution % ------------------- %
% FOV: Relative to figure (Normalized) e.g. [1 1];

plot_par.res = 1./plot_par.dim(1:2); % Axis e.g. voxel resolution

% Axis: Grid % ------------------- %

% Loop X Y and Z
% X, Y and Z == csi-space. Position in figure starts at point zero (0). 
% Resolution equals the nr of CSI-voxels to plot in the row (y in 
% image)and column (x in image) dimension of the figure. 
for kk = 1:numel(plot_par.res)
    plot_par.range{kk} = ...
    0 : plot_par.res(kk) : ...
       (plot_par.res(kk) * plot_par.dim(kk)-plot_par.res(kk));
end
% Caluclate a grid for each axis in the figure representing the voxels in
% CSI data.
[x,y] = meshgrid(plot_par.range{1},plot_par.range{2});
plot_par.grid.x = (x); plot_par.grid.y = flipud(y);

% 1D CORRECTION: Transpose x/col and y/row
% Creats Nx1 vs 1xN lists of the axis grid coordinates.
if plot_par.data_dim == 1, plot_par.grid.y = y'; plot_par.grid.x = x'; end

% Figure: Position % ------------------- %
if ~isequal(sum(init_pos),0)
    % 1. Use old position vector of previously opened CSI_plot2D figure.
    set(fh, 'Position', init_pos);
else
    % 1. Default figure size and screen pixel size
    def_sz = 640; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
    % 2. Ratio to define figure height to def_size
    fig_sz = [def_sz def_sz.*(scr_sz(2)/scr_sz(1))];
    % 4. Position of figure.
    fig_ps = [40 scr_sz(2)-(1.15*fig_sz(2))];
    % 5. Apply
    set(fh, 'Position', [fig_ps fig_sz]);
end

% Save figure object
plot_par.fh = fh;

% SNAP 2 PLOT % --------------------------------------- % DEV Experimental
% Add snapping of 2D panel to main plot figure.
% panel_2D_followPlot2D_initiate(); panel_2D_followPlot2D();

% --- Executed by CSI_plot2D_initiate and other to get plot-data.
function plot_par = CSI_2D_getData(plot_par, gui, data_volume)
% Adds data with correct settings used for plotting to plot_par field/
%
% Data:    display slice-index (plotindex), slice data (data2D), data unit
%          (dat_unit), slice-struct.index.data

% Data: Slice-Index: % ------------------- %

% Index legenda:
% Dim 1 = time/freq, Dim2/Dim3 = 2D roaster, DimN = panel/"slices";

% Get Panel2D GUI-object.
panelobj = findobj('Tag','CSIpanel_2D_DataToDisplay'); 
if ~isempty(panelobj)
    % If it exists get guidata panel 2D.
    pan2D_gui = guidata(panelobj);
    % Get values from panel
    plot_par.plotindex = cell(1,size(pan2D_gui.sliders,2));
    for kk = 1:size(pan2D_gui.sliders,2)
        plot_par.plotindex{kk} = get(pan2D_gui.sliders{kk}, 'Value');
    end
else
    % Only possible if no extra dimensions are present
    % E.g. data is 1D! No panel displayed!
    plot_par.plotindex = {1}; % Slice and other indices set to 1.
end

% Data: Slice % ------------------- %
            
% Get data-unit to plot and extract data from array.
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};
plot_par.data_unit = data_unit;

% Data save to field - Used for plotting.
plot_par.data2D = data_volume(:,:,:,plot_par.plotindex{:});

% --- Ececuted to return colorscaling settings for plotting.
function plot_par = CSI_2D_getPlotSettings_ColorScaling(plot_par, data_volume)
% Input-fields: plot_par...
%                    data_unit, plotindex.
%               Data total volume (3D)
% Output-fields: plot_par.clrs, plot_par.clrs_data_range

% Scaling: Plot Color % ------------- %

csiobj = findobj('Tag','CSIgui_main'); gui = guidata(csiobj);

% Returns colors gradient and related data-values, a range set within the
% given limits. This is used to color the plot of each voxel in the
% displayed CSI slice relative to the limits of the slice.
% E.g. visualise data amplitude using colors allowing individual voxel
% y-axis scaling!
%

% Color scaling by SLICE, VOLUME or STATIC.
vol_data = CSI_getUnit(data_volume, plot_par.data_unit);
scaleby = CSI_2D_Scaling_Color_Get(gui);
switch scaleby
    case 'vol' % Scale by volume
        % Plot color scaled for all voxels in the *volume*
        max_per_voxel = max(vol_data,[],1);
        data_ylimits_color = [min(max_per_voxel(:)), max(max_per_voxel(:))]; 
    case 'sli' % Scale by slice
        % Calc limits in slice
        max_per_voxel = max(vol_data(:,:,:,plot_par.plotindex{:}),[],1);
        data_ylimits_color = [min(max_per_voxel(:)) max(max_per_voxel(:))];           
    case 'sta' % Static color
        % Color scaling limited to 1 color: static line color
        data_ylimits_color = NaN;
end

% Get plot colors range for different max-limits.
if ~isnan(data_ylimits_color)
    % Check limits agrees with rules: lim(1) < lim(2)
    if data_ylimits_color(2) <= data_ylimits_color(1),...
            data_ylimits_color(2) = data_ylimits_color(1)+1; 
    end
    [plot_par.clrs, plot_par.clrs_data_range] = ...
        CSI_2D_Scaling_calc_ColorOfPlots(data_ylimits_color);
else
    % Set static line color.
    tmp = vol_data(:,:,:,plot_par.plotindex{:}); 
    plot_par.clrs = gui.colors.lines1; 
    plot_par.clrs_data_range = max(tmp(:)); 
end

% --- % Executed by CSI_plot2D_initiate: get plot2D settings
function plot_par = CSI_2D_getPlotSettings(plot_par, gui, data_volume)
% Add to structure plot_par the following plot-settings and plot-data
% fields:
% 
% Scaling: plot color, axis-y and x limit (by volume/static/voxel);


% Scaling: Plot Color % ------------- %
% Create - plot_par.clrs, plot_par.clrs_data_range
plot_par = CSI_2D_getPlotSettings_ColorScaling(plot_par, data_volume);


% Scaling: Axis Y Limit % ------------------- %
            
% Get data-unit to plot and extract data from array.
% Vol data is only used if volume color scaling to ylimits is on.
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};
plot_par.data_unit = data_unit;

% Y-limit of axis scaling by VOXEL, SLICE OR VOLUME.
scaleby = CSI_2D_Scaling_Axis_Get(gui);

% If voxel, set scale to zero, slice to one and volume to two. For slice
% and volume, already calculate the limits. 

switch scaleby 
    case 'vox', plot_par.scale = 0; 
    case 'sli', plot_par.scale = 1; 
        tmp_data = CSI_getUnit(plot_par.data2D, plot_par.data_unit);
        plot_par.axScale_ylimit = [min(tmp_data(:)) max(tmp_data(:))]; 
    case 'vol', plot_par.scale = 2; 
        tmp_dataVol = CSI_getUnit(data_volume, plot_par.data_unit);
        plot_par.axScale_ylimit = [min(tmp_dataVol(:)) max(tmp_dataVol(:))];
end

% Axis: x-axis data % ------------------- %

% Find possible frequenty axis.
if isfield(plot_par.xaxis, 'ppm'), plot_par.xaxisdata = plot_par.xaxis.ppm; 
else,                              plot_par.xaxisdata = plot_par.xaxis.none; 
end

% Display xlimit
plot_par.xlimit = plot_par.xaxis.xlimit;

% --- % Executed by CSI_plot2D_initiate: plot images
function plot_par = CSI_2D_plotImages(plot_par, csiguiObj)
% GUI of 2D figure

            % --------------% IMAGES TO PLOT % --------------%
                  % %%% MRI MATCHED TO CSI PLOT %%% %

% Get MRI image data and set to CSI space to plot.
if isappdata(csiguiObj, 'conv')
    conv = getappdata(csiguiObj, 'conv');
    try 
    
        % Get image @ this slice: uses CSI2MRI
        plot_par.img = MRI_matchSlices(csiguiObj);  
        plot_par.img = plot_par.img(:,:,plot_par.plotindex{1});
       
        % Plot IMG at CSI-coordinates! %%%%
        if ~isfield(plot_par, 'imax')
            plot_par.fh.Units = 'Pixels';
            ax_sz = plot_par.fh.Position(3:4);
            plot_par.imax = ...
                axes('parent',plot_par.fh,...
                'Unit','pixels','InnerPosition',[1 1 ax_sz], 'Color', 'blue');
        end
        plot_par.imax.Units = 'normalized';
        
        % hold(plot_par.imax, 'on'); 
        if (sum(plot_par.img(:)) == 0) % Image is only zeroes.
            imagesc(plot_par.imax,plot_par.img); 
            colormap(plot_par.imax,gray(255)); 
            set(plot_par.imax,'Color', 'Black'); alpha(plot_par.imax,0); 
        else
            % Image plotting
            % Imscale as it plots over the entire figure and does not
            % imply any border issues as with imshow-function.
            % hold(plot_par.imax, 'on'); 
           imagesc(plot_par.imax,plot_par.img); 
           %colormap(gray(2^15)); set(imax, 'Color', 'Black'); alpha(1);
           
            % Image Contrast.
            if isfield(conv, 'contrast')
                if (conv.contrast(1) >= conv.contrast(2))
                    caxis(plot_par.imax, ...
                        [conv.contrast(2)-1 conv.contrast(2)]);
                else
                    caxis(plot_par.imax, conv.contrast);
                end
                caxis(plot_par.imax, conv.contrast);
            else
                if min(plot_par.img(:)) >= max(plot_par.img(:)) 
                else
                    caxis(plot_par.imax,...
                        [min(plot_par.img(:))-1 0.5*max(plot_par.img(:))]);
                end
            end
            colormap(plot_par.imax, gray(255));
        end
        
        
        
        uistack(plot_par.imax,'bottom')
        plot_par.colors.main = [0 0 0];
        
    catch err
        CSI_Log({[err.stack(1).name ':'], 'Line: '},...
                       {err.message, err.stack(1).line});
    end        
end % End "if converted image data"

% --- % Executed by CSI_plot2D_initiate: plot voxels
function CSI_2D_plotVoxels(plot_par, gui)
% Plot the axis and its data in plot_par.fh; Input struct created from
% functions: plot2D_getSettings and plot2D_plotImages.
%
% If the figure is still opened
% e.g. axis-cell "ax" is present in input-struct plot_par, dont replot the
% axis nor the grid.
%
% plot_par.:
% fh, dim, grid.y, grid.x, res, data2D, clrs, clrs_data_range, xaxisdata,
% colors-structuct, xlimit, scale;

% Axes linewidth: One point == 1/72 inch.
axlw_pt = 1;  

% Replot the entire figure or use open figure.
if isfield(plot_par, 'ax'), createFig = 0; else, createFig = 1; end

% Plot csi voxel per axis in plot_par.grid
for ci = 1:plot_par.dim(1)                  % Col loop.
    for ri = 1:plot_par.dim(2)              % Row loop.
        
        % AXIS DETAILS % -------------- %
        
        % X and Y position in the figure of the axis to plot.
        x   = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
        % Position of axis to plot
        pos = [x y plot_par.res(1) plot_par.res(2)];
        % Create axis with pos(3,4) size at pos(1,2) position
        if ~ishandle(plot_par.fh), return; end
        if createFig
            plot_par.ax{ri,ci} = axes('parent',plot_par.fh,'position',pos);
        end

        % VOXEL DATA % ----------------- %
        plot_data_vox = plot_par.data2D(:,ci,ri);
        plot_data_vox = CSI_getUnit(plot_data_vox, plot_par.data_unit);
        
        % COLOR SCALING % -------------- %
        % This is NOT the Ylimit eventually used for axis-display.
        % See below plot();
        ylimit = [min(plot_data_vox(:)) max(plot_data_vox(:))];
        % Relative to maximum Y-data.
        [~, clr_ind] = min(abs(plot_par.clrs_data_range-ylimit(2)));
        plot_color = plot_par.clrs(clr_ind,:); % See before ri/ci for-loops.
  
       
        %%%% PLOT CSI data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(plot_par.ax{ri,ci}, plot_par.xaxisdata, plot_data_vox,...
                                'color', plot_color, 'LineWidth', 1.25);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % YAXIS SCALING % -------------- %
        
        % Scale Y axis by limits of voxel, slice or volume.
        if     plot_par.scale == 0                                          % Voxel - No action required.
        elseif plot_par.scale == 1 || plot_par.scale == 2                   % Slice or Volume
            ylimit = plot_par.axScale_ylimit;                               % Use pre-calculated limit.
        end
            
        % Largest absolute limit: min/max.
        if     abs(ylimit(2)) > abs(ylimit(1)), ylimfac = abs(ylimit(2));
        elseif abs(ylimit(1)) > abs(ylimit(2)), ylimfac = abs(ylimit(1));
        else
            ylimfac = 1;
        end
       
        % Create visualised Ylimit. Centered around zero!
        ylimit = [-1.05*ylimfac  1.05.*ylimfac];
        % Safety for the Ylimits;
        if ylimit(2) <= ylimit(1), ylimit(2) = ylimit(1)+1; end
        % Set Ylimit
        ylim(plot_par.ax{ri,ci}, ylimit);

        % X LIMIT % -------------------- %
        % Use limits of ppm axis(1) and (end) plus reverse the x-axis
        % direction to correctly display MRSI data.
        xlim(plot_par.ax{ri,ci}, plot_par.xlimit(1:2)); 
        set(plot_par.ax{ri,ci}, 'xdir', 'reverse');
        
        % AXIS COSMETICS % -------------- %
        set(plot_par.ax{ri,ci},...
               'Color', 'None',...
               'XColor', plot_par.colors.main,'YColor', plot_par.colors.main,...
               'LineWidth', axlw_pt, 'Xtick',[], 'Ytick',[],...
               'TickLength',[0 0.00001],'Box', 'off');             
        
        % Set action when clicked on this axis (ci, ri);
        set(get(plot_par.ax{ri,ci},'Children'),'HitTest', 'Off');
        set(plot_par.ax{ri,ci},'ButtonDownFcn',@CSI_2D_voxel_select);
        
        % Update guidata of axis, every axis gets the plot_parameters send
        % with it, thus its own handle too! 
        
        plot_par.ax{ri,ci}.UserData = [ri, ci];
        % guidata(plot_par.ax{ri,ci}, plot_par);
    end
end

guidata(plot_par.fh, plot_par);

% PLOT GRID: Overlay of the axis. % ----------------------------------- %
% Why: Due to a axis linewidth error in Matlab, the ancor point of position
% of the axis is not the bottom left point but the middle of the axis
% linewidth. This create descripancies in the plot layout due to
% pixel/normalized/inches conversion on the screen - especially when saving
% figures. Using seperate grids, this error is not visible.
if createFig
    CSI_2D_grid(plot_par.fh,  plot_par.fh.Position(3:4), ...
                plot_par.dim, plot_par.range, gui.colors.grid);
end

% Bring figure to front.
figure(plot_par.fh);

% --- Executes on button press in button_CSI_setFigure_ratio.
function button_CSI_setFigure_ratio_Callback(~, ~, gui)
CSI_2D_setFigure_ratio(gui);
 
function CSI_2D_setFigure_ratio(gui)
% Calculate ratio of screen and voxels; set figure ratio to voxel ratio.

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if ~isfield(csi,'ori'), return; end
obj2D = findobj('Tag', 'CSIgui_plot2D');
if isempty(obj2D), return; end

% Calculate ratio difference
res = csi.ori.res; scrsz = get(0,'screensize'); 
ratio_vox = res(1)./res(2); ratio_pc = scrsz(3)./scrsz(4); % X/W div Y/H
pos = obj2D.Position;

% Figure to monitor ratio
pos(3) = pos(4).*ratio_pc; obj2D.Position = pos; % Figure to screen ratio.

% Figure to voxel ratio
pos(3) = pos(4).*ratio_vox; obj2D.Position = pos;


% PLOT 2D FUNCTIONS % -------------------------------------------------- %

% --- Calculate and add the voxel grid in 2D-plots
function CSI_2D_grid(target, target_sz, dim, range, grid_clr)

% Get figure size for normalization of grid thickness
w = target_sz(1); h = target_sz(2);
% Define using vertical line the thickness of horizontal ones 
wv = 3; wh = wv; wv = wv./w; wh = wh./h; 
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

% --- Close the 2D MRSI plot figure gui
function CSI_close2D(hObj, ~)
% Custom close request function of the 2D CSI plot figure.

% Close the data to display panel.
panelobj = findobj('Tag','CSIpanel_2D_DataToDisplay');
if ~isempty(panelobj), delete(panelobj); end

% Close the 2D CSI figure.
delete(hObj);

% --- Executes when user clicks a voxel in the CSIgui 2D plot.
function CSI_2D_voxel_select(hObject, evt, ~)
% Runs when clicked on a voxel in a 2D CSI plot.
%
% Get clicked index of the voxel, higlight the voxel and remove any
% highlighting.

% Get clicked mouse button.
mouse_button = evt.Button;
if mouse_button == 1 % Left
    private_instance = 0 ; private_tag = '';
elseif mouse_button == 2 % Middle
    private_instance = 1 ; private_tag = datestr(now,'HHMMSS');
elseif mouse_button == 3 % Right
    private_instance = 0 ; private_tag = datestr(now,'HHMMSS');
end

% CLICKED VOXEL % --------------------- %

% Get data from axis object.
% This gui-obj includes slice info and data plus 2D-plot grid information
data2D = guidata(hObject.Parent); 

% Find which axis is clicked
% Use the input object (axis) and compare with object handles in the 
% spatially dependent cell array to find data index.
data2D.click_index = hObject.UserData; 
% WAY FASTER - In axis-user data field added the coordinates

% tic
% eq_bool = ...
%     cellfun(@isequal,data2D.ax,repmat({hObject},size(data2D.ax)));
% [r,c] = find(eq_bool == 1);
% data2D.click_index = [r, c];
% toc

% HIGHLIGHT VOXEL % --------------------- %

% Position and sizes
ind = data2D.click_index;

% Find any previous highlights
ax_color = extractField(data2D.ax,'Color');
no_color = ~cellfun(@isequal, ax_color, repmat({'none'},size(ax_color)));
[r,c] = find(no_color == 1);
% Remove this highlight
if isempty(r) % Store default background data
    data2D.ax_def_bgcolor = data2D.ax{1,1}.Color;
else          % Replace highlight with default
    data2D.ax{r,c}.Color = data2D.ax_def_bgcolor;
end
% Highlight clicked voxel (add transparent background)
data2D.ax{ind(1),ind(2)}.Color = [data2D.colors.hilight1 0.4];


% Save any update to data
guidata(hObject,data2D);

% Get data for 1D plot
% "data1D" also contains the data of data2D;
data1D = CSI_2D_voxel_selected_getData(data2D);

% Send update to plotCSIvoxel
CSI_1D_initiateGUI(data1D, private_instance, private_tag);

% --- Executed for selecting multiple voxel, if enabled in 2D plot.
function CSI_2D_voxel_selectMultiple(hObj, ~, ~)
% Every selected voxel is highlighted or the highlighting is removed.
% Specific for creating selected voxels.
%
% Used by merge_Voxels. Not intended for main 2D Plot functions.

% Clicked axes object and its data
gui = guidata(hObj); 

            % ----------- % Clicked index % ------------- %

% Find current tab number
tab_nr = ...
    cellfun(@isequal,gui.tabh, repmat({hObj.Parent},size(gui.tabh)));
tab_nr = find(tab_nr == 1);

% Find clicked axis index
% Compare ax-obj with stored ax-objects

gui.click_index = [hObj.UserData tab_nr];
r = gui.click_index (1); c = gui.click_index(2);


            % ----------- % Highlight voxel % ----------- %

% Position and sizes
ind = gui.click_index;

% If highlighted, remove it.
if ~ischar(gui.ax{ind(1),ind(2),tab_nr}.Color)
    gui.ax{r,c,tab_nr}.Color = 'none';
else
    % Highlight clicked voxel (add transparent background)
    gui.ax{ind(1),ind(2),tab_nr}.Color = [gui.highlight 0.2];
end

% Save any update to data
guidata(hObj,gui);

% --- % Executed to get data of clicked voxel
function data2D = CSI_2D_voxel_selected_getData(data2D)
% Get data of the clicked voxel
%
% Requires structure with fields:
%   clicked_index, slice.data, data_unit, xaxis-structure
%
% Returns structure with:
%   voxel.index, voxel.original, voxel.unit, axis, axis.unit;
%
% Output "data2D" applicable for CSIgui_1D or processing.


% CSIgui-1D: 2D to 1D data transfer % ---------------------- %
% ---------------------------------------------------------- %

% 2D to 1D: voxel index in CSIgui 2D
ind = data2D.click_index; % This is row/col == y/x;

% 2D to 1D: save the full voxel index (row,col,slice,etc..)
data2D.voxel.index = [flip(ind,2) cell2mat(data2D.plotindex)];

% 2D to 1D: voxel data @ index from CSIgui 2D
ydata_1D = data2D.data2D(:,ind(2),ind(1));
data2D.voxel.original = ydata_1D; % Add spectrum to 1D-appdata struct;

% 2D to 1D: store correct data unit for display
data2D.voxel.unit = data2D.data_unit;

% Add the frequentie struct containing x-axis data for the 1D data.
data2D.axis = data2D.xaxis;

% 2D to 1D: x-axis to 1D-appdata struct; find unit to display
if isfield(data2D.xaxis,'ppm')    
    % Save the ppm axis unit to appdata
    data2D.axis.unit = 'ppm';
else
    
    % Check for new axis data in CSIgui main appdata
    CSImain_obj = findobj('Tag','CSIgui_main');
    if isempty(CSImain_obj), return; end
    
    % Get CSI data
    csi = getappdata(CSImain_obj, 'csi');
    % Find if ppm axis is available
    if isfield(csi.xaxis,'ppm')
        data2D.axis = csi.xaxis;
        data2D.axis.unit = 'ppm';  % Save the axis unit to appdata
    else
        data2D.axis.unit = 'none'; % Save the axis unit to appdata
    end
     
end

% --- % Executed by CSI_2D_setFigure.
function panel_2D_followPlot2D_initiate(hobj)
% Initiate the panel2D to follow CSIgui2D plot figure if this screen is 
% moved.
%
% Requires some javascript to work.

% Master object
if nargin ~= 1
    hobj = findobj('Tag', 'CSIgui_plot2D', 'Type', 'Figure');
    if isempty(hobj), return; end
end

try
    % Just chill for a moment
    pause(0.05); % Wait for the figure construction complete.

    % Get JavaFrame: warnings
    jFig = get(hobj, 'JavaFrame'); 
    % Get Windowclient
    jWindow = jFig.fHG2Client.getWindow; 
    % Prevent memory leak
    jbh = handle(jWindow,'CallbackProperties');
    % Set new callback
    set(jbh,'ComponentMovedCallback',@panel_2D_followPlot2D);
catch err
    CSI_Log({'Error in followPlot2D_initiate'},{err.Message});
    return;
end

% --- % Executed if user moves plot2D window
function panel_2D_followPlot2D(~,~)
% Moves the 2D panel with sliders to the 2D plot figure automatically.
% Executed when user moves the main 2D plot window.

try 
    % Master
    mgui = findobj('Tag', 'CSIgui_plot2D', 'Type', 'Figure');
    % Slave
    sgui = findobj('Tag', 'CSIpanel_2D_DataToDisplay', 'Type', 'Figure');

    % Master and slave position
    mpos = mgui.Position; spos = sgui.Position;

    % Calculate new position slave
    spos_new = spos; spos_new(1:2) = mpos(1:2); 
    spos_new(1) = spos_new(1) + mpos(3);
    % Set new position of slave
    sgui.Position = spos_new;
    
    figure(sgui);
catch err
    CSI_Log({'Error in followPlot2D'},{err.message});
    return;
end

% PLOT 2D SCALING % ---------------------------------------------------- %

% --- Plot a colorbar for 2D plot information
function CSI_2D_Scaling_plotColorbar(hObj,~)
% Open up a colorbar representing color limits of the 2D CSI plot.
%
% Requires full csi-data acces.

try 
% Get guidata
gui = guidata(hObj);

if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');


% colors = plot_par.clrs;
% values = plot_par.clrs_data_range;

% --- Data unit
data_unit = get(gui.popup_plotUnit,'String');
data_unit = data_unit{get(gui.popup_plotUnit,'Value')};
plot_par.data_unit = data_unit;

% --- Plot index

% Get Panel2D GUI-object.
panelobj = findobj('Tag','CSIpanel_2D_DataToDisplay'); 
if ~isempty(panelobj)
    % If it exists get guidata panel 2D.
    pan2D_gui = guidata(panelobj);
    % Get values from panel
    plot_par.plotindex = cell(1,size(pan2D_gui.sliders,2));
    for kk = 1:size(pan2D_gui.sliders,2)
        plot_par.plotindex{kk} = get(pan2D_gui.sliders{kk}, 'Value');
    end
else
    % Only possible if no extra dimensions are present
    % E.g. data is 1D! No panel displayed!
    plot_par.plotindex = {1}; % Slice and other indices set to 1.
end

% Get plot-color scaling range
plot_par = CSI_2D_getPlotSettings_ColorScaling(plot_par, csi.data.raw);



% --- Open colorbar with color-options from plot_par

% Set input
colors = plot_par.clrs; values = plot_par.clrs_data_range;

colorbarQ(colors, values)

catch err
    CSI_Log({[err.stack(1).name ':'], 'Line: '},...
                       {err.message, err.stack(1).line});
     
end

% --- Executed by plot2D and others to get current color scaling state
function scaleby = CSI_2D_Scaling_Color_Get(gui)
% Returns what to color scale by; slice(SLI), volume(VOL) or static(STA)

% Enabled menu option
sli_scaling = get(gui.menubar.MRSI.ColorScale.Slice,  'Checked');
vol_scaling = get(gui.menubar.MRSI.ColorScale.Volume, 'Checked');
sta_scaling = get(gui.menubar.MRSI.ColorScale.Static, 'Checked');

% Set option which is checked to output.
if     strcmp(sli_scaling,'on'), scaleby = 'sli'; 
elseif strcmp(vol_scaling,'on'), scaleby = 'vol'; 
elseif strcmp(sta_scaling,'on'), scaleby = 'sta'; 
end

% --- Executed by plot2D and others to get current axis scaling state
function scaleby = CSI_2D_Scaling_Axis_Get(gui)
% Returns what to axis scale by; voxel(VOX), slice(SLI), volume(VOL)

% Enabled menu option
vox_scaling = get(gui.menubar.MRSI.AxisScale.Voxel,  'Checked');
sli_scaling = get(gui.menubar.MRSI.AxisScale.Slice,  'Checked');
vol_scaling = get(gui.menubar.MRSI.AxisScale.Volume, 'Checked');

% Set option which is checked to output.
if     strcmp(sli_scaling,'on'), scaleby = 'sli'; 
elseif strcmp(vol_scaling,'on'), scaleby = 'vol'; 
elseif strcmp(vox_scaling,'on'), scaleby = 'vox'; 
end

% --- Executed by user in menubar: CSI > Scaling Color
function CSI_2D_Scaling_Color_Set(hObject, evt)
% Change the color scaling setting of the 2D plot to slice or volume 
% maxima or turn it to static colors e.g. only one.
%
% Called by the menubar, to ensure only one option is checked.
% See CSI_2Dplot for processing of set option.

% Get guidata
gui = guidata(hObject);

% Get matlab year for correct field assessment.
matyr = version('-Release');
matyr = matyr(1:end-1); matyr = str2double(matyr);
if matyr >= 2017
    % Get clicked
    evt_click = evt.Source.Text; 
else
    % Get clicked
    evt_click = evt.Source.Label; 
end


% 2. Check what to change -> Scale by slice or volume.
switch evt_click
    case 'Slice'
        set(gui.menubar.MRSI.ColorScale.Slice, 'Checked', 'on');
        set(gui.menubar.MRSI.ColorScale.Volume,'Checked', 'off');
        set(gui.menubar.MRSI.ColorScale.Static,'Checked', 'off');
    case 'Volume'
        set(gui.menubar.MRSI.ColorScale.Slice, 'Checked', 'off');
        set(gui.menubar.MRSI.ColorScale.Volume,'Checked', 'on');
        set(gui.menubar.MRSI.ColorScale.Static,'Checked', 'off')
    case 'Static'
        set(gui.menubar.MRSI.ColorScale.Slice, 'Checked', 'off');
        set(gui.menubar.MRSI.ColorScale.Volume,'Checked', 'off');
        set(gui.menubar.MRSI.ColorScale.Static,'Checked', 'on');
end


% 3. Update GUI data.
guidata(hObject,gui);
    
% --- Executed by user in menubar: CSI > Scaling Axis
function CSI_2D_Scaling_Axis_Set(hObject, evt)
% Change the y-axis limit scaling setting of the 2D plot to
% on voxel, slice or volume maxima.
% 
% Called by the menubar, to ensure only one option is checked.
% See CSI_2Dplot for processing of set option.

% Get guidata
gui = guidata(hObject);

% Get matlab year for correct field assessment.
matyr = version('-Release');
matyr = matyr(1:end-1); matyr = str2double(matyr);
if matyr >= 2017
    % Get clicked
    evt_click = evt.Source.Text; 
else
    % Get clicked
    evt_click = evt.Source.Label; 
end

switch evt_click
    case 'Voxel'
        set(gui.menubar.MRSI.AxisScale.Voxel, 'Checked', 'on');
        set(gui.menubar.MRSI.AxisScale.Slice, 'Checked', 'off');
        set(gui.menubar.MRSI.AxisScale.Volume,'Checked', 'off');
    case 'Slice'
        set(gui.menubar.MRSI.AxisScale.Voxel, 'Checked', 'off');
        set(gui.menubar.MRSI.AxisScale.Slice, 'Checked', 'on');
        set(gui.menubar.MRSI.AxisScale.Volume,'Checked', 'off');
    case 'Volume'
        set(gui.menubar.MRSI.AxisScale.Voxel, 'Checked', 'off');
        set(gui.menubar.MRSI.AxisScale.Slice, 'Checked', 'off');
        set(gui.menubar.MRSI.AxisScale.Volume,'Checked', 'on');
end
    
% 3. Update GUI data.
guidata(hObject,gui);

% --- Executed by CSI_plot2D; y-axis color scaling.
function [clrs, clrs_range] = CSI_2D_Scaling_calc_ColorOfPlots(data_limits)
% Given limits of the data in a CSI slice:
%       a vector of size (1 x #colors) of values within the CSI-slice limit
%       is returned.
%       an array of size (#colors x 3) with a color per row.
% 
% To get a plot color associated with the plot-data per voxel
% relative to the slice its data e.g. y-limits:
%
% Example: [m, clr_ind] = min(abs(clrs_data_range-upper_limit_vox(2)));
%            plot_color = clrs(clr_ind,:); % See before ri/ci for-loops.

% Create plot-color of the line depending on the y-axis maxima in the 
% plotted CSI-data-slice.
clrs = jet(96);
% Yellow and blue index in jet-array
yind1 = find(ismember(clrs, [1 1 0],'rows') == 1);
yind2 = find(ismember(clrs, [0 0 1],'rows') == 1);
% Blue to Yellow and take only N colors from this list.
clrs = clrs(yind2+5:yind1+5,:);  
clrs = clrs(1:3:size(clrs,1),:);

% Assign value, within data_limits, to each color.
% Lower limit : difference limits / #colors = step : upper limit - step.
% The minus step at upper limit creates a range with equal size to colors.
clrs_range = ...
    data_limits(1):...                             
    abs(diff(data_limits))/size(clrs,1):...       
    data_limits(2)- abs(diff(data_limits))/size(clrs,1) ; 
if isempty(clrs_range), clrs_range = data_limits(1); end

% --- Executed by parameters button callback and several other scripts
function CSI_2D_Scaling_calc_xaxis(hObject, ~, auto)
% Creates struct xaxis which contains all frequenty information and
% corresponding axis for display.
%
% If auto equals 1, no user answer is required; 
% works only for SPAR/SDAT files or if specific frequency info is 
% available. Without this data, only arbitrary units are returned.

% If no automatic option is given, turn if off.
if nargin <= 2, auto = 0; end

% Get time axis.
gui = guidata(hObject); 

% Require csi: data, samples, csisl
csi = getappdata(gui.CSIgui_main,'csi'); 
if isempty(csi), return; end

% Xaxis structure
xaxis = struct;

% Store previous static data
if isfield(csi, 'xaxis')
    foi = {'BW','nucleus','trans','gyro','tesla','xlimit', 'shift'};
    for kk = 1:size(foi,2)
        tmp = foi{kk};
        if isfield(csi.xaxis, tmp)
            xaxis.(tmp) = csi.xaxis.(tmp);
        end
    end
end

% READ HEADER % ------------------------------------------ %
% Works for SPAR file only.
if strcmp(csi.ext ,'spar') || strcmp(csi.ext ,'sdat') 
    
    % Get from header file: BW, nucleus, transmit freq
    xaxis.BW      = csi.list.sample_frequency;
    xaxis.nucleus = csi.list.nucleus;
    xaxis.trans   = csi.list.synthesizer_frequency;
    
    % Gyromagnetic constant in MHz to calc magnet strength
    switch xaxis.nucleus 
        case '1H',  xaxis.gyro = 42.57747892*10^6;
        case '2H',  xaxis.gyro = 6.536*10^6;
        case '31P', xaxis.gyro = 17.235*10^6;
        case '23Na',xaxis.gyro = 11.262*10^6;
        case '19F', xaxis.gryo = 40.052*10^6;
        otherwise,  xaxis.gyro = 42.57747892*10^6;
    end
    
    % Calc magnet strength
    xaxis.tesla = xaxis.trans./xaxis.gyro;
end

% USER INPUT % ------------------------------------------ %
% If auto == 1, this is skipped

if auto == 0
    % Request required values
    qry = { 'Nucleus (1H): ', 'Magnet strength(T): ', 'Bandwidth (Hz): ',...
            'PPM Shift (Num): '} ;

    % Create possible answers if exist: previous input or header file.
    if isfield(xaxis,'nucleus'),an{1} = xaxis.nucleus; 
    else, an{1}= '31P';  
    end    
    if isfield(xaxis,'tesla'),  an{2} = xaxis.tesla;   
    else, an{2}= '7T';   
    end
    if isfield(xaxis,'BW'),     an{3} = xaxis.BW;      
    else, an{3}= '8000'; 
    end
    if isfield(xaxis,'shift'),  an{4} = xaxis.shift;   
    else, an{4}= '0';    
    end

    % Display UI and get user input
    inp = getUserInput(qry, an);
    if isempty(inp)
        CSI_Log({'Skipped setting parameters.'},{''}); return; 
    end

    % Set user input.
    xaxis.nucleus = inp{1};
    xaxis.tesla   = str2double(inp{2}(1));
    xaxis.BW      = str2double(inp{3});
    xaxis.shift   = str2double(inp{4});

    switch inp{1} % Gyromagnetic constant in MHz
        case '1H',  xaxis.gyro = 42.57747892*10^6;
        case '2H',  xaxis.gyro = 6.536*10^6;
        case '31P', xaxis.gyro = 17.235*10^6;
        case '23Na',xaxis.gyro = 11.262*10^6;
        case '19F', xaxis.gryo = 40.052*10^6;
        otherwise,  xaxis.gyro = 42.57747892*10^6;
    end

    % Process input
    xaxis.trans = xaxis.gyro.*xaxis.tesla;

else
    if ~isfield(xaxis,'shift') % Only if no shift is present!
        xaxis.shift = 0;
    end
end

% CALC AXES % -------------------------------------------- %

% Calculate axis-parameters
ind = find(strcmp(csi.data.labels,'sec') == 1);
if ~isempty(ind), xaxis.N = csi.data.dim(ind);
else,             xaxis.N = csi.data.dim(1);
end

if isfield(xaxis, 'trans')
    % Calculate axis-parameters

    % One PPM
    xaxis.img  = xaxis.trans.*1e-6;
    % Dwelltime
    xaxis.dwellTime = 1./xaxis.BW;
    % Total acquisition time
    xaxis.acqTime = xaxis.dwellTime*xaxis.N;
    
    % Unitless N
    xaxis.none = (0 : (xaxis.N)-1);
    % N -> Time             
    xaxis.time = (0 : (xaxis.N)-1) .* xaxis.dwellTime; 
    % N -> HZ
    xaxis.freq = ( -xaxis.N/2 : (xaxis.N/2) - 1) ./ xaxis.acqTime; 
    % HZ -> PPM
    xaxis.ppm  =  (xaxis.freq / xaxis.img) - xaxis.shift; 
    
    % Display xlimit
    if ~isfield(xaxis,'xlimit')
        xaxis.xlimit = [xaxis.ppm(1) xaxis.ppm(end)];
    end
else
    xaxis.none = (0 : (xaxis.N)-1);
    if ~isfield(xaxis,'xlimit')
        xaxis.xlimit = [xaxis.none(1) xaxis.none(end)];
    end
end

% If PPM available, set new limit
if isfield(xaxis,'ppm')
    xaxis.xlimit = [xaxis.ppm(1) xaxis.ppm(end)];
end

% CLEAN UP % -------------------------------------------- %

% Add CSI-output.
csi.xaxis = xaxis; % Contains frequency details of MRS data

% Save in appdata
setappdata(gui.CSIgui_main,'csi',csi);

% --- Executed to manually set various 2D-plot scaling & display options
function CSI_2D_Scaling_Options(gui)
% Set display options: x-axis control, y-axis control, data color
% scaling 
%
% Uses CSI_Scaling2D_Axis; CSI_Scaling2D_Color


% Get CSI app data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

% USER INPUT % -------------------------------------- %

% Y-axis scaling options % -------- %

% Get current y-axis scaling options: full name
current_clr = CSI_2D_Scaling_Color_Get(gui);
switch current_clr
    case 'vol', current_clr = 'Volume';
    case 'sli', current_clr = 'Slice';
    case 'sta', current_clr = 'Static';
end
current_axs = CSI_2D_Scaling_Axis_Get(gui);
switch current_axs
    case 'vox', current_axs = 'Voxel';
    case 'sli', current_axs = 'Slice';
    case 'vol', current_axs = 'Volume';
end
% Sort with current setting at top.
inp_clr = {'Voxel', 'Volume', 'Static'};
ind = contains(inp_clr, current_clr);inp_clr = {inp_clr{ind} inp_clr{~ind}};
inp_axs = {'Voxel','Slice','Volume'};
ind = contains(inp_axs, current_axs);inp_axs = {inp_axs{ind} inp_axs{~ind}};

% Ask user.
yaxis_ans = getUserInput_Popup({'Color scaling:','Y-axis scaling:'},...
                          {inp_clr,inp_axs});
if isempty(yaxis_ans) 
    CSI_Log({'Skipped setting 2D plot options.'},{''}); return; 
end  
                      
% X-axis scaling options % -------- %

% Get current x-axis settings
xlimit = csi.xaxis.xlimit;                     

% Ask user
xaxis_ans = getUserInput({'X-axis display range: (ppm or unitless)'},...
                         {xlimit});                     
% Return if user pressed skip.
if isempty(xaxis_ans)
    CSI_Log({'Skipped setting 2D plot options.'},{''}); return; 
end  

% PROCESS INPUT % -------------------------------------- %

% Process y-axis scaling answer
evt.Source.Text = yaxis_ans{1}; evt.Source.Label = yaxis_ans{1};
CSI_2D_Scaling_Color_Set(gui.CSIgui_main, evt);
evt.Source.Text = yaxis_ans{2}; evt.Source.Label = yaxis_ans{2};
CSI_2D_Scaling_Axis_Set(gui.CSIgui_main, evt);

% Process x-axis range answer
csi.xaxis.xlimit = str2double(strsplit(xaxis_ans{1},' '));

% Figure ratio
button_CSI_setFigure_ratio_Callback([], [], gui);

% Save appdata.
setappdata(gui.CSIgui_main,'csi',csi); 

% --- Executes on button press in button_CSI_DisplayOptions.
function button_CSI_DisplayOptions_Callback(~, ~, gui)
% Set all scaling options for 2D plot using a GUI
CSI_2D_Scaling_Options(gui);


% PLOT 2D PANEL % ----------------------------------------------------- %


% --- Panel for controlling which 2D CSI slice is displayed.
function panel_2D_DataSliders(~, ~, gui)
% Opens up a figure with sliders to control slice display in the CSIgui 2D
% plot figure.
%
% Will open according to color theme of CSIgui.

% CHECK FIG INSTANCES % -------------------------------- %

if ~isappdata(gui.CSIgui_main,'csi'), return; end

% CSIpanel_2D_DataToDisplay
% If panel exist, do not proceed again.
panelobj = findobj('Tag','CSIpanel_2D_DataToDisplay');
if ~isempty(panelobj), close(panelobj); end

% PREP DATA % ------------------------------------------ %

% Get CSI data-structure
csi = getappdata(gui.CSIgui_main,'csi');
% Analyze data size
dim = csi.data.dim; % To exclude index dimensions equal to 1.

% Safety - No "Slice" or other dimensions available - Dont open & return;
if numel(dim) <= 3, return; end 
% Exclude time: X/row and Y/col.
dim = dim(4:end); 

% PREP FIGURE % -------------------------------------- %

% Figure: Coloring
% Get coloring from CSIgui
if isfield(gui, 'colors')
    clr_bg = gui.colors.main; clr_butt = gui.colors.text_main;
else
    clr_bg = [0 0 0]; clr_butt = [1 0 0];
end

% Figure: CSI_displayDataPanel
fh = figure('Tag', 'CSIpanel_2D_DataToDisplay','Name', 'CSIgui display panel',...
            'Color', clr_bg,'MenuBar','none','Toolbar', 'none', ...
            'resize', 'off','NumberTitle', 'off');                        

% Set size of controls
w = 320; h = 20; sz_slider = [w*0.85 h]; sz_txt = [w*0.15 h];
% Size of figure;
hfig = (numel(dim)*h)+(2*h); wfig = w+(0.10*w);
% Get position of figure CSIgui_2Dplot -> Add width of this figure to its
% position to plot it next to the figure.
set(fh, 'position', [770 683 wfig hfig]);

% PREP SLIDERS % -------------------------------------- %

for sli = 1:numel(dim)
    % Slider XY
    sldX = (wfig-w)/2; sldY = hfig-(h*sli+h);
    % Slider min, max, val and step
    sMax=dim(sli); sMin=1; sVal = 1;
    slStep = ([1 1]./(dim(sli)-1));
    if sum(~isfinite(slStep)) >=1 , slStep = [1 1]; end
    
    % Create slider.
    sgui.sliders{sli} = uicontrol('Style', 'Slider', 'Parent', fh, ...
        'Units', 'pixels','tag', sprintf('dimslider_%i',sli), ...
        'position', [sldX sldY  sz_slider(1) sz_slider(2)], ...
        'min', sMin, 'max', sMax,'value', sVal,...
        'SliderStep', slStep, ...
        'BackgroundColor', clr_bg, 'ForegroundColor', clr_butt,...
        'Callback', @panel_2D_sliders);
    
    % Text of slider nDim/MaxDim;
    txtX = ((wfig-w)/2)+sz_slider(1); txtY = hfig-(h*sli+h);
    sgui.texts{sli} = uicontrol('Style', 'Text', 'Parent', fh, ...
        'Units', 'pixels','tag', sprintf('dimtxt_%i',sli), ...
        'String', sprintf('%i/%i',sVal,dim(sli)),...
        'position', [txtX txtY sz_txt(1) sz_txt(2)], ...
        'ForegroundColor', clr_butt, 'BackgroundColor', clr_bg);
    
end
% Update figure with its sgui - subgui.
guidata(fh, sgui);




% --- Executes when user changes slider in CSIpanel_2D_DataToDisplay.
function panel_2D_sliders(hObject, eventdata)
% Get new value from slider and store rounded. Update text-element in gui
% using to nr found in tag of the slider e.g. index to handle of texts{ind}.

% Round value of slider and set it as value.
val = round(get(hObject, 'Value')); 
if val <= get(hObject, 'Max') && val >= get(hObject,'Min')
    set(hObject, 'Value', val);
else
    set(hObject, 'Value', get(hObject,'Min'));
    val = get(hObject,'Min');
end

% Update txt of displayDataPanel
gui = guidata(hObject); ind = get(hObject,'Tag'); ind = ind(end); %TxtNr
str_old = get(gui.texts{str2double(ind)}, 'String');
str_sls = strfind(str_old, '/'); str_old(1:str_sls-1) = []; 
set(gui.texts{str2double(ind)}, 'String', [num2str(val) str_old]);

% Update plot!
CSI_2D_initiate2D(hObject,eventdata,gui)


% DATA RETRIEVE FUNCTIONS % ------------------------------------------- %


% --- Get CSI domain; frequency or time
function domain = CSI_getDomain(gui)

% Get both domain states
tstate = gui.menubar.MRSI.domain.time.Checked;
fstate = gui.menubar.MRSI.domain.frequency.Checked;

% Set output argument domain
if strcmp(tstate,'on'),     domain = 'time';
elseif strcmp(fstate,'on'), domain = 'freq';
else,                       domain = 'none';
end

% --- Set CSI domain; frequency or time
function CSI_setDomain(hObj,evt)
% Set the data domain of MRSI data to time or frequency e.g. spectrum or
% FID.

% Get gui appdata.
gui = guidata(hObj);

% Get matlab year for correct field assessment.
matyr = version('-Release');
matyr = matyr(1:end-1); matyr = str2double(matyr);
if matyr >= 2017
    % Get clicked
    evt_click = evt.Source.Text; 
else
    % Get clicked
    evt_click = evt.Source.Label; 
end

% Set clicked property to on in menubar
switch evt_click
    case 'Frequency'
        gui.menubar.MRSI.domain.frequency.Checked = 'on';
        gui.menubar.MRSI.domain.time.Checked = 'off';
    case 'Time'
        gui.menubar.MRSI.domain.time.Checked = 'on';
        gui.menubar.MRSI.domain.frequency.Checked = 'off';
end

% --- Get data in unit; real, imaginary, absolute phase
function data = CSI_getUnit(data, unit_str)
% Returns the correct data format e.g. real, imaginary, magnitude or phase
% of input data.

switch unit_str
    case 'Real',      data = real(data);
    case 'Imaginary', data = imag(data);
    case 'Magnitude', data = abs(data);
    case 'Phase',     data = angle(data);
    otherwise,        data = abs(data);
end



% --- Get data at peak of interest
function [doi, doi_axis, range] = CSI_getDataAtPeak(spec, xaxis, range)
% Input 1: Spec is the full spectrum;
% Input 2: xaxis can be empty and a arbitrary range will be used otherwise 
%          a ppm-axis for the spec data is expected.
% Input 3: range can be empty and the user will be asked to enter a range
%          otherwise the index range is requested (low to high - not ppm)

% Get peak of interest from user if absent
if nargin <= 2
    % Get range from user.
    range = CSI_getPeakOfInterest(xaxis); 
    
    % If no range is found
    if isempty(range)
        CSI_Log({'No data in given range.'},{'Returning.'});
        doi = NaN; doi_axis = NaN; range = NaN; return;
    end
end

% Get data at range % ------------------------------- %

% Convert first get dimensions of data matrix to linear vector per cell.
sz = size(spec); 
if numel(sz) >=2, nDimC = num2cell(sz(2:end));else, nDimC = {1}; end
nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

% Cut data with given range
doi = spec(range(1):range(2),nDimC{:}); 
if isfield(xaxis,'ppm'), doi_axis = xaxis.ppm(:,range(1):range(2)); 
else,                    doi_axis = range(1):range(2);
end

% --- Get peak of interest from user
function range = CSI_getPeakOfInterest(xaxis, poi_tag)
% Returns a peak of interest from the user using the xaxis field of CSIgui.
% 
% Output: unitless index range, not the ppm range.
% This output can be used as an index in xaxis.unitless or ppm;

if nargin == 1, poi_tag = ''; end

if isfield(xaxis, 'ppm')
    ppm = 1;
    unit_str = '(ppm):'; unit_ans = [xaxis.ppm(1) xaxis.ppm(end)];
else
    ppm = 0;
    unit_str = '(Unitless):'; unit_ans = [1 xaxis.none(end)+1];
end
unit_str = ['Peak range of interest ' unit_str ' ' poi_tag];

% Get range from user
uans = getUserInput({unit_str},{unit_ans});
if isempty(uans), range = []; return; end  

% Set ranges in double
poi = str2double(strsplit(uans{1}, ' ')); poi = sort(poi);

% If PPM available convert to unitless.
if ppm
    % If xaxis is not zero - use the xaxis input
    CSI_Log({'Converting ppm-range for POI: '},{poi});
    [~, poi(1)] = min(abs(xaxis.ppm - poi(1)));
    [~, poi(2)] = min(abs(xaxis.ppm - poi(2)));
    CSI_Log({'Converted range-index: '},{poi});
end

% Sort the range from low to high
range = poi; 

% --- Convert FID/ECHO poi to FID/ECHO poi 
function [doi, doi_axis, doi_range] = CSI_getDataAtPeak_Stored(range, gui)
% If data is split up into FID and echoes:
% Returned CSI_getDataAtPeak range from the active type will be converted
% to the stored type: e.g. fid2echo or vice versa.
%
% 1. If no ppm-scale is present, the sample ratio difference is used.
% 2. No ppm-shift is taken into account if no ppm-scale is available.


% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Active and stored: fid or echo
active = csi.data.split.type;
switch active
    case 'FID', stored = 'echo';
    case 'Echo',stored = 'fid';
end

% Get the xaxis structure of the stored data set
data = csi.data.split.(stored); xaxis_stored = csi.data.split.xaxis.(stored);

% Use limits of xaxis 
poi = [range(1) range(size(range,2))]; 
if isfield(xaxis_stored,'ppm') && isfield(csi.xaxis,'ppm')
    % Find PPM range in axis of the stored one.
    
    [~, range_converted(1)] = min(abs(xaxis_stored.ppm - poi(1)));
    [~, range_converted(2)] = min(abs(xaxis_stored.ppm - poi(2)));
else
    % Use ratio of diff FID & echo.

    % Ratio 
    active_sz = size(csi.data.raw,1); 
    stored_sz = size(csi.data.split.(stored),1); 
    ratio = stored_sz/active_sz; 
    
    % Range in echo data
	range_converted = poi.*round(ratio);

end

% Get peak of interest data using converted range
[doi, doi_axis, doi_range] = ...
    CSI_getDataAtPeak(data, xaxis_stored, range_converted);



% CSIgui-1D: DISPLAY & PANEL % ----------------------------------------- %


% --- Executed after CSIgui_2D voxel click to initate 1D plot GUI
function CSI_1D_initiateGUI(data1D, varargin)
% Initiates the 1D spectrum plot figure to plot data from the 
% 2D CSIgui plot figure.
%
% Requires data1D; @CSI_2D_voxel_selected_getData(data2D)
%           the CSIgui 1D application data containing all important
%           data for plotting and processing.
%
%
%  data1D fields
%  voxel:       voxel information
%  data:        spectrum/fid actual plot data
%  axis:        x-axis data for plotting
%
% varargin{1}:  Start a private (1) or coupled(0, default) instance of
%               CSIgui 1D.
% varargin{2}:  Tag to show to user in CSIgui-1D

if     size(varargin,2) == 1
    instance.private = varargin{1};  tag = ''; 
elseif size(varargin,2) == 2
    instance.private = varargin{1};  tag = varargin{2};
else
    instance.private = 0;            tag = '';
end

if instance.private, instance.tag = datestr(now,'HHMMSSFFF');
else,                instance.tag = '';
end

% CSIgui-1D: Create figure % ------------------------------- %
% ---------------------------------------------------------- %

% Close existing window: only if not a private session
if ~instance.private
    CSI_1Dobj = findobj('Type','Figure','Tag','CSIgui_plot1D',...
                        'Name', 'CSIgui - plot 1D');  
    for kk = 1:size(CSI_1Dobj,1)
        tmp = getappdata(CSI_1Dobj(kk),'instance');
        if strcmp(tmp.tag,instance.tag)
            delete(CSI_1Dobj(kk)); break;
        end    
    end
end

% Create plot figure: CSIgui_plot1D.
CSI_1Dobj = figure('Tag', 'CSIgui_plot1D','Name', 'CSIgui - plot 1D',...
               'numbertitle', 'off');
data1D.tag = tag; 

% Add custom close function
set(CSI_1Dobj, 'CloseRequestFcn', @CSI_1D_closeGUI); 

% Set figure position; keep original size.
scrsz = get(0, 'screensize'); figsz = get(CSI_1Dobj,'Position');
fig_loc = [ceil(scrsz(3)/2 - (figsz(3)/2)) ceil(scrsz(4)/2 - (figsz(4)/2))];
set(CSI_1Dobj,'Position',[fig_loc(1) fig_loc(2) figsz(3) figsz(4)]);

% Create axis
data1D.axes = axes('parent',CSI_1Dobj);

% Create instance text for instance-tag
uicontrol(CSI_1Dobj, 'Style','Text', 'String', ['Instance: ' instance.tag],...
    'Unit','Pixels','Position',[1 1 125 15],...
    'Fontsize', 8, 'FontWeight', 'Bold',...
    'HorizontalAlignment','Left');


% Create text for text-tag
uicontrol(CSI_1Dobj, 'Style','Text', 'String', ['Tag: ' data1D.tag],...
    'Unit','Pixels','Position',[126 1 125 15],...
    'Fontsize', 8, 'FontWeight', 'Bold',...
    'HorizontalAlignment','Left');

% Create panel button
uicontrol(CSI_1Dobj, 'Style','pushbutton','String','Panel',...
    'Unit','Pixels','Position',[figsz(3)-50 figsz(4)-20 50 20],...
    'Fontsize', 8,'FontWeight','Bold','HorizontalAlignment', 'Center',...
    'Callback',@CSI_1D_movePanelToGUI,...
    'TooltipString',...
    'Show the 1D-options panel linked to this CSIgui 1D instance.',...
    'Unit','normalized'); % Revert to normalized units!

% CSIgui-1D: Clean up % ----------------------------------- %
% ---------------------------------------------------------- %

% Save appdata
setappdata(CSI_1Dobj, 'data1D' , data1D);

% Save instance data
setappdata(CSI_1Dobj, 'instance', instance);

% CSIgui-1D: PLOT + PANEL % -------------------------------- %
% ---------------------------------------------------------- %

% Launch 1D plot-func
CSI_1D_displayData(CSI_1Dobj);

% Launch 1D control panel
try   
     CSI_1D_panel(instance); % Execute
     CSI_1D_movePanelToGUI(CSI_1Dobj); % Position
catch err
    CSI_Log({'Error: CSI 1D Panel could not open'},{err.message});
end

% --- Executed by panel-button in CSIgui_1D; panel to figure snap
function CSI_1D_movePanelToGUI(hObj, ~)
% Snap the panel to the current figure.
% Pa = panel; obj1D = figure;

% Get CSIgui 1D figure object
if ~strcmp(hObj.Type,'figure'), obj1D = hObj.Parent; else, obj1D = hObj; end

% Get its panel window
objPa = panel_1D_getInstanceData(hObj,'CSIpanel_1D'); 
if ~ishandle(objPa), return; end

% Get size and position 
posPa = objPa.Position; pos1D = obj1D.Position;

% Calculate new position
szPa = posPa(3:4); sz1D = pos1D(3:4); ps1D = pos1D(1:2);
newYX = (sz1D + ps1D); 
newYX(1) = newYX(1) + 5; newYX(2) = newYX(2) - sz1D(2);

% Set new position of panel
objPa.Position = [newYX szPa];

figure(objPa); % Figure to top.

% --- Executed to close the 1D Plot GUI
function CSI_1D_closeGUI(hObj, ~)
% Close both the 1D figure plot and the 1D buttons panel gui.

instance = getappdata(hObj,'instance');

% Close panel
panelObj = findobj('Tag','CSIpanel_1D');
for kk = 1:size(panelObj,1)
    tmp = getappdata(panelObj(kk),'instance');
    
    if strcmp(instance.tag, tmp.tag), delete(panelObj(kk)); end
end

% Close 1D plot figure.
delete(hObj);

% --- Executed to display data in CSIgui_1D; 1D plot GUI
function CSI_1D_displayData(CSI_1D_obj)
% Plots the 1D data to the axis in CSI_1D_obj (figure obj);

% Get 1D application object and actual data
[~, appdata1D, data] = CSI_1D_getData(CSI_1D_obj);

% CSIgui-1D: get plotdata % ------------------------- %

% 1. DATA (YAXIS)
% Convert to correct display type
plot_y = CSI_getUnit(data, appdata1D.voxel.unit);

% 2. AXIS (XAXIS)
if isfield(appdata1D.axis,appdata1D.axis.unit)
    plot_x = appdata1D.axis.(appdata1D.axis.unit);
else
    plot_x = appdata1D.axis.none;
end
if strcmp(appdata1D.axis.unit,'time'), revX = 0; else, revX = 1; end


% CSIgui-1D: plot data % ------------------------- %
plot(appdata1D.axes, plot_x, plot_y);
% ------------------------------------------------ %

% CSIgui-1D: plot cosmetics % -------------------- %

% X axis
if isequal(appdata1D.axis.unit ,'ppm')
    appdata1D.axes.XLim = appdata1D.axis.xlimit;
end
if revX, appdata1D.axes.XDir = 'reverse'; end

% Y axis
limy = [min(plot_y) max(plot_y)].*1.05;
if limy(1) >= limy(2), limy(1) = limy(1)-1; end
ylim(appdata1D.axes, limy);

% Title
title_str = CSI_1D_display_createTitle(appdata1D.voxel.index);
title(appdata1D.axes, title_str);

% --- Executed by CSI_1D_display to create voxel specific title string
function title_str = CSI_1D_display_createTitle(index)
% 2D to 1D: create plot title includes index in MRS data
% X index =  second index in csi.data.raw = col;
% Y index =  third index in csi.data.raw = row;
title_str = ['Col:' num2str(index(1)) ' ' ...   % Col
             'Row:' num2str(index(2)) '    ' ]; % Row
for kk = 3:size(index,2)
    title_str = [title_str sprintf('Z%i:%i  ',kk, index(kk))];
end

% --- Executed at start of CSIpanel_1D functions: data & objects
function [CSI_1DObj, appdata1D, data] = CSI_1D_getData(CSI_1DObj)
% Returns the 1D figure object, the 1D application data and the correct
% data for processing.

% Get instance
% instance = getappdata(CSI_1DObj, 'instance');

% Get 1D-appdata.
appdata1D = getappdata(CSI_1DObj, 'data1D');

% Get actual plot-data - original e.g. possibly complex.
if isfield(appdata1D.voxel, 'processed')
    % Take edited data if available.
    data = appdata1D.voxel.processed;
else
    % Else use original.
    data = appdata1D.voxel.original;
end

% --- Panel with additional functions to enable 1D data editting.
% All available function in 1D_panel found with panel_1D_[Func name];
function CSI_1D_panel(instance)
% Create a small window with specific functions applicable to the 1D
% spectrum. Processed data is stored in this 1D plot figure. Use
% Replace_Voxel e.g. SaveToMain to store the data back into CSIgui-main.
%
% Varargin - if private instance, the private tag is required.
% 


% CHECK EXISTANCE % ---------------------------------- %

% If panel exist, do not proceed again.
panelobj = findobj('Tag','CSIpanel_1D');
if instance.private == 0 
    for kk = 1:size(panelobj, 1)
       instanceOpened  = getappdata(panelobj(kk), 'instance'); 
       if instanceOpened.private == 0, delete(panelobj(kk)); end
    end
end

% CREATE FIGURE % ------------------------------------ %

% Get coloring from CSIgui
CSIgui_obj = findobj('type','figure','tag','CSIgui_main');
gui = guidata(CSIgui_obj);
if isfield(gui, 'colors')
    clr_bg = gui.colors.main; clr_butt = gui.colors.text_main;
else, clr_bg = [0 0 0]; clr_butt = [0.94 0.94 0.94];
end

% Figure: CSI_displayDataPanel
fh_1D = figure('Tag', 'CSIpanel_1D', ...
            'Name', 'CSIgui 1D Control Panel',...
            'Color', clr_bg, 'MenuBar','none', 'Toolbar', 'none', ...
            'resize', 'off','numbertitle', 'off'); 

% Get PC screensize for panel size
scrsz = get(0,'Screensize'); 
w = 180; h = 320;
% Set panel figure size: linked to the 1D plot!
set(fh_1D, 'position', [ceil(scrsz(3).*22/32) ceil(scrsz(4).*2/4) w h]);

% Save panel instance data
setappdata(fh_1D,'instance',instance);

% ADD BUTTONS % ------------------------------------ %

% Add buttons
bw = w-20; bh = 20; % Normalised w and h of button

% Buttons, handles and their info description.
bName = {'Phasing','FFT','iFFT','Apodization', 'Zero Fill', 'SNR',...
         'Linewidth','Baseline','Data Display','Replace Voxel', 'Export'};
bCall = {@panel_1D_PhasingMethod,...
         @panel_1D_FFT, ...
         @panel_1D_iFFT, ...
         @panel_1D_Apodization, ...
         @panel_1D_ZeroFill, ...
         @panel_1D_SNR,  ...
         @panel_1D_FWHM, ...
         @panel_1D_Baseline, ...
         @panel_1D_DataDisplay, ...
         @panel_1D_SaveToMain, ...
         @panel_1D_Export};
bInfo = {'Apply phasing corrections. Multiple methods available',...
         'Forward fourier to frequency domain: FID to spectrum.',...
         'Backwards fourier to time domain: Spectrum to FID.',...
         'Apply apodization: requires time domain data.',...
         'Apply zero filling: requires time domain data.',...
         'Calculate Signal-to-Noise-ratio of the 1D spectrum.',...
         'Calculate the full width at half maximum of a peak.',...
         'Apply a baseline correction to the spectrum.',...
         'Change display settings of the 1D plot.',...
         'Replace voxel data in CSIgui with edited 1D data.',...
         'Export spectrum to file.'};
         
% Add buttons     
nButtons = size(bName,2);
for bi = 1:nButtons
    uicontrol(fh_1D, 'Style','pushbutton','String', bName{1,bi},...
         'Tag', ['button_1D_' bName{1,bi}], 'Callback',bCall{bi},...
         'Position',[(w-bw)/2 ((h-(bi*bh))-bh)-((bh/4)*(bi-1)) bw bh],...
         'ForegroundColor', clr_butt,'BackgroundColor', clr_bg,...
         'ToolTipString',bInfo{bi});
end


% Create instance text
uicontrol(fh_1D, 'Style','Text', 'String', ['Instance: ' instance.tag],...
    'Unit','Pixels','Position',[0 0 125 15],...
    'Fontsize', 8, 'FontWeight', 'Bold',...
    'HorizontalAlignment','Left',...
    'ForegroundColor', clr_butt,'BackgroundColor', clr_bg);


% CSIgui-1D: FUNCTIONS % ----------------------------------------------- %

% --- Executes to get the correct instance of panel data and CSIgui1D
function obj1D  = panel_1D_getInstanceData(hObj,tag)

if nargin == 1, tag = 'CSIgui_plot1D'; end

% Get instance from object
if ~strcmp(hObj.Type,'figure')
    instance = getappdata(hObj.Parent,'instance');
else
    instance = getappdata(hObj,'instance');
end

% Get figures
obj1D_act = findobj('Tag',tag); % Figure handles
if isempty(obj1D_act), obj1D = NaN; return; end

% Get correct figure matching instance
for kk = 1:size(obj1D_act,1)
    instanceTmp = getappdata(obj1D_act(kk),'instance');
    if strcmp(instanceTmp.tag, instance.tag)
        obj1D = obj1D_act(kk); break;
    end
end

% No open CSI 1D figure nor matching found
if ~exist('obj1D','var') || isempty(obj1D), obj1D = NaN; return; 
end

% --- Executes when user presses "Phasing" in panel_1D
function panel_1D_PhasingMethod(hObj, ~)
% Gives the user choice to correct zeroth and first order phases manually,
% automatic (zeroth only) or apply manual phasing to all.


uans = getUserInput_Popup({'Choose desired phasing method: '},...
                    {{'Manual', 'Automatic', 'Phase All'}});
if isempty(uans), return; end  

switch uans{1}
    case 'Manual'
        panel_1D_PhaseCorrection_Manual(hObj);
    case 'Automatic'
        panel_1D_PhaseCorrection_Auto(hObj);
    case 'Phase All'
        panel_1D_PhaseCorrection_ApplyToAll(hObj);
end

% --- Executes if user chooses in "Phasing" -> "Manual"
function panel_1D_PhaseCorrection_Manual(hObj, ~)
% Apply manual phase correction to a 1D voxel.

% PROCESSING % --------------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);

if isfield(appdata1D.axis, 'ppm')
    xax = appdata1D.axis.ppm;
else
    xax = appdata1D.axis.none;
end

% 3. Send to PhaseCorrectionGUI
[data_phased, zero_ph, first_ph] = csi_PhaseCorrectionGUI(data,xax);

% Add 1D data structure to the 1D figure.
appdata1D.voxel.processed = data_phased;
appdata1D.phasing.zero = zero_ph; appdata1D.phasing.first = first_ph;

% CLEAN UP % ----------------------------------------- %

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D:'},{'Applied manual phase correction.'});

% --- Executes if user chooses in "Phasing" -> "Phase All"
function panel_1D_PhaseCorrection_ApplyToAll(hObj,~)
% See if there is any phase correction data. Else call this function first
% and then apply it to all voxel in the data set. Allow ability to exclude
% certain voxels.

% INITIATE % ------------------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D figure object, data_1D appdata and data array
[CSI_1D_obj, appdata1D] = CSI_1D_getData(obj1D);

% 2D/3D appdata
CSImain_obj = findobj('Tag','CSIgui_main');
if isempty(CSImain_obj)
    fprintf('Error: CSIgui appears to be closed. Returning.');
    return; % CSIgui is closed!
end
csi = getappdata(CSImain_obj, 'csi');

% FUNCTION PROCESS STARTS % ----------------------------- %

% 1. Check for available phasing corrections data
if ~isfield(appdata1D, 'phasing')
    % Manual Phasing required
    panel_1D_PhaseCorrection_Manual(hObj, []);
    % 1D appdata
    appdata1D = getappdata(CSI_1D_obj, 'data1D');
end

% 2. Apply phasing to all voxels.
magn = abs(csi.data.raw); pha = angle(csi.data.raw);

% Add zero order phase correction
pha_new = pha + appdata1D.phasing.zero;
% Add first order phase correction
pha_new = pha_new + appdata1D.phasing.first;
% Create complex data
csi.data.raw = complex(magn.*cos(pha_new), magn.*sin(pha_new));

% Save the appdata
setappdata(CSImain_obj, 'csi', csi);

% Replot the 2D plot!
CSI_2D_initiate2D();

% Update info to user.
CSI_Log({'CSIgui-1D:'},...
    {'Applied the phase correction to all voxel in the MRSI dataset.'});

% --- Executes if user chooses in "Phasing" -> "Automatic"
function panel_1D_PhaseCorrection_Auto(hObj, ~)
% Apply auto first order phase correction.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);

% USER INPUT % ---------------------------------- %

% Get peak of interest.
range = CSI_getPeakOfInterest(appdata1D.axis);

uans = getUserInput_Popup({'Method'},{{1,2}});
if isempty(uans), return; end
ph_meth = str2double(uans{1});

% APPLY CORRECTION % ----------------------------- %

% POI from user.
if length(range) > 1, poi = range(1):range(2); end

% Apply autophasing function
data_phased = csi_autoZeroPhase(data, poi, ph_meth, 0);

% PLOT & SAVE % ---------------------------------- %

% Add 1D data structure to the 1D figure.
appdata1D.voxel.processed = data_phased;

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D:'},...
    {'Applied automatic zeroth order phase correction.'});

% --- Executes when user presse "Zero Fill" in panel_1D
function panel_1D_ZeroFill(hObj, ~)

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);

% USER INPUT % --------------------------------- %

% Get zerofill target sample size + if equal for all dimensions
uans = getUserInput({'Requested sample size:'}, {size(data,1),1});
if isempty(uans), return; end

% Integer of target sample size
target_sz = str2double(uans{1});
if target_sz < size(data,1), return; end

% Get zerofill direction (pre/post/both)
uans = getUserInput_Popup({'Append direction:'},{{'Post','Pre','Both'}});
if isempty(uans), return; end
dir = uans{1}; % Direction for zerofilling of fids

% PROCESSING % --------------------------------- %

% Zero fill
appdata1D.voxel.processed = csi_zeroFill(data, target_sz, dir);

% UPDATE XAXIS LENGTHS % ---------------------- %
axs = appdata1D.axis;

% Update sample size
axs.N = target_sz;

% Update all axis: time, ppm, freq and none
if isfield(appdata1D.axis,'time')
    axs.time = linspace(axs.time(1),axs.time(end),target_sz);
    axs.ppm  = linspace(axs.ppm(1),axs.ppm(end),target_sz);
    axs.freq = linspace(axs.freq(1),axs.freq(end),target_sz);
end
axs.none = 1:target_sz;

% Save axis update
appdata1D.axis = axs;

% PLOT & SAVE % -------------------------------- %

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D: Applied zerofilling, target sample size:'},...
               {target_sz});

% --- Executes when user presses "SNR" in panel_1D
function panel_1D_SNR(hObj, ~)
% Calculate SNR for each peak or an individual peak in the current
% spectrum.


% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, appdata1D, data] = CSI_1D_getData(obj1D);

% ----- % Userinput
% All or single peak and snr-method

% Automatic or manual calculation.
uans = getUserInput_Buttons('Calculate SNR:',{'Automatic', 'Individual'});
if isempty(uans), return; end

% SNR method
method = getUserInput_Popup({'Signal unit: '},{{'Real','Magnitude'}});
if isempty(method), return; end
switch method{1}, case 'Real', method = 1; case 'Magnitude', method = 0; end

% ----- % Calculate SNR  

% Offset in y-axis used for plotting at peak positons
yst = diff(appdata1D.axes.YLim)/50; 

% Automatic or for single peak
if strcmp(uans,'Automatic') % All peaks
    
    % ------ % Find peaks
    
    % Positions of peaks
    [peak_pos, ~] = csi_findPeaks(data);
    
    % Convert peak position to ppm
    if isfield(appdata1D.axis,'ppm')
        peak_pos_plot = appdata1D.axis.ppm(peak_pos); 
    else                       
        peak_pos_plot = peak_pos;
    end
    
    % Plot peak locations
    hold(appdata1D.axes,'on'); 
    plot(appdata1D.axes, peak_pos_plot, real(data(peak_pos)), 'or');

    % ------ % Calculate SNR
    
    % Calculate SNR for found peaks
    np = size(data(peak_pos),1); snr = NaN(1,np);
    
    for kk = 1:np
        range = peak_pos(kk) + [-5 5]; % Range of the peak
        snr(kk) = csi_SNR(data, 50, method, range);
        % Plot SNR as text at peak location
        text(appdata1D.axes, peak_pos_plot(kk), data(peak_pos(kk))+yst, ...
            sprintf('%3.1f',snr(kk)),'FontSize',8,'FontWeight','Bold');
    end

else % Individual peak
    
    % ----- % Calculate SNR
    
    % Get range from user (easy to use the dataAtPeak function)
    range = CSI_getPeakOfInterest(appdata1D.axis);
    if isempty(range), return; end
    
    % Calculate SNR
    snr = csi_SNR(data, 50, method, range); % Output SNR
    
    % ----- % Display SNR 
    
    % At max vlaue in range
    tmp = range(1):range(2); [val, ind] = max(real(data(tmp,:))); 
    ind = tmp(ind);
    
    % Convert position to ppm values
    if isfield(appdata1D.axis,'ppm')
        ind = appdata1D.axis.ppm(ind); 
    end
    
    % Plot marker at peak location
    hold(appdata1D.axes,'on'); plot(appdata1D.axes, ind, real(val), 'or');
    % Plot SNR as text at peak location
    text(appdata1D.axes, ind, val+yst, ...
       sprintf('%3.1f',snr),'FontSize',8,'FontWeight','Bold');
end
hold(appdata1D.axes,'off'); % Turn off hold

% Show SNR to user.
CSI_Log({'CSIgui-1D: Voxel -','SNR:'},{appdata1D.voxel.index, snr});

% --- Executes when user presses "Linewidth" in panel_1D
function panel_1D_FWHM(hObj, ~)
% Calculate linewidth at FWHM at a certain peak of for all peaks found
% automatically.

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, appdata1D, data] = CSI_1D_getData(obj1D);

% Automatic or manual calculation?
uans = getUserInput_Buttons('Calculate linewidth at FHWM:',...
                           {'Automatic', 'Individual'});
if isempty(uans), return; end

% Prep xaxis % -------- %

% Use unitless or ppm
if isfield(appdata1D.axis, 'ppm')
    ax = appdata1D.axis.ppm;
else                         
    ax = appdata1D.axis.none;
end

% Prep peaks % -------- %

% Process individual peak or loop each found peak.
if strcmp(uans,'Individual')
    
    % Get range from user (Unitless range given)
    range = CSI_getPeakOfInterest(appdata1D.axis);
    if isempty(range), return; end
    
    % Maximum in this range
    [~,mi] = max(real(data(range(1):range(2))));
    % Set max as peak centre.
    range
    mi
    peak_pos = mi+range(1)
    
elseif strcmp(uans,'Automatic')
    % Get peak positions automatically (Maximum of peak!)
    peak_pos = csi_findPeaks(data);
end
    
% Peak width 
peak_width = ceil(appdata1D.axis.N/100);
if peak_width < 3, peak_width = 3; end


% Loop each peak. %  ------ %

sz = size(peak_pos); lw = NaN(sz); lwv = NaN(sz(1),2); lwp = NaN(sz(1),2);
yst = diff(appdata1D.axes.YLim)/50; xst = diff(appdata1D.axes.XLim)/100;
for pp = 1:size(peak_pos,1)

    % Set range by adding and subtracting peak width from peak position.
    range = peak_pos(pp); range = range + [-peak_width peak_width];

    if range(1) <0, range(1) = 1; end
    % Calculate FWHM % ---- %
    [lw(pp), lwv(pp,:), lwp(pp,:)] = ...
        csi_LineWidth(data, ax, range(1):range(2), 1);
    
    % Plot data
    hold(appdata1D.axes,'on'); % Hold 1D plot figure
    plot(appdata1D.axes, lwp(pp,:), real(lwv(pp,:)),'or');    % Plot marker
    text(appdata1D.axes, lwp(pp,1)-xst, lwv(pp,1)+yst, ...    % Plot text
       sprintf('%2.3f',lw(pp)),'FontSize',8,'FontWeight','Bold');
    
end
hold(appdata1D.axes,'off');

% Update in UI
CSI_Log({'Lowest position to highest, FWHM: '},{lw'});

% --- Executes when user presses "Baseline" in panel_1D
function panel_1D_Baseline(hObj,~)
% Apply a baseline correction to the spectrum.

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);

% Baseline correction % ----------------- %

% ... , bline, blevel, best
data_new = csi_baseline(data, 1, 0);

% Store in 1D app structure
appdata1D.voxel.processed = data_new;

% Plot and Save % ----------------------- %

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D: Applied baseline correction.'},{''});

% --- Executes when user presses "Replace" in panel_1D.
function panel_1D_SaveToMain(hObj, ~)
% Replace the original voxel data in the main CSIgui appdata with the
% edited 1D data shown.


% PREP DATA % ------------------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, appdata1D, data4replace] = CSI_1D_getData(obj1D);

% Check if CSIgui is still opened.
CSIgui_main_obj = findobj('Tag','CSIgui_main');
if isempty(CSIgui_main_obj)
    fprintf('Error: CSIgui appears to be closed. Returning.');
    return;
end % CSIgui is closed!

% REPLACE VOXEL % --------------------------------------- %

% Get voxel location in 3D data set
if sum(appdata1D.voxel.index) == 0, return; end 
ind = num2cell(appdata1D.voxel.index); % C/R/D4/D5 etc

% Get csi-main data struct.
if isappdata(CSIgui_main_obj, 'csi')
    csi = getappdata(CSIgui_main_obj, 'csi');
else
    CSI_Log({'Error: No main CSI data present in CSIgui.'},{''});
    return;
end

% Safety if sample size changes
if size(csi.data.raw,1) < size(data4replace,1)
    CSI_Log({'CSIgui-1D: Sample size changes detected.',...
       'Use multivoxel functions or export 1D spectrum to file.'},...
      {'Substitution into the dimensional data set is not possible.',''});
    return
end

% Replace the data in the voxel with the data from the 1D plot.
csi.data.raw(:,ind{:}) = data4replace;

% SAVE AND PLOT % --------------------------------------- %

% Save the appdata
setappdata(CSIgui_main_obj, 'csi', csi);
% Replot the 2D plot!
CSI_2D_initiate2D();

% --- Executes when user presses "FFT" in panel_1D.
function panel_1D_FFT(hObj, ~)
% Time to frequency domain.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end

% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data2fft] = CSI_1D_getData(obj1D);

% FFT % ---------------------------------------- %

% Create a N x 1 vector!
if size(data2fft,2) > size(data2fft,1), data2fft = data2fft'; end

% Apply 1D fourier transform
appdata1D.voxel.processed = csi_fft(data2fft);

% PLOT AND SAVE % ------------------------------ %

% Save the complex result to the 1D-plot figure's appdata.
setappdata(CSI_1D_obj, 'data1D', appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj)

% Show info to user
CSI_Log({'CSIgui-1D:'},{'Applied FFT.'});

% --- Executes when user presses "iFFT" in panel_1D.
function panel_1D_iFFT(hObj, ~)
% Frequency to time domain.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data2ifft] = CSI_1D_getData(obj1D);

% iFFT % --------------------------------------- %

% Create a N x 1 vector.
if size(data2ifft,2) > size(data2ifft,1), data2ifft = data2ifft'; end

% Apply 1D fourier transform
appdata1D.voxel.processed = csi_ifft(data2ifft);

% PLOT AND SAVE % ------------------------------ %

% Save it to the 1D-plot figure's appdata.
setappdata(CSI_1D_obj, 'data1D', appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj);

% Show info to user
CSI_Log({'CSIgui-1D:'},{'Applied iFFT.'});

% --- Executes when user presses "Apodization" in panel_1D.
function panel_1D_Apodization(hObj,~)
% Apply apodization using different filter windows.
% Gaussian, exponential, hamming, hann, flattop, blackman, sinebell.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);


% USER INPUT % --------------------------------- %

% Get filter to apply.
uanstype = getUserInput_Popup({'Filter Method'},...
                          {{'Gaussian', 'Hamming','Hann', 'Exponential',...
                            'Blackman', 'Flattop'}});
if isempty(uanstype), return; end

% Additional options for specific filters
switch uanstype{1}
    case 'Gaussian'
        uansopts = getUserInput(...
            {'Standard deviation e.g. lentgh of FWHM (Samples): '},{5});
    case 'Exponential'
        uansopts = getUserInput(...
            {'Exponential decay strength (Samples): '},{0.5});
end

% Set options if applicable
if exist('uansopts', 'var'), opts = str2double(strsplit(uansopts{1},' '));
else, uansopts = {''}; opts = 0;
end

% APODIZATION % ---------------------------------- %
% Create window and apply filter
data = CSI_filterSpectra(data, uanstype{1}, opts);

% Apply window.
appdata1D.voxel.processed = data;

% hann, hamming, blackman, flattop, gaussian, exponentional decay
% Res enhancement: gauss, sinebell
% Line broadening - Exponential like

% PLOT AND SAVE % ------------------------------ %

% Save it to the 1D-plot figure's appdata.
setappdata(CSI_1D_obj, 'data1D', appdata1D);

% Display data.
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({...
    ['CSIgui-1D: Filter ' uanstype{1} ' applied. Possible opts:']},...
    {uansopts{1}});

% --- Executes when user presses "Data display" in panel_1D;
function panel_1D_DataDisplay(hObj, ~)
% Change specific display settings of the 1D data display.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D] = CSI_1D_getData(obj1D);


% USER INPUT % -------------------------------- %
uans = getUserInput_Popup({'X-axis type: ','Plot unit: '},...
                          {{'PPM', 'Frequency','Time','None'},...
                           {'Real', 'Magnitude','Phase','Imaginary'}});
if isempty(uans), return; end

% User input handling
axis_unit = lower(uans{1}); plot_unit = uans{2};

% Edit some x axis labels...
if strcmp(axis_unit,'frequency'), axis_unit = 'freq'; end 

% APPLY SETTINGS % ---------------------------- %
appdata1D.voxel.unit = plot_unit;
appdata1D.axis.unit = axis_unit;
          
% SAVE AND DISPLAY % -------------------------- %

% Save the new appdata
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot
CSI_1D_displayData(CSI_1D_obj);

% --- Executes when user presses "Export" in panel_1D.
function panel_1D_Export(hObj, ~)
% Export shown spectrum or fid to file
% Return if no CSIgui plot 1D is opened.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, ~, data2export] = CSI_1D_getData(obj1D);

% SAVE TO FILE % ------------------------------ %

% Get file path and extension from user
[fn, fp, fi] = uiputfile(...
            {'*.txt', 'Text file';'*.SDAT','SDAT file';...
            '*.mat','MAT-file'},'Save 1D data as...');
if fi == 0, return; end

% Get extension.
[~, fn ,ext] = fileparts(fn);
% Create full path to filename.
fpn = [fp fn ext];

% Operate accoring extension
switch lower(ext)
    case '.txt'         % TXT % ------------------------------ %
        
        % Create complex and real part columns
        dataR = real(data2export); dataI = imag(data2export);
        
        % Data array to save
        dataA = cat(2,dataR, dataI);
        
        % Write array to file
        fid = fopen(fpn,'wt');
        for ri = 1:size(dataA,1)
            fprintf(fid,'%32.32f %32.32f\n', dataA(ri,:)); 
            % fprintf(fid,'%16.16f %16.16f\n', dataA(ri,:)); 
        end
        fclose(fid);
    case '.sdat'        % SDAT % ---------------------------- %
        
        % Ask user if freq -> time convert is required.
        qst = 'SDAT requires time-domain data. Convert data? (y/n)';
        sparQ1 = 'Sample frequency: ';
        ua = getUserInput({qst,sparQ1},{'n', '8000'});
        % Convert if requested.
        if isempty(ua), return; end
        
        if strcmp(ua{1},'y')
            % Convert to time domain by using iFFT
            data2export = csi_ifft(data2export);
            CSI_Log({'CSIgui-1D: Converted data to time domain'},...
                           {'using iFFT before saving to file.'});
        end
        CSI_Log({'CSIgui-1D: SPAR requires manual updating',...
                        'CSIgui-1D: SPAR tested with Matlab and 3D CSI'},...
                       {'to work with specific applications.',''});
        % Save.
        csi_writeSDAT(fpn , data2export, 'fid');
        
        % Add required parameters to SPAR file
        % 3D CSI: samples, sample_frequency
        % JMRUI:  untested for 1D data
        % Matlab: enabled in csi_writeSDAT - requires sizes per dimension.
        csi_editSPAR([fpn(1:end-4) 'SPAR'], ...
            {'samples' num2str(size(data2export,1)) ...
             'sample_frequency' ua{2}});
    case '.mat'
        data_1D = data2export; save(fpn, 'data_1D');
end

% Show info to user
CSI_Log({'CSIgui-1D:'},{'Saved 1D data to file.'});

% ---------------------------------------------------------------------- %

% --- Called everywhere to update CSIinfo listbox
function CSI_Log(varargin)
% Expects: infolabel, infoval as input
% Updates listbox in main CSIgui displaying csi-data information.
% Requires the info to be send as input
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
if nargin == 2
    infolabel = varargin{1}; infoval  = varargin{2};

    % 1. PROCESS USER INPUT
    new_info = cell(1,size(infolabel,2));
    for kk = 1:size(infolabel,2)
       val = infoval{kk}; 
       if ~ischar(val), val = num2str(val); end % Convert to a double.
       new_info{kk} = sprintf('%s - %s %s',datestr(rem(now,1),'HH:MM:SS'), ...
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
    % Always shows the ACTUAL array size to catch errors in csi.data.dim which
    % also saves the dimensions equal to one. 

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
% str_new = cat(1, stand_info, new_info, str_save);
str_new = {stand_info{:} new_info{:} str_save{:}};

try 
    % Print to info-box.
    set(gui.listbox_CSIinfo, 'Value', 1, 'String', str_new);
catch 
    % Warn if something goes wrong.
    warning('Something in CSI_Log went wrong. No list box update!');
end
drawnow;

% --- Executed by X in gui to delete LOG-line
function CSI_Log_deleteLine(hObj)
% Delete a line from LOG list box.

% Get GUIdata from CSIgui_main;
gui = guidata(hObj);

% Get selected line value and str-size
str = gui.listbox_CSIinfo.String; val = gui.listbox_CSIinfo.Value;
sz = size(str,1);

% Delete line and update LOG
if val == 1
    gui.listbox_CSIinfo.Value = 1;
elseif val < sz-1
    gui.listbox_CSIinfo.Value = val-1;
elseif val == sz
    gui.listbox_CSIinfo.Value = sz-1;
end
if val <= sz, str(val) = []; gui.listbox_CSIinfo.String = str; end

% ---------------------------------------------------------------------- %

% --- View noise component executed in menu
function CSI_viewNoise(hObject, ~, ~)
% Gui-data struct.
gui = guidata(hObject);


% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if isfield(csi.data, 'noise')
    
    if isfield(csi.data, 'backup_toview_noise')
        % Return to original data
        csi.data.raw = csi.data.backup_toview_noise;
        % Remove the backup field.
        csi.data = rmfield(csi.data,'backup_toview_noise');
        % Message
        msg = 'Replaced noise data by backup of original.';
        % Set correct menu-label
        gui.menubar.MRSI.noise.Label = 'Import Noise..';
    else    
        % Create backup raw
        csi.data.backup_toview_noise = csi.data.raw;
        % Replace raw with noise.
        csi.data.raw = csi.data.noise;
        % Info for user
        msg = 'Noise data inserted. Backup of original created.';
        % Set correct menu-label
        gui.menubar.MRSI.noise.Label = 'Revert to Data..';
    end
    
    % Save appdata.
    setappdata(gui.CSIgui_main, 'csi', csi);
    
    %Show user.
    CSI_Log({msg},{''}); 
else
    CSI_Log({'No noise data present.'},{''}); return;
end
   
% ---------------------------------------------------------------------- %

% --- Executed via menubar to save data to file
function CSI_saveData(hObject, ~, ~)
% Gui-data struct.
gui = guidata(hObject);

% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Find default path to open
if isfield(csi.data, 'filepath'), fp = csi.data.filepath;
else, fp = [];
end
           
% UI for file from user
[fn, fp] = uiputfile({'*.sdat','Sdat file'; '*.txt','Text file'; ...
                  '*.mat','MATLAB file'},'Save MRSI data to file...',fp); 
if fn == 0, return; end
% Analyze file extension.
[~,fn,ext] = fileparts(fn);

switch ext
    case '.sdat'
         CSI_saveSDAT(hObject, [], [], fp, fn, ext)
    case '.txt'
         CSI_saveTXT(hObject, [], [], fp, fn, ext)
    case '.mat'
         CSI_saveMAT(hObject, [], [], fp, fn, ext)
end

% --- Save CSIgui main CSI data as SDAT.
function CSI_saveSDAT(hObject, ~, ~, varargin)
% Save CSI data stored as SDAT file.
% varargin{1:3} expected as filepath, filename and extension.

% Gui-data struct.
gui = guidata(hObject);

% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');


% USERINPUT % ---------------------------------------- %

% File location and name if not as input (varargin)
if nargin <= 3
    % UI for file from user
    [fn, fp] = uiputfile('*.SDAT','Save as SDAT'); 
    if fn == 0, return; end
    [~,fn,ext] = fileparts(fn);
else
    fp = varargin{1}; fn = varargin{2}; ext = varargin{3};
end

% Create full filepath for SDAT.
if strcmp(fp(end), '\'), fpn_sdat = [fp fn ext];
else,                    fpn_sdat = [fp '\' fn ext];
end

% Inverse fourier for correct domain e.g. time.
uans = getUserInput_Buttons(...
    {'File format requires data in time domain. Apply iFFT?'},...
    {'Yes','No'});
if isempty(uans), return; end


% PREP SDAT DATA % ----------------------------------------- %
% Required to work with 3D CSI R2 Proton loading:
% nr_of_slices_for_multislice
% samples
% sample_frequency
% rows
% nr_phase_encoding_profiles
% nr_of_phase_encoding_profiles_ky
% slice_distance                % Slice thickness.
% slice_thickness               % FOV in slice direction.
% phase_encoding_fov            % FOV in phase encoding e.h Kx direction.
% nr_of_slices_for_multislice

% Get data
data = csi.data.raw; labels = csi.data.labels;

% Convert by iFFT if requested.
if strcmp(uans,'Yes'), data = csi_ifft(data); end

% SAVE SDAT  % --------------------------------------------- %

status = csi_writeSDAT(fpn_sdat , data, 'fid', labels);
if (status == numel(data)*2), sdat_result = 'Success! SDAT file created.';
else, sdat_result = 'Fail. Data corruption possible.';
end

% Update info to user.
CSI_Log({'Filepath:', 'SDAT write status:'},...
               {fpn_sdat, sdat_result});

% EDIT SPAR % --------------------------------------------- %  
CSI_Log({'Initiate updating SPAR parameters.'},...
               {'Possible input required.'});

% SPAR info edit struct;
spans = struct; % SPAR answers struct;

% AVAILABLE INFO % ---------------------------------------- %

% Size of data in X, Y and Z.
spans.samples = size(data,1); 
spans.rows = size(data,2).*size(data,3);
spans.nr_phase_encoding_profiles = size(data,2);
spans.nr_of_phase_encoding_profiles_ky = size(data,3);
spans.nr_of_slices_for_multislice = size(data,4);
% Nr of averages.
aver_ind = csi_findDimLabel(csi.data.labels, {'aver'});
if ~isnan(aver_ind)
    spans.averages = size(data,aver_ind);
end

% AVAILABLITY SPECIFIC INFO % -------------- %

% Required fields for proper spar file handling in other software:
req_freq = {'nucleus', 'synthesizer_frequency', 'sample_frequency'};
req_ori  = {'phase_encoding_fov', 'slice_distance', 'slice_thickness'};

% SPAR: Frequency
if isfield(csi.xaxis, 'ppm')
    spans.sample_frequency = csi.xaxis.BW; 
    spans.synthesizer_frequency = csi.xaxis.trans;
    spans.nucleus = csi.xaxis.nucleus;
end

% SPAR: Data dimensions
if isfield(csi, 'ori')
    spans.slice_distance = csi.ori.res(3); 
    spans.slice_thickness = csi.ori.fov(3);
    spans.phase_encoding_fov = csi.ori.fov(1);
end

% USERINPUT % -------------------------------------- %
% Create remaining questionaire for user.
req_freq_bool = NaN(1,size(req_freq,2));
for fri = 1:size(req_freq,2) 
    if ~isfield(spans, req_freq{fri}), req_freq_bool(fri) = fri; end
end
req_freq_bool(isnan(req_freq_bool)) = [];

req_ori_bool = NaN(1,size(req_ori,2));
for or = 1:size(req_ori,2)
    if ~isfield(spans, req_ori{or}),req_ori_bool(or) = or; end
end
req_ori_bool(isnan(req_ori_bool)) = [];
% Extract quests of interest for user
spqst = cat(2,req_freq(req_freq_bool), req_ori(req_ori_bool));

% Get userinput and process.
if ~isempty(spqst)
    uans = getUserInput(spqst,repmat({'1'},1,size(spqst,2)));
    if isempty(uans), return; end
    % Add answers to spans struct;
    for qi = 1:size(spqst,2), spans.(spqst{qi}) = uans{qi};
    end
end

% CREATE SPAR EDIT INPUT % ----------------------------- %
% Create correct input for editting SPAR file function
sparfn  = fieldnames(spans); % All parameters to edit in SPAR
spar_edit_input = cell(1,size(sparfn,1).*2); % Empty cell for SPAR function
sparinp = getFields(spans, sparfn); % All input parameters gathered
for kk = 1:size(sparfn,1)
    spar_edit_input{(kk*2)-1} = sparfn{kk};
    spar_edit_input{(kk*2)} = sparinp{kk};
end

% UPDATE SPAR FILE % -------------------------------------- %
csi_editSPAR([fpn_sdat(1:end-5) '.SPAR'], spar_edit_input);

% Update info to user.
CSI_Log({'SPAR parameters updated successfully.'},{''});

% --- Save CSIgui main CSI data as TEXT.
function CSI_saveTXT(hObject, ~, ~, varargin)
% varargin{1:3} expected as filepath, filename and extension.

% Get guidata
gui = guidata(hObject);

% Check for CSI data
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

if nargin <= 3
    % Default filepath to start UI in.
    if isfield(csi.data,'filepath'), fp_def = csi.data.filepath; 
    else, fp_def = []; 
    end

    % Get filepath destination
    % Get file path and extension from user
    [fn, fp, fi] = ...
        uiputfile({'*.txt', 'Text file'},'Save MRSI data...',fp_def);
    if fi == 0, return; end
    [~, fn, ext] = fileparts(fn);
else
    fp = varargin{1}; fn = varargin{2}; ext = varargin{3};
end
if ~strcmp(fp(end),'\'), fp = [fp '\']; end

% Get data to export
data2exp = csi.data.raw;

% Display information to user.
CSI_Log({'Writing to text file. Please be patient.'},{''});

% Write text file
csi_writeText(data2exp, [fp fn ext]);

% % Convert to a two column array
% resh2exp = reshape(data2exp,[], 1); 
% % Create complex and real part columns
% dataR = real(resh2exp); dataI = imag(resh2exp);
% % Create 2 column array with real and imaginary part.
% dataA = cat(2,dataR, dataI);
%
% % Write array to file
% fid = fopen([fp fn ext],'wt');
% for ri = 1:size(dataA,1)
%     fprintf(fid,'%32.32f %32.32f\n', dataA(ri,:)); 
% end
% fclose(fid);
% 
% % Also write away a size to reconstruct the data.
% fid = fopen([fp fn '_ArraySize' ext],'wt');
% fprintf(fid,'%32.32f ', size(csi.data.raw)); 
% fclose(fid);

% Display information to user.
CSI_Log({'MRSI data and array size written to text file:'}, {[fp fn]});
        
           
           
           
% --- Save CSIgui main CSI data as MAT.
function CSI_saveMAT(hObject, ~, ~,varargin)
% varargin{1:3} expected as filepath, filename and extension.

% Get guidata
gui = guidata(hObject);

                 % ------- % Save CSI data % ------- %

% Check for CSI data       
if isappdata(gui.CSIgui_main,'csi')
    csi = getappdata(gui.CSIgui_main,'csi'); 

    % Get file path from user if not as input
    if nargin <= 3
        if isfield(csi.data,'filepath'), fp = csi.data.filepath;
        else, fp = [];
        end
           
        % Get file path and extension from user
        [fn, fp, fi] = uiputfile({'*.mat', 'Matlab file'},...
            'Save MRSI data...',fp);
        if fi == 0, return; end
        [~,fn, ext] = fileparts(fn);
    else
        fp = varargin{1}; fn = varargin{2}; ext = varargin{3};
    end
    % Get MRSI data
    csigui = csi.data;

    % Add log
    log = get(gui.listbox_CSIinfo,'String');
    csigui.log = char(log);

    % Add xaxis structure
    csigui.xaxis = csi.xaxis;
    
    % Add Ori
    if isfield(csi,'ori'), csigui.ori = csi.ori; end
end

                  % ------- % Save CONV data % ------- %

% Check for converted data
if isappdata(gui.CSIgui_main,'conv')
    csigui.conv = getappdata(gui.CSIgui_main,'conv');
end

                  % ------- % Save MRI data % ------- %

% Check for mri data
if isappdata(gui.CSIgui_main,'mri')    
    csigui.mri = getappdata(gui.CSIgui_main,'mri');
end


                   % ------- % Save file  % ------- %

% Save it.
save([fp fn ext], 'csigui','-v7.3');

% Done.
CSI_Log({'Saved MRSI data to MAT-file:'},{[fp fn ext]});

% --- Save CSIgui 2D plot to file       
function CSI_saveFig(hObject, ~, ~)
% Save the 2D CSIgui plot to file supporting multiple formats;
% jpg, png, eps, tiff, fig.


%%%%% SOME INFO %%%%%
% Okay how to do it for now:
% 1. Export options --> resolution (in dpi or screen size) 
        % (redun for eps) use painters inst opengl
% 2. Export spectra only or merged w image
        % Uhm.. how? Get rid of that object in figure . print
% 3. Export image seperatly at MRI > Export > Images
        % Mat + PNG etc. ---> Save current displayed slice or all images as
        % PNG :)
% export_fig([fp fn], '-transparent','-png','-nocrop','-m8', fig_obj);


% Get guidata
gui = guidata(hObject);

% Check for CSI data
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

% Get filepath destination
% Get file path and extension from user
[fn, fp, fi] = uiputfile({'*.png','Portable network graphic';...
                          '*.jpg', 'JPEG image';'*.eps','EPS file';...
                          '*.tiff','TIFF image';'*.bmp','Bitmap file';...
                          '*.fig','MATLAB Figure'},...
                          'Save CSIgui 2D figure...');
if fi == 0, return; end % Return if canceled.
% Analyse userinput
[~, fn, ext] = fileparts(fn);

% Set correct extension format type
switch ext
    case '.jpg',    ftype = 'jpeg';
    case '.eps',    ftype = 'epsc';
    case '.tiff',   ftype = 'tiffn';
    otherwise,      ftype = ext(2:end);
end

% Set correct renderer.
switch ext
    case 'epsc', render = '-painters';
    otherwise,   render = '-opengl';
end

uans = getUserInput_Popup({'Specify figure(s) to save:','Resolution (DPI):'},...
                         {{'Current figure', 'All'},...
                           cat(2,num2cell(0:200:1200),'Custom')});
if isempty(uans), CSI_Log({'Skipped exporting figures.'},{''}); return; 
end  

% Process waht to save
switch uans{1}
    case 'Current figure', save_all = 0;
    case 'All', save_all = 1;
end
% Process resolution
if strcmp(uans{2},'Custom')
        uans{2} = getUserInput({'Custom DPI scaling:'},{''});
end
dpi = str2double(uans{2});


% FIGURE OBJECT % ----------------------------------------- %

% Check if 2D Plot is active
fig_obj =  findobj('type','figure','tag','CSIgui_plot2D');
% If 2D Plot not active, open it.
if isempty(fig_obj),panel_2D_DataSliders([],[],gui); CSI_2D_initiate2D(); end

% Get 2D-Plot and 2D-Panel figure
fig_obj = findobj('type','figure','tag','CSIgui_plot2D');
pan_obj = findobj('type','figure','tag','CSIpanel_2D_DataToDisplay');
if ~isempty(pan_obj)
    pan_gui = guidata(pan_obj);
else
    
    % 4. Save the figure to image file
    fig_obj.InvertHardcopy = 'off';
    print(fig_obj,[fp fn sprintf('_%i_%i',[1 1 1]) ext],...
          ['-d' ftype], ['-r' num2str(dpi)], render);
      
%     export_fig([fp fn sprintf('_%i_%i',[1 1 1]) '_B'], ...
%         '-transparent',['-' ftype],'-nocrop', '-m2', fig_obj);
    return;
end

% FIGURE DATA % ------------------------------------------- %

if save_all
    % Data slice dimensions
    dim = csi.data.dim(4:end); 
    % All figure indexes to cover
    indArray = slice2rowIndex(num2cell(dim));
else
    % Data slice dimensions
    dim = csi.data.dim(4:end);
    % Get index of current figure - Loop each slider
    indArray = NaN(1,numel(dim));
    for sli = 1:size(pan_gui.sliders,2)
        indArray(sli) = pan_gui.sliders{sli}.Value;
    end
end


% LOOP & SAVE % -------------------------------------------- %
% Loop each index & Save the figure
for indi = 1:size(indArray,1)
    % 1. Set index at the slider
    for sli = 1:size(pan_gui.sliders,2)
        pan_gui.sliders{sli}.Value = indArray(indi,sli);
    end
    
    % 2. Replot 2D CSI plot.
    CSI_2D_initiate2D();
    
    % 3. Get new 2D CSI plot figure object
    fig_obj = findobj('type','figure','tag','CSIgui_plot2D');
    
    % 4. Save the figure to image file
    fig_obj.InvertHardcopy = 'off';
    print(fig_obj,[fp fn sprintf('_%i_%i',indArray(indi,:)) ext],...
          ['-d' ftype], ['-r' num2str(dpi)], render);
      
%     export_fig([fp fn sprintf('_%i_%i',indArray(indi,:)) '_B'], ...
%         '-transparent','-png','-nocrop', '-m2', fig_obj);
end

% --- Save CSIinfo listbox log as text-file
function CSI_saveLog(hObject, ~, ~)
% Get GUI data
gui = guidata(hObject);

% Get string of listbox
lb_str = get(gui.listbox_CSIinfo,'String');

% Save as text-file
[fn,fp] = uiputfile( {'*.txt','Text-file'},'Save log-file');

% Check if user pressed cancel, closed the window or something else.
if isequal(fn, 0) || isequal(fp, 0), return; end

% Write file.
fid = fopen([fp fn],'wt');
fprintf(fid, '%s\n', lb_str{:});
fclose(fid);


% CSI Backup System % -------------------------------------------------- %
% ---------------------------------------------------------------------- %

% --- Executes on button press in button_setBackup.
function button_setBackup_Callback(~, ~, gui)
% Store the csi.data field into csi.backup
% csi.backup.(bu time).data = csi.data;
% csi.backup.(bu time).tag  = info_str;
CSI_backupSet(gui, 'User backup.');

% --- Executes on button press in button_getBackup.
function button_getBackup_Callback(~, ~, gui)
% Open up a menu and select backup of interest to revert to.
% Requires fields shown below:
% csi.backup.(bu time).data = csi.data;
% csi.backup.(bu time).tag  = info_str;
CSI_backupGet(gui);

% --- Create a backup of mrsi data
function CSI_backupSet(gui, info_str)
% Create a backup of the csi appdata. Info_str is an option to show a
% specific reason for backing up. Displayed in the listbox of CSIinfo.
%
% E.g. Store the csi.data field into csi.backup
%      csi.backup.(time of backup).data = csi.data;
%      csi.backup.(time of backup).tag  = info str;

% If not info_str given.
if nargin == 1, info_str = ''; end

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Safety pause.
pause(0.15);

% Create backup details.
% Backup time
bup_time = sprintf('%s',datestr(rem(now,1),'HH:MM:SS'));
% Backup tag
bup_tag = info_str;
% Backup fieldname
bup_fn = ['backup_' strrep(bup_time,':','_')];

% Create a new backup.
csi.backup.(bup_fn).data = csi.data;
csi.backup.(bup_fn).tag  = bup_tag;

% Check if there was a "lastBackup" field and remove.
if isfield(csi.backup,'lastBackup')
    csi.backup = rmfield(csi.backup,'lastBackup'); 
end


% Store and display
setappdata(gui.CSIgui_main, 'csi', csi);

% Show user.
CSI_Log({'Created a backup of CSI app-data.'},{info_str});

% Check with backupManager if memory is not full.
CSI_backupManager(gui);

% --- Revert to a backup of mrsi data.
function CSI_backupGet(gui, askUser)
% Overwrite the csi.data field with the backup field.
% csi.data = csi.backup;

if nargin == 1, askUser = 1; end

% Check for CSI data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Check for backup data
if ~isfield(csi,'backup')
    CSI_Log({'No MRSI backup available.'},{''}); return;
end

% AVAILABLE BACKUP CHOICE % -------------------- %

% Get list of available backup fields
bup_list = fieldnames(csi.backup);

% Remove lastBackup if included
loc = contains(bup_list,'lastBackup'); bup_list(loc) = [];

if askUser == 1 % Ask user which backup to get...
    % Get tags of backups
    tag = cell(size(bup_list,1),1);
    for fi = 1:size(bup_list, 1)
        tag{fi} = csi.backup.(bup_list{fi}).tag;
    end

    % flip the tag and backup list
    bup_list = flipud(bup_list); tag = flipud(tag);

    % Create backup list for display
    bup_disp = cellfun(@cat, repmat({2},size(tag)), ...
                    bup_list, repmat({' '}, size(tag)),tag, 'UniformOutput',0);
    % Replace underscores            
    bup_disp = strrep(bup_disp, '_', ':'); 
    col_ind = strfind(bup_disp,':');  col_ind = col_ind{1}(1);
    for kk = 1:size(bup_disp,1), bup_disp{kk}(col_ind) = ' '; end

    % Display to user to select backup of interest
    uans = getUserInput_Popup({'Available backups: '},{bup_disp});
    if isempty(uans), return; end

    % Process user choice: backup of interest
    boi_ind = strcmp(uans, bup_disp); boi = bup_list{boi_ind};
    tagoi = tag{boi_ind};

else % NO user input required.
    
    % This is for undo/ctrl+z
    if isfield(csi.backup,'lastBackup') % Is there a previous backup used.
        prev_boi = csi.backup.lastBackup; 
        % Get index of last backup
        ind = contains(bup_list, prev_boi); ind = find(ind,1);
        % Get backup before last backup.
        if ind>1
            boi = bup_list{ind-1}; tagoi = csi.backup.(boi).tag;
        else
            CSI_Log({'No previous backup available.'},...
                           {'Returning.'});
            return; 
        end
    else  % Get latest set backup;
        boi_ind = size(bup_list,1);
        boi = bup_list{boi_ind}; tagoi = csi.backup.(boi).tag;
    end
end


% SET CHOSEN BACKUP % --------------------------- %
csi.data = csi.backup.(boi).data;
csi.backup.lastBackup = boi;

% SAVE AND CLEAN UP % --------------------------- %
% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);
% Show user.
CSI_Log(...
    {['Reverted to backup ' strrep(boi(8:end),'_',':') ]},{tagoi});

% --- Revert to previous backup: ctrl+z or menubar
function CSIgui_Undo(hObj,~,~)
% This function launches CSI_backupGet without any required user input. The
% CSI_backupGet will use lastBackup field in csi if available and will get
% the latest backup created or a backup before the previously restored 
% backup.

% Get current gui data
gui = guidata(hObj);
% Get the latest backup: setting zero to get latest backup without option
% to choose between all backups
CSI_backupGet(gui,0);

% --- Executes by CSIbackupSet
function CSI_backupManager(gui)
% To manage data memory, if more then N backups are present, the last
% created is deleted.

N = 15; % Max number of backups.

% Check for CSI data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Check for backup data
if ~isfield(csi, 'backup'), return; end

% Check for #backups
bup_list = fieldnames(csi.backup);
bup_sz = size(bup_list,1);

if bup_sz > N % Delete oldest backup
    csi.backup = rmfield(csi.backup,bup_list{1});
end

% Save appdata.
setappdata(gui.CSIgui_main, 'csi', csi); % Done.


% IMG BUTTONS % -------------------------------------------------------- %
% ---------------------------------------------------------------------- %

% ---------------------------------------------------------------------- %
% Some IMG coodinate functions are located at COORDINATES paragraph which
% also containes the CSI coordinates function(s).
% ---------------------------------------------------------------------- %

% --- Executes on button press in button_plotIMG.
function button_plotIMG_Callback(hObj, ~, ~)

% Get GUI data
gui = guidata(hObj);

% Plot MRI data as is
if isappdata(gui.CSIgui_main, 'mri')
    mri = getappdata(gui.CSIgui_main,'mri'); 
    % Get all possible image type fields
    fns = fieldnames(mri.data);
    % Plot each image type with display3D
    for kk =1:size(fns)
        im = mri.data.(fns{kk});
        display3D(im,'tag',fns{kk},...
            'limit',[200 max(im(:)).*0.02]); 
        
    end
end

% Plot converted MRI data.
if isappdata(gui.CSIgui_main, 'conv')
    conv = getappdata(gui.CSIgui_main,'conv'); 
    display3D(conv.data,'tag','Converted');
end

% --- Executes on button press in button_MRI_setContrast.
function button_MRI_setContrast_Callback(hObject, eventdata, gui)
% Set the contrast for plotting the image behind the CSI data

% Get data.
if ~isappdata(gui.CSIgui_main, 'mri'), return; end
mri = getappdata(gui.CSIgui_main, 'mri');

% Get current contrast or calculate.
if isfield(mri, 'contrast')
    disp_cont =  mri.contrast; 
else
    fn = fieldnames(mri.data);
    disp_cont = [min(mri.data.(fn{1})(:)) max(mri.data.(fn{1})(:))];
end

% Get new contrast from user.
uans = getUserInput({'Contrast (min/max):'},{disp_cont});
if isempty(uans), return; end

% Define the new input contrast.
contr_new = strsplit(uans{1},' ');
contr_new = str2double(contr_new);

% Save - done.
mri.contrast = contr_new; setappdata(gui.CSIgui_main, 'mri',mri);
if isappdata(gui.CSIgui_main, 'conv')
    conv = getappdata(gui.CSIgui_main, 'conv');
    conv.contrast = mri.contrast;
    setappdata(gui.CSIgui_main, 'conv', conv);
end


% Close CSIgui % ------------------------------------------------------- %


% --- Executes when user attempts to close CSIgui_main from menu.
function CSIgui_main_CloseFromMenu(hObject, ~)
gui = guidata(hObject); CSIgui_main_CloseRequestFcn(gui.CSIgui_main);

% --- Executes when user attempts to close CSIgui_main.
function CSIgui_main_CloseRequestFcn(hObject, ~, ~)

% Close the data to display panel: 2D
obj_panel2D = findobj('Tag','CSIpanel_2D_DataToDisplay');
if ~isempty(obj_panel2D), delete(obj_panel2D); end

% Close the 2D data figure
obj_fig2D = findobj('Tag', 'CSIgui_plot2D');
if ~isempty(obj_fig2D), delete(obj_fig2D); end

% Close the edit 1D data panel figure
obj_panel1D = findobj('Tag', 'CSIpanel_1D');
if ~isempty(obj_panel1D), delete(obj_panel1D); end

% Close the 1D data figure
obj_fig1D = findobj('Tag', 'CSIgui_plot1D');
if ~isempty(obj_fig1D), delete(obj_fig1D); end

% Close CSIgui_main.
delete(hObject);


% GUI COLOR THEME % --------------------------------------------------- %


function setGUIcolor_bymenu(hObj, evt, ~)
% User selected theme color in menu bar. Change settings and run
% setColorgui function to change colors accordingly.

% Get GUI data handle.
gui = guidata(hObj);

% Get event click info
matyr = version('-Release'); matyr = matyr(1:end-1);
switch matyr
    case '2017', evt_click = evt.Source.Text;       
    otherwise,   evt_click = evt.Source.Label;    
end

% Did user change/click to day or night theme?
switch evt_click
    case 'Day'
        gui.menubar.View.theme.day.Checked = 'on';
        gui.menubar.View.theme.night.Checked = 'off';
        gui.menubar.View.theme.custom.Checked = 'off';
    case 'Night'
        gui.menubar.View.theme.day.Checked = 'off';
        gui.menubar.View.theme.night.Checked = 'on';
        gui.menubar.View.theme.custom.Checked = 'off';
    case 'Custom'
        gui.menubar.View.theme.day.Checked = 'off';
        gui.menubar.View.theme.night.Checked = 'off';
        gui.menubar.View.theme.custom.Checked = 'on'; 
end
guidata(hObj, gui); % Update gui data structure

% Apply color scheme
setGUIcolor(hObj);

function setGUIcolor(hObject)

% Change colors of CSIgui between day and night theme.
gui = guidata(hObject);

% Check which theme to be set
day    = gui.menubar.View.theme.day.Checked;
night  = gui.menubar.View.theme.night.Checked;
custom = gui.menubar.View.theme.custom.Checked;
if strcmp(day,'on')
    % GUI color palet: Day theme
    colors.main       = [0.94 0.94 0.94];    % Backgrounds
    colors.text_main  = [0.5 0.5 0.5];       % Used for button text
    colors.text_title = [0 0 0];             % Used for titles
    colors.hilight1   = [0.64 0.15 0.15];    % Highlighted text
    colors.hilight2   = [0.2 0 0.2];         % Highlighted background
    colors.lines1     = [0 0 1];             % Static line color
    colors.lines2     = [1 0 0];             % Second line color
    colors.grid       = [0.4 0.4 0.4];       % Grid color 2D plot
elseif strcmp(night,'on')
    % GUI color palet: Night theme
    colors.main       = [0.1 0.1 0.1];             % Backgrounds
    colors.text_main  = [0.94 0.94 0.94];    % Used for button text
    colors.text_title = [0.502 0.502 0.502]; % Used for descriptions/titles
    colors.hilight1   = [0.8 0 0];             % Highlight text
    colors.hilight2   = [0.2 0 0.2];         % Highlight bg
    colors.lines1     = [1 1 1];             % Static line color
    colors.lines2     = [1 0 0];             % Second line color
    colors.grid       = [0.4 0 0];           % Grid color 2D plot
elseif strcmp(custom, 'on')
   
   % Read colors from file
   colors = setGUIcolor_custom_read;

end

% Special object to apply highlight colors
obj_hilight1 = {'txt_infoCSI',  'txt_infoIMG', 'button_ws'};   % Text hili
obj_hilight2 = {'button_plotCSI',  'button_plotIMG'};          % BG hili

% Get al fieldnames in gui-data struct;
fn = fieldnames(gui);

% Loop each object in gui data structure
for fi = 1:size(fn,1)
    
   
    % Check if its a object handle.
    if sum(ishandle(gui.(fn{fi}))) == 1 
        
        % If its a figure object
        if strcmp((gui.(fn{fi}).Type),'figure')         % Figure
            % Set backgroundcolor of figure
            gui.(fn{fi}).Color = colors.main;
            
        % If its a uicontrol object
        elseif strcmp((gui.(fn{fi}).Type),'uicontrol')  % GUI control
            
            if sum(strcmp(gui.(fn{fi}).Tag, obj_hilight1)) == 1
                % SPECIAL UICONTROL OBJECT 1
                % Highlight foreground
                gui.(fn{fi}).BackgroundColor = colors.main;
                gui.(fn{fi}).ForegroundColor = colors.hilight1;
                
            elseif sum(strcmp(gui.(fn{fi}).Tag, obj_hilight2)) == 1
                % SPECIAL UICONTROL OBJECT 2
                % Highlight background
                gui.(fn{fi}).BackgroundColor = colors.hilight2;
                gui.(fn{fi}).ForegroundColor = colors.text_main;
                
            else
                
                % DEFAULT UICONTROL OBJECT
                if strcmp(gui.(fn{fi}).Style, 'pushbutton') || ...
                    strcmp(gui.(fn{fi}).Style, 'listbox')   || ...
                        strcmp(gui.(fn{fi}).Style, 'popupmenu') 
                    % Buttons, listbox and popupmenu objects
                    gui.(fn{fi}).BackgroundColor = colors.main;
                    gui.(fn{fi}).ForegroundColor = colors.text_main;
                        
                else
                    % Text and other objects
                    gui.(fn{fi}).BackgroundColor = colors.main;
                    gui.(fn{fi}).ForegroundColor = colors.text_title;
                end
            end

        end
    end
end

% Save used colors
gui.colors = colors;

% Update gui handle.
guidata(hObject, gui);

function setGUIcolor_custom(hObj,~,~)
% Allow changing of each color in CSIgui

gui = guidata(hObj);
clrs = setGUIcolor_custom_read();
clrs_fn = fieldnames(clrs);

                % ------- % Create Window % -------- %

% Create figure
fh = figure('Tag', 'CSIgui_Custom_Colors','Name', 'Custom Colors',...
            'Color', gui.colors.main,...
            'Toolbar', 'None', 'MenuBar', 'None','NumberTitle', 'Off');                   

txt_sz = [150 15];        
        
% 1. Default figure size and screen pixel size
def_sz = 240; scr_sz = get(0, 'screensize'); scr_sz(1:2) = [];
% 2. Ratio to define figure height to def_size
fig_sz = [def_sz ((txt_sz(2)*2).*(size(clrs_fn,1)+1))+(txt_sz(2)*2)];
% 4. Position of figure.
scr_middle = scr_sz/2-(fig_sz/2); 
% 5. Apply
set(fh, 'Position', [scr_middle fig_sz]);

% Loop each color: create change color ui line
xpos = (fig_sz(1)-200)/2; bh = cell(1,size(clrs_fn,1));
for kk = 1:size(clrs_fn,1)
    
    % Step in y in window of each line.
    line_step = (txt_sz(2)+(kk*txt_sz(2)))+((kk-1)*txt_sz(2));
    
    % 1. Title
    uicontrol('Style','text','String',['Color for "' clrs_fn{kk} '":'], ...
        'BackgroundColor',gui.colors.main,...
        'ForegroundColor',gui.colors.text_main,...
        'HorizontalAlignment','Left',...
        'Position',[xpos fig_sz(2)-line_step txt_sz(1) txt_sz(2)]);
    
    % 2. Color "block" button
    bh{kk} = uicontrol('Style','pushbutton','String','',...
        'Tag',clrs_fn{kk},...
        'BackgroundColor',clrs.(clrs_fn{kk}),...
        'ForegroundColor',clrs.text_main,...
        'Position',[xpos+txt_sz(1) fig_sz(2)-line_step 50 txt_sz(2)],...
        'Callback', @setGUIcolor_custom_pickColor);

end

% Done button
uicontrol('Style','pushbutton','String','Done',...
    'BackgroundColor',gui.colors.main,...
    'ForegroundColor',gui.colors.text_main,...
    'Position',[xpos txt_sz(2) 50 txt_sz(2)],...
    'Callback', @setGuiColor_custom_Done);

% Store button object in figure gui data
fig_gui.bh = bh; guidata(fh,fig_gui);


function colors = setGUIcolor_custom_read()

% Load custom color set from file
fid = fopen('custom_colors.txt');
nn = 1;
while ~feof(fid), lines{nn} = fgetl(fid) ; nn = nn+1; end
fclose(fid);

for kk = 1:size(lines,2)
    if strcmp(lines{kk},'\\'), break; end 
    tmp = strsplit(lines{kk}, '=');
    fn = strip(tmp{1}, ' '); val = str2double(strsplit(tmp{2},' '));
    colors.(fn) = val;
end

function setGUIcolor_custom_pickColor(hObj,evt,~)

% Guidata of color pick window
gui = guidata(hObj);

% Source of pushed button
src = evt.Source.Tag;    

% Create a dialog to pick colors
clr = uisetcolor('Select color');

% Find click button index and set picked color.
button_tags = extractField(gui.bh,'Tag');
ind = strcmp(button_tags, src);
gui.bh{ind}.BackgroundColor = clr;

function setGuiColor_custom_Done(hObj,~, ~)
% Read the colors of all color-buttons in the Custom Color window and saves
% it to the custom_color text file.

% Gui data in colorpick window
gui = guidata(hObj);

             % -------- %  Extract new colors % -------- %

% Get each button color
new_clrs = extractField(gui.bh,'BackgroundColor');
tags = extractField(gui.bh,'Tag');

              % -------- %  Store new colors % -------- %
          
% Read custom colors file              
fid = fopen('custom_colors.txt'); kk = 1;
while ~feof(fid), lines{kk} = fgetl(fid); kk = kk +1; end
fclose(fid);

% Replace colors in correct line
for kk =1:size(tags,2)
    tmp_color = ...
        strjoin(cellfun(@num2str,num2cell(new_clrs{kk}),'Uniform',0),' ');
    lines{kk} = [tags{kk} ' =' tmp_color];
end

% Store new data to file
fid = fopen('custom_colors.txt','wt');
for kk = 1:size(lines,2)
    fprintf(fid, '%s\n', lines{kk});
end
fclose(fid);

% Close color pick window
delete(gui.bh{1}.Parent);

csigui_obj = findobj('Tag','CSIgui_main');
gui = guidata(csigui_obj);
gui.menubar.View.theme.day.Checked    = 'Off';
gui.menubar.View.theme.night.Checked  = 'Off';
gui.menubar.View.theme.custom.Checked = 'On';
setGUIcolor(csigui_obj);



%%% MERGE VOXELS % ------------------------------------------------------ %

% --- Executed if user presses merge voxels button
function CSI_MergeVoxels_Initiate(hObj, gui)
% Get the data and create plot_par structure with all settings and data for
% plotting the MergeVoxels figure: image + grid + click-axis only.

% Prepare data % ------------------- %

% Get CSI data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get colors
plot_par.colors = gui.colors;

% Get x-axis
plot_par.xaxis = csi.xaxis;

% Data: dimensions % ------------------- %

% Axis size to fit figure
plot_par.dim      = size(csi.data.raw);   % Data dimensions
plot_par.dim(1)   = [];                   % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 

% Correction 1D % ------------------- %
% Set in plot_par the 1D dimension to the second index e.g. column.
if plot_par.data_dim == 1
    % Add a Y dimension.
    plot_par.dim(2) = 1; plot_par.data_dim = 2;  
end

% ------- % Create Figure for plotting 2D
plot_par = CSI_2D_setFigure(plot_par, 0, 'CSIgui - Merge Voxels');

% ------- % Get data for plotting
plot_par = CSI_2D_getData(plot_par, gui,  csi.data.raw);

% ------- % Get and set all other plot options
plot_par = CSI_2D_getPlotSettings(plot_par, gui, csi.data.raw);

% ------- % Plot data
CSI_MergeVoxels_plotVoxels(hObj, gui, plot_par);

% --- Executes on button press in button_CSI_MergeVoxels.
function button_CSI_MergeVoxels_Callback(hObj, ~, gui)
% Initiate MergeVoxels figure to allow selecting multiple voxels to average
% using specific options.

% Initiate a figure for CSI_mergeVoxels
CSI_MergeVoxels_Initiate(hObj, gui);

% --- Executes by CSI_MergeVoxels_Initiate
function CSI_MergeVoxels_plotVoxels(hObj, gui, plot_par)
% Plot voxels, axis only, each slice in a tab.       

% Add tab-group
plot_par.tabg = uitabgroup(plot_par.fh);

% --------- % Prepare images % --------- %
              
img = MRI_matchSlices(hObj);
if sum(size(img)) == 1, img = []; end

              % --------- % Loop all slices % --------- %

save_data = struct;
for sli = 1:plot_par.dim(3) 
    
    loadBar(sli./plot_par.dim(3), 'Creating figure..');
    
    % Create a tab for this slice in in tabgroup   
    save_data.tabh{sli} = uitab(plot_par.tabg, 'Title',int2str(sli),...
        'BackgroundColor', plot_par.colors.main,...
        'ForegroundColor',gui.colors.text_title);

    % Images: Plot % ------- %

    if isappdata(gui.CSIgui_main, 'conv') && ~isempty(img)

        conv_data = getappdata(gui.CSIgui_main,'conv'); % Conv data
        img2plot = img(:,:,sli);                        % Get image to plot.
        
        save([datestr(now,'HHMMSSFFF') '.mat' ],'img');
        
        % Create axis for image
        hold on; 
        imax = axes('parent',save_data.tabh{sli},...
            'Position',[0 0 1 1], 'Color', 'None');

        % Plot Images
        if (sum(img2plot(:)) == 0) % Image is only zeroes.
            colormap(gray(255)); set(imax,'Color', 'Black'); alpha(1); 
        else
            % Image plotting:
            % Imscale as it plots over the entire figure and does not
            % imply any border issues as with imshow-function.
            imagesc(img2plot, 'parent', imax); 

            % Image Contrast.
            if isfield(conv_data, 'contrast')
                caxis(imax, conv_data.contrast);
            else
                caxis(imax,[min(img2plot(:)) max(img2plot(:))*0.5]);
            end
            colormap(gray(255));
        end 

        plot_par.colors.main = [0 0 0];

    end


    % Loop: each axis/voxel % -------- %
    axlw_pt = 1;
    for ci = 1:plot_par.dim(1)                  % Col loop.
        for ri = 1:plot_par.dim(2)              % Row loop.

            % Axis: create % -------- %

            % X and Y position in the figure of the axis to plot.
            x   = plot_par.grid.x(ri,ci); y = plot_par.grid.y(ri,ci);
            % Position of axis to plot
            pos = [x y plot_par.res(1) plot_par.res(2)];
            % Create axis with pos(3,4) size at pos(1,2) position
            if ishandle(plot_par.fh)
                save_data.ax{ri,ci,sli} =...
                    axes('parent',save_data.tabh{sli},'position',pos);
            else, return; 
            end
            
            set(save_data.ax{ri,ci,sli},...
               'Color', 'None',...
               'XColor', plot_par.colors.main,'YColor', plot_par.colors.main,...
               'LineWidth', axlw_pt, 'Xtick',[], 'Ytick',[],...
               'TickLength',[0 0.00001],'Box', 'off');     

            % Set action when clicked on this axis (ci, ri);
            set(get(save_data.ax{ri,ci,sli},'Children'),'HitTest', 'Off');
            set(save_data.ax{ri,ci,sli},...
                'ButtonDownFcn',@CSI_2D_voxel_selectMultiple);

            save_data.ax{ri,ci,sli}.UserData = [ri ci sli];
            
        end
    end


    % Grid: plot % -------- %

    CSI_2D_grid(save_data.tabh{sli},plot_par.fh.Position(3:4),...
        plot_par.dim, plot_par.range, plot_par.colors.grid);


    % Buttons: create % --------- %

    % Align selected voxels
    uicontrol('Parent',save_data.tabh{sli},'Style','pushbutton','String', 'Merge',...
              'BackgroundColor',plot_par.colors.main,...
              'ForegroundColor',gui.colors.text_main,...
              'Position',[5 5 50 15],...
              'Callback',@CSI_MergeVoxels_Button_Merge);                  
    % Save selected voxel data or figure
    uicontrol('Parent',save_data.tabh{sli},'Style','pushbutton','String', 'Save',...
              'BackgroundColor',plot_par.colors.main,...
              'ForegroundColor',gui.colors.text_main,...
              'Position',[5 20 50 15],...
              'Callback',@CSI_MergeVoxels_Button_Save);

end % ------ % End of Slice Loop % ------ %

save_data.highlight = gui.colors.hilight1;
% Store data to figure.
guidata(plot_par.fh, save_data);

loadBar(NaN);

% --- Executed by functions related to CSI_MergeVoxels
function selected = CSI_MergeVoxels_getSelected(hObj)
% Find all selected voxels in MergeVoxels figure and return data of each
% selected voxel.
%
% Output: selected-struct
%                       .data  = selected voxels
%                       .index = index of selected voxels in csi.data.raw;

% ---- % Fig Obj and gui-data
fh = hObj.Parent; gui = guidata(fh);

% ---- % Find clicked voxels
% Color of axis is not "none".
clr = extractField(gui.ax,'Color');
ind = cellfun(@strcmp,clr,repmat({'none'},size(clr)));
ii = find(ind == 0); [ri, ci, si] = ind2sub(size(ind),ii); % Swap x/y/r/c

% ---- % Get CSI data
csigui_obj = findobj('Tag','CSIgui_main');
csi = getappdata(csigui_obj,'csi'); 
if isempty(csi), CSI_Log({'No CSI data loaded!'},{''}); return; end

% ---- % Loop each selected voxel and get data
selected.voxels  = NaN(size(csi.data.raw,1), size(ri,2));
selected.index = NaN(size(ri,2),3);
for vx = 1:size(ri,1)
    selected.voxels(:,vx)  = csi.data.raw(:,ci(vx),ri(vx),si(vx));
    selected.index(vx,:) = [ci(vx) ri(vx) si(vx)];
end


% --- % Add additional data

% Add frequency parameters.
selected.xaxis = csi.xaxis;

% Add Log
gui = guidata(csigui_obj); selected.log = gui.listbox_CSIinfo.String;

% --- % Executes on button press: Align @ CSI_MergeVoxel figure
function CSI_MergeVoxels_Button_Merge(hObj,~)

% --- % Get the selected voxels data
sel = CSI_MergeVoxels_getSelected(hObj);

% --- % Userinput

quest = {'SNR Filtering','Aligning','Weighting'};
uans = getUserInput_Radio(quest,ones(1,size(quest,2)));
if isempty(uans), CSI_Log({'Skipped merging voxels.'},{''}); return; 
end  


% --- % Apply options
if sum(uans) > 0 
    
    % --- % Get POI
    [sel.poi, sel.pex] = CSI_MergeVoxels_POI(sel.xaxis);
    
    % --- % SNR Filtering
    if uans(1), sel = CSI_MergeVoxels_SNR_Filter(sel);  end

    % --- % Aligning
    if uans(2), sel = CSI_MergeVoxels_Align(sel); end

    % --- % Weighting
    if uans(3), sel = CSI_MergeVoxels_Weighted(sel); end

else
    % No Options
end

% --- % Average & display
CSI_MergeVoxels_Average(sel);

% --- % Executes on button press: Save @ CSI_MergeVoxel figure
function CSI_MergeVoxels_Button_Save(hObj,~)

% ---- % File destination and type: figure or data
[fn,fp,idx] = uiputfile(...
            {'*.mat',       'MATLAB File (*.mat)';...
             '*.fig',       'MATLAB Figure (*.fig)';...
             '*.mat;*.fig', 'MATLAB File (*.mat) & Figure (*.fig)'},...
             'Save File');
if isequal(fn,0) && isequal(fp,0), return; end

% Split name, path and extension
[fp,fn] = fileparts([fp fn]);

if     idx ==1 % Save Data
    CSI_MergeVoxels_SaveData(hObj, fp, fn); 
elseif idx ==2 % Save Figure
    CSI_MergeVoxels_SaveFig(hObj,  fp, fn); 
elseif idx ==3 % Save both
    CSI_MergeVoxels_SaveData(hObj, fp, fn);
    CSI_MergeVoxels_SaveFig(hObj,  fp, fn); 
end

% --- % Executed by CSI_MergeVoxels_Button_Save
function CSI_MergeVoxels_SaveData(hObj, fp, fn)
% Get the selected voxel data and save to file.

% --- % Get selected voxels.
selected = CSI_MergeVoxels_getSelected(hObj);

% --- % Store data
save([fp '\' fn '.mat'], 'selected');

% --- % LOG
CSI_Log({'Saved selected voxel data to:'},{[fp '\' fn '.mat']})

% --- % Executed by CSI_MergeVoxels_Button_Save
function CSI_MergeVoxels_SaveFig(hObj, fp, fn)

% ---- % Figure object
fh = hObj.Parent; 

% ---- % Save figure

matver = version('-release'); matyr = str2double(matver(1:end-1));
% Compact for MATLAB 2015 and later versions only.
if (matyr <= 2014), savefig(fh, [fp '\' fn '.fig']);
else,               savefig(fh, [fp '\' fn '.fig'], 'compact');
end

% LOG
CSI_Log({'Saved selected voxel figure to:'},{[fp '\' fn '.fig']});

% --- % Executes is user requests merge with options none
function CSI_MergeVoxels_Average(sel)
% Requires data structure data with fields:
% xaxis, data.
%
% Returns average of the given voxels in CSIgui-1D plot

% Average all voxels
sel.avg = mean(sel.voxels,2);

% Display the merged voxel: CSIgui-1D
data1D.voxel.original = sel.avg; data1D.data = sel.avg;  
data1D.axis = sel.xaxis; data1D.voxel.index = [0 0 0];
data1D.voxel.unit = 'Real';
if isfield(sel.xaxis,'ppm'), data1D.axis.unit = 'ppm';
else,                        data1D.axis.unit = 'none';
end

CSI_1D_initiateGUI(data1D, 1, sprintf('#Voxels %i',size(sel.index,1)));

% --- % Peak of interest for merging and filter criteria
function [poi, pex] = CSI_MergeVoxels_POI(xaxis)
% Get POI and possibly PEX for merging voxels
%
% Returns: sel.poi and sel.pex 
%          Ranges in the voxel data: index.
%
% If no peak of exclusion is used, pex is set to NaN;

% ---- %  Userinput
% Get the index of the peak of interest and peak of exclusion.

tags = {'for alignment', 'for exclusion'}; range = cell(1,2);

for kk = 1:2, range{kk} = CSI_getPeakOfInterest(xaxis,tags{kk}); end

if isempty(range{1}), return; end % User presses skip.
if isempty(range{2}), range{2} = NaN; end

% ---- % Set output
poi = range{1}; pex = range{2};

% --- % Executed to filter merging voxels by SNR
function sel = CSI_MergeVoxels_SNR_Filter(sel)
% Include and/or exclude voxels using SNR.

voxels = sel.voxels; pex_on = ~isnan(sel.pex); index = sel.index;

% ---- % Userinput
% SNR Limits

% SNR Limits to include or exclude
tags = {'Alignment: Exclude SNR < Limit','Exclusion: exclude SNR > Limit'}; 
defans = [1,12]; uans = cell(1,2);
if pex_on, lims = 2; else, lims = 1; end
for kk = 1:lims
    uans{kk} = getUserInput({['Peak SNR threshold for ' tags{kk}]},...
                            {defans(kk)});
end
if isempty(uans{1})
        CSI_Log({'Skipped SMR filtering merging voxels.'},{''}); 
        return; 
end  

% SNR limit for POI
poi_snr_lim = str2double(uans{1});                  
% SNR limit for PEX
if pex_on, pex_snr_lim = str2double(uans{2}); else, pex_snr_lim = 0; end

% ---- % Phasing and SNR

% SNR of POI 
poi_snr = CSI_MergeVoxels_SNR(sel, 'alignment');         

% SNR of PEX
if pex_on
    seltmp = sel; seltmp.poi = sel.pex;
    pex_snr = CSI_MergeVoxels_SNR(seltmp, 'exlusion'); 
end

% ---- % Select Voxels by SNR

% For POI
poi_incl = (poi_snr >= poi_snr_lim); poi_excl = (poi_snr <  poi_snr_lim);

% For PEX
if pex_on, pex_excl = (pex_snr >  pex_snr_lim); else, pex_excl = 0; end

% Merge SNR incl/excl from POI and PEX
vox_incl = poi_incl; vox_excl = poi_excl;
if pex_on, vox_incl(pex_excl == 1) = 0; vox_excl(pex_excl == 1) = 1; end

% Voxels remaining check.
if sum(vox_incl) == 0
    CSI_Log({'Voxel selection criteria to strict.'},...
            {'No Voxels to average.'}); return;
end

% Get all voxels above POI snr limit and below PEX snr limit
voi     = voxels(:,logical(vox_incl));
voi_snr = poi_snr(vox_incl); voi_ind = index(vox_incl,:);


CSI_Log({'-----------------------','Total #Voxels for merging:',...
         'Used SNR Limits: ',...
         '#Voxels after filtering:',...
         '#Excluded by SNR of POI/PEX:','-----------------------'},...
         {'',size(voxels,2),...
         [num2str(poi_snr_lim) ' | ' num2str(pex_snr_lim)],...
         sum(vox_incl), ...
         [num2str(sum(poi_incl)) ' | ' num2str(sum(pex_excl))],''});

% ---- % Set Output
sel.voxels = voi; sel.snr = voi_snr; sel.index = voi_ind;

% --- % Executed to calculate SNR and apply phasing
function snr = CSI_MergeVoxels_SNR(sel, poi_tag)
% Calculate SNR and apply phasing prior to calculations if necessary.
%
% poi_tag = 'Name for peak'

if nargin == 1, poi_tag = ''; end

voxels = sel.voxels; % Get voxel data.
poi = sel.poi;       % Peak of interest index

% ---- % Userinput
% Phasing

uans = getUserInput_Popup({['Apply phasing to ' poi_tag ' peak?']},...
                          {{'No','Yes'}});
if isempty(uans{1})
        CSI_Log({'Skipped SNR filtering merging voxels.'},{''}); 
        return; 
end  
if strcmp(uans{1},'No'), poi_phase = 0; else, poi_phase = 1; end
   

% ---- % Phasing

if poi_phase 
    
    % Voxel data to cell format
    sz = size(voxels); 
    cell_layout = arrayfun(@ones,...
        ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
    voxels_cell = mat2cell(voxels, sz(1), cell_layout{:});                    
    
    % Apply zeroth order phase correction focused on POI
    poi_vox_4snr = cellfun(@csi_autoZeroPhase,...
        voxels_cell, repmat({poi}, size(voxels_cell)),... % Peak range
                     repmat({2},   size(voxels_cell)),... % Method
                     repmat({0},   size(voxels_cell)),... % Display
                     'UniformOutput', 0);
    poi_vox_4snr = cell2mat(poi_vox_4snr);
else
    poi_vox_4snr = voxels;
end      

% ---- % SNR
snr = csi_SNR(poi_vox_4snr, 50, 1, sel.poi); % Calc SNR for POI   

% --- % Executed to weight voxels to SNR
function sel = CSI_MergeVoxels_Weighted(sel)
% Apply weighting to each voxel before averaging
% Weighting is calculated from maximum SNR at a POI from all voxels.
 
sel.snr = csi_SNR(sel.voxels, 50, 1, sel.poi); % Calculate SNR

% ---- % Weighting

% Weighting values by SNR
snr_max = max(sel.snr); sel.weights = sel.snr ./ snr_max;

% Apply weighting
for kk = 1:size(sel.voxels,2)
    sel.voxels(:,kk)  = sel.voxels(:,kk) .* sel.weights(kk);
end

% --- % Executed by CSI_MergeVoxels_Button_Align
function sel = CSI_MergeVoxels_Align(sel)
% Align voxels before averaging.
%
% Sel-structure requires field voxels and peak of interest poi.

% Undress the input structure:
voxels = sel.voxels; poi = sel.poi; 

                   % ---- % Aligning voxels % ---- %
                    
% Find the index of each voxels POI maximum
[~ , voi_poi_max_pos] = max( real( voxels(poi(1):poi(2),:)), [], 1);
% Correct position to full sample size (e.g. not in poi_ind)
voi_max_pos = voi_poi_max_pos + poi(1)-1;

% Shift to center of given range.
voi_shift_to = poi(1)+round((diff(poi)/2)); 
% Number of samples to shift each voxel
voi_shifts = voi_shift_to - voi_max_pos;

% Shift each voxel by its voi_shift amount
voi_cell = mat2cell(voxels, size(voxels,1), ones(1,size(voxels,2)));
voi_aligned_cell = ...
    cellfun(@circshift, voi_cell, num2cell(voi_shifts), 'uniform',0);
voi_aligned = cell2mat(voi_aligned_cell);

% ---- % Set Output
sel.voxels = voi_aligned;


% Additional plots

nsli = unique(sel.index(:,3));
nfig = numel(nsli);

for kk = 1:nfig
    
    % N voxels for this figure.
    nvox_perfig{kk} = find(sel.index(:,3) == nsli(kk));
    vox_ind_all = nvox_perfig{kk};
    
    fh = figure();
    nvox = size(nvox_perfig{kk},1);
    for vi = 1:nvox
        subplot(nvox,1,vi);
        
        vox_ind = vox_ind_all(vi);
        
        if isfield(sel.xaxis,'ppm'), xax = sel.xaxis.ppm;
        else,                        xax = sel.xaxis.none;
        end

        % Plot aligned
        plA = plot(xax, real(voi_aligned(:,vox_ind)),'-r');hold on;
        % Plot voxel
        plV = plot(xax, real(voxels(:,vox_ind)),'--b'); 

        
        % Axis cosmetics
        

        axh = plV.Parent;
        plV.Parent.XDir = 'reverse';
        
        axh.XLim = [-16 8];
    end
    
    str = ['Slice ' num2str(nsli(kk)) ...
           ' | Avg shift: ' num2str(mean(voi_shifts(vox_ind_all)))];
       
    uicontrol(fh, 'Style','text','Units','Normalized',...
        'Position',[0.35 0.925 0.3 0.05],...
        'string', str);
        % 'ForegroundColor','White','BackgroundColor','Black',...

    
end


% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %

                            % MISCELLANEOUS %

% --- Load bar
function loadBar(perc, info)
% Opens a load bar if not existant and sets the load-bar percentage to
% perc. If perc equals NaN, the bar is closed.
%
% Input info, a string, is optional and can be used to set text in the
% window of the loadbar. 
% 
% Mainly used for launching of CSIgui

% If no info string is given
if nargin == 1, info = 'Busy...'; end

% Normalize percentage
if ~isnan(perc) && perc > 1, perc = perc/100; end

% Check for bar-figure
barObj = findobj('Type','Figure','Tag','loadBar');

                % ------- % Create loadBar % ------- %
if isempty(barObj)
    % Create figure
    figh = figure('Tag', 'loadBar','Name', ['CSIgui: ' info],...
           'Color', 'Black','Toolbar', 'None', 'MenuBar', 'None',...
           'NumberTitle', 'Off','Resize','off');  
      
    % Set figure position
    w = 480; h = 20;
    scrsz = get(0,'Screensize'); 
    figpos = round(scrsz(3:4)./2) - ([w h]./2);
    set(figh,'Position', [figpos w h]);
    
    % GUI data
    bgui = guidata(figh); bgui.fig = figh;
    
    % Set a (text) bar
    bgui.bar = uicontrol('style','text', ...
                'units','normalized','position',[0 0 0 1],...
                'BackgroundColor', [0.9 0.0 0.0]);
            
    % Save guidata
    guidata(bgui.fig, bgui);
else
    % Get loadbar guidata
    bgui = guidata(barObj);
end

                % ------- % UPDATING loadBAR % ------- %

% Close if NaN
if isnan(perc), delete(bgui.fig); return; end

% Set bar progress
bgui.bar.Position = [0 0 perc 1];

% Set info string
bgui.fig.Name = ['CSIgui: ' info];

% Flush
drawnow;                            
                            
                            
% --- Executes on selection change in listbox_CSIinfo.
function listbox_CSIinfo_Callback(hObject, eventdata, gui)
% --- Executes during object creation, after setting all properties.
function listbox_CSIinfo_CreateFcn(hObject, eventdata, gui)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in listbox_MRIinfo.
function listbox_MRIinfo_Callback(hObject, eventdata, gui)
% --- Executes during object creation, after setting all properties.
function listbox_MRIinfo_CreateFcn(hObject, eventdata, gui)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popup_plotUnit.
function popup_plotUnit_Callback(hObject, eventdata, gui)
% --- Executes during object creation, after setting all properties.
function popup_plotUnit_CreateFcn(hObject, eventdata, gui)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popup_plotIMG.
function popup_plotIMG_Callback(hObject, eventdata, gui)
% --- Executes during object creation, after setting all properties.
function popup_plotIMG_CreateFcn(hObject, eventdata, gui)
% hObject    handle to popup_plotIMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------- %
            %   Everything below here is new OR beta!  %
                        % Probably not finished %
% --------------------------------------------------------------------- %


%%% ----------- % ------------------------------------------------------ %

% --- Executes on button press in button_TestSomething.
function button_TestSomething_Callback(hObj, eventdata, gui)
% warndlg('Watch out! Developers testing button. Panic!');

% Get csi data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get image data
if ~isappdata(gui.CSIgui_main, 'conv'), return; end
conv = getappdata(gui.CSIgui_main, 'conv');

mri = getappdata(gui.CSIgui_main,'mri');

% Tra to Sagital % ---------------------- %
% A (top) L (right) to R (top) A(right)
% Current: kx rl = x = 2;  To  rl to z = 4;
%          ky ap = y = 3;      ap to x = 2; 
%          kz fh = z = 4;      fh to y = 3;


plane = [3 4]; % Plane for MRS data e.g. time/x/y/z...
% For no-time index data use plane-1;



% ROTATE MRS
[csi.data.raw, permv] = CSI_rotate(csi.data.raw, plane, 1);



% ROTATE MRS RELATED
% Apply rotation to labels
csi.data.labels(plane) =  csi.data.labels(fliplr(plane));
% And dimensions
csi.data.dim(plane) = csi.data.dim(fliplr(plane));



% ROTATE IMG
conv.data = CSI_rotate(conv.data, plane-1, 1);



% ROTATE SPATIAL INFORMATION  MRS

% Orientation info.
ori = csi.ori;
ori.res(plane-1)       = ori.res(fliplr(plane-1));       
ori.offcenter(plane-1) = ori.offcenter(fliplr(plane-1));  
ori.dim(plane-1)       = ori.dim(fliplr(plane-1));       

% Coordinates of CSI data.
% Adds fields: coordinate vector (vector) & coordinate limits (limit) &
% volume limit (limit_vol).
csi.ori = csi_coordinates(ori,'center', ori.vox_cor, ori.fft_cor); 

% Meshgrid
[csi.ori.mesh.x, csi.ori.mesh.y, csi.ori.mesh.z] = ...
    meshgrid(csi.ori.vector{1} , csi.ori.vector{2}, csi.ori.vector{3});


% ROTATE SPATIAL INFORMATION CONVERTED IMG

% Swap resolutions.
% conv.res = mri.ori.res;
conv.res(plane-1) = conv.res(fliplr(plane-1)); % Initial... May change!

% % Calculate a resolution fitting the CSI space such that there is a integer
% % amount of image pixels fitted in each CSI direction of space.
% res_fit = csi.ori.res ./ conv.res;      % #MRpix / CSIpix
% res_rem = res_fit - floor(res_fit);     % Pixel change
% res_new = csi.ori.res ./ floor(res_fit);% New MRpix resolution 
% 
% % New resolution for each direction
% conv.res = res_new;

% Volume limits of CSI but with half a voxel distance for voxel limits e.g.
% a total voxel (MRI res) vs the volume.
conv.fov       = csi.ori.fov;     % Does not change regards to CSI
conv.lim_vol   = csi.ori.lim_vol; % Volume of MRSI grid
conv.lim(:,1)  = conv.lim_vol(:,1) + (0.5.*conv.res)'; % Voxel limits
conv.lim(:,2)  = conv.lim_vol(:,2) - (0.5.*conv.res)'; % Used for coords

% Range of volume/grid of MRSI for MRI
for kk = 1:size(conv.lim,1)    
    % Number of pixels in kk-direction
    N = (conv.fov(kk)./conv.res(kk));
    % Grid vector defining volume
    conv.vec{kk} = linspace(conv.lim(kk,1), conv.lim(kk,2), N);
end

% Image grid in MRSI space 
% This grid lays in the CSI-space e.g. within limits of CSI FOV but is
% sampled in x, y and z as close to the resolution of the image as possible
[x,y,z] = meshgrid(conv.vec{1} ,conv.vec{2}, conv.vec{3});
conv.mesh.x = x; conv.mesh.y = y; conv.mesh.z = z; 

% Store data
setappdata(gui.CSIgui_main, 'csi',csi);
setappdata(gui.CSIgui_main, 'conv',conv);
CSI_Log({'Done testing.'},{''});


