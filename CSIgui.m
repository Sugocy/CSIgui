function varargout = CSIgui(varargin)
% Spectroscopy GUI; initial purpose to merge MRI and MRS data in Matlab.
% See Zonedo.org or github for latest releases. Help file is present in the
% "files" directory or via help in the GUI.
%
% Possible labels for input arguments:
%            'data','list','csi','spec','image', 'mrs', 'labels'
%            {filepath}, {filepathi}
% Input:
% CSIgui(datafield, label);
%
% ------------------------------------------------------------------------
% 2023/08 - CSIgui v2.2 ELH.
% Tools for quantitative MR imaging and spectroscopy for the improvement of
% therapy evaluation in oncology. doi(.org/) 10.33540/52
%
% Quincy van Houtum, PhD. quincyvanhoutum@gmail.com

% Last Modified by GUIDE v2.5 02-May-2024 14:31:12

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
id = strrep(strjoin(string(datevec(datetime("now"))),''),'.','');
gui.ID = sprintf('%s', id);
% Define CSIgui version here.
gui.version = '2.0';


%% Add menu-bar

loadBar(0.1, 'Creating menubar...');
gui = CSIgui_menubar_Set(gui);

loadBar(0.3, 'Executing JavaScript...');
gui = CSIgui_menubar_javaScript(gui);

% Update gui structure
guidata(hObject, gui);

%% Add Color-scheme to GUI
% Run GUI coloring theme function
loadBar(0.4, 'Creating GUI...');
setGUIcolor(hObject); % This stores the right color palet in gui-handle
gui = guidata(hObject);


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
                            'mrs', 'filepath', 'filepathi', 'fp','fpi',...
                            'twix'};
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


% Close menubar
loadBar(1, 'Launching GUI...'); pause(0.0005); loadBar(NaN);

% --- Executes at launch
function gui = CSIgui_menubar_Set(gui)
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
    'Accelerator', 'z', 'Callback', @CSIgui_menubar_backupUndo);

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

% 2. MRSI > Color Scaling (2D) > Scale by window
gui.menubar.MRSI.ColorScale.ScalebyWindow = ...
    uimenu(gui.menubar.MRSI.ColorScale.main, 'Label', 'Scale by Window',...
    'Check', 'Off','Separator','on',...
    'Callback',@CSI_2D_Scaling_Color_ScaleByWindow, 'Enable', 'On');

% 2. MRSI > Color Scaling (2D) > Colorbar
gui.menubar.MRSI.ColorScale.Colorbar = ...
    uimenu(gui.menubar.MRSI.ColorScale.main, 'Label', 'Colorbar',...
    'Check', 'Off','Separator','off',...
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
        
% 2.MRSI > Axis Scaling (2D) >  Scale by window
gui.menubar.MRSI.AxisScale.ScalebyWindow = ...
    uimenu(gui.menubar.MRSI.AxisScale.main, 'Label', 'Scale by Window',...
    'Check', 'Off','Separator','on',...
    'Callback',@CSI_2D_Scaling_Axis_ScaleByWindow, 'Enable', 'On');

% 2.MRSI > Axis Scaling (2D) >  Grid
gui.menubar.MRSI.AxisScale.Grid = ...
    uimenu(gui.menubar.MRSI.AxisScale.main, 'Label', 'Grid',...
    'Check', 'Off','Separator','off',...
    'Callback',@CSI_2D_Scaling_GridVoxels, 'Enable', 'On');



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
    'Callback',@CSI_Noise_ViewManager);        


% 2. MRSI > SliderfollowPlot
gui.menubar.MRSI.snapWindow.main = ...
    uimenu(gui.menubar.MRSI.main,'Label', 'Snap Window','Enable', 'on',...
    'Check', 'off', 'Callback', @panel_2D_followPlot2D_menubar );        


% 2.MRSI > Save > ...
gui.menubar.MRSI.export.main = ...
    uimenu(gui.menubar.MRSI.main,'Separator','on','Label', 'Export',...
    'Enable', 'on');        

% 2.MRSI > Save > DATA

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
    'Callback', @button_MRI_PlotIMG_Callback);

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
    
gui.menubar.About.github = ...
    uimenu(gui.menubar.About.main,'Label', 'Github', 'Separator','off',...
        'Callback', @openGithub, 'Enable', 'On');
    

end

% --- Executes at launch
function gui = CSIgui_menubar_javaScript(gui)
% Sets tooltips to menu-bar entries.

% ----------- %%% Some Jave Menu ToolTip String Magic %%% -------------- %
try
    % Flush java.
    drawnow; pause(0.015); 
    % Surpress obsolete function message
    warning('off', 'all');
    
    % Get figure java frame
    jFrame = get(handle(gui.CSIgui_main),'JavaFrame');
    % Get menubar component
    jMenuBar = jFrame.fHG2Client.getMenuBar; 
    
    % Get CSI menubar; java starts at zero, 3-1 = 2;
    jMRSI = jMenuBar.getComponent(1); % CSI menu 
    
    % Check if right component
    if strcmp(jMRSI.getText,'MRSI')
        % Click once to open the menu 
        jMRSI.doClick;  pause(0.015); 
        
        menu_name = {'Color Scaling', 'Axis Scaling','Noise',...
            'Set Domain', 'Snap Window'};
        menu_tip = {...
     ['Set spectra color scaling in 2D MRSI plot relative to,',... 
      'maximum per slice, volume or as static (e.g. single color), '...
      'from the total spectrum or within the displayed x-axis window. '],...
      ['Scale the spectra y-axis in 2D MRSI plot to the maximum '...
      'per voxel, slice or volume, '...
      'from the total spectrum or within the displayed x-axis window. '...
      'Enable/Disable the grid within a voxel in the 2D MRSI plot.'],...
     ['View noise data from list/data file. Use the show '...
      'CSI-button or plot in menu > CSI to display. ' ...
      'To revert back to the original data, click here again.'], ...
      'Set data domain of MRSI data to frequency or time.',...
      'DEV: Enable/disable slider-window to snap to the CSI-2D plot-window.'};
        
        nMenus = size(fieldnames(gui.menubar.MRSI),1)-1; % Remove main CSI
        for kk = 0:nMenus-1 % Start indexing @ zero - Matlab vs Java.
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
    
    warning('on', 'all'); % Enable warnings again
catch err
    fprintf('Jave Menubar setToolTipText error: %s', err.message);
    save(join([string(datetime('now','format', 'HHmmSS')) '_JavaError.mat'],''),  'err');
end

% --- Outputs from this function are returned to the command line.
function varargout = CSIgui_OutputFcn(~, ~, gui)
% varargout  cell array for returning output args (see VARARGOUT);
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


% Open and Parse files % ---------------------------------------------- %
% --------------------------------------------------------------------- %

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
                'Raw MRS files Philips (*.list, *.data)';...
     '*.dcm;*.DCM;*.par;*.rec;*.PAR;*.REC',...
                'Image files (*.dcm, *.par, *.rec)';...
     '*.spar;*.SPAR;*.sdat;*.SDAT',...
                'MRS files Philips (*.sdat, *.spar)';...
     '*.dat;*.DAT;',...
                'Raw MRS file Siemens (*.dat)';...
     '*.txt;*.TXT;',...
                'Text files (*.txt)'},...
     'Select a file', fp, 'MultiSelect', 'on'); 
    % Canceled file selection
    if fi == 0, return; end 
    

    % MULTIPLE FILES
    if iscell(fn)

        [~,fn,ext] = cellfun(@(x,y) fileparts([x y]),...
            repmat({fp},1,size(fn,2)),fn,'UniformOutput', 0);
        ext = cellfun(@lower,ext, 'UniformOutput', 0);
        if size(unique(ext),2) > 1
            % Different extensions selected, unsupported for now.
            CSI_Log({'Multiple extensions detected, aborting.'},{''});
        else
            ext_tmp = ext{1}; ext = 'multi';            
        end

    % SINGLE FILE
    else
        % Get file-extension for further processing, see below.
        [fp, fn, ext] = fileparts([fp, fn]); ext = lower(ext);
    end

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
    set(gui.listbox_CSIinfo, 'String', {''}, 'value', 1); % Clear info in listbox.    
end

% Update LOG
CSI_Log({['Loading ' strrep(ext,'.','') ' file.']}, {'Please wait.'});


% Parse data per extension.
if strcmp(ext,'.dcm') || strcmp(ext,'.ima')                    % DICOM

    % Parse dicom file
    success = parse_dicom(fp, [fn ext], gui);
    
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
    % Data domain set to 2 = frequency
    domain = 3; 
    
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
    
    % Data domain set to 1 = none 2 = frequency 3 = time
    domain = 3; 
    
elseif strcmpi(ext, '.txt')                                    % TEXT
    
    % Parse the text file
    success = parse_text(fp, fn, gui);
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);
    
    % Data domain set to 2 = frequency
    domain = 2; 
    
elseif strcmpi(ext,'userinput')                                % USER
    
    % Parse the userinput
    success = parse_userinput(gui);
    
    % Data domain set to 2 = frequency
    domain = 1; 
    
elseif strcmpi(ext,'.dat')                                     % DAT-file
    
    % Parse the dat-file
    success = parse_datfile(fp, fn, gui);
   
    % Calculate xaxis data struct
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);
    
    % Frequency domain - kspace
    domain = 2; 
    
elseif strcmpi(ext,'multi')
    
    % Parse multiple files
    success = parse_multifile(fp, fn, ext_tmp, gui);

    % Calculate xaxis if applicable
    CSI_2D_Scaling_calc_xaxis(hObj,[],1);

    % Set domain
    domain = 1; % Is it though?

else                                                           % ERROR 
    % Update LOG and return
    CSI_Log({'Warning. File format not supported. '},{ext});
    warning('The selected file format is not supported.'); return;
end


% Update LOG
if success
    CSI_Log({['Loading ' strrep(ext,'.','') ' file succeeded.']},{''});
else
    CSI_Log({['Loading ' strrep(ext,'.','') ' file failed.']},{''});
end

% Delete any user input from app-memory
delete_UserInput(hObj, gui);

% Fresh gui data
gui = guidata(hObj);    

                           % --- GUI updating --- %
                           
% Set domain - MRS data
if exist('domain', 'var') && exist('evt','var')
    if     domain == 2, evt.Source.Text = 'Frequency';        
    elseif domain == 3, evt.Source.Text = 'Time';   
    end
    gui.popup_domain.Value = domain;  guidata(hObj,gui); 
    CSI_setDomain(hObj,evt);
end

% Fresh gui data
gui = guidata(hObj);  

% Update visuals              
maxL = 40; % Max name-length
if isappdata(gui.CSIgui_main,'csi') 
    csi = getappdata(gui.CSIgui_main,'csi');
    if isfield(csi, 'filename')
        name = csi.filename;                             % Get filename
        if size(name,2) > maxL, name = name(1:maxL); end % Max length name
        set(gui.txt_fnCSI, 'String', name);              % Set in GUI
        gui.txt_fnCSI.TooltipString = csi.filename;
    end
end
if isappdata(gui.CSIgui_main,'mri')
    mri = getappdata(gui.CSIgui_main,'mri');
    if isfield(mri, 'filename')
        name = mri.filename;                             % Get filename
        if size(name,2) > maxL, name = name(1:maxL); end % Max length name
        set(gui.txt_fnIMG, 'String', name);              % Set in GUI
        gui.txt_fnIMG.TooltipString = mri.filename;      % Full name ttstr
        
    end
end

% --- Parse DAT file
function success = parse_datfile(fp, fn, gui)
% Load dat-file from Siemens platform into memory using the mapVBVD
% function from Philipp Ehses
%
% This function is basic as possible, loading the data into memory and
% setting the dimension labels. 
% Other data-operations and header-nfo will be parsed in separate functions

% Twix-loading
twix = mapVBVD([fp '\' fn '.dat']);

% Header
csi.twix = twix.hdr;

% Filename
csi.filepath = fp;

% Read data dimension and labels
dims = num2cell(twix.image.dataSize);
dims_rng = cellfun(@(x) 1:x, dims,'uniform', 0); % Range

% Dimension labels from header
dims_txt = twix.image.dataDims;
dims_txt = dims_txt(twix.image.dataSize>1);

% Correct data labels from twix-file.
% Labels in dat-file
labels_dat_file      = {'col','cha', 'lin','par','seg', 'ave', 'set'}; 
% Library for other name
labels_match_library = {'fid','chan','ky','kz','kx', 'aver', 'aver'};    
dims_txt_corrected = cell(1,size(dims_txt,2));
for ti = 1:size(dims_txt,2)        
    ind = strcmpi(labels_dat_file, dims_txt{ti});
    if sum(ind) > 0
        dims_txt_corrected{ti} = labels_match_library{ind};
    end
end
csi.data.labels = dims_txt_corrected;

% MRS-data
tmp = squeeze(twix.image(dims_rng{:}));

% Free up memory
clear('twix')

% Data sizes.
nfo = whos('tmp'); 
matlab_nfo = memory;

if nfo.bytes > matlab_nfo.MaxPossibleArrayBytes
    % Array is too large for matlab-array-memory
    fprintf('CSIgui: data requires much memory. Loading as single.\n')
    csi.data.raw = tmp; 
else
    % Array will fit matlab-max-array-memory
    csi.data.raw = double(tmp);    
end
clear('tmp')
csi.data.dim = size(csi.data.raw);

success = 0;
if csi.data.dim(1) > 0, success = 1; end

% Meta-data
csi.ext = 'dat'; csi.filename = fn; % Save filename

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi);  

% Close file ID
fclose('all');

% --- Parse LIST/DATA file
function success = parse_listdata(fp, fn, gui)
% Load list & data file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

[fp, fn, ~] = fileparts([fp '\' fn]); data_nfo = dir([fp '\' fn '.data']);
if ((data_nfo.bytes/1024^2) > 500)
    
    % Load list/data file into memory: use other script for larger files to
    % ease memory usage.
    [csi.data, csi.list] = csi_loadData_largeFiles([fp '\' fn], 0);
    
else

    % Load list/data file into memory.
    [csi.data, csi.list] = csi_loadData([fp '\' fn], 0);

end

% Restructure noise
if isfield(csi.data,'noise')
    tmp = csi.data.noise;
    csi.data = rmfield(csi.data, 'noise');
    csi.data.noise.raw = tmp;
    csi.data.noise.dim = size(csi.data.noise.raw);
    csi.data.noise.labels = cellfun(@char,...
        num2cell(65:65+numel(csi.data.noise.dim)-1), 'Uniform', 0);
    csi.data.noise.labels{1} = 'sec'; 
end

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
% Set domain
gui.popup_domain.Value = 2;


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
% Add dimension labels - Create labels 
csi.data.labels = csi.list.dim_labels;
% Set domain
gui.popup_domain.Value = 2;

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi);  

% --- Parse TEXT file
function success = parse_text(fp, fn, gui)
% Load text file, parse the data and store it into CSIgui appdata.
%
% Created: csi-struct.

% Load textfile
[csi.data.raw] = csi_readText([fp '\' fn '.txt']);

% Check if loading succeeded --> ELSE TRY JMRUI LOAD-TEXT FUNC
if isnan(csi.data.raw)
    CSI_Log({'Initial loading of text-file failed.'},...
            {'Reading data as a JMRUI-exported text-file.'}); 
    csi.data.raw = csi_readText_jmrui([fp '\' fn '.txt']);
    
    if isnan(csi.data.raw)
        success = 0; 
        CSI_Log({'Initial and JMRUI-method loading of text-file failed.'},...
                {''}); 
        return; 
    else
        success = 1; 
    end
   
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
old_format = 0; new_format = 0;
if strcmp(mat_cont(1).name,'csigui'),  old_format = 1;
elseif strcmp(mat_cont(1).name,'csi'), new_format = 1;
end

if (old_format + new_format) == 0
    success = 0;
    CSI_Log(...
    {'Incorrect mat-file. Use by CSIgui generated mat-file.',...
     'Expected fields:'},...
    {'Required structure: csigui (old) or csi (new).',...
     'raw, dim, filepath and name, noise, split, conv, mri,'}); return;
else, success = 1; 
end

% Integrety verified, read matfile
if old_format
    inp = load([fp '\' fn '.mat'], 'csigui');
    csigui = inp.csigui; clear inp;
elseif new_format
    inp = load([fp '\' fn '.mat'], 'csi');
    csi = inp.csi; clear inp;
end

              % ------- % Process MRI struct % ------- %

if old_format            

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
    
                  % ------- % Process CSI HDRNFO % ------- %
                  
    if isfield(csigui,'twix')
        csi.twix = csigui.twix; csigui = rmfield(csigui,'twix');
    end
    if isfield(csigui,'list')
        csi.list = csigui.list; csigui = rmfield(csigui,'list');
    end
    
                % ------- % Process CSI VOXMASK % ------- %

    if isfield(csigui,'voxelmask')
        csi.voxelmask = csigui.voxelmask;
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

elseif new_format
    
    if isfield(csi, 'conv')
        setappdata(gui.CSIgui_main,'conv',csi.conv); 
        csi = rmfield(csi,'conv');
    end

    if isfield(csi, 'mri')
        setappdata(gui.CSIgui_main,'mri', csi.mri); 
        csi = rmfield(csi,'mri');
    end

    if isfield(csi, 'log')
        CSI_Log({'Copying log from mat-file.'},{'See below.'})
        log_data = cellstr(csi.log);
        % Exclude first two lines
        log_data = log_data(3:end);
        
        % Add to log listbox
        gui.listbox_CSIinfo.String = ...
        cat(1, gui.listbox_CSIinfo.String,{''} ,log_data);
    
        % Remove log-field
        csi = rmfield(csi,'log');
    end

    % Set extensions
    csi.ext = 'mat';    
end

% Save CSI data in app-data
setappdata(gui.CSIgui_main,'csi',csi); 

% Clear memory
clear csigui
      
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
        csi.ext = '.data'; csi.filename = gui.inp.data.filename; % File nfo

        % Save CSI data in app-data
        setappdata(gui.CSIgui_main,'csi',csi);   
        
        button_CSI_ReorderDim_Auto_Callback(gui.CSIgui_main);
        gui = guidata(gui.CSIgui_main);
        

    catch err
        CSI_Log({'Failed processing user input!'},{err.message});   
        succes = 0; return;
    end

    % Update process info for user.
    CSI_Log({'User input processed: '},{userInfo});

    % Calculate xaxis from available data
    CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1); 
    gui = guidata(gui.CSIgui_main);
    csi = getappdata(gui.CSIgui_main,'csi'); 
    
    % Parse succes
    succes = 1;
end

% Label LIST % ----------------------------------------------------- %
% If list-struct is given (from csi_loadList/csi_loadData)
if isfield(gui.inp,'list')
    userInfo = 'List-struct loaded.';

    % Set input.
    csi.list = gui.inp.list; csi.filename = gui.inp.list.filename;
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
    gui = guidata(gui.CSIgui_main);
    csi = getappdata(gui.CSIgui_main,'csi'); 
    
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

% Label TWIX % ---------------------------------------------------- %
if isfield(gui.inp, 'twix')
     userInfo = 'Saved twix-header.';
     
     % Save twix data to CSIgui appdata
     csi.twix = gui.inp.twix;
     
     % Save CSI data in app-data
     setappdata(gui.CSIgui_main,'csi',csi);  
     
     % Update processed info to user.
     CSI_Log({'User input processed: '},{userInfo});
     
     % Pase succes
     succes = 1;
end

% --- Open multiple files
function success = parse_multifile(fp, fn, ext, gui)
% Open up multiple files to be loaded into memory. 
%
% Expects the same file extensions and resulting data requires the same
% dimensions. The header of the first file will be used.

% First we need to correct the order as Matlab does not save the order of
% the files as clicked... 
order = 0;
while size(unique(order),2) ~= size(fn,2)

    numbers = cell(1,size(fn,2));
    for nf = 1:size(fn,2)
        if nf == 1
            tmp = nf:size(fn,2);
        else
            tmp = [nf:size(fn,2) 1:nf-1];
        end
        numbers{nf} = strsplit(num2str(tmp));
    end
    uans = getUserInput_Popup(fn, numbers,...
           [gui.colors.main; gui.colors.text_main], 'Set load order.');
    if isempty(uans), order = 1:size(fn,2); break; end
    order = str2double(uans);
    
end
[~, sort_index] = sort(order);
fn = fn(sort_index);

switch ext
    case '.dat' % Spectroscopy - Siemens

        % Loop all files
        for kk = 1:size(fn,2)
            
            % Twix-loading
            twix = mapVBVD([fp '\' fn{kk} '.dat']);
            
            % Header
            if kk == 1
                csi.twix = twix.hdr;
            end

            % Read data dimension and labels
            dims = num2cell(twix.image.dataSize);
            dims_rng = cellfun(@(x) 1:x, dims,'uniform', 0); % Range
            
            % Dimension labels from header
            dims_txt = twix.image.dataDims;
            dims_txt = dims_txt(twix.image.dataSize>1);
            
            % Correct data labels from twix-file.
            % Labels in dat-file
            labels_dat_file      = {'col','cha', 'lin','par','seg', 'ave'}; 
            % Library for other name
            labels_match_library = {'fid','chan','ky','kz','kx', 'aver'};    
            dims_txt_corrected = cell(1,size(dims_txt,2));
            for ti = 1:size(dims_txt,2)        
                ind = strcmpi(labels_dat_file, dims_txt{ti});
                if sum(ind) > 0
                    dims_txt_corrected{ti} = labels_match_library{ind};
                end
            end
            
            % Store labels, only first iteration
            if kk == 1, csi.data.labels = [dims_txt_corrected 'file']; end
            
            % MRS-data
            tmp = squeeze(twix.image(dims_rng{:}));
            
            % Free up memory
            clear('twix')

            % Data sizes.
            nfo = whos('tmp'); 
            matlab_nfo = memory;
            
            if nfo.bytes > matlab_nfo.MaxPossibleArrayBytes
                % Array is too large for matlab-array-memory
                CSI_Log({'CSIgui: data requires much memory.'},...
                        {'Loading as single.'});
                tmp = single(tmp);
            else
                % Array will fit matlab-max-array-memory
                tmp = double(tmp);    
            end
        
            if kk == 1
                data = tmp; data_sz = size(data);
                file_dim = numel(size(tmp))+1;                
            else
                if (data_sz == size(tmp))
                    data = cat(file_dim,data,tmp);
                else
                    CSI_Log({'Aborting, data dimensions between'},...
                            {'files did no agree. Aborted loading.'});
                    success = 0;
                    return;
                end
            end
            clear('tmp');

            % Close file ID
            fclose('all');
        end

        % Store to CSI memory
        csi.data.raw = data;
        clear('data');

        % Data dimensions
        csi.data.dim = size(csi.data.raw);
                
        success = 0;
        if csi.data.dim(1) > 0, success = 1; end
        
        % Meta-data
        csi.ext = 'dat'; csi.filename = strjoin(fn, ' | '); % Save filename
        
        % Filename
        csi.filepath = fp;

        % Save CSI data in app-data
        setappdata(gui.CSIgui_main,'csi',csi);  


    case '.ima' % Image - DICOM - Siemens
    CSI_Log({'Loading multiple files for file type not implemented!'},...
            {''});
    case '.dcm' % Image - DICOM - Other
    CSI_Log({'Loading multiple files for file type not implemented!'},...
            {''});
    case '.data' % Spectroscopy - Philips
    CSI_Log({'Loading multiple files for file type not implemented!'},...
            {''});        
    case 'spar' % Spectroscopy - Philips
    CSI_Log({'Loading multiple files for file type not implemented!'},...
            {''});
    case 'par' % Image - Philips
    CSI_Log({'Loading multiple files for file type not implemented!'},...
            {''});
    
    case '.mat' % Spectroscopy - CSIgui
    CSI_Log({'Reading all mat-files.'},...
            {'Please be patient, load-time dependent on file-sizes.'});
        
        % Load in every mat-file
        % Get data.raw catenate - if sizes are equal!
        % Save one file its conv/mri/csi-structs 

        % Loop all files
        mri = cell(1,numel(fn)); conv = cell(1,numel(fn));
        csi_dump = cell(1,numel(fn));
        for kk = 1:size(fn,2)

            % Check mat-file integrity   
            mat_cont = whos('-file',[fp '\' fn{kk} '.mat']); 
            old_format = 0; new_format = 0;
            if strcmp(mat_cont(1).name,'csigui'),  old_format = 1;
            elseif strcmp(mat_cont(1).name,'csi'), new_format = 1;
            end
            
            if (old_format + new_format) == 0
                success = 0;
                CSI_Log(...
                {'Incorrect mat-file. Use by CSIgui generated mat-file.',...
                 'Expected fields:'},...
                {'Required structure: csigui (old) or csi (new).',...
                 'raw, dim, filepath and name, noise, split, conv, mri,'}); 
                return;
            else, success = 1; 
            end
            
            % Integrity verified, read matfile
            if old_format
                inp = load([fp '\' fn{kk} '.mat'], 'csigui');
                csigui = inp.csigui; clear inp;
            elseif new_format
                inp = load([fp '\' fn{kk} '.mat'], 'csi');
                csi = inp.csi; clear inp;
            end            
                                      
            if old_format            
                
                % ------- % Process MRI struct % ------- %  
                % Find MRI struct and store             
                if isfield(csigui,'mri')
                    mri{kk} = csigui.mri; csigui = rmfield(csigui,'mri');
                end

                % ------- % Process CONV struct % ------- %               
                % Find conv struct and store
                if isfield(csigui,'conv')
                    conv{kk} = csigui.conv; csigui = rmfield(csigui,'conv');
                end
                
                % ------- % Process CSI HDRNFO % ------- %
                if isfield(csigui,'twix')
                    if kk == 1, twix = cell(1,numel(fn)); end
                    twix{kk} = csigui.twix; csigui = rmfield(csigui,'twix');
                end
                if isfield(csigui,'list')
                    if kk == 1, list = cell(1,numel(fn)); end
                    list{kk} = csigui.list; csigui = rmfield(csigui,'list');
                end
                
                % ------- % Process CSI VOXMASK % ------- %        
                if isfield(csigui,'voxelmask')
                    csi.voxelmask = csigui.voxelmask;
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
                    csi.xaxis = csigui.xaxis; 
                    csi.data = rmfield(csi.data,'xaxis');
                end
                
                if isfield(csigui, 'ori')
                    csi.ori = csigui.ori; 
                    csi.data = rmfield(csi.data,'ori');
                end

                % ------- % Process LOG % ------- %                
                % Resubmit log if available
                if isfield(csigui,'log')
                    if kk == 1, log_data = cell(1, numel(fn)); end
                    log_data_tmp = cellstr(csigui.log);
                    % Exclude first two lines
                    log_data{kk} = log_data_tmp(3:end);                                       
                
                    % Remove log-field
                    csi.data = rmfield(csi.data,'log');
                end
            
            elseif new_format
                
                if isfield(csi, 'conv')
                    conv{kk} = csi.conv; csi = rmfield(csi,'conv');
                end
            
                if isfield(csi, 'mri')
                    mri{kk} = csi.mri; csi = rmfield(csi,'mri');
                end
            
                if isfield(csi, 'log')
                    if kk == 1, log_data = cell(1, numel(fn)); end
                    log_data_tmp = cellstr(csi.log);
                    % Exclude first two lines
                    log_data{kk} = log_data_tmp(3:end);   
                
                    % Remove log-field
                    csi = rmfield(csi,'log');
                end
            
                % Set extensions
                csi.ext = 'mat';    
            end
            
            % Store app-data csi-struct of this file
            csi_dump{kk} = csi;

        end

        % Ask the user which file to use for default header info.
        % Save only data.raw from other files
        qry = {'Which file-nfo and header should be stored and used?'};
        uans_foi = getUserInput_Popup(qry ,{fn});
        if isempty(uans_foi)
            CSI_Log({'Loading canceled.'},{''}); return; 
        end
        foi_ind = contains(fn, uans_foi);

        % Check size-compatibility
        data_sz = cellfun(...
            @size, extractField(csi_dump, 'data.raw'), 'uniform', 0);
        data_sz_check = cellfun(...
            @isequal, data_sz, repmat(data_sz(foi_ind),size(data_sz)));
        if sum(data_sz_check) ~= numel(fn)
            % File of doom - wrong size
            fod_ind = find(data_sz_check == 0);
            
            % Remove FODs from memory
            csi_dump(fod_ind) = []; mri(fod_ind) = [];
            conv(fod_ind) = []; log_data(fod_ind) = [];
            
            % Show user which files are removed
            tmp = repmat({'Removed:'}, 1, numel(fod_ind));
            qry = ...
            [{'The files below have incompatible data-dimensions.'},...
            tmp(:)'];
            CSI_Log({'% -------------------------------------- %'}, {''});
            CSI_Log(qry, {'', fn{fod_ind}});
            CSI_Log({'% -------------------------------------- %'}, {''});

            fn(fod_ind) = [];
            foi_ind = contains(fn, uans_foi);
        end

        % Concatenate data and add to struct of interest
        CSI_Log({'Combining mat-files and storing it into app-data.'},{''});
        csi = csi_dump{foi_ind};
        nDim = numel(csi.data.labels);
        data = extractField(csi_dump, 'data.raw');                     
        csi.data.raw = cat(nDim + 1, data{:});
        label_name = 'file';
        if sum(contains(csi.data.labels, 'file')) >= 1, label_name = 'inp';
        end
        csi.data.labels{end+1} = label_name;
        csi.data.dim = size(csi.data.raw);
        csi.filepath = fp;
        csi.filename = strjoin(fn, ' | ');
        clear csi_dump;

        % Save CSI data in app-data
        setappdata(gui.CSIgui_main,'csi',csi); 

        % Save Conv data in app-data
        if ~isempty(conv{foi_ind})
            setappdata(gui.CSIgui_main,'conv', conv{foi_ind});
            clear conv;
        end

        % Save MRI data in app-data
        if ~isempty(mri{foi_ind})
            setappdata(gui.CSIgui_main,'mri', mri{foi_ind}); 
            clear mri;
        end

        CSI_Log({'Copying logs from mat-files.'},{'See below.'})
        % Add to logs listbox
        log = cell(1,numel(fn)*2);
        for kk = 1:numel(fn)
            ind = (kk*2) - 1;
            log{ind+1} = log_data{kk};
            log{ind} = {'', '% -------------------------------- %',...
                        fn{kk},...
                            '% -------------------------------- %', ''}';                     
        end

        gui.listbox_CSIinfo.String = ...
                cat(1, gui.listbox_CSIinfo.String,{''} ,cat(1, log{:}));
end % End of switch-case loop

% --- Parse DICOM file
function success = parse_dicom(fp, fn, gui)
% Analyse dicom file, parse the data and store it into CSIgui appdata.
%
% Created: mri-struct


% \\ Is it Philips or Siemens?
[~, fn, ext] = fileparts(fn); ext = lower(ext);
switch ext
    case '.dcm'

        % Read and order dicom data
        [m,i,g] = dicomread7T({[fp '\' fn]}); mri = struct;
        if isnan(m)
            CSI_Log({'Loading DICOM file failed.'},{''}); 
            success = 0; return; 
        else
            success = 1;
        end

        % Order dicom data
        [mri.data, mri.par] = dicomorder3(m, i); 
        
        mri.examinfo = g; mri.ext = 'dcm'; 
        mri.filename = fn; mri.filepath = fp;
    
    case '.ima'
               
        % Question user to load in the single file or load group of dicom
        % files into memory.
        qry = 'Load single file or all files from image directory: ';
        def = {'All','Single'};
        uans = getUserInput_Popup({qry},{def}, [], 'Load DICOM');
        if isempty(uans), success = 0; return; end
        if strcmpi(uans{1},'all')
            % get all files.
            files = dicomreadSiemens_getSeries(fp,fn);
        else
            files = {[fn ext]};
        end
        
        try
            % Load dicom files
            [img, nfo] = dicomreadSiemens(fp, files);  
            mri = struct;
            mri.data.M = img; mri.par = nfo;
            success = 1;
        catch err
            fprintf('%s\n',err.message);
            CSI_Log({'Loading DICOM file failed.'},{''}); 
            success = 0;
            return;
        end
        
        try
            % Process nfo
            examinfo = dicomreadSiemens_sortParameters(nfo);
            % Store data            
            mri.examinfo = examinfo; mri.ext = 'ima';
            mri.filename = fn; mri.filepath = fp;
        catch err
            fprintf('%s\n',err.message);
            % Exam info loading failed.
            mri.examinfo = NaN; mri.ext = 'ima';
            mri.filename = fn; mri.filepath = fp;
            CSI_Log({'Analysing DICOM header failed.'},{'Images loaded.'}); 
        end
               
        % Siemens dicom fields of interest.
        % nfo{1}.SliceThickness
        % nfo{1}.PatientPosition
        % nfo{1}.ImagePositionPatient
        % nfo{1}.ImageOrientationPatient
        % nfo{1}.SliceLocation
        % nfo{1}.SlicePosition_PCS      
end
        


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


% Res RL - INDEX 1
if isfield(nfo,'VoxRL_mm') 
   csi.ori.res(1) = nfo.VoxRL_mm;     
   csi.ori.offcenter(1) = [nfo.RL_L__mm]; 
   
elseif isfield(nfo,'ACQVoxelSizeRL_mm')
   csi.ori.res(1) = nfo.ACQVoxelSizeRL_mm;     
   csi.ori.offcenter(1) = [nfo.RL_L__mm]; 
   
elseif isfield(nfo, 'VoxelSizeRL_mm') 
   csi.ori.res(1) = nfo.VoxelSizeRL_mm;
   csi.ori.offcenter(1) = [nfo.RL_L__mm];
end

% Res AP - INDEX 2
if isfield(nfo,'VoxAP_mm') 
    csi.ori.res(2) = nfo.VoxAP_mm;
    csi.ori.offcenter(2) = [nfo.SliceOffc_AP_P__mm]; 
elseif isfield(nfo,'ACQAP_mm')
   csi.ori.res(2) = nfo.ACQAP_mm;     
   csi.ori.offcenter(2) = [nfo.SliceOffc_AP_P__mm];
end

% Res FH - INDEX 3
if isfield(nfo,'VoxFH_mm')
   csi.ori.res(3) = nfo.VoxFH_mm;
   csi.ori.offcenter(3) = [nfo.FH_H__mm]; 
elseif isfield(nfo,'ACQAP_mm')
   csi.ori.res(3) = nfo.ACQFH_mm;     
   csi.ori.offcenter(3) = [nfo.FH_H__mm];
end       


    
csi.nfo = nfo;

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
% DimN_step                                 Resolution
% SliceOffc_AP_P__mm RL_L__mm FH_H__mm'     Offcenter
       
% Copy some frequency parameters

% Covert nucleus 
% OFFCENTER
if isfield(nfo,'si_lr_off_center')
    csi.ori.offcenter(1) = nfo.si_lr_off_center;
end
if isfield(nfo,'si_ap_off_center')
    csi.ori.offcenter(2) = nfo.si_ap_off_center;
end
if isfield(nfo,'si_cc_off_center')
    csi.ori.offcenter(3) = nfo.si_cc_off_center;
end


% RESOLUTION
if isfield(nfo,'slice_distance')
    csi.ori.res(3) = nfo.slice_distance;
end
if isfield(nfo, 'dim2_step'), csi.ori.res(1) = nfo.dim2_step; end
if isfield(nfo, 'dim3_step'), csi.ori.res(2) = nfo.dim3_step; end


% Save nucleus and BW to axis info
csi.xaxis.nucleus = nfo.nucleus;
csi.xaxis.BW = nfo.sample_frequency;

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
            
% --- Executes on button press in button_CSI_ReadInfo2.
function button_CSI_ReadInfo2_Callback(~, ~, gui)
button_CSI_ReadInfo_Callback([],[],gui);


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
          
% --- Execute by user click in menubar
function openGithub(~,~,~)
web('https://github.com/Sugocy/CSIgui', '-browser');

    
% getUserInput EDITS GUI % --------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI to get user input
function userInput = getUserInput(question, defans, clrs, fig_title)
% Opens up a menu for user to input data.
% question: cell with each required user input
% defans  : default answers.
%
% If user canceled, userInput is returned empty.


% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2 || isempty(clrs)
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

if nargin < 4, fig_title = 'CSIgui - UserInput';
else, fig_title = ['CSIgui - ' fig_title];
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name',fig_title,'Tag', 'CSIgui_UI_edits'); 
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


% GetUserInput: TICKS GUI % -------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI to get user input via tick-boxes
function userInput = getUserInput_Tick(tick_title, tick_input, fig_title, tick_def, clrs)
% Create a simple gui with dropdown e.g. tick marks to quickly choose
% between two options

% If no color input is present, use default dark theme or get from main 
% app.
if nargin < 5
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

if nargin < 3, fig_title = 'CSIgui - UserInput'; 
else,          fig_title = ['CSIgui - ' fig_title];
end

if nargin < 4, tick_def = zeros(size(tick_title)); end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name',fig_title,'Tag', 'CSIgui_UI_tick'); 
subdat = guidata(fig_userinp); axis off;

% Nr of dropdown menu's
Ninput = size(tick_title,2);

% Total size of figure.
dw = 10;  dy = 10; w = 320+2*dw; 
% Define values for EDIT and TEXT sizes.
sz_tick  = [(w-(2*dw))/2.5 25]; sz_text = [w-2*dw 20]; 
sz_button = [(w-2*dw)/2 20];

% Define height
h = Ninput * (sz_tick(2)*2) + (sz_button(2)*2) + 2*dy;

% Size tune parameters and set figure size and window-position
scrsz = get(0,'Screensize'); 
figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

% Set each title text and popup menu
for kk = 1:Ninput
    %Tick menu title
    subdat.h_title{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Text','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_tick(2)*kk + sz_text(2)*(kk-1)) ...
                    sz_text(1) sz_text(2)],...
        'Tag', 'h_popup','String',tick_title{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left',...
        'TooltipString', tick_title{kk});
    
    % Tick 1  
    subdat.h_tick{kk,1} = ...
        uicontrol(fig_userinp, 'Style', 'checkbox','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_tick(2)*kk + sz_text(2)*(kk))...
                    sz_tick(1) sz_tick(2)],...
        'Tag', sprintf('%i_%i',kk,0),'String',tick_input{kk}{1},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left', 'value', tick_def(kk),...
        'Callback', @getUserInput_Tick_setToggle);

    % Tick 2   
    subdat.h_tick{kk,2} = ...
        uicontrol(fig_userinp, 'Style', 'checkbox', 'Unit', 'Pixels',...
        'Position',[dw+sz_tick(1) h-dy-(sz_tick(2)*kk + sz_text(2)*(kk))...
                    sz_tick(1) sz_tick(2)],...
        'Tag', sprintf('%i_%i',kk,1),'String',tick_input{kk}{2},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left', 'value', abs(tick_def(kk)-1),...
        'Callback', @getUserInput_Tick_setToggle);
    
end

% Set button CONTINUE
subdat.h_buttonContinue = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw        dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_CONTINUE','String','Continue','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Tick_setOutput);

% Set button SKIP
subdat.h_buttonSkip = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw+sz_button(1)  dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_SKIP','String','Skip','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @getUserInput_Tick_setOutput);

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

% -. Exected by ticking a box - allows toggle.
function getUserInput_Tick_setToggle(hobj,~,~)
% Toggle the ticks of the query of interest. Only one can be ticked per
% title. 

% Which tick was tagged: first is qry-number, second is tick-number
tick_tag = str2double(strsplit(hobj.Tag, '_'));

% Gui data
gui = guidata(hobj.Parent);

% Set other tick-box of this query to opposite value
gui.h_tick{tick_tag(1), abs(tick_tag(2)-1)+1}.Value = abs(hobj.Value-1);

% -. Executed by button in getUserInput_Tick to save answer
function getUserInput_Tick_setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

%userInputSource = get(hObject,'Style');
userInputString = get(hObject,'String'); % What did get us here? 

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','CSIgui_UI_tick');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

if strcmp(userInputString, 'Skip')
    val = [];
else
    % Get strings input from popup menu
    val = NaN(1,size(subdat.h_tick,2));
    for kk = 1:size(subdat.h_tick,1)
        val(kk) = subdat.h_tick{kk}.Value;        
    end
end

% Set userinput value and Update storage
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);



% getUserInput: POPUP GUI % -------------------------------------------- %
% ---------------------------------------------------------------------- %

% -. Opens UI with popup menu's to get user input
function userInput = getUserInput_Popup(popup_title, popup_input, clrs, fig_title)
% Create a simple gui with dropdown e.g. popup menu to choose certain 
% options.
%
% Colors can be kept empty to set figure title.

% If no color input is present, use default dark theme or get from main 
% app.
if nargin == 2 || isempty(clrs)
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

if nargin < 4, fig_title = 'CSIgui - UserInput'; 
else,          fig_title = ['CSIgui - ' fig_title];
end

% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name',fig_title,'Tag', 'CSIgui_UI_popup'); 
subdat = guidata(fig_userinp); axis off;

% Nr of dropdown menu's
Ninput = size(popup_title,2);

% Total size of figure.
dw = 10;  dy = 10; w = 320+2*dw; 
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
        'HorizontalAlignment', 'Left',...
        'TooltipString', popup_title{kk});
    
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
function userInput = getUserInput_Buttons(qst, button_text, clrs, fig_title)
% Create a simple gui with two buttons to choose between two options.
% qst:          One strin/char array for the question.
% Button text:  Two string/char arrays for both buttons.

% If no color input is present, use default dark theme or get from main 
% app.
if nargin < 3 || isempty(clrs)
    csi_obj = findobj('type','figure','Tag','CSIgui_main');
    if isempty(csi_obj)
        clrs = [0.000 0.000 0.000; 0.941 0.941 0.941]; 
    else
        csi_gui = guidata(csi_obj);
        clrs = [csi_gui.colors.main; csi_gui.colors.text_main];
    end
end

if nargin < 3, fig_title = 'CSIgui - UserInput';
else,          fig_title = ['CSIgui - ' fig_title];
end
% --------------------- % Set GUI details

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', clrs(1,:),...
                     'Name',fig_title ,'Tag', 'CSIgui_UI_button'); 
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
% --- Executed by button in getUserInput_Button to save answer
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
if nargin < 3
    obj = findobj('Type', 'Figure','Tag', 'CSIgui_main');
    gui = guidata(obj);
end
if nargin < 4, backup = gui.checkbox_backup.Value; end

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
[data_averaged, nsa_ind] = ...
    CSI_average(csi.data.raw,nsa_ind, csi.data.labels);

% Checksum
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
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before apodization.'); end

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
function [data_averaged, index] = CSI_average(data, index, labels)
% Average data over specific index. If the index not given or empty 
% empty, the user is asked to enter an index e.g. dimension to average 
% over.

if nargin < 3, labels = strsplit(num2str(1:numel(size(data)))); end

if nargin == 1 || isempty(index) % If not or an empty index is given
    uans = getUserInput_Popup({'Index to average:'},{labels});    
    if isempty(index), data_averaged = []; return; end
    index = csi_findDimLabel(csi.data.labels,uans);
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
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % ------------------- %

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'time')
    CSI_Log({'MRS data is in time domain; '},...
    {'change to frequency domain (k-space) to apply FFT to spatial time domain.'});
    return;
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before spatial FFT.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% GET OPTION INPUT % ------------------- %

quest = {'Specify the spatial FFT method, use automatic if unknown:'};
defan = {{'Automatic', 'Circular shift', 'FFT shift'}};
uans = getUserInput_Popup(quest, defan, [], 'Spatial FFT');
if isempty(uans), return; end

% Set option
switch uans{1}
    case 'Automatic', shift_opt = 2;
    case 'Circular shift', shift_opt = 0;
    case 'FFT shift', shift_opt = 1;
end

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

% Set domain
CSI_setDomain(gui.CSIgui_main,[],'time');

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
% Apply forward fourier to transform the MRSI data from the spatial(voxels)
% time(FID)domain to the spatial(voxels) frequency(spectrum) domain.
% 
% Uses csi_fft();
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % ------------------------------- %

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'freq')
    CSI_Log({'MRS data is in frequency domain; '},...
            {'change it to time domain to apply forward FFT.'});
    qry = {'MRS-data is in frequency domain, apply FFT or abort?'};
    uans = ...
    getUserInput_Buttons(qry,{'Continue', 'Abort'});
    if isempty(uans) || strcmp(uans, 'Abort'); return; end
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before forward FFT.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% FFT % -------------------------------------------- %

% Get user input
uans = getUserInput_Popup(...
           {'Correct for N term?',...
            'Correct for onesided fft in spectroscopy?',...
            'Shift before and after FFT? (Echo)'},...
           {{'No','Yes'},{'Yes','No'},{'No','Yes'}},[],'FFT');
if isempty(uans), CSI_Log({'Skipped FFT.'},{''}); return; end

if strcmpi(uans{1},'yes'), correct_N = 1; else, correct_N = 0; end
if strcmpi(uans{2},'yes'), onesided = 1;  else, onesided = 0; end
if strcmpi(uans{3},'yes'), dbl_shift = 1; else, dbl_shift = 0; end

% Apply FFT
csi.data.raw = csi_fft(csi.data.raw, correct_N, onesided, dbl_shift);

% CLEANUP % ---------------------------------------- %

% Set domain
CSI_setDomain(gui.CSIgui_main,[],'freq');

% Save appdata
setappdata(gui.CSIgui_main, 'csi', csi);
% After update info
CSI_Log({'Forward FFT applied.'},{''});

% --- Executes on button press in button_CSI_iFFT.
function button_CSI_iFFT_Callback(~, ~, gui, backup)
% Apply inverse fourier to transform MRS data from the frequency domain to
% the time domain.
%
% Uses csi_ifft();
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % ------------------- %

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'time')
    CSI_Log({'MRS data is in time domain; '},...
            {'change to frequency domain to apply backward FFT.'});
    qry = {'MRS-data is in time domain, apply iFFT or abort?'};
    uans = ...
    getUserInput_Buttons(qry,{'Continue', 'Abort'});
    if isempty(uans) || strcmp(uans, 'Abort'); return; end
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before inverse FFT.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USER DATA % -------------------------- %

uans = getUserInput_Popup(...
           {'Correct for onesided fft in spectroscopy?'},...
           {{'Yes','No'}}, [], 'iFFT');
if isempty(uans), CSI_Log({'Skipped FFT.'},{''}); return; end
if strcmpi(uans{1},'yes'), onesided = 1;  else, onesided = 0; end

% iFFT % ------------------------------- %

% Apply iFFT
csi.data.raw = csi_ifft(csi.data.raw, onesided);

% CLEAN UP % --------------------------- %

% Set domain
CSI_setDomain(gui.CSIgui_main,[],'time');

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
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % --------------------------------- %

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'time')
    CSI_Log({'MRS data is in time domain; '},...
     {'change to frequency domain (k-space) to apply 2D/3D apodization.'});
    return;
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before k-space apodization.'); end

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
                           {'No','Yes'}}, [] ,'K-space Apodization');
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
if nargin < 4, backup = gui.checkbox_backup.Value; end



% BACKUP + APPDATA % --------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before 1D apodization.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% DOMAIN % ------------------------------------- %

% Check data domain
domain = CSI_getDomain(gui); doFFT = 0;
if ~strcmp(domain, 'time')
    CSI_Log({'MRS data is possibly in frequency domain; '},...
            {'Apodization requires time domain data.'});
    qry = {'MRS-data is possibly in frequency domain:'};
    def = {{'Ignore', 'Convert', 'Abort'}};
    uans = getUserInput_Popup(qry,def,[],'Apodization');
    if isempty(uans) || strcmp(uans{1}, 'Abort'); return; end
    if strcmp(uans{1}, 'Convert'), doFFT = 1; end
end 

% USER INPUT % ---------------------------------- %

% Get filter to apply.
quest = {'Filter algorithm:', 'FID or Echo:', 'Flip window:'};
defan = {{'Gaussian', 'Hamming','Hann', 'Exponential','Blackman',...
          'Flattop', 'Sinebell'}, {'FID','Echo'},{'No','Yes'}};
uanstype = getUserInput_Popup(quest,defan,[], 'Apodization');
if isempty(uanstype), CSI_Log({'Skipped apodization.'},{''}); return; end

% Get options for specific algorithms
switch uanstype{1}
    case 'Gaussian'
        % ((1/(1/BW x nSamples) * (ans/nSamples)) / pi) == (bw/ans)/pi 
        % Hz = bw / (N*pi) && N = (Hz * Pi)/BW
        if isfield(csi.xaxis,'BW')
            uopts = getUserInput({'Apodization factor: (Hz)'},{20});    
            if isempty(uopts)
                CSI_Log({'Skipped apodization.'},{''}); return; 
            end
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

% Flip window
if strcmp(uanstype{3},'Yes'), opts(3) = 1; else, opts(3) = 0; end

% DO FFT % ---------------------------------------- %
if doFFT, csi.data.raw = csi_ifft(csi.data.raw); end

% APPLY FILTER % ---------------------------------- %

% Create window and apply filter
[csi.data.raw, win] = CSI_filterSpectra(csi.data.raw, uanstype{1}, opts);


% DO FFT % ---------------------------------------- %
if doFFT, csi.data.raw = csi_fft(csi.data.raw); end


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
uans = getUserInput_Popup({'Normalize by: ','', 'Use component: '},...
    {{'Maximum','Peak','Noise','Value'},...
     {'per voxel','in volume'},...
     {'Real','Imaginary'}},[],'Normalize');
if isempty(uans), CSI_Log({'Skipped normalize data.'},{''}); return; 
end  

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get data-array
data = csi.data.raw;

% Split units: use real for normalization
switch uans{3}
    case 'Real'
        dataForNorm = CSI_getUnit(data,'Real');
    case 'Imaginary'
        dataForNorm = CSI_getUnit(data,'Imaginary');
end

% dataR = CSI_getUnit(data,'Real'); dataI = CSI_getUnit(data,'Imaginary');

% Data as cell index layout
sz = size(data); 
cell_layout = ...
arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);

% Create cell of data for normalization.
datacForNorm = mat2cell(dataForNorm, sz(1), cell_layout{:});

switch uans{1}
    case 'Peak'
        
        % Get peak of interest
        [~, ~, peak_range] = CSI_getDataAtPeak(data, csi.xaxis);
        if isnan(peak_range), return; end
                
        switch uans{2}
                        
            case 'per voxel'
                % Normalize to a specific peak

                % Normalization factor
                doi = dataForNorm(peak_range(1):peak_range(2),...
                    cell_layout{:}); 
                val = max(doi,[],1);

                % Normalize
                data = data ./ val;
                
            case 'in volume'
                % Normalize to maximum peak in the volume
               
                if numel(sz) >=2, nDimC = num2cell(sz(2:end));
                else,             nDimC = {1};
                end
                nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

                % Normalization factor
                doi = dataForNorm(peak_range(1):peak_range(2), nDimC{:}); 
                val = max(doi,[],1); [val, pos] = max(val(:));

                % Normalize
                data = data ./ val;
                
                sub = ind2subQ(size(max(doi,[],1)), pos);
                CSI_Log({'Index of max-peak in volume:'},...
                    {sprintf(' %i |',sub{2:end})});

        end
    case 'Maximum'
        
        switch uans{2}
            case 'per voxel'
                % Normalize to maximum of each voxel
                
                % Get normalization factor
                val = max(dataForNorm,[],1); 
                
                % Normalize
                data = data./val;
               
            case 'in volume'
                
                % Normalize to maximum in volume
                val = max(dataForNorm,[],1); val = max(val(:));
                
                ind = find(dataForNorm == val);
                sub = ind2subQ(size(dataForNorm), ind);
                
                % Normalize
                data = data ./ val;
    
                CSI_Log({'Index of maximum in volume: '},...
                    {sprintf(' %i |',sub{:})});
        end
        
    case 'Noise'
        
        % Per voxel
        switch uans{2}
            case 'per voxel'
                % Normalize to noise per voxel

                % Noise mask
                uans = getUserInput({'Size of noise mask: ',...
                    'Save noise/voxel to file (y/n) : '},{'50','y'});
                if isempty(uans), return; end
                msksz = str2double(uans{1}); savenoise = uans{2};
               
                % Noise = absolute std of last N samples
                noiseV = cellfun(@(x) abs(std(x(end-msksz+1:end))), ...
                    datacForNorm,'Uniform',0);
                noiseV = cell2mat(noiseV);

                % Normalize
                data = data ./ noiseV;
                
                % Save noise if requested.
                if strcmp(savenoise,'y')
                    currdir = cd; filt = {'*.mat','MATLAB File (*.mat)'};
                    [fn, fp] = uiputfile(filt,'Save Noise Data',currdir);
                    save([fp '\' fn],'noiseV');
                end
        
            case 'in volume'
                
                % Volume
                CSI_Log({'Normalize to noise in volume unavailable!'},{''})
        end
        
    case 'Value'
        
        % Get value from user
        tuans = getUserInput({'Value for normalization:'},{'1'});
        if isempty(tuans), CSI_Log({'Skipped normalization.'},{''}); 
            return; 
        end
        val = str2double(tuans{1});
        
        data = data ./ val;
        
        CSI_Log({'Normalized data with user-value: '},{val});
end

% Store data.
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

if size(opts,2) > 2
    if opts(3) == 1
        window = flip(window,1);
    end
end

% Rearrange data to cell: every FID in a cell
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
% Uses: csi_lineWidth();

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Userinput % --------------------------------- %
uans = getUserInput_Popup({'Save line width data (.txt): ',...
                           'Use FWHM method: (BETA)'},...
                         {{'No','Yes'},  {'Intersect','Local minima'}},...
                         [], 'FWHM');
if isempty(uans)
    CSI_Log({'Skipped linewidth calculations.'},{''}); return; 
end


% Save data boolean
switch uans{1}, case 'Yes', dataSave = 1; case 'No', dataSave = 0; end

% Method
lwmeth = uans{2};

% Get data at peak of interest
[doi, ~, range] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);



% Calculate fwhm via prefered method
if strcmp('Local minima', lwmeth) % Local minima method
    % doi:      data of peak of interest only
    % tmp:      full data array
    % range:    x-axis range of the peak of interest
    
    % Maximum in this range
    [~,mi] = max(real(doi),[],1);

    % Set the maximum value as peak centre.
    peak_pos = mi+range(1);  % And correct doi-size vs full axis

    % Set 1/2(peak width)
    peak_width = ceil(csi.xaxis.N/100); % One percent of #samples
    if peak_width < 3, peak_width = 3; end 
    % If 1% of #samples is less than 3, 1/2*(peak width) set to 3.
    
    % Peak range used to find FWHM
    peak_pLow = num2cell(peak_pos - peak_width);
    peak_pHig = num2cell(peak_pos + peak_width);
    % Peak of interest range in x-axis
    poi_range_perVox = cellfun(@(x,y) x:y, peak_pLow, peak_pHig,'uniform',0); 

    % Input cell data
    sz = size(csi.data.raw); 
    cell_layout = arrayfun(...
        @ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
    tmp = mat2cell(csi.data.raw, sz(1), cell_layout{:}); 

    % Input axis
    ax_inp = repmat({csi.xaxis.ppm}, size(tmp));
    % Set plotting off - Overload matlab otherwise!
    plot_off = repmat({0},size(tmp)); % plot_off{1,2,4,1} = 1;


    % Calculate line width % --------------------------------- %
    % Input: data/x-axis/peak range/plot
    linewidth = cellfun(@csi_lineWidth, tmp, ax_inp, poi_range_perVox,...
        plot_off, 'Uniform',0);

elseif strcmp('Intersect', lwmeth)  % Intersect method
    
    % Log
    CSI_Log({'Calculating linewidth, please be patient.'},{''});
    
    % Rearrange data to cell 
    sz = size(csi.data.raw); 
    cell_layout = arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),...
        'UniformOutput',0);
    data = mat2cell(csi.data.raw, sz(1), cell_layout{:}); 
    
    % Axis input
    if isfield(csi.xaxis,'ppm')
        inp_xaxis = repmat({csi.xaxis.ppm}, size(data));
    else
        inp_xaxis = repmat({csi.xaxis.none}, size(data));
    end
    inp_range = repmat({range}, size(data));
    
    % Run script for intersect fwhm
    linewidth = cellfun(@csi_linewidth_intersect,...
        data, inp_xaxis, inp_range, 'Uniform',0);
    
end

% Show statistics nfo

linewidth(cellfun(@isnan,linewidth)) = {0};
linewidth = cellfun(@double,linewidth);
stats = csi_statistics_of_volume(linewidth);
CSI_Log({'Linewidth statistics ------------------------------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});

% \\ Display Data
CSI_dataAs_Initiate(linewidth, 'Linewidth', gui, csi.data.labels);

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

% --- Executes on button press in button_CSI_ISIS.
function button_CSI_ISIS_Callback(~, ~, gui)
disp('This function is in development!');

% Set a backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before ISIS recon.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

uans = getUserInput_Popup({'Cycle dimension: '},{csi.data.labels(2:end)},...
    [], 'ISIS');
if isempty(uans), CSI_Log({'Skipped ISIS'},{''}); return;  end

% Add/Subtr scheme of ISIS data
udimstr = uans{1}; udim = csi_findDimLabel(csi.data.labels, {udimstr});
udimsz = size(csi.data.raw, udim); soi = ones(1, udimsz);

% Go nuts and add/subtract according to chosen scheme over the dimension
% called udimstr at index udim.

% Data to ISIS recon
data = csi.data.raw; sz = size(data); 

% Set chosen cycli index (udim) to the first index. 
pvec = 1:numel(sz); pvec(udim) = []; pvec = [udim pvec];
data = permute(data,pvec);

% New data output size for empty variable
sznew = size(data); sznew(1) = 1; 


% Add and subtract using scheme of interest (SOI) at data-index udim.
outp = zeros(sznew); 
misc_index = arrayfun(@(x) 1:x, sznew(2:end), 'uniform', 0);
for kk = 1:sz(udim)
    if soi(kk)  % Add
       outp = outp + data(kk,misc_index{:});
       
    else        % Subtract
       outp = outp - data(kk,misc_index{:});
    end  
end
 
% Permute to original index
data = ipermute(outp,pvec);

% Replace in appdata
csi.data.raw = data; csi.data.dim = size(data);

% Save data
setappdata(gui.CSIgui_main, 'csi',csi);

% Update info
schemes_str = repmat('+', 1, sz(udim));
CSI_Log({'Applied ISIS add/subtract scheme to data:'},{schemes_str} );

% --- Executes on button press in button_CSI_setParameters.
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
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before reordering.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% PERMUTE % ---------------------------- %
% Get order from user
uans = getUserInput(...
    {'New order of MRS data dimensions?: '},{1:size(csi.data.labels,2)});
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

% Create new labels.
empty_lab = cellfun(@isempty, new_label) == 1;
if sum(empty_lab) ~= 0
    % String characters to double
    labelAsNum = cellfun(@double, new_label, 'UniformOutput', false);
    % Number of characters
    nChar = cellfun(@size,new_label, repmat({2},1,size(new_label,2)));
    % Index of single charactes
    strt = max(cell2mat(labelAsNum(nChar == 1)))+1;
    % If no single character label found
    if isempty(strt), strt = 65; end

    % Create new label names for empty indexes.
    ind = find(empty_lab == 1);
    for kk = 1:size(ind,2)
        new_label(ind(kk)) = {char(strt+kk-1)}; 
    end
end
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
    sz = size(csi.data.raw);
    [~, sort_ind] = sort(sz(2:end), 'descend');
    new_order = [1 sort_ind+1];
else
   dim_order = 1:numel(size(csi.data.raw));
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

if nargin < 4, backup = gui.checkbox_backup.Value; end


% BACKUP + APPDATA % ------------------------------- %
% Create backup
if backup, CSI_backupSet(gui, 'Before applying bins.'); end

% Check if csi app data exists
if ~isappdata(gui.CSIgui_main, 'csi'), return; end


% User Input % -------- %
% Dimension to bin over and number of bins?
uans = getUserInput({'Dimension to bin: ','Number of bins: '},{1,2});
if isempty(uans),CSI_Log({'Skipped binning.'},{''}); return; end
dim = str2double(uans{1}); nbin = str2double(uans{2});


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
bin = NaN([sz_temp(:)', nbin]);
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
function button_CSI_Sum_Callback(hobj, ~, gui, backup)
% Summate MRSI data over a specific dimensions

% BACKUP + APPDATA % ------------------- %
if nargin < 3, gui = guidata(hobj); end
if nargin < 4, backup = gui.checkbox_backup.Value; end
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
function button_T1_MRS_Callback(~, ~, gui)
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
function button_CSI_T2_Callback(~, ~, gui)
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
uans = getUserInput_Popup({'Plot Images: '},{{'Yes','No'}},...
    [], 'T2-map');
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
        {{'Yes', 'No'}}, [] , 'T2-map');
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
            doi_echo = CSI_getDataAtPeak_Stored(doi_xaxis);
            
            % ECHO: Convert echo data and calculate maximum
            data_echo = CSI_getUnit(doi_echo, data_unit);
            max_echo = max(data_echo,[],1);
            
        case 'Echo'
            
            % ECHO: Convert doi e.g. data at peak of interest
            data_echo = CSI_getUnit(doi, data_unit);
            max_echo = max(data_echo,[],1); 
            
            % FID: get peak in FID data
            doi_fid = CSI_getDataAtPeak_Stored(doi_xaxis);
            
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
str_time = string(datetime('now','Format','yyMMdd_HHmmss'));
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
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % ------------------- %

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'freq')
    CSI_Log({'MRS data is in frequency domain; '},...
{'Converting to time domain and back, use the backup (ctrl+z) to undo.'});
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before zero filling.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USER INPUT % ------------------------- %

uans = getUserInput({'Requested #samples after zero filling: '},...
                    {size(csi.data.raw,1).*2}, [], 'Zerofill');
if isempty(uans), CSI_Log({'Skipped zero filling.'},{''}) ; return; end

% Convert to length after zero filling.
N = str2double(uans{1});

% Get zerofill direction from FID/Echo e.g post/both
uans = getUserInput_Popup({'FID or Echo: (Post/Both)'},{{'FID','Echo'}},...
    [], 'Zerofill');
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

% Update info
CSI_Log({'Applied zero filling. Sample size increased to'},{N});

% --- Executes on button press in button_CSI_AutoPhase.
function button_CSI_AutoPhase_Callback(~, ~, gui,backup)
% Apply zero order phase correction to all voxels in the data set.
if nargin < 4, backup = gui.checkbox_backup.Value; end


% USER INPUT % ---------------------------------- %

% Zeroth or First order corrections:
uans = getUserInput_Popup({'Phase corrections for: '},...
                          {{'Zeroth-order', 'First-order', 'Both'}},[],...
                          'Phase Corrections');
if isempty(uans), CSI_Log({'Aborted phase-corrections.'},{''}); return; end


switch uans{1}
    case 'Zeroth-order'
        CSI_AutoPhase_Zero([],[], gui, backup);
    case 'First-order'
        CSI_AutoPhase_First([], [], gui, backup);
    case 'Both'
        CSI_AutoPhase_Zero([],[], gui, backup);
        CSI_AutoPhase_First([],[], gui, backup);
end

% --- Executes via button_CSI_AutoPhase.
function CSI_AutoPhase_First(~, ~, gui, backup)
% Automatic first-order phasing
if nargin < 4, backup = gui.checkbox_backup.Value; end

% BACKUP + APPDATA % ------------------------------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before first-order auto-phase corrections.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% USER INPUT % ------------------------------------------------------- %


uans = getUserInput_Tick(...
    {'Apply additional zeroth-order corrections afterwards: ',...
     'Accuracy: ', 'Parallel Computing:'},...
    {{'Yes', 'No'},...
     {'Default','Custom'}, {'Yes', 'No'}},...
     'First-order phase corrections',[1, 1, 0]);
if isempty(uans)
    CSI_Log({'Aborted first-order phase corrections.'},{''});
    return; 
end

if uans(2), acc = 4e-2; 
else
    uans_acc = getUserInput(...
        {'First-order phase correction accuracy:'},{'3e-2'});
    if isempty(uans_acc)
        CSI_Log({'Aborted first-order phase corrections.'},{''});
        return;
    end
    acc = str2double(uans_acc{1}); 
end
do_zero = uans(1);
do_parallel = uans(3);


% CORRECT DATA-DOMAIN (I/II) % ----------------------------------------- %

% Correct data domain (I/II)
domain = CSI_getDomain(gui);  

% Correct domain if possible
if strcmp(domain, 'time')
    CSI_Log({[ 'MRS data is in time domain; Autophase method requires ',...
   'frequency domain data. Converting to time domain and back, use backup',...
   '(ctrl+z) to undo.'] },{''});
    % FFT from time to frequency domain.
    csi.data.raw = csi_fft(csi.data.raw); 
end

% AUTO-PHASING % ------------------------------------------------------ %

% MRSI data to cell format
sz = size(csi.data.raw); 
cell_layout = arrayfun(@ones,...
    ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
cell_mrsi = mat2cell(csi.data.raw, sz(1), cell_layout{:});

CSI_Log({'Calculating first-order phase corrections.'},{''});

% Use cell-fun
if ~do_parallel
    
    tic
    % Apply auto zerophase to each cell
    [cell_mrsi, ~, ~] = ...
    cellfun(@csi_autoFirstPhase, ...
                cell_mrsi, ...                              % data
                repmat({do_zero}, size(cell_mrsi)),...      % method
                repmat({acc},      size(cell_mrsi)),...     % accuracy
                repmat({0},size(cell_mrsi)),...             % display
                'UniformOutput', 0);  
    dt = toc;
    
    % MRSI data to array        
    csi.data.raw = cell2mat(cell_mrsi);
end

% Use parallel computing
if do_parallel
    
    % Reshape data to {nDim x nChan} x nVox
    sz = size(cell_mrsi);
    cell_mrsi = reshape(cell_mrsi,[],1);

    p = gcp('nocreate'); 
    if isempty(p)
        nCores = feature('numcores');
        parpool('Processes', nCores);
    end

    tic
    parfor vi = 1:size(cell_mrsi,1)
        cell_mrsi{vi} = csi_autoFirstPhase(cell_mrsi{vi},...
                                           do_zero, acc, 0);
    end
    dt = toc;

    % Undo list-shape
    cell_mrsi = reshape(cell_mrsi, sz);
    
    % MRSI data to array        
    csi.data.raw = cell2mat(cell_mrsi);
end

% CORRECT DATA-DOMAIN (II/II) % ---------------------------------------- %

% Convert to starting data domain (II/II)
if strcmp(domain,'time'), csi.data.raw = csi_ifft(csi.data.raw); end

% SAVE AND CLEANUP % -------------------------------------------------- %

% Save
setappdata(gui.CSIgui_main,'csi',csi);

% Update info to user.
if do_zero, uans = {'Yes'}; else, uans = {'No'}; end
CSI_Log({'Applied automatic first-order phase corrections. Accuracy: '},{acc});
CSI_Log({'Included additional zeroth-order phase corrections:'}, uans);
CSI_Log({'Elapsed time:'},{dt});

% --- Executes via button_CSI_AutoPhase.
function CSI_AutoPhase_Zero(~, ~, gui, backup)
% Automatic zeroth-order phasing

% BACKUP + APPDATA % ---------------------------- %

% Create backup
if backup, CSI_backupSet(gui, 'Before zeroth-order auto-phase correction.'); end

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');


% USER INPUT % ---------------------------------- %

% Get peak of interest.
[poi] = CSI_getPeakOfInterest(csi.xaxis, 'Automatic zeroth-order phasing');
if isempty(poi), return; end

% Get method of auto-phasing
qry = {'Automatic zeroth-order phase correction method:'}; 
def = {'Match real part to maximum absolute signal.',...
       'Maximize real part of signal.',...
       'Correct for TE delay in FID by shift.'}; 
uans = getUserInput_Popup(qry, {def}, [], 'Phase corrections');   
if isempty(uans), CSI_Log({'Skipped zero-phasing.'},{''}) ; return; end

switch uans{1}
    case 'Correct for TE delay in FID by shift.',       phase_method = 3;
    case 'Match real part to maximum absolute signal.', phase_method = 2;
    case 'Maximize real part of signal.',               phase_method = 1;
end


% APPLY CORRECTION % ----------------------------- %

% POI from user.
if length(poi) > 1, poi = poi(1):poi(2); end

% Correct data domain (I/II)
domain = CSI_getDomain(gui);  % THIS DOES NOT UPDATE

% Correct domain if possible
if (phase_method == 3) && strcmp(domain,'freq')
    CSI_Log({['MRS data is in frequency domain; Autophase method requires '...
     'time domain data. Converting to time domain and back, use backup'...
     '(ctrl+z) to undo. ']},{''});
    csi.data.raw = csi_ifft(csi.data.raw); 
elseif strcmp(domain,'time') && (phase_method ~= 3)
    CSI_Log({[ 'MRS data is in time domain; Autophase method requires ',...
     'frequency domain data. Converting to time domain and back, use backup',...
     '(ctrl+z) to undo.'] },{''});
         
    csi.data.raw = csi_fft(csi.data.raw); 
end

% MRSI data to cell format
sz = size(csi.data.raw); 
cell_layout = arrayfun(@ones,...
    ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
cell_mrsi = mat2cell(csi.data.raw, sz(1), cell_layout{:});

if phase_method ~= 3
% Apply auto zerophase to each cell
[cell_mrsi_phased, ~] = ...
cellfun(@csi_autoZeroPhase, ...
            cell_mrsi, ...                              % data
            repmat({poi},    size(cell_mrsi)),...       % range
            repmat({phase_method},size(cell_mrsi)),...  % method
            repmat({0},      size(cell_mrsi)),...       % plot
            'UniformOutput', 0);            
else
    
    data = cell_mrsi;
    [~,ind] = cellfun(@max,data,'uniform', 0);
%     sz = size(data); sz(1) = size(data{1},1);    
    cell_mrsi_phased = cellfun(@CSI_shift_phasing,data,ind,'uniform',0); 

end

% MRSI data to array        
csi.data.raw = cell2mat(cell_mrsi_phased);

% Write to csi-struct
% csi.data.raw = array_mrsi;

% Convert to starting data domain (II/II)
if (phase_method == 3) && strcmp(domain,'freq')
    csi.data.raw = csi_fft(csi.data.raw); 
elseif strcmp(domain,'time') && (phase_method ~= 3) 
    csi.data.raw = csi_ifft(csi.data.raw); 
end

% CLEAN UP % ------------------------------------ %

% Save
setappdata(gui.CSIgui_main,'csi',csi);

% Update info to user.
CSI_Log({'Applied automatic zero-order phase corrections:'}, uans);

% --- Executes on button press phase-correction by shift of fid.
function new = CSI_shift_phasing(data, shift_index)
new = zeros(size(data,1),1);
new(1:(size(data,1) - shift_index+1))= data(shift_index:end); 

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

% --- Rotate ND-plane executed by CSI_Rotate_Callback
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
function button_CSI_Multiply_Callback(hObj, ~, gui, backup)
% Multiply the MRS data

if nargin < 4, backup = gui.checkbox_backup.Value; end

if ~exist('gui', 'var'), gui = guidata(hObj); end

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
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before dividing'); end


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

% --- Executes on button press in button_CSI_Normalize.
function button_CSI_Normalize_Callback(~ , ~, gui, backup)
% Normalize data to specific peak maximum: multiple methods available. See
% CSI_Normalize();
if nargin < 4, backup = gui.checkbox_backup.Value; end

% Create backup
if backup, CSI_backupSet(gui, 'Before normalization.'); end

% Normalize
CSI_Normalize(gui);

% --- Executes on button press in button_CSI_Noise_openFile.
function button_CSI_Noise_openFile_Callback(~, ~, gui)
% Load a separate data-file containing noise data.
%
%
% Noise is loaded into csi.data.noise with field raw, dim, labels and any
% possible header.


% USERINPUT % ---------------------------------------------------------- %

% Get default path if available
if isappdata(gui.CSIgui_main,'csi')
    csi = getappdata(gui.CSIgui_main,'csi');
    if isfield(csi,'filepath'), fp = csi.filepath; else, fp = []; end
else, CSI_Log({'No MRS data in memory.'},{'Load MRS-data first.'}); return; 
end
    
% Select-file UI.
[fn, fp, fi] = uigetfile({'*.*','All Files (*.*)';...
 '*.list;*.data;*.LIST;*.DATA',...
            'Raw MRS files Philips (*.list, *.data)';...
 '*.spar;*.SPAR;*.sdat;*.SDAT',...
            'MRS files Philips (*.sdat, *.spar)';...
 '*.dat;*.DAT;',...
            'Raw MRS file Siemens (*.dat)';...
 '*.txt;*.TXT;',...
            'Text files (*.txt)'},...
 'Select a noise-data file', fp, 'MultiSelect', 'off'); 
% Canceled file selection
if fi == 0, return; end 

% Get file-extension for further processing.
[fp, fn, ext] = fileparts([fp, fn]); ext = lower(ext);


% Load Data % ---------------------------------------------------------- %
switch ext
    case '.list'
        CSI_Log({'Noise from list-file not implemented.'},...
            {'Noise-data is usually present in main CSI data-file.'});
        return;
    case '.data'
        CSI_Log({'Noise from data-file not implemented.'},...
            {'Noise-data is usually present in main CSI data-file.'});
        return;
    case '.txt'
        CSI_Log({'Noise from text-file not implemented.'},{''});
    case '.dat'
        % Load data
        twix = mapVBVD([fp '\' fn ext]);
        
        % Save hdr
        csi.data.noise.twix = twix.hdr;

        % Read data dimension and labels
        dims = num2cell(twix.image.dataSize);
        dims_rng = cellfun(@(x) 1:x, dims,'uniform', 0); % Range
        
        % Dimension labels from header
        dims_txt = twix.image.dataDims;
        dims_txt = dims_txt(twix.image.dataSize>1);
        
        % Correct data labels from twix-file.
        % Labels in dat-file
        labels_dat_file      = {'col','cha', 'lin','par','seg', 'ave'}; 
        % Library for other name
        labels_match_library = {'fid','chan','ky','kz','kx', 'aver'};    
        dims_txt_corrected = cell(1,size(dims_txt,2));
        for ti = 1:size(dims_txt,2)        
            ind = strcmpi(labels_dat_file, dims_txt{ti});
            if sum(ind) > 0
                dims_txt_corrected{ti} = labels_match_library{ind};
            end
        end
        csi.data.noise.labels = dims_txt_corrected;

        % MRS-data
        tmp = squeeze(twix.image(dims_rng{:}));
        
        % Free up memory
        clear('twix');

        % Data memory size
        nfo = whos('tmp'); 
        matlab_nfo = memory;
        
        if nfo.bytes > matlab_nfo.MaxPossibleArrayBytes
            % Array is too large for matlab-array-memory
            fprintf('CSIgui: data requires much memory. Loading as single.\n')
            csi.data.noise.raw = tmp; 
        else
            % Array will fit matlab-max-array-memory
            csi.data.noise.raw = double(tmp);    
        end
        clear('tmp');
        
        % Noise dimensions
        csi.data.noise.dim = size(csi.data.noise.raw);
end


% Clean up % ----------------------------------------------------------- %
setappdata(gui.CSIgui_main,'csi', csi);

CSI_Log({'Loaded noise-data into memory:'},{[fn ext]});
CSI_Log({'Use noise-options in the MRSI menubar to view.'},{''});

% --- Executes on button press in button_CSI_PCA_Denoising.
function button_CSI_PCA_Denoising_Callback(hobj, ~, gui, backup)
% Apply PCA denoising according to:
% PCA denoising, Froeling et. al. doi.org/10.1002/mrm.28654
% Denoising, Veraart et. al. 10.1016/j.neuroimage.2016.08.016

if nargin < 4, backup = gui.checkbox_backup.Value; end

% Create backup
if backup, CSI_backupSet(gui, 'Before PCA Denoising.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');


% User Input % ------------------------------------------------------- %
psz_def = 3:2:19; psz_def([1 2]) = psz_def([2 1]);
qry = {'Apply noise decorrelation: ', ...
       'Calculate noise-covariance matrix for decorrelation using:', ...
       'PCA patch-size:', 'PCA SVD-method:'};
opt = {{'Yes','No'},  {'Measurement', 'Data'}, {psz_def},...
       {'Fast', 'Default'}};
uans = getUserInput_Popup(qry, opt, [], 'PCA-Denoising');
if isempty(uans), CSI_Log({'PCA Denoising cancelled'},{''}); return; end

% //Process user input % 
    
% uans{1} - Noise decorrelation
if strcmp(uans{1},'Yes'), do_NoiseDecorrelation = 1; 
else, do_NoiseDecorrelation = 0; end

% uans{2} - nCov using noise measurement or data
if strcmp(uans{2},'Measurement'), do_NoiseData = 1; 
else, do_NoiseData = 0; end 

% uans{3} - Patch size
patch_size = str2double(uans{3});

% uans{4} - SVD method for PCA denoising
if strcmp(uans{4},'Default'), do_svd = 0; else, do_svd = 1; end 


% Channel Index
chan_ind = csi_findDimLabel(csi.data.labels,{'chan'});

% Dimensions
sz = size(csi.data.raw); 

% Prepare Noise Data % ----------------------------------------------- %
if do_NoiseData
    if isfield(csi.data,'noise')
        % Run noise-prepare fcn
        [csi, ~, gui] = CSI_Noise_Prepare(hobj, gui);
        if ~isstruct(csi), CSI_Log({'PCA:'}, {'Aborted.'}); return; end

        % Message to user
        CSI_Log({'PCA:'},{'Noise-data processed and stored.'});  
    else
        do_NoiseData = 0; 
        CSI_Log({'PCA: No noise-data present,'},...
                {'using a noise-mask per voxel instead.'});
    end
end

% Decorrelate Noise % ------------------------------------------------ %
if do_NoiseDecorrelation 
    if isfield(csi.data, 'noise') && do_NoiseData
        
        % Get noise-cov using noise data
        noise_cov = CSI_NoiseCov_usingMeasurement(csi.data.noise.raw,...
                                              csi.data.noise.labels);

        % Copy it for nVoxels
        non_spat_ind = ...
            csi_findDimLabel(csi.data.labels,...
            {'chan','cha','fid','sec', 'avg', 'aver', 'nsa'});
        non_spat_ind(isnan(non_spat_ind)) = []; 
        non_spat_ind(non_spat_ind > numel(sz)) = [];
        sz(non_spat_ind) = [];

        % Covariance matrix
        noise_cov = repmat({noise_cov}, sz);
    
        % NFO Update
        CSI_Log({'PCA: Noise-data used to calculate'},...
                    {'noise-covariance matrix for noise decorrelation.'});
    else
        % Noise covariance matrix required for noise-decorrelation
        % calculated using noise in data itself.
        noise_cov = CSI_NoiseCov_usingData(csi.data.raw, chan_ind);
        CSI_Log({'PCA: Voxel-data used to calculate'},...
                {'noise-covariance matrices for noise decorrelation.'});
    end

    % Decorrelate signal (Cholensky)
    [csi.data.raw, ~] = ...
        csi_decorrelate_noise(csi.data.raw, chan_ind, noise_cov);
    CSI_Log({'PCA:'},{'Decorrelated noise via Cholesky-decomposition'});
end


% Reshape data % ---------------------------------------------------- %


% Cell-ify to: {nS x Kx x Ky x Kz x nChan} x Other
[csi.data.raw, permv] = ...
    CSI_PCA_Denoising_Reshape(csi.data.raw, csi.data.labels, chan_ind);
sz_cell_matrix = size(csi.data.raw);

% Reshape to list
csi.data.raw = reshape(csi.data.raw,[],1);

% Number of volumes to PCA-denoise
nVolumes = size(csi.data.raw,1);

% PCA Denoising % ---------------------------------------------------- %

% NFO update for user
CSI_Log({'PCA: Applying PCA-denoising. Patch size:'},{patch_size});
CSI_Log({'PCA: SVD-Method,'}, uans(4));
CSI_Log({'PCA:'},{'Starting a parallel pool of workers.'}); 

% Principle Component Analysis - Denoising
dt = 0;
for vi = 1:nVolumes
    tic
    csi.data.raw{vi} = ...
        csi_pca_denoising(csi.data.raw{vi}, chan_ind, patch_size, do_svd);
    dt = dt + toc;
end

% Processing-time output
 CSI_Log({'PCA: Denoising applied. Duration:'} , {dt} );

% Revert Reshape % --------------------------------------------------- %

% Undo cell array to list
csi.data.raw = reshape(csi.data.raw, sz_cell_matrix);
% Undo array to cell
csi.data.raw = cell2mat(csi.data.raw);
% Undo permute of dimensions
csi.data.raw = ipermute(csi.data.raw, permv);

% Clean Up % --------------------------------------------------------- %
setappdata(gui.CSIgui_main, 'csi', csi);

% Output
CSI_Log({'PCA: Clean up successful.'},{'Finished.'});

% --- Executes by CSI_PCA_Denoising
function [data, permv] = CSI_PCA_Denoising_Reshape(data, labels, chan_ind)
% Reshape data for PCA denoising such that dimensions represent:
% {nS x Kx x Ky x Kz x nChan} x Other Dimensions..
%
% Prepare data-structure
% PCA-denoising assumes a volume with nChannels and no other dimensions.
% Here the other dimensions are split.

% Size of volume
sz = size(data);

% Create [nS x Kx x Ky x Kz x nChan] x Other Dims...
spat_labels =  {'kx','ky', 'kz', 'x', 'y', 'z'};
spat_dim = csi_findDimLabel(labels, spat_labels);
spat_dim(isnan(spat_dim)) = [];

dim_order = 1:numel(sz);
% Set the new order with kspace dimensions first
permv = [1 spat_dim chan_ind];
% Find excluded dimensions
excl_dim_order = ismember(dim_order, permv);
% Append new_order.
permv = [permv dim_order(excl_dim_order == 0)];

% Reorder using new_order.
data   = permute(data, permv);

% Rearrange data to cell 
sz = size(data); cell_layout = {1};
if numel(sz) > 5
    cell_layout = arrayfun(@ones, ones(1,size(sz(6:end),2)),sz(6:end),...
        'UniformOutput',0);
end   
if numel(sz) < 5
    cell_layout = num2cell(sz);
else
    cell_layout = [num2cell(sz(1:5)) cell_layout];
end
% Array 2 Cell: {nS x Kx x Ky x Kz x nChan} x Other
data = mat2cell(data, cell_layout{:}); 

% --- Executes on button press in button_CSI_Fieldmap.
function button_CSI_Fieldmap_Callback(~, ~, gui)
% Calculate and show fieldmaps as used by Roemer weighting method.

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Preperation of data % ---------------------------------------------- %

% Channel Index
chan_ind = csi_findDimLabel(csi.data.labels,{'chan'});

% iFFT to FID
if strcmp(CSI_getDomain(gui), 'freq')
    CSI_Log({'Fieldmaps:'},{'Converting data to time domain (iFFT).'});
    fid = csi_ifft(csi.data.raw);
else
    fid = csi.data.raw;
end

% Sensitivity Maps % -------------------------------------------------- %
[sens_maps, permv, szr]  = CSI_Sensitivity_Maps(fid, chan_ind);
CSI_Log({'Fieldmaps'},{'calculated.'});

% Reshape the data
szr(1) = 1; szr(2) = 1;  sens_maps = reshape(sens_maps, szr); 

% Undo cell
sens_maps = cell2mat(sens_maps);

% Create undo-permute-permute-vector
nindex = numel(permv); restore_permv = NaN(1,nindex);
for kk = 1:nindex, restore_permv(kk) = find(kk == permv); end

% Permute the data
sens_maps = permute(sens_maps, restore_permv);

% Display % ------------------------------------------------------------ %

% User input
qry = {'Display data:', 'Data Type', 'Save data'};
opt = {{'Yes','No'},  {'Real','Imaginary', 'Absolute'}, {'No','Yes'}};
uans = getUserInput_Popup(qry, opt, [], 'Display Fieldmaps');
if isempty(uans), CSI_Log({'Skipped field map display.'},{''}); return; end

if strcmp(uans{1}, 'Yes')
    % Output
    CSI_Log({'Displaying fieldmaps.'},{''});

    if strcmp(uans{2}, 'Real')
        CSI_dataAs_Initiate(real(sens_maps), 'Fieldmaps (Real)', gui,...
        csi.data.labels)
    elseif strcmp(uans{2}, 'Imaginary')
        CSI_dataAs_Initiate(imag(sens_maps), 'Fieldmaps (Imaginary)', gui,...
        csi.data.labels)
    elseif strcmp(uans{2}, 'Absolute')
        CSI_dataAs_Initiate(abs(sens_maps), 'Fieldmaps (Abs)', gui,...
        csi.data.labels)
    end
end

% Save data
if strcmp(uans{3}, 'Yes')
    if isfield(csi,'filepath'), fp = csi.filepath;
    else, fp = [];
    end
       
    % Get file path and extension from user
    [fn, fp, fi] = uiputfile({'*.mat', 'Matlab file'},...
                             'Save fieldmap data...',fp);
    if fi == 0, return; end
    [~,fn, ext] = fileparts(fn);

    % Save file
    save([fp '\' fn ext], 'sens_maps');
end

% --- Opens on button press of [i]nformation app-data
function CSI_openDataDirectory(gui)
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

if isfield(csi,'filepath')
    winopen(csi.filepath)
end

% --- Executes on button press in button_CSI_Single.
function button_CSI_Single_Callback(~, ~, gui)
% Conver csi.data.raw to single.

if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

if isa(csi.data.raw, 'double')
    csi.data.raw = single(csi.data.raw); msgA = 'single'; msgB = 'double';
else
    csi.data.raw = double(csi.data.raw); msgB = 'single'; msgA = 'double';
end

CSI_Log({['CSI data converted to ' msgA '-precision.']},...
        {['Press again to set to ' msgB '-precision.']});

setappdata(gui.CSIgui_main, 'csi', csi);

% --- Executes on button press in button_CSI_FlipOrientation.
function button_CSI_FlipOrientation_Callback(hobj, ~, gui)
% Rotate CSI and MRI space around an axis and look at different
% orientations.

warning('CSI_FlipOrientation_Callback is in beta and not functional yet.');

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before rotation'); end
gui = guidata(hobj);

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');



% --------------------------------------------------------- %

% Check for image-data
doMRI = 0;
if isappdata(gui.CSIgui_main, 'conv'), doMRI = 1; end


% Rotate z-slices into plane, with y-slices as result
% x = 1 y = 2 z = 3
rot = [3 2 1]; fldim = 3;

rotimg = [3 2 1];

% Spatial dimensions
spat_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
spat_dim(isnan(spat_dim)) = [];
if isempty(spat_dim), spat_dim = [2 3 4]; end

permv = 1:numel(size(csi.data.raw));
permv(spat_dim) =  spat_dim(rot);

% \\ CSI-data % --------------------- %

% Permute CSI-data
csi.data.raw = permute(csi.data.raw, permv);
if fldim > 0, csi.data.raw = flip(csi.data.raw, fldim); end

% Permute resolution
csi.ori.res = csi.ori.res(rot);

% Permute labels
% csi.data.labels(spat_dim) = csi.data.labels(spat_dim(rot));
csi.data.dim = size(csi.data.raw);


% Permute mesh-matrix AND x/y/z swapped according to rot.
x = csi.ori.mesh.x;
csi.ori.mesh.x = csi.ori.mesh.y; csi.ori.mesh.y = x;

mesh_str = fieldnames(csi.ori.mesh); mesh_str_rot = mesh_str(rot);
mesh = struct; 
for mi = 1:size(mesh_str,1)
    mesh.(mesh_str_rot{mi}) = permute(csi.ori.mesh.(mesh_str{mi}), rot);
end
csi.ori.mesh = mesh; clear('mesh');


% \\ MRI-data % --------------------- %

if doMRI

    
    % Interpolate    
    % button_MRI_Interpolate_Callback([], [], gui);
    
    gui = guidata(gui.CSIgui_main);

    % Retrieve data
    conv = getappdata(gui.CSIgui_main, 'conv');
    
       
    % Permute image-data
    conv.data = permute(conv.data, rotimg);
    conv.data = flip(conv.data, fldim);

    % Permute resolution
    conv.res = conv.res(rotimg);
    

    % Permute mesh-matrix AND x/y/z swapped according to rot.
    % 
    % Or recalculate mesh-matrix from volume vectors..?
    mesh_str = fieldnames(conv.mesh); mesh_str_rot = mesh_str(rotimg);
    mesh = struct;
    for mi = 1:size(mesh_str,1)
        mesh.(mesh_str_rot{mi}) = permute(conv.mesh.(mesh_str{mi}), rotimg);
    end
    conv.mesh = mesh; clear('mesh');


    % mri = getappdata(gui.CSIgui_main, 'mri')


% % Converted resolution equals original MR image resolution. (3D)
% % However, to fit correctly in csi grid, the resolution is changed (below).
% conv.res = mri.ori.res; % Initial... May change!!
% 
% % Calculate a resolution fitting the CSI space such that there is a integer
% % amount of image pixels fitted in each CSI direction of space.
% res_fit = csi.ori.res ./ conv.res;      % #MRpix / CSIpix
% res_rem = res_fit - floor(res_fit);     % Pixel change
% res_new = csi.ori.res ./ floor(res_fit);% New MRpix resolution 
% 
% % New resolution for each direction
% conv.res = res_new;
% 
% % Volume limits of CSI but with half a voxel distance for voxel limits e.g.
% % a total voxel (MRI res) vs the volume.
% conv.fov       = csi.ori.fov;     % Does not change regards to CSI
% conv.lim_vol   = csi.ori.lim_vol; % Volume of MRSI grid
% conv.lim(:,1)  = conv.lim_vol(:,1) + (0.5.*conv.res)'; % Voxel limits
% conv.lim(:,2)  = conv.lim_vol(:,2) - (0.5.*conv.res)'; % Used for coords
% 
% % Range of volume/grid of MRSI for MRI
% for kk = 1:size(conv.lim,1)    
%     % Number of pixels in kk-direction
%     N = (conv.fov(kk)./conv.res(kk));
% 
%     % Grid vector defining volume
%     conv.vec{kk} = linspace(conv.lim(kk,1), conv.lim(kk,2), N);
% end
% 
% % Image grid in MRSI space 
% % This grid lays in the CSI-space e.g. within limits of CSI FOV but is
% % sampled in x, y and z as close to the resolution of the image as possible
% [x,y,z] = meshgrid(conv.vec{2} ,conv.vec{1}, conv.vec{3});
% conv.mesh.x = x; conv.mesh.y = y; conv.mesh.z = z; 
% 
% % Interp values @ CSI space % --------------------------------- %
% % Interp3 is used to convert images to csi-grid in conv-struct.
% 
% % Check for stack availability - these images need seperate inclusion!
% if isfield(mri.ori,'stack_of_interest')
%     stack_ind = mri.ori.stack_of_interest;
% else
%     stack_ind = 1; 
% end
% 
% if strcmp(mri.ext ,'par') || strcmp(mri.ext ,'ima')   
%     image_convert_data = mri.data.(imtoi)(:,:,:);
% else
%     image_convert_data = mri.data.(imtoi)(:,:,stack_ind,:);
% end
% 
% % Interp
% conv.data = interp3(mri.ori.mesh.x, mri.ori.mesh.y, mri.ori.mesh.z,...      % Original MRI coordinates
%                     squeeze(image_convert_data),...                         % Original MRI values
%                     conv.mesh.x, conv.mesh.y, conv.mesh.z,'Linear',0);      % Requested coordinates
% conv.dim  = size(conv.data);


    msg = 'MRI-data included in transformation.';
else
    msg = 'MRI-data not available.';
end

% --------------------------------------------------------- %

% Store data %
setappdata(gui.CSIgui_main,'csi', csi);
if doMRI, setappdata(gui.CSIgui_main,'conv', conv); end

% Update x-axis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show nfo
CSI_Log({'Changed orientation of data.'},{msg});

% --- Executes on button press in button_MRI_Interpolate.
function button_MRI_Interpolate_Callback(~, ~, gui)
%  Interpolate the slice direction of the converted data set.

% Get appdata
if ~isappdata(gui.CSIgui_main, 'conv'), return; end
conv = getappdata(gui.CSIgui_main, 'conv');

% --------------------------------------------------------- %

% Get userinput: Resolution
uans = getInput({'edit'}, {'New resolution for slice-dimensions:'}, ...
    {5}, 'Interpolate MRI');
if isempty(uans), return; end
res = str2double(uans{1});

% Interpolate    
N = round((conv.lim(3,2) - conv.lim(3,1)) ./ res);
[y,x,z] = ndgrid(conv.vec{1}, conv.vec{2}, ...
    linspace(conv.lim(3,1),conv.lim(3,2), N));
img = interp3(conv.mesh.x, conv.mesh.y, conv.mesh.z, conv.data, x, y, z,...
              "spline");

conv.mesh.x = x; conv.mesh.y = y; conv.mesh.z = z; 
conv.data = img; conv.res(3) = res;

% --------------------------------------------------------- %

% Store data %
setappdata(gui.CSIgui_main,'conv', conv);

% Show nfo
CSI_Log({'Interpolated the converted MRI-data. Resolution:'},{res});

% --- Executes on button press in button_CSI_FlipSpace.
function button_CSI_FlipSpace_Callback(~, ~, gui)
% Flip data 180 degrees in all spatial directions.

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before flipping data.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Flip over k-space
ind = csi_findDimLabel(csi.data.labels,{'ky','kz'});
for kk = 1:size(ind,2)
    csi.data.raw = flip(csi.data.raw,ind(kk));
end

% Store data
setappdata(gui.CSIgui_main,'csi', csi);

% Log
CSI_Log({'Corrected orientation by flipping dimensions: '},...
        {strjoin(csi.data.labels(ind),' | ')});

% --- Executes on button press in button_CSI_Interpolate.
function button_CSI_Interpolate_Callback(~, ~, gui)
% Interpolate data over the spatial dimensions

% BACKUP + APPDATA % ------------------------------- %
% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before interpolation.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% k-space indexes
kspace = csi_findDimLabel(csi.data.labels,{'kx','ky','kz'});

% Userinput % --------------------------------------------- %

% Default 
if isfield(csi,'twix')
    % Get interpolation matrix from dat-file header.
    
    % How to figure out what the phase/read direction in data is?
    % csigui.twix.Meas....
    % FinalMatrixSizePhase FinalMatrixSizeRead FinalMatrixSizeSlice
    
    % From twix-header par = kz; lin = ky; seg = kx;
    if isfield(csi.twix.Meas,'RawSeg')
        def_ans(1) = csi.twix.Meas.RawSeg;
    end
    if isfield(csi.twix.Meas,'RawLin')
        def_ans(2) = csi.twix.Meas.RawLin;
    end
    if isfield(csi.twix.Meas,'RawPar')
        def_ans(3) = csi.twix.Meas.RawPar;
    end
else
    % Set as double of current matrix size
    def_ans = csi.data.dim(kspace) .* 2;
end

% Ask user
uans = getUserInput({'New interpolation dimensions:'},{def_ans});
if isempty(uans),CSI_Log({'Aborted interpolation.'},{''}); return; end
new_dim = str2double(strsplit(uans{1}, ' '));

uans = getUserInput_Popup({'Interpolation method: '},...
    {{'Spline', 'Cubic', 'Linear', 'Nearest', 'Makima'}},...
    [], 'Spatial Interpolation');
if isempty(uans{1}),CSI_Log({'Aborted interpolation.'},{''}); return; end
meth = uans{1};

% Interpolate --------------------------------------------- %
% Data will be automatically rearranged if kspace indexes are not after the
% sample-index (1).
CSI_Log({'Interpolation started, please be patient.'},{''});

csi.data.raw = csi_interpolate(csi.data.raw, new_dim, kspace, meth);
csi.data.dim = size(csi.data.raw);

% Store data       
setappdata(gui.CSIgui_main,'csi', csi);

% Log
CSI_Log({'Applied interpolation over all spatial dimensions: '},...
        {strjoin(csi.data.labels(kspace),' | ')});

% --- Executes on button press in button_CSI_Maps.
function button_CSI_Maps_Callback(~, ~, gui)
% Calculate specific maps to visualize data
%
% Includes: SNR, FWHM, peak, linewidth, and maximum maps.

% Check appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end

% Ask user what to map
qry = 'Choose a map to calculate: ';
opts = {'SNR','Linewidth','Maximum','Peak', 'Ratio'};

maptype = getUserInput_Popup({qry},{opts}, [], 'Mapping');
if isempty(maptype), return; end

switch maptype{1}
    case 'SNR'
        button_CSI_SNR_Callback([], [], gui);
    case 'Linewidth'
        button_CSI_Linewidth_Callback([], [], gui);
    case 'Maximum'
        button_CSI_MaxValue_Callback([], [], gui);
    case 'Peak'
        button_CSI_Peak_Map_Callback([], [], gui);
    case 'Ratio'
        button_CSI_Peak_Ratio_Map_Callback([], [], gui);
end

% --- Executes via map-selection to calculate peak-ratio maps
function button_CSI_Peak_Ratio_Map_Callback(~, ~, gui)
% Calculate a ratio-map from two given peaks and display it.
% Check if csi appdata is present

% CSI data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get peak-data for two-peaks of interest.
map = cell(1,2); doi_range = cell(1,2);
for kk = 1:2
    % Get peak of interest
    [doi, ~, doi_range{kk}] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
    doir = real(doi);
        
    % Maximum positions and values: e.g. a map
    map{kk} = max(doir, [], 1);    
end

% Calculate ratio.
nmap = map{1}./map{2};
% Normalize map
% nfac =  max(map(:)); nmap = map./nfac;

     
% Filter by SNR
nmap = CSI_dataAs_SNRfilter(nmap, 'peak-ratio', gui, doi_range{1});

% \\ Display Data
CSI_dataAs_Initiate(nmap, 'Peak-Ratio', gui);



% Coil Combination Methods % ---------------------------------------- %
% ------------------------------------------------------------------- %


% --- Executes on button press in button_CSI_Combine.
function button_CSI_Combine_Callback(hobj, ~, gui)
% Combine coil data with a variety of different combination methods: WSVD,
% Roemer, Manual, SNR weighted and more.
%
% This function acts only as a gateway to the script of the other
% functions.

% BACKUP + APPDATA % --------------------- %

% Get app data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end

% Get user input for coil-combination method.
methods = {'WSVD', 'Roemer', 'SNR Weighted', 'Manual', 'WSVD - Static',...
           'Sum-of-Squares'};
uans = getUserInput_Popup({'Coil weighting method:'},{methods},...
                          [], 'Coil Weighting');
if isempty(uans), return; end


switch uans{1}
    case 'Roemer'
        CSI_Combine_Roemer(hobj)
    case 'WSVD'
        CSI_Log({'WSVD Source:'},{'http://dx.doi.org/10.1002/mrm.22230'});
        CSI_Combine_WSVD(hobj);
    case 'SNR Weighted'
        CSI_Combine_SNRweighted(hobj)
    case 'Manual'
        CSI_Combine_Manual(hobj)
    case 'WSVD - Static'
        CSI_Combine_WSVD_StaticWeights(hobj)
    case 'Sum-of-Squares'
        CSI_Combine_SumOfSquares(hobj)
end

% --- Executeds on Roemer-method selection via CSI_Combine
% --- Executes if user selects Roemer in Combine
function CSI_Combine_Roemer(hobj,~)
% Combine data using the Roemer et al method described in:
% https://doi.org/10.1002/mrm.1910160203 - The NRM Phases Array 1990
%
% This method' quality improves with separate noise-measurements.
%
% Required Variables:
% FID of data, noise-covariance matrix, sensitivity maps.
%
% Optional: Noise decorrelation, PCA-denoising, noise covariance matrix
% from noise-measurements or data itself.
%
% NB. If noise decorrelation and PCA-denoising are applied, the Roemer
% algorithm can be calculated using the identity matrix instead of the
% noise covariance matrix.


% Get GUI-handle
gui = guidata(hobj); 

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before combining channels (Roemer).'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% User Input % ------------------------------------------------------- %

% Define questions and options
qry = {'Apply noise decorrelation: ', ...
       'Calculate noise-covariance matrix using:', ...
       'Apply PCA denoising:',...
       'Noise-covariance for Roemer:'};
opt = {{'Yes','No'},  {'Measurement', 'Data'}, {'Yes','No'},...
       {'Default','Cholesky', 'Identity-matrix'}};

% Ask user
uans = getUserInput_Popup(qry,opt, [], 'Roemer');
if isempty(uans), CSI_Log({'Roemer:'},{'Aborted.'}); return; end

% Process user input % 
    
% uans{1} - Noise decorrelation
if strcmp(uans{1},'Yes'), do_NoiseDecorrelation = 1; 
else, do_NoiseDecorrelation = 0; end

% uans{2} - nCov using noise measurement or data
if strcmp(uans{2},'Measurement'), do_NoiseData = 1; 
else, do_NoiseData = 0; end 

% uans{3} - PCA denoising
if strcmp(uans{3},'Yes'), do_PCA = 1; else, do_PCA = 0; end

% uans{4} - Noise Covariance Matrix
if strcmpi(uans{4},'Default')
    do_nCovDef = 1; 
elseif strcmpi(uans{4}, 'Cholesky')
    do_nCovDef = 2; 
    do_NoiseDecorrelation = 1;
    CSI_Log({'Roemer: Using decorrelated noise-cov for coil combinations'},...
        {'enabled noise decorrelation by default.'});
elseif strcmpi(uans{4}, 'Identity-matrix')
    do_nCovDef = 0;
end


% -------------------------------------------------------------------- %
% Prepare data ------------------------------------------------------- %
CSI_Log({'Roemer:'},{'Starting preparations and calculations.'});

% Channel Index
chan_ind = csi_findDimLabel(csi.data.labels,{'chan'});

% Dimensions
sz = size(csi.data.raw); 

% Prepare Noise Data % ----------------------------------------------- %
if do_NoiseData
    if isfield(csi.data,'noise')
        % Run noise-prepare fcn
        [csi, ~, gui] = CSI_Noise_Prepare(hobj, gui);
        if ~isstruct(csi), CSI_Log({'Roemer:'}, {'Aborted.'}); return; end

        % Message to user
        CSI_Log({'Roemer:'},{'Noise-data processed and stored.'});  
    else
        do_NoiseData = 0; 
        CSI_Log({'Roemer: No noise-data present,'},...
                {'using a noise-mask per voxel instead.'});
    end
end

% Decorrelate Noise % ------------------------------------------------ %
if do_NoiseDecorrelation 
    if isfield(csi.data, 'noise') && do_NoiseData
        
        % Get noise-cov using noise data
        noise_cov = CSI_NoiseCov_usingMeasurement(csi.data.noise.raw,...
                                              csi.data.noise.labels);

        % Copy it for nVoxels
        non_spat_ind = ...
            csi_findDimLabel(csi.data.labels,...
            {'chan','cha','fid','sec', 'avg', 'aver', 'nsa'});
        non_spat_ind(isnan(non_spat_ind)) = []; 
        non_spat_ind(non_spat_ind > numel(sz)) = [];
        sz(non_spat_ind) = [];

        if isempty(sz), sz = 1; end

        % Covariance matrix
        noise_cov = repmat({noise_cov}, sz);
    
        % NFO Update
        msg = ...
        'Used noise-data to calculate noise-cov matrix for decorrelation.';
    else
        % Noise covariance matrix required for noise-decorrelation
        % calculated using noise in data itself.
        noise_cov = CSI_NoiseCov_usingData(csi.data.raw, chan_ind);
        msg = ...
        'Used voxel-data to calculate noise-cov matrix for decorrelation.';
    end
    CSI_Log({'Roemer:'},{msg});

    % Decorrelate signal (Cholensky)
    [csi.data.raw, noise_cov_chol] = ...
        csi_decorrelate_noise(csi.data.raw, chan_ind, noise_cov);
    CSI_Log({'Roemer:'},{'Decorrelated noise via Cholesky-decomposition'});
end


% PCA Denoising % ---------------------------------------------------- %
if do_PCA    
    % Run PCA
    button_CSI_PCA_Denoising_Callback(hobj, [], gui, 0);
    gui = guidata(hobj);
    CSI_Log({'Roemer:'},{'PCA Denoising applied, continuing.'});
end


% Noise Cov Matrix % ------------------------------------------------- %
% For Roemer method, if noise decorrelation and PCA, the ID-matrix can 
% be used.
if do_nCovDef == 1 % Use default noise-covariance matrix

    % Use noise-covariance matrix from data or noise-measurement
    % NB. noise-cov is present if noise decorrelation is applied.
    if ~do_NoiseDecorrelation 
        if do_NoiseData
            noise_cov = CSI_NoiseCov_usingMeasurement(...
                csi.data.noise.raw, csi.data.noise.labels);
            noise_cov = repmat({noise_cov}, sz);
            CSI_Log({'Roemer: Using noise-cov matrix from noise-data'},...
            {'for coil combinations.'});
        else
            noise_cov = CSI_NoiseCov_usingData(csi.data.raw, chan_ind);
            CSI_Log({'Roemer: Using noise-cov from voxel-data'},...
                {'for coil combinations'});
        end
    end

elseif do_nCovDef == 0  % Use identity matrix
    nVox = prod(csi.data.dim(2:end));
    noise_cov = ...
        repmat({diag(ones(csi.data.dim(chan_ind),1))},nVox,1); 
    CSI_Log({'Roemer:'},{'Using identity matrix for coil-combinations.'});

elseif do_nCovDef == 2
    % Use cholesky decomposed noise covariance matrix.
    noise_cov = noise_cov_chol;
    CSI_Log({'Roemer:'},...
    {'Using noise-cov matrix from decorrelation for coil-combinations.'});
end


% FID of data as cell-list {nS x nChan} x nVox
[csi.data.raw, permv, szr] = csi_combine_reshape(csi.data.raw, chan_ind);

% iFFT to FID
CSI_Log({'Roemer:'},{'Converting data to time domain (iFFT).'});
fid = csi_ifft(csi.data.raw);

% Sensitivity Maps
sens_maps = CSI_Sensitivity_Maps(fid, chan_ind);
CSI_Log({'Roemer:'},{'Generated sensitivity maps.'});


% Reference channel % ------------------------------------------------ %
refChannel = CSI_Combine_Roemer_refChannel(fid);
CSI_Log({'Roemer: Reference channel'},{refChannel});

% Roemer Method % ---------------------------------------------------- %
CSI_Log({'Roemer:'},{'Calculating weights for each channel.'});

% Loop each voxel
for vi = 1:size(fid,1)
    % Sensitivity for each channel
    S = sens_maps{vi}'; 
    % Noise Cov Matrix or ID-matrix if (pca)-denoised.
    N = noise_cov{vi};
    % Moore-Penrose Pseudo-Inverse (uses SVD)
    U = pinv(sqrt(S'*pinv(N)*S))*S'*pinv(N); % The magic is here.
    
    % Make an arbitrary value to scale to, get back into same magnitude.
    % rescale = norm(U) * (U / norm(U(refChannel)));
    % fid{vi} = ((U .* rescale) * fid{vi}' )';

    % Weight the FID    
    fid{vi} = (U * fid{vi}')';
end
CSI_Log({'Roemer:'},{'Combined all channels using calculated weights.'});

% Reshape to original and store data
fid = csi_combine_reshape_revert(cell2mat(fid'), permv, szr);
CSI_Log({'Roemer:'},{'Converting data to frequency-domain (FFT).'});
csi.data.raw = csi_fft(fid);
csi.data.dim = size(csi.data.raw);

% Clean Up % --------------------------------------------------------- %
setappdata(gui.CSIgui_main, 'csi', csi);

% Update some calculations
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Output
CSI_Log({'Roemer: Done.'},{''});

% --- Reference channel for channel weighting-calculations in combine
function refChannel = CSI_Combine_Roemer_refChannel(fid)
% Find the channel with the most signal and return the corresponding
% channel index.
%
% This function simply looks at the maximum value of the fid for every
% voxel and every channel.

% Real values
fid_real = cellfun(@real, fid, 'UniformOutput', false);

% Maximum per channel per voxel
mx_val_perchan_pervox = cellfun(@max, fid_real, 'UniformOutput', false);

% Maximum value per voxel with its channel-number
[mx_val_pervox, mx_val_pervox_chanind] = ...
    cellfun(@max, mx_val_perchan_pervox, 'UniformOutput', false);

% Maximum from all voxels with its voxel-number
[~, mx_val_voi] = max(cell2mat(mx_val_pervox));

% Reference channel
% Thus: max(fid_real{mx_val_voi}(:,refChannel)) == mx_val
refChannel = mx_val_pervox_chanind{mx_val_voi};

% --- Sensitivity maps for Roemer-combination method
function [sens_maps, permv, szr] = CSI_Sensitivity_Maps(fid, chan_ind)
% Calculate sensitivity maps as the mean value of the FID for the 2 to 5th
% sample.
%
% Returned are the maps and the permute vector and matrix size for
% rearranging the data such that the channel dimensions is on the second
% index.

% Data to cell
if ~iscell(fid)
    [fid, permv, szr] = csi_combine_reshape(fid, chan_ind);
end

% Calculate maps
sens_maps = cellfun(@(x) mean(x(2:5,:),1), fid, 'UniformOutput',0);

% --- Noise Covariance matrix using noise-data
function noise_cov = CSI_NoiseCov_usingMeasurement(noise, labels)
% Return the noise covariance matrix (nCov) using the noise measurement.
% Option to average or concatenate averages.

% Noise average index
chan_ind_avg = csi_findDimLabel(labels,{'avg', 'aver','nsa'});
chan_ind_avg(isnan(chan_ind_avg)) = [];

% Handle averages in noise-data (if present)
if ~isempty(chan_ind_avg) && ...
        size(noise,chan_ind_avg) > 1
    
    % Get userinput request to average of concatenate
    uans_noise = getUserInput_Popup(...
        {'Noise covariance matrix calculations, concatenate averages:'},...
        {{'Yes', 'No'}}, [], 'Noise Covariance');
    if isempty(uans_noise), noise_cov = NaN; return; end
    
    % Concatenate or average
    if strcmp('Yes', uans_noise{1}) % Concatenate noise nsa
        % Take noise
        sz = size(noise); nind = chan_ind_avg;
        
        % Permute avg-ind to second index.
        permv = 1:numel(sz); permv([2 nind]) = permv([nind 2]);
        noise = permute(noise, permv);
        
        % Reshape to catenate
        sz =  size(noise); nsz = num2cell(sz(3:end));
        noise = reshape(noise, sz(1) * sz(2),1 , nsz{:});
        
        % Permute
        noise = permute(noise, permv);

    else % Average noise nsa
        noise = mean(noise, chan_ind_avg);
    end
end

% Covariance matrix
noise_cov = cov(noise);

% --- Noise Covariance matrix using data
function noise_cov = CSI_NoiseCov_usingData(spec, chan_ind, noise_mask)
% Returns a noise covariance matrix with respect to each channel for 
% every voxel using the data itself.
%
% output noise_cov = {nChan x nChan} x Spatial Dimensions ...
%
% !! If no separate noise-data is available, use this function
% !! If no noise_mask given - it will be calculated as 2x(1/6)% of the
%                              start and end of the spectrum.


dim = size(spec);
if nargin < 3
    nS = dim(1); half_nm_size = round(nS./6);
    noise_mask = [1:half_nm_size nS - half_nm_size];
end

% Index for cutting data
cut_ind = arrayfun(@(x) 1:x, dim,'UniformOutput',0);
cut_ind{1} = noise_mask;
noise_data = spec(cut_ind{:});

% Reshape channel index: {nS x nChan} x spatial dimensions...
if isempty(chan_ind), chan_ind = numel(dim)+1; end

% Cell-ify and reshape data to {nS x nChan} x nVox
[noise_data, ~, szr] = csi_combine_reshape(noise_data, chan_ind);

% Noise covariance Matrix
noise_cov = cellfun(@cov, noise_data,'UniformOutput',0);
noise_cov = cellfun(@diag, noise_cov,'UniformOutput',0);
noise_cov = cellfun(@diag, noise_cov,'UniformOutput',0);

% Reshape to {noise_cov} x Spatial Dimensions
if numel(szr) >= 3 && (numel(dim)-2 > 2)
    noise_cov = reshape(noise_cov, szr(3:end)); 
end

% --- Executes if SNR weighted combinations is selected
function CSI_Combine_SNRweighted(hobj, ~)
% Combine channels using the SNR to weight the spectra from each coil.


% If exists, close the combine options menu
if exist('hobj', 'var')
    gui_combCoil = guidata(hobj);close(gui_combCoil.fig);
end

% Get GUI and object: CSIgui_main
CSIgui_obj = findobj('Tag', 'CSIgui_main'); gui = guidata(CSIgui_obj);

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before combining channels (SNR).'); end


% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');


% User Input -------------------------------------------------------- %

% Get user input for noise component

% Get peak of interest.
poi = CSI_getPeakOfInterest(csi.xaxis, 'SNR Weighted Combination');
if isempty(poi), return; end

% SNR Noise mask
uans = getUserInput({'Size of SNR noise mask: ','Exclude channels: '}, ...
    {round(csi.data.dim(1)/12), ''});
if isempty(uans)
    CSI_Log({'Skipped SNR Weighted Combination.'},{''}); return; 
end

% SNR noise mask
noiseMask = str2double(uans{1});

if ~isempty(uans{2})
    chan_excl = str2double(strsplit(uans{2}, ' '));
else
    chan_excl = [];
end

% Channel index % ----------------------------------------------------- %
ind_cha = csi_findDimLabel(csi.data.labels,{'chan', 'cha'});
ind_cha(isnan(ind_cha)) = [];

nChan = size(csi.data.raw,ind_cha);
ch_incl = 1:nChan;
ch_incl(chan_excl) = [];
nChan = size(ch_incl,2);

% Index-vector to exclude data with
sz_vec = arrayfun(@(x) 1:x, size(csi.data.raw), 'UniformOutput', 0);
sz_vec{ind_cha} = ch_incl;    

% Calculate SNR & Weights --------------------------------------------- %
snr = csi_SNR(csi.data.raw, noiseMask, 0, poi(1):poi(2));

% Maximum over channels and divide by max
snr_max_chan = max(snr,[], ind_cha);
snr_norm = snr./snr_max_chan;

% Apply weighting
csi.data.raw = csi.data.raw .* snr_norm;

% Mean
csi.data.raw = mean(csi.data.raw(sz_vec{:}),ind_cha) .* nChan;
csi.data.dim = size(csi.data.raw);

% --------------------------------------------------------- %

% Store data %
setappdata(gui.CSIgui_main,'csi', csi);

% Update x-axis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show nfo
CSI_Log({'Combined channels via SNR weighting, noise mask: '},{noiseMask});

% --- Executes if Static WSVD combinations is selected
function CSI_Combine_WSVD_StaticWeights(~, ~)
% 1. Calculate weights of one non-spatial dimensions for each voxel.
% 2. Apply to other non-spatial dimensions for each voxel.
% Uses WSVD to calculate weights.

% GUI prepwork % ------------------------------------------- %

% Get GUI and object: CSIgui_main
CSIgui_obj = findobj('Tag', 'CSIgui_main'); gui = guidata(CSIgui_obj);

% Create a backup of the current data set.
backup = gui.checkbox_backup.Value; 
if backup,CSI_backupSet(gui, 'Before combing channels (WSVD - Static).');
end
gui = guidata(CSIgui_obj);

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

% User input % ----------------------------------------------- %

% Get user input for noise component
uans = getUserInput({'Use the noise prescans? (y/n):','Noise mask size:',...
                     'Exclude channels:'},...
                    {'n', round(csi.data.dim(1).*0.1),''}); 
if isempty(uans), CSI_Log({'Aborted WSVD.'},{''}) ; return; end
                

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

% Get user non-spatial dimension of interest input % ------------------- %
spat_dim = csi_findDimLabel(csi.data.labels, {'kx','ky','kz', 'x','y','z'});
spat_dim(isnan(spat_dim)) = [];
lab = csi.data.labels; lab([1 spat_dim]) = [];

uans_dim = ...
    getUserInput_Popup({'Dimension used to calculate weights: '},{lab},...
    [], 'WSVD - Static');
if isempty(uans_dim), CSI_Log({'Aborted WSVD.'},{''}); return; end
dim_of_int = csi_findDimLabel(csi.data.labels,uans_dim);

% Value ...
uans_dim_val = getUserInput_Popup({'Index to use for weights: '},...
    {num2cell(1:csi.data.dim(dim_of_int))}, [], 'WSVD - Static');
if isempty(uans_dim_val), CSI_Log({'Aborted WSVD.'},{''}); return; end
dim_of_int_val = str2double(uans_dim_val{1});

% Number of samples and channels.
ndimf = size(csi.data.raw,1);

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


% Data of interest % ------------------------------------------------- %

% Here, the data needs to be cut for the dimension of interest.
dim = csi.data.dim;

% Cell with vector per dimension to access data in csi.data.raw.
sz_vec = arrayfun(@(x) 1:x, dim, 'UniformOutput', 0);
sz_vec{dim_of_int} = dim_of_int_val;

% Calculate weights using THIS data set
data = csi.data.raw(sz_vec{:});

% Reshaping (1/2) % -------------------------------------------------- %

% Reshape data to {nDim x nChan} x nVox
data = csi_combine_reshape(data, ind_cha);

% WSVD % ------------------------------------------------------------- %

% Number of voxels
nvox = size(data,1);

% Containers for WSVD
comb = struct; comb.data = zeros(dim(1),nvox); 
comb.W = zeros(size(ch_incl,2), nvox); 
% comb.qual = zeros(nvox,1); comb.ampl = zeros(nvox, size(ch_incl,2)); 

% WSVD loop. 
% Apply for every indices excluding the channel index: e.g. every voxel.
for vi = 1:nvox
    % Get spectrum and exclude channels given by user.
    tmp_spec = data{vi,1}; tmp_spec = tmp_spec(:,ch_incl);
    
    % Create noise coVariance matrix.
    if mask == 1 % Use mask
        
        noiseMask = ndimf-mask_size+1:ndimf;
        % noiseCov = cov(tmp_spec(noiseMask,:));
        noiseCov=diag(diag(cov(tmp_spec(noiseMask,:))));
           
    else         % Use noise pre-scans
        
        % FFT of noise prescans. (all channels)
        noise_spec = csi_fft(csi.data.noise.raw);
        noiseMask = 1:ndimf; noiseCov = cov(noise_spec(:,ch_incl));              
           
    end

    % WSVD algorithm
    % Data, quality, coil amplitude and weights.
    [comb.data(:,vi), ~, ~, comb.W(:,vi)] = ...
        wsvd(tmp_spec, noiseMask', 'noiseCov', noiseCov);
end
% comb.W will be used to combine the entire data-set. Previous calculations
% have shown that there is a small rounding-error (10e-12) for data from
% the wsvd-script and manually applying the weights. Therefor the
% calculated combination above is discarded.

% Sizes and output container
sz = size(csi.data.raw); sz(dim_of_int) = 1;  
sz_outp = size(csi.data.raw); sz_outp(ind_cha) = 1;
outp = NaN(sz_outp);
for di = 1:csi.data.dim(dim_of_int)
    
    % Grab data of interest at this iteration    
    sz_vec = arrayfun(@(x) 1:x, sz, 'UniformOutput', 0);
    sz_vec{dim_of_int} = di;    
    doi = csi.data.raw(sz_vec{:});  
 
    % Convert to {dt x nChan} x nVoxels
    [doi, permv, szr] = csi_combine_reshape(doi, ind_cha);

    % Loop each voxel
    doi_out = NaN(size(csi.data.raw,1), size(doi,1));
    for vi = 1:size(doi,1)
           % Get spectrum and exclude channels given by user.
            spec = doi{vi,1}; 
            spec = spec(:,ch_incl);           
        
            % WSVD algorithm
            % Data, quality, coil amplitude and weights.
            doi_out(:,vi) = spec * comb.W(:,vi);
    end

    % Convert to original matrix-shape and size.
    ind4outp = sz_vec; ind4outp{ind_cha} = 1; ind4outp{dim_of_int} = di;
    outp(ind4outp{:}) = ...
        csi_combine_reshape_revert(doi_out, permv, szr);
end

% Set output
csi.data.raw = outp;
csi.data.dim = size(csi.data.raw);

% Store data %
setappdata(gui.CSIgui_main,'csi', csi);

% Update x-axis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show nfo
CSI_Log({'WSVD: data combined by applying weights from'},...
    {sprintf('%s(%i)', csi.data.labels{dim_of_int}, dim_of_int_val)});

% --- Executes if manual coil combination is selected.
function CSI_Combine_Manual(hobj, ~)
% Comine channels using a manual method. Include specific channels and
% summate the combined channels or calculate the mean.
%
% Requires prior user input. (String).
% userInp{1} = 'Channels to exclude'
% userInp{2} = 'Summation only'
%
% The above will be in this function --> See to do list Combine update.


% GUI prepwork % ------------------------------------------------------ %
 
% If exists, close the combine options menu
if exist('hobj', 'var')
    gui_combCoil = guidata(hobj);close(gui_combCoil.fig);
end

% Get GUI and object: CSIgui_main
CSIgui_obj = findobj('Tag', 'CSIgui_main'); gui = guidata(CSIgui_obj);

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
% Get CSI data
csi = getappdata(gui.CSIgui_main, 'csi');



% User input % ------------------------------------------------------ %

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

% Data prepwork % ------------------------------------------------------ %

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
csi.data.dim = size(csi.data.raw);

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

% --- Executes if WSVD coil combination is selected.
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
 
% Get GUI and object: CSIgui_main
gui = guidata(hobj);

% Create a backup of the current data set.
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before combining channels (WSVD).'); end
gui = guidata(hobj);

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Data preperations % -------------------------------------------------- %

% Is there channel/coil data
ind_cha = csi_findDimLabel(csi.data.labels, {'chan', 'cha'});
ind_cha(isnan(ind_cha)) = [];
if isempty(ind_cha)
    CSI_Log({'WSVD: Aborted.'},{'No "chan" or "cha" dimension found.'})
    return; 
end

% Number of samples and channels.
% ndimf = size(csi.data.raw,1);
nchan = size(csi.data.raw,ind_cha);

% Are there averages available:
ind_avg = csi_findDimLabel(csi.data.labels, {'aver', 'avg', 'nsa'});
ind_avg(isnan(ind_avg)) = [];
if isempty(ind_avg), avg_opt = {'No'}; else, avg_opt = {'Yes','No'}; end

% User Input % ------------------------------------------------------- %

% Define questions and options
qry = {'Calculate noise-covariance matrix using:', ...
       'Apply PCA denoising:','Average data:', 'Exclude channels:'};
opt = {{'Measurement', 'Data', 'ID-Matrix'}, {'Yes','No'}, avg_opt, ...
       {'No', 'Yes'}};

% Ask user
uans = getUserInput_Popup(qry, opt, [], 'WSVD');
if isempty(uans), CSI_Log({'WSVD:'},{'Aborted.'}); return; end

% Process user input % 

% uans{1} - nCov using noise measurement or data or ID-mat
if strcmp(uans{1},'Measurement'),       do_NoiseData = 1; 
elseif strcmp(uans{1}, 'ID-Matrix'),    do_NoiseData = 2; 
else,                                   do_NoiseData = 0; 
end 

% uans{2} - PCA denoising
if strcmp(uans{2},'Yes'), do_PCA = 1; else, do_PCA = 0; end

% uans{3} - Average data
if strcmp(uans{3},'No'), do_avg = 0; else, do_avg = 1; end 

% uans{4} - Exclude channels
if strcmp(uans{4},'No'), do_excl = 0; else, do_excl = 1; end 
              
% Get user input for displaying results
qry = {'Quality maps: ', 'Amplitudes table: ', 'Weights table: '};
opt = repmat({{'Yes', 'No'}}, size(qry));
dans = ones(1,size(qry,2));  
uans_disp = getUserInput_Tick(qry,opt,dans); 
if isempty(uans_disp), CSI_Log({'WSVD:'},{'Aborted.'}) ; return; end


% Exclusion of channels % ------------------------------------------- %

% Convert user input.
if do_excl
    % Get channels to exclude from user
    uans = getUserInput({'Exclude channels:'}, {''},[], 'WSVD');
    if isempty(uans), CSI_Log({'WSDV: Aborted.'},{''}); end

    % Channels to exclude as double
    ch_excl = str2double(strsplit(uans{1},' ')); 
    if isnan(ch_excl), ch_excl = []; end

    % Create included channels array
    ch_incl = 1:size(csi.data.raw, ind_cha); ch_incl(ch_excl) = [];
else
    ch_incl = 1:nchan;
end

% Prepare Noise Data % ----------------------------------------------- %
if do_NoiseData == 1
    if isfield(csi.data,'noise')
        % Run noise-prepare fcn
        [csi, ~, gui] = CSI_Noise_Prepare(hobj, gui);
        if ~isstruct(csi), CSI_Log({'WSVD:'}, {'Aborted.'}); return; end

        % Message to user
        CSI_Log({'WSVD:'},{'Noise-data processed and stored.'});  
    else
        do_NoiseData = 0; 
        CSI_Log({'WSVD: No noise-data present,'},...
                {'using a noise-mask per voxel instead.'});
    end
end



% Averaging % ------------------------------------------------------ %

% If data not averaged - ask user if to do so.
if do_avg        
    csi.data.raw = mean(csi.data.raw, ind_avg);  % Average
    CSI_Log({'WSVD: Averaged data over index:'},csi.data.labels(ind_avg));    
end


% PCA Denoising % ---------------------------------------------------- %
if do_PCA    
    % Run PCA
    button_CSI_PCA_Denoising_Callback(hobj, [], gui, 0);
    gui = guidata(hobj);
    CSI_Log({'WSVD:'},{'PCA Denoising applied, continuing.'});
end

% Noise Covariance matrix % ------------------------------------------ %

if do_NoiseData == 1 % Use noise-data
    % Get noise-cov using noise data
    noise_cov = CSI_NoiseCov_usingMeasurement(csi.data.noise.raw(:,ch_incl),...
                                              csi.data.noise.labels);
    if isnan(noise_cov), CSI_Log({'WSVD:'},{'Aborted.'}); return; end
    
    % Repeat cov-mat for volume
    noise_cov = repmat({noise_cov}, size(csi.data.raw));     
    msg = 'Used noise-data to calculate the noise-cov matrix.';

elseif do_NoiseData == 2 % Use identity matrix
    noise_cov = diag(ones(csi.data.dim(ind_cha),1));
    noise_cov = repmat({noise_cov},size(csi.data.raw)); 
    msg = 'Used identity matrix as a substitute for the noise-cov matrix.';

else % Use voxel-data
    noise_cov = CSI_NoiseCov_usingData(csi.data.raw, ind_cha);
    msg = 'Used voxel-data to calculate the noise-cov matrix.';
end
% Update LOG
CSI_Log({'WSVD:'},{msg});

% Reshaping (1/2) % -------------------------------------------------- %

% Reorder data setting the channel data to the second dimension index.
% Data format: [dt x channels x rest ... ]

% Reshape data to {nDim x nChan} x nVox
[tmp_cell, permv, szr] = csi_combine_reshape(csi.data.raw, ind_cha);


% WSVD % ------------------------------------------------------------ %

% Number of voxels to combine
nvox     = size(tmp_cell, 1);

% Containers for WSVD output
comb = struct; 
comb.data = zeros(szr(1),nvox); 
comb.Q = zeros(nvox,1); 
comb.W = zeros(size(ch_incl,2), nvox); 
comb.A = zeros(nvox, size(ch_incl,2)); 


% WSVD loop. 
% Apply for every indices excluding the channel index: e.g. every voxel.
for vi = 1:nvox
          
    % Noise covariance matrix
    N = noise_cov{vi};

    % WSVD algorithm
    [comb.data(:,vi), comb.Q(vi),comb.W(:,vi),comb.A(vi,:)] = ...
        WSVD2(tmp_cell{vi}(:,ch_incl), N);
end


% Reshaping (2/2) % ------------------------------------------------- %
% Reshape the cell array back to the an array: dt x chan x residual.
% Reorder to original order as before reshape (1/2).
csi.data.raw = csi_combine_reshape_revert(comb.data, permv, szr);

% Update csi.data.dim data.
csi.data.dim = size(csi.data.raw); 

% Display resulting statical data % --------------------------------- %

if uans_disp(1)
    % One quality value per voxel
    qual = reshape(comb.Q, [1 szr(3:end)]); 

    % Plot as map
    CSI_dataAsTabs(gui, qual, 'WSVD Quality',csi.data.labels, [0 1]);
end

if uans_disp(2)
    % N amplitudes (number of channels) for each voxel
    amp = reshape(comb.A, [size(comb.A,2) szr(3:end)]);
    ndims = numel(size(amp));
    amp = permute(amp , [ndims+1 2:ndims 1]); 

    % Plot as table
    CSI_dataAsTable(real(amp), 'WSVD Amplitudes real-values')
    CSI_dataAsTable(imag(amp), 'WSVD Amplitudes imag-values')
end

if uans_disp(3)
    % N weights (number of channels) for each voxel
    W = reshape(comb.W', [size(comb.W,1) szr(3:end)]);
    ndims = numel(size(W));
    W = permute(W , [ndims+1 2:ndims 1]); 
    % Plot as table
    CSI_dataAsTable(real(W), 'WSVD Weights real-values')
    CSI_dataAsTable(imag(W), 'WSVD Weights imag-values')
end

% Save data % ------------------------------------------------------- %

% Set CSI app data
setappdata(gui.CSIgui_main, 'csi', csi);

% Recalculate xaxis properties
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show statistics nfo
stats = csi_statistics_of_volume(comb.Q);
CSI_Log({'WSVD Quality statistics: --------------------------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});

% Display info to user
CSI_Log({'WSVD; Channels are combined.','WSVD; Included channels:'},...
        {'',ch_incl});   
 
% --- Executes if Sum-of-Squares coil combination is selected
function CSI_Combine_SumOfSquares(hobj, ~)
% Combine the coils using the square root sum of squares method.
%
% combine image data from array coils is a pixel-by-pixel sum of coil
% signals, with each signal weighted by the individual coil sensi-
% tivity at the location of the pixel. A pragmatic alternative com-
% bines the images from the coils as the square root of the sum of
% squares (SOS), which can reduce the signal-to-noise ratio (SNR)
% and introduce bias.

% GUI prepwork % ------------------------------------------- %
 
% Get GUI and object: CSIgui_main
gui = guidata(hobj);

% Create a backup of the current data set.
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before combining channels (SOS).'); end
gui = guidata(hobj);

% Return if no CSI data present.
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Data preperations % -------------------------------------------------- %

% Is there channel/coil data
ind_cha = csi_findDimLabel(csi.data.labels, {'chan', 'cha'});
ind_cha(isnan(ind_cha)) = [];
if isempty(ind_cha)
    CSI_Log({'SOS: Aborted.'},{'No "chan" or "cha" dimension found.'})
    return; 
end

% N-Channels.
nchan = size(csi.data.raw,ind_cha);

% Are there averages available:
ind_avg = csi_findDimLabel(csi.data.labels, {'aver', 'avg', 'nsa'});
ind_avg(isnan(ind_avg)) = [];
if isempty(ind_avg), avg_opt = {'No'}; else, avg_opt = {'Yes','No'}; end

% Exclusion of channels % ------------------------------------------- %

% Convert user input.
do_excl = 1;
if do_excl
    % Get channels to exclude from user
    uans = getUserInput({'Exclude channels:'}, {''},[], 'SOS');
    if isempty(uans), CSI_Log({'SOS: Aborted.'},{''}); end

    % Channels to exclude as double
    ch_excl = str2double(strsplit(uans{1},' ')); 
    if isnan(ch_excl), ch_excl = []; end

    % Create included channels array
    ch_incl = 1:size(csi.data.raw, ind_cha); ch_incl(ch_excl) = [];
else
    ch_incl = 1:nchan;
end


% Reshaping (1/2) % -------------------------------------------------- %

% Reorder data setting the channel data to the second dimension index.
% Data format: [dt x channels x rest ... ]

% Reshape data to {nDim x nChan} x nVox
[tmp_cell, permv, szr] = csi_combine_reshape(csi.data.raw, ind_cha);


% SOS % --------------------------------------------------------------- %
nvox = size(tmp_cell,1);

% tmp_cell = ...
% cellfun(@rssq, tmp_cell, repmat({2},size(tmp_cell)), 'UniformOutput', 0);

% Do some magic here...
for vi = 1:nvox
   vdata = tmp_cell{vi}(:,ch_incl);

   % Including the noise-cov matrix would make this Roemers method.
   % tmp_cell{vi} = sqrt(sum(vdata' * Noise-Cov * conj(vdata),2));
   tmp_cell{vi} = sqrt(sum(vdata .* conj(vdata),2));

end

% Reshaping (2/2) % ------------------------------------------------- %
% Reshape the cell array back to the array: dt x chan x residual.
% Reorder to original order as before reshape (1/2).
szc = size(tmp_cell); szs = size(tmp_cell{1});
csi.data.raw = csi_combine_reshape_revert(...
    reshape(cell2mat(tmp_cell), szs(1), szc(1)), permv, szr);

% Update csi.data.dim data.
csi.data.dim = size(csi.data.raw); 


% Save data % ------------------------------------------------------- %

% Set CSI app data
setappdata(gui.CSIgui_main, 'csi', csi);

% Recalculate xaxis properties
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Update user
CSI_Log({'SOS:'},{'Finished.'});


% ------------------------------------------------------------------ %
% ------------------------------------------------------------------ %


% --- Executes on button press in button_CSI_MaxValue.
function button_CSI_MaxValue_Callback(~, ~, gui)
% Calculate the max for each voxel per slice (Dim 2,3, and 4). Data is
% shown as graph per slice with higher index dimensions represented as a
% line for each voxel or as a table per slice including every higher 
% index dimension per slice-table, seperated by its index.

% Get option from user: per slice or in 3D.
uans = getUserInput_Popup({'Select display of maxima option: '},...
                         {{'Maximum per slice','Maximum in 3D'}}, [], ...
                         'Maximum Maps');                     
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
uans = getUserInput_Popup({'Select peak:'},{{'Yes','No'}}, ...
                          [], 'Maximum per slice');
if isempty(uans)
    CSI_Log({'Skipped max/slice.'},{''}) ; return; 
end

% User requested peak selection
switch uans{1}
    case 'Yes', selectPeak = 1;
    case 'No',  selectPeak = 0;
end

% FID & Echo detection
combine_fe = 0;
if isfield(csi.data,'split')
    % From user: Calculate using FID and echoes?
    uans = getUserInput_Popup({'Use both FID and Echo data?'},...
        {{'Yes', 'No'}}, [], 'Maximum per Slice');
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
    doi = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
else
    doi   = csi.data.raw;
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
CSI_dataAs_Initiate(data_max, 'Maximum per slice', gui, csi.data.labels);

% Show statistics nfo
stats = csi_statistics_of_volume((data_max));
CSI_Log({'Maxima statistics ------------------------------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});
          
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


% Calculate using noise mask
SNR = csi_SNR(data, round(size(data,1)./12), 1, [1 size(data,1)]);

% Get maximum for each fid/spec entry: dim Dt x .. x .. x .. etc.
data = (max(data,[],1)); 
% Remove first index only by permuting to the end.
data = permute(data,[2:numel(size(data)) 1]);
% Normalize
data = data./max(data(:));

% SNR Filter
data(SNR < 10) = 0;

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
%     sc.MarkerFaceAlpha =  0.2;
    
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
range_none = CSI_getPeakOfInterest(csi.xaxis, 'Calculate SNR');
if isempty(range_none), return; end

% SNR Noise mask
uans = getUserInput({'Size of SNR noise mask: '},{round(csi.data.dim(1)/12)});
if isempty(uans), CSI_Log({'Skipped SNR calculations.'},{''}) ; return; end

% Convert user input
mask_size = str2double(uans{1});

% SNR-method (real/magnitude) and display method (table or graph)
uans = getUserInput_Popup({'SNR Signal unit: '},...
                         {{'Real', 'Magnitude'}}, [], 'SNR');
if isempty(uans), CSI_Log({'Skipped SNR calculations.'},{''}) ; return; end
switch uans{1} % SNR method
    case 'Real', SNRmethod = 1; case 'Magnitude', SNRmethod = 0; 
end

              % --------------- % SNR % --------------- %

% Display Info %
CSI_Log({'Calculating SNR per voxel, please be patient.'},{''});

% Calculate using noise mask
SNR_all = csi_SNR(csi.data.raw, mask_size, SNRmethod, range_none);
% Convert NaNs to zero
SNR_all(isnan(SNR_all)) = 0; 

% Show statistics nfo
stats = csi_statistics_of_volume(SNR_all);
CSI_Log({'SNR statistics: Full Volume ---------------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});


           % --------------- % DISPLAY DATA % --------------- %

% \\ Display Data
CSI_dataAs_Initiate(SNR_all, 'SNR', gui, csi.data.labels);



% --- Get 2D-CSI plot slider data
function [pan_gui, gui2D] = CSI_2D_getDataSliders(gui)
% FIGURE OBJECT % ------------------------------------------------------ %
% Get 2D-plot figure its object, apply chosen options and save as image.

% Check if 2D Plot is active (open).
fig_obj =  findobj('type','figure','tag','CSIgui_plot2D');
% If 2D Plot not active, open it.
if isempty(fig_obj)
    panel_2D_DataSliders([],[],gui); CSI_2D_initiate2D(); 
end

% Get 2D-Plot object and 2D-Panel figure
fig_obj = findobj('type','figure','tag','CSIgui_plot2D');
pan_obj = findobj('type','figure','tag','CSIpanel_2D_DataToDisplay');
if ~isempty(pan_obj)
    pan_gui = guidata(pan_obj);
    gui2D = guidata(fig_obj);
else
    % Open slider panel
    panel_2D_DataSliders([],[],gui);
    % Get the object
    pan_obj = findobj('type','figure','tag','CSIpanel_2D_DataToDisplay');
    pan_gui = guidata(pan_obj);
    % Get plot index for 2D plot figure gui data.
    gui2D = guidata(fig_obj); gui2D.plotindex
    
    % Set sliders from values in current 2D-plot figure.
    for sli = 1:size(pan_gui.sliders,2)
        pan_gui.sliders{sli}.Value = gui2D.plotindex{sli};
        pan_gui.texts{sli}.String = sprintf('%i/%i', gui2D.plotindex{sli},...
            gui2D.dim(sli+2));
    end
end

% --- Executes on button press in button_ws: export appdata to workspace.
function button_ws_Callback(~, ~, gui)
% Export appdata from CSIgui to workspace; mainly used for debugging.

appdat = getappdata(gui.CSIgui_main);
assignin('base', 'appdat', appdat);

% --- Executes on button press in button_Info.
function button_Info_Callback(~, ~, gui)

% ---- % Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Update LOG
log = get(gui.listbox_CSIinfo,'String');
csi.log = char(log);

% Send to workspace
assignin('base', 'csi', csi);
fprintf('Variable "csi" accessible from workspace.\n');

% ---- % Get MRI data-structure
if isappdata(gui.CSIgui_main,'mri') 
    mri = getappdata(gui.CSIgui_main,'mri');
    % Send to workspace
    assignin('base', 'mri', mri);
    fprintf('Variable "mri" accessible from workspace.\n');
end

% ---- % Get CONV data-structure
if isappdata(gui.CSIgui_main,'conv') 
    conv = getappdata(gui.CSIgui_main,'conv');
    % Send to workspace
    assignin('base', 'conv', conv);
    fprintf('Variable "conv" accessible from workspace.\n');
end

% Open data-directory
CSI_openDataDirectory(gui)

% Open viewer
fprintf('Opening csi-struct in NFO-viewer.\n'); NFOviewer(csi);



% --- Executes on button press in button_CSI_Peak_Map.
function button_CSI_Peak_Map_Callback(~, ~, gui)
% Plot a map of a specific peak maximum including images and voxel grid.


% INITIATE % --- %

% Check if csi appdata is present
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get peak of interest
[doi, ~, ~] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
doir = real(doi); % doim = imag(doi);

% MAP % --- %

% Maximum positions and values: e.g. Map
map = max(doir, [], 1);
% Normalize map
% nfac =  max(map(:)); nmap = map./nfac;

CSI_dataAs_Initiate(map, 'Peak', gui, csi.data.labels);

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

% --- Executes on button press in button_CSI_PhaseShift.
function button_CSI_PhaseShift_Callback(~, ~, gui, backup)
% Add function csi_frequencyShift
if nargin < 4, backup = gui.checkbox_backup.Value; end

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

% --- Executes on button press in button_CSI_SpatialShift_kspace.
function button_CSI_SpatialShift_kspace_Callback(~, ~, gui, backup)
% Shift CSI k-space to spatially shift the volume a number of voxels in the
% prefered direction.

if nargin < 4, backup = gui.checkbox_backup.Value; end
% Create backup
if backup, CSI_backupSet(gui, 'Before voxel shift'); end

% Check data domain
domain = CSI_getDomain(gui);
if strcmp(domain, 'time')
    CSI_Log({'MRS data is in time domain; '},...
            {'Convert it to frequency domain to use this shift method.'});
    return;
end 


% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% USER INPUT % ------------------- %

% Get shift from user
uans = getUserInput(...
    {'Requires k-space domain CSI data. Number of voxels to shift: '},...
    {[0.5 0.5 0.5]});
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


% Shift the data
[csi.data.raw, ~] = csi_voxelshift(csi.data.raw, shift, k_dim);


% CLEAN UP % -------------------- %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% LOG
CSI_Log({'Phase change applied to spatially shift the CSI volume. Shifted by:'},...
        {strjoin(strsplit(uans{1},' '),' | ')});
                
% --- Executes on button press in button_CSI_FidOrEcho.
function button_CSI_FidOrEcho_Callback(~, ~, gui)
% Split data to FID and Echo data, specially design for use with (A)MESING
% data. Echoes and FIDs require different processing resulting in different
% sample size; (zerofilling to match bandwidth)

% Get CSI data-structure
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if ~isfield(csi.data,'split')
    
                    % ------- % Cut FID/ECHO % ----------%
    
    % Create backup
    backup = gui.checkbox_backup.Value; 
    if backup, CSI_backupSet(gui, 'Before FID/Echo split'); end                  
                    
    % Userinput % ------------ %

    % Get FID/Echo dimension
    uansDim = getUserInput_Popup(...
    {'Which dimension represents the FID and echoes?'},{csi.data.labels},...
    [], 'FID/ECHO');
    if isempty(uansDim), CSI_Log({'Skipped FID & Echo split.'},{''}); return; 
    end
    
    % FID/Echo dimensions and size
    fidecho_dim = strcmp(csi.data.labels, uansDim{1});
    dim_sz = csi.data.dim(fidecho_dim); 

    % Get FID index @ FID/Echo dimension
    uansInd = getUserInput_Popup({'Which index in this dimension is the FID?'},...
                          {{1:dim_sz}}, [], 'FID/ECHO' );
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
    log_msg = 'FID and echo data swapped. In memory: ';
end


              % ------------ % CLEAN UP % ------------ %

% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);

% Update LOG.
CSI_Log({log_msg}, {csi.data.split.type});

% --- Executes on button press in button_Log_DeleteLine.
function button_Log_DeleteLine_Callback(hObj, ~, ~)
CSI_Log_deleteLine(hObj);

% --- Executes on button press in button_CSI_AutoProcessing
function button_CSI_AutoProcessing_Callback(~, ~, gui)
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
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui,'Initiated auto-processing.'); end

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
dw = 5; dh = 4; szt = [320 18]; szl = round(dh/2); % DEFINES FIGURE SIZE

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
    'BackgroundColor',gui.colors.main,...
    'ForegroundColor',gui.colors.hilight2,...
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


% --- Executes by map-scripts to start any visualiziation of data
function CSI_dataAs_Initiate(data, data_tag, gui, labels)
% After calculating some maps or anything 3D, and one wants to display it.
% Call this function. It will ask the user the proper info; display type,
% filter by SNR and more.


if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if nargin < 4    
    labels = csi.data.labels;
end

% Display type from user
uans = getUserInput_Popup({'Display type: ','Apply excluding voxel mask:'},...
                          {{'Map','Graph','Table','Histogram'}, {'Yes', 'No'}},...
                          [], data_tag );
if isempty(uans)
    CSI_Log({sprintf('%s mapping skipped.', data_tag)},{''}); 
    return; 
end
dataDisp = lower(uans{1}); % Display method
switch uans{2}, case 'Yes', doMask = 1; otherwise, doMask = 0; end

% If voxel mask is requested.
if doMask
    if ~isfield(csi, 'voxelmask')
        button_CSI_VoxelMask_Callback([], [], gui)
    end
    csi = getappdata(gui.CSIgui_main, 'csi');
    % User can quit/cancel voxel-mask creation - this catches that error
    % and continues without creating masked-data
    if isfield(csi, 'voxelmask')
        mask = csi.voxelmask;  
    
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
    else
        CSI_Log({'Voxel-Mask was not loaded into memory.'}, ...
                {'Continuing with full data-volume without masking.'})
    end
end

% Show statistics nfo
stats = csi_statistics_of_volume(data);
CSI_Log({[ data_tag ' statistics (VOI) ---------------------- %'],...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: ', 'Applied mask:'},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max),...
          uans{2}});

% Display Info %
CSI_Log({'Starting data display:'}, {dataDisp});

% Switch to data-display type.
switch dataDisp
    case 'table'      % Table % ----- %    
        % Send to tableData function, to show each slice as a table.
        % If higher dimensions are available, the data will be 
        % concentonated per slice. Graph may be better suited.
        CSI_dataAsTable(data, data_tag)
    
    case 'graph'      % Graph % ----- %
        % Send to displayData function, to show each slice as a graph.
        % If higher dimensions represent only 1 value/voxel, a dot is 
        % shown. Adviced is using table in those cases.
        CSI_dataAsGraph(data, gui, data_tag);

        % Optional older function
        % CSI_dataAsGraph2(data, data_text, datatag, gui, img)

    case 'map'
        % Send data to dataAsTabs and create maps in a tabbed figure for 
        % all slices (or only the current slice plotted).
        
        % Prepare dataAsTabs input
        [data, color_scale, sloi] = CSI_dataAsTabs_Prepare(data, gui);
        if isnan(data), return; end
        % Plot as tab
        CSI_dataAsTabs(gui, data, data_tag, labels, color_scale, sloi);
    
        % Show statistics nfo
        stats = csi_statistics_of_volume(data);
        % Show statistics nfo
        CSI_Log({[ data_tag ' statistics (Plotted) ------------------- %'],...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: ', 'Applied mask:'},...
         {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max),...
          uans{2}});

    case 'histogram'
        % Plot data as histogram - simple display.    
        fig = figure(); ax = axes(fig);
        mx = max(data(:), [], 'omitnan'); md = median(data(:),'omitnan');  
        tmp = data(:);
        histogram(ax, tmp(~isnan(tmp)), round(100+(mx./md)) );
        title([data_tag ' Histogram']); xlabel([data_tag ' Bins']);
end

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

toolbar_create(figure_object)
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

% Store all other indexes..
if numel(plot_par.dim) > 2, nDimCtmp = num2cell(plot_par.dim(3:end));
else, nDimCtmp = {1};
end
nDimCtmp = cellfun(@(x) 1:x, nDimCtmp,'UniformOutput',0);
plot_par.select_all_dim = nDimCtmp;

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
    if numel(plot_par.dim) >3
        % Dirty quick edit QvH
        sz = size(data);
        nDimC = num2cell(sz(5:end));
        
        %nDimC = num2cell(plot_par.dim(4:end));
    else, nDimC = {1};
    end
    % To linear vector per cell.
    nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

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
             mean(plot_data(:),'omitnan'), min(plot_data(:)), max(plot_data(:))),...
             'Color', [0.6 0.6 0.6],'FontSize',8,'FontWeight','Bold',...
             'HorizontalAlignment', 'Right','VerticalAlignment', 'Top');    
        
        % Cosmetics % -------------------- %
        set(plot_par.ax{ri,ci},'Color', 'None',...
                       'XColor', [0.4 0 0],'YColor', [0.4 0 0],...
                       'LineWidth', 1.75, 'Xtick',[], 'Ytick',[]);
        
        end % row loop
    end % column loop
end % slice loop

%tgui.fig.ToolBar = 'figure';
toolbar_create(plot_par.fh)

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

    % Convert double-array to cell-array and then to string
    tableData_str = cellfun(@sprintf, repmat({'%f'},...
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

% Tabgroup
tabgroup = src.Parent.Parent.Children(2);
% Number of tabs
nTabs = size(tabgroup.Children);
for tabi = 1:nTabs
    % Tab-obj
    tab = tabgroup.Children(tabi);
    % Children in tab: contains axes
    tab_childs = tab.Children;
    % Number of children
    nChildren = size(tab_childs,1);
    for kk = 1:nChildren
        if strcmp(tab_childs(kk).Type,'axes')
            tab_childs(kk).Color = [tab_childs(kk).Color 0.33];
        end
    end
end

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
     
    pause(0.75); % Make sure display is updated.       

    check = 0; n = 0;
    while check ~= 1
        tabg = figobj.Children(ind);
        if isequal(tabg.SelectedTab, tabs(kk)), check = 1; end
        n = n + 1; if n == 100, break; end
    end

    pause(0.75);

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

% \\ USER INPUT: new color range
cur_clr_rng = plot_par.clr_rng;
uans = getUserInput({'New color scale range:'}, {num2str(cur_clr_rng)});
if isempty(uans), return; end
new_clr_rng = str2double(strsplit(uans{1}));

% Calculate new color-range table
clr_map = plot_par.clr_map;
clr_val = linspace(new_clr_rng(1),new_clr_rng(2),size(clr_map,1));      

% Store into gui-data-obj
tgui.plot_par.clr_val = clr_val; tgui.plot_par.clr_rng = new_clr_rng;            
tgui.plot_par.clr_map = clr_map;

% \\ ADJUST: axis colors
ax = plot_par.ax(:);
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        % Axis data value
        vox_val = ax{axi}.Children.YData;
        % New color index: this index relates to clr_map
        [~, clr_ind]= min(abs(clr_val-vox_val)); 
        % Adjust color of axis background
        ax{axi}.Color = [clr_map(clr_ind,:) 0.33];
        % Adjust color of data-point in axis // Line-point
        ax{axi}.Children.Color = [clr_map(clr_ind,:) plot_par.alpha];
    end    
end
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
uans = getUserInput({'New transparency value:'}, {num2str(cur_alpha_val)});
if isempty(uans), return; end
new_alpha = str2double((uans{1}));
tgui.plot_par.alpha = new_alpha;

% Store into gui-data-obj
clr_val = tgui.plot_par.clr_val;
clr_map = tgui.plot_par.clr_map;

% \\ ADJUST: axis alpha value
ax = plot_par.ax(:);
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        % Axis data value
        vox_val = ax{axi}.Children.YData;
        % New color index: this index relates to clr_map
        [~, clr_ind]= min(abs(clr_val-vox_val)); 
        % Adjust color of axis background
        ax{axi}.Color = [clr_map(clr_ind,:) new_alpha];
        % Adjust color of data-point in axis // Line-point
        ax{axi}.Children.Color = [clr_map(clr_ind,:) new_alpha];
    end    
end
guidata(figobj, tgui);

% --- Executes on press of toolbar's convert to table button
function toolbar_MapToTable(src,~)
% Get data from maps and display in table

% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

ax = tgui.plot_par.ax(:);
data = NaN(size(ax));
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        data(axi) = ax{axi}.Children.YData;
    end    
end
data = reshape(data, size(tgui.plot_par.ax));
sz = size(data);
data = permute(data, [numel(sz)+1 1:numel(sz)]);

CSI_dataAsTable(data,tgui.fig.Name);

function toolbar_Statistics(src,~)
% do something
% \\ GET: object handles for figure and tab-group
figobj = src.Parent.Parent;
tgui = guidata(figobj);

% Which data to visualize
uans = getUserInput_Popup({'Calculate stastics for:'},...
                         {{'Current', 'All'}});
if isempty(uans), return; end


% Get data
ax = tgui.plot_par.ax(:);
data = NaN(size(ax));
for axi = 1:size(ax,1)
    if ~isempty(ax{axi})
        data(axi) = ax{axi}.Children.YData;
    end    
end
data = reshape(data, size(tgui.plot_par.ax));
sz = size(data);
data = permute(data, [numel(sz)+1 1:numel(sz)]);

% Tab of interest
if strcmp(uans{1}, 'Current')
    toi = cell2mat(cellfun(@isequal,tgui.tabh, ...
        repmat({tgui.tabg.SelectedTab}, size(tgui.tabh)), 'uniform', 0));
    toi = tgui.tabh{toi == 1};
    tab_dims = strsplit(toi.Title, ' ');
    tab_dims = num2cell(str2double(tab_dims));
    data = squeeze(data(:,:,:,tab_dims{:}));
end

% Calculate and show results.
stats = csi_statistics_of_volume(data);
Statistics_Viewer(stats);



% --- Executes by CSI_dataAs_Initiate to prepare data for tabbed-maps
function [data, color_scale, sloi] = CSI_dataAsTabs_Prepare(data, gui)
% This function will get userinput, prep the data and return required
% variables to plot data-array in a tabbed-figure.

% USERINPUT: What data to show?
qry = {'Data range to show: ','Color scale range: '};
def = {{'All', 'Current Slice'}, {'Min to Max', 'Histogram optimized'}};
uans = getUserInput_Popup(qry, def, [], 'Data Display - Maps');                                                 
if isempty(uans), data = nan; color_scale = nan; sloi = nan; return; end

% \\ Get part of data-array to visualize
if strcmpi(uans{1},'Current Slice')
    % \\ Get current plotted slice index from CSIgui 2D-plot
    % \\ Get curent slice data from SNR_all array.

    % Get the panel-slider object and 2d-plot object
    [~, gui2D] = CSI_2D_getDataSliders(gui);

    % Plotted slice index of interest
    sloi = gui2D.plotindex{1};

    % Get data of interest.
    data_dim = size(data);    
    data_dim_range = arrayfun(@(x) 1:x,data_dim, 'uniform',0);
    data_dim(4) = 1; data_dim_range{4} = sloi;

    % 2. Get data of interest from SNR_all.
    data_cut = data(data_dim_range{:});
    permvec = 1:numel(data_dim);
    permvec(4) = []; permvec(end+1) = 4;
    data_cut = permute(data_cut, permvec);

    % Replace SNR-all with SNR-cut
    data = data_cut;
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
function fh_all = CSI_dataAsTabs(gui, data, tag, labels, color_range, sloi)
% Data is shown as a color map and also a numeric representation in the
% figure window. 
%
% gui           = CSIgui its gui-object requiring only gui.colors.
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
    CSI_Log({'Aborted. Function dataAsTab requires 2D or 3D data',...
             'Reshape the data as such to enable maps.'},{'',''}); return; 
end



% Number of tab-windows
nwindows = 1; index_range = {1};
if numel(dim) > 4
    nwindows = dim(5); 
    if numel(dim) > 5
        % Data ranges
        index_range = cellfun(@(x) 1:x, ...
            num2cell(dim(6:end)), 'uniform', 0); 
    end
end


% Restore data
data_main = data; tag_main = tag;
if numel(size(labels)) > 4
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
fh = CSI_dataAsTabs_create_figure(tag, gui.colors.main);


            % -------- % Figure: Create tabs % ------------------------- %
           
% Create a tab per slice including a voxel-grid overlay and
% create plot_par            
CSI_dataAsTabs_create_griddedTabs(fh, data, gui.colors) ;       

            % -------- % Figure: Create image axis % ------------------- %                


% Get images
plot_img = 0;
if isappdata(gui.CSIgui_main,'conv')
    conv = getappdata(gui.CSIgui_main,'conv');
    img = MRI_matchSlices(gui.CSIgui_main); 
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
    loadBar(tabi./plot_par.tabs_total , 'Plotting data...');

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
        if isnan(sloi)
            sli_img = tab_index_cell{1}; % Slice from data.
        else
            % Slice of interest - applicable when Current Slice is chosen
            % for displaying data. Correct image-slice is taken.
            sli_img = sloi; 
        end
        
       
        img2plot = img(:,:,sli_img); 

        if size(img2plot,1) > 1 % Safety against NaN-values
        imagesc(img2plot, 'parent', tgui.himg{tab_index_cell{:}}); 


        % Image Contrast.
        if isfield(conv, 'contrast')
            clim(tgui.himg{tab_index_cell{:}}, conv.contrast);
        else
            contrast_min = min(img2plot(:));
            contrast_max = max(img2plot(:))*0.75;
            if contrast_max <= contrast_min
                contrast_max = contrast_min +1; 
            end 
           

            v = version('-release'); v(end) = []; v = str2double(v);
            if v < 2023
                caxis(tgui.himg{tab_index_cell{:}},...
                    [contrast_min contrast_max]);
            else
                clim(tgui.himg{tab_index_cell{:}},...
                    [contrast_min contrast_max]);
            end
        end
        colormap(tgui.himg{tab_index_cell{:}}, gray(255));
        
        end
    end


    for xi = 1:plot_par.dim(1)                  % Col loop. X
        for yi = 1:plot_par.dim(2)              % Row loop. Y

            % VOXEL VALUE
            vox_val = data(:,xi,yi,sli,1,tab_index_for_data{:});
            
            % MAP COLOR
            if isnan(vox_val)
                clr = [0 0 0];
            else
                % Get color from value2color-table
                [~, clr_ind]= min(abs(clr_val-vox_val));
                clr = clr_map(clr_ind,:);
            end
            
            
            
            % Add value in center of voxel
            % // Set limit of yaxis
            ylimit = [0 2*vox_val]; ylimit(isnan(ylimit)) =  1;
            if ylimit(2) <= ylimit(1), ylimit(2) = ylimit(1)+1; end
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.YLim = ylimit;
            
            % // Set limit of xaxis
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.XLim = [0 2];
            
            % // plot value
            plot(tgui.plot_par.ax{xi,yi,tab_index_cell{:}},...
                1,vox_val, '.', 'MarkerSize', 0.5, ...
                                'Color', [clr plot_par.alpha])
            
            % // Turn off all axis cosmetic;
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.YTick = [];
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.XTick = [];
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.YLabel = [];
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.XLabel = [];
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.Box = 'off';
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.XColor = ...
                plot_par.colors.main;
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.YColor = ...
                plot_par.colors.main;
            
            
            % Set background color of axis according to colormap            
            tgui.plot_par.ax{xi,yi,tab_index_cell{:}}.Color = ...
                [clr 0.33];            
            
        end 
    end
    


end
loadBar(NaN);

% Save tgui-data to gui itself 
% NB. Should use setappdata if data is involved!
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
        'ForegroundColor', plot_par.colors.text_title);

    % Plot voxel grid
    % Input: target figure, target figure size, data dimensions, range and color.
    CSI_2D_grid(tgui.tabh{tmp_plotindex_tabs_cell{:}},...
        fh.Position(3:4), plot_par.dim, ...
        plot_par.range, plot_par.colors.text_title);

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
ax = cell(plot_par.dim);

% Loop each tab of figure
for sli = 1:plot_par.tabs_total                % Sli loop/tabs loop
    loadBar(sli./plot_par.tabs_total , 'Adding voxel-axis...');


    tab_index = plot_par.tabs_index_table(sli,:);
    tab_index_cell = num2cell(tab_index);
    
    % Plot csi voxel per axis in plot_par.grid
    for xi = 1:plot_par.dim(1)                  % Col loop. X
        for yi = 1:plot_par.dim(2)              % Row loop. Y

            % AXIS DETAILS 
            % X and Y position in the figure of the axis to plot.
            x   = plot_par.grid.x(yi,xi); y = plot_par.grid.y(yi,xi);
            % Position of axis to plot
            pos = [x y plot_par.res(1) plot_par.res(2)];
            % Create axis with pos(3,4) size at pos(1,2) position
            if ~ishandle(tgui.tabh{tab_index_cell{:}}), return; end
            ax{xi,yi,tab_index_cell{:}} = axes('parent',tgui.tabh{tab_index_cell{:}},'position',pos);           

            % AXIS COSMETICS
            set(ax{xi,yi,tab_index_cell{:}},...
                   'Color','None',...
                   'XColor', plot_par.colors.main,...
                   'YColor', plot_par.colors.main,...
                   'LineWidth', 1.7, 'Xtick',[], 'Ytick', [],...
                   'TickLength',[0 0.00001], 'Box', 'off');             

        end 
    end 
end
tgui.plot_par.ax = ax; guidata(fh, tgui);
loadBar(NaN);

% --- Executes by dataAs-scripts to filter calculated data
function data = CSI_dataAs_SNRfilter(data, tag, gui, doi_range)
% Apply an SNR filter to the data-volume that needs to be displayed.
% SNR is calculated for a peak given by doi_range, if not given, the user
% will be prompted with peak-selection. The main CSI data in memory will be
% used to calculate the SNR and the filter is applied on data.

        % --------------- % SNR FILTER % --------------- %

% GET DATA 
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if nargin < 4 || isempty(doi_range)
    % Get peak of interest
    [~, ~, doi_range] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);
end

% \\ USER INPUT
filter_snr = getUserInput(...
    {'Minimum value for SNR Filtering: [0 = off]', 'SNR window:'},...
    {'0', round(csi.data.dim(1)./10)});
if isempty(filter_snr)
    CSI_Log({['Skipped ' tag ' display.']},{''}) ; return; 
end

if str2double(filter_snr{1}) ~= 0
    
    snr_limit = str2double(filter_snr{1});
    snr_mask = str2double(filter_snr{2});
    
    % Display Info %
    CSI_Log({['Calculating SNR per voxel, '... 
              'and filter data using minimum SNR limit: ']},...
            {snr_limit});

    % Noise mask
    mask_size = snr_mask;

    % Dimensions of array to calculate SNR.
    snr_dims = cellfun(@(x) 1:x, num2cell(csi.data.dim),'uniform', 0);
    snr_dims(1) = [];
    snr_dims = [{1:csi.data.dim(1)}, snr_dims]; 
    
    % Calculate using noise mask
    SNR_all = csi_SNR(csi.data.raw(snr_dims{:}), mask_size, 1, doi_range);
    
    % Convert NaNs to zero
    SNR_all(isnan(SNR_all)) = 0; 
    
    % Filter boolean 
    SNR_bool = SNR_all < snr_limit;
    
    % Filter data
    data(SNR_bool) = NaN;

end




% COORDINATES % -------------------------------------------------------- %
% MRI and CSI

% --- Executes on button press in button_CSI_setCoordinates.
function button_CSI_setCoordinates_Callback(~, ~, gui)
% Calculate coordinates of the CSI-volume:
%
% 1. Gather information from header and ask for user intervention if 
%    required
% 2. All nfo and userinput is used to calculate coordinates for each voxel.
%
% Required: fov, resolution, offcenter and volume dimensions.

% Return if no csi data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Get possibly available ori-struct
if isfield(csi,'ori'), ori = csi.ori; else, ori = struct; end

% SPATIAL INDEX % ----------- %
space_str = {'kx','ky','kz','x','y','z'};
space_dim = csi_findDimLabel(csi.data.labels, space_str);  
space_dim(isnan(space_dim)) = []; 

% Ask user if no kx/ky/kz are found
if isempty(space_dim)
    uans = getUserInput({'Spatial dimensions MRS data? (x/y/z)'},...
                        {[2 3 4]});
    if isempty(uans), CSI_Log({'Skipped MRS coordinates.'},{''}); return; 
    end
    % Convert to double
    space_dim = str2double(strsplit(uans{1},' '));
end
    

% HEADER DATA % ----------------- %
% SPAR-INFO
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

% HEADER DATA % ----------------- %
% TWIX-INFO  (dat-file)
if isfield(csi,'twix')
    
    % Image orientation
    if csi.twix.Meas.VoiNormalSag,     ori.image_orientation = 'Sag';
    elseif csi.twix.Meas.VoiNormalCor, ori.image_orientation = 'Cor';
    elseif csi.twix.Meas.VoiNormalTra, ori.image_orientation = 'Tra';
    end

    % Offcentre
    % Also located in Config.VoI_Position_Cor/Sag/Tra
    % Sagital means slices in LR dir
    if isempty(csi.twix.Meas.VoiPositionSag), ori.offcenter(2) = 0;
    else, ori.offcenter(2) = csi.twix.Meas.VoiPositionSag;
    end
    % Coronal means slices in AP dir
    if isempty(csi.twix.Meas.VoiPositionCor), ori.offcenter(1) = 0;
    else, ori.offcenter(1) = csi.twix.Meas.VoiPositionCor;
    end
    % Transverse means slices in FH dir
    if isempty(csi.twix.Meas.VoiPositionTra), ori.offcenter(3) = 0;
    else, ori.offcenter(3) = csi.twix.Meas.VoiPositionTra;
    end
    
    
    % FOV
    if isfield(csi.twix.Config,'VoI_RoFOV')
        ori.fov(1) = csi.twix.Config.VoI_RoFOV; % Read out fov
        ori.fov(2) = csi.twix.Config.VoI_PeFOV; % Phased enc fov
        ori.fov(3) = csi.twix.Config.VoI_SliceThickness; % Slice-dir FOV
    elseif isfield(csi.twix.MeasYaps.sSpecPara,'sVoI')                
        ori.fov(1) = csi.twix.MeasYaps.sSpecPara.sVoI.dReadoutFOV;
        ori.fov(2) = csi.twix.MeasYaps.sSpecPara.sVoI.dPhaseFOV;
        ori.fov(3) = csi.twix.MeasYaps.sSpecPara.sVoI.dThickness;
    end

    % Dimensions: [AP LR FH]
    ori.dim(1) = csi.twix.Meas.ReadResolution;
    ori.dim(2) = csi.twix.Meas.PhaseEncodingLines;
    ori.dim(3) = csi.twix.Meas.SliceResolution;

    % RES:  [AP LR FH] - [Y X Z]
    ori.res(2)= ori.fov(1)./csi.data.dim(space_dim(1)); % LR = x
    ori.res(1)= ori.fov(2)./csi.data.dim(space_dim(2)); % AP = y
    ori.res(3)= ori.fov(3)./csi.data.dim(space_dim(3)); % FH = z
    
    % Orientation (HFS/FFS etc.)
    % hdr.Meas.Matrix Meas.PatientPosition
    % twix.Meas.VoiNormalTra or Cor/Sag

    % After multi-file analysis: Data needs flip in YZ direction
    if ~isfield(ori, 'flipCorrection')
        ori.flipCorrection = 1;
    end

    % After multi-file analysis: Data needs voxel shift correction
    if ~isfield(ori, 'voxShiftCorrection')
        % If READ i.e. LR-dimension is ODD...
        if mod(csi.twix.Meas.ReadResolution,2) ~= 0 %  ODD
            ori.voxShiftCorrection = 1;
        else
            ori.voxShiftCorrection = 0;
        end
    end
end



% ----------------------------------------------------------------------- %
% USERINPUT  - SPATIAL % ------------------------------------------------ %

                    % RESOLUTION  +   OFFCENTER %

% Default answers
if isfield(ori,'res'), dans{1} = ori.res; 
else, dans{1} = '20 20 20';
end
if isfield(ori,'offcenter'), dans{2} = ori.offcenter;
else, dans{2} = '0 0 0';
end
if isfield(ori,'shift'), dans{3} = ori.shift;
else, dans{3} = '0.5 0.5 0.5';
end

% Flip correction for dat-file
if isfield(ori, 'flipCorrection'), dans{4} = ori.flipCorrection;
else, dans{4} = 0;
end

% Voxel flip correction for dat-file
if isfield(ori, 'voxShiftCorrection'), dans{5} = ori.voxShiftCorrection;
else, dans{5} = 0;
end


% Resolution and Offset userinput
qry = {'Resolution (AP LR FH) [mm]:', 'Offcenter (AP LR FH) [mm]: ', ...
    'Voxel shift (AP LR FH) [-]: ', ...
    'Flip-fix for dat-file [0 = No, 1 = Yes, -1 = Applied]: ',...
    'Shift-fix for dat-file [0 = No, 1 = Yes, -1 = Applied]: '};
uans = getUserInput(qry, dans);
if isempty(uans)
    CSI_Log({'Skipped calculating coordinates.'},{''}); 
    return; 
end


% ----------------------------------------------------------------------- %
% Calculate Coordinates % ----------------------------------------------- %

% Corrections % ----------------- %                      
vox_cor = 1; fft_cor = 0; % Default settings.

% Prep misc. parameters
ori.res       = str2double(strsplit(uans{1},' '));
ori.offcenter = str2double(strsplit(uans{2},' '));
ori.shift     = str2double(strsplit(uans{3},' '));
ori.flipCorrection = str2double(strsplit(uans{4}));
ori.voxShiftCorrection = str2double(strsplit(uans{5}));

% Dimensions of data  [AP LR FH]
ori.dim       = csi.data.dim(space_dim); % in [X Y Z] == [C R S]
% AP to 2nd dim, LR to first. 
ori.dim([2 1]) = ori.dim([1 2]); % in [Y X Z] == [R C Z]

% Correction factors saved.
ori.vox_cor = vox_cor; ori.fft_cor = fft_cor;       

% Coordinates of CSI data. % ------------------------- %
% Adds fields: coordinate vector (vector) & coordinate limits (limit) &
% volume limit (limit_vol).

% Default voxel shift
% AP [negetive is image down] LR [positive is image to left] FH
def_shift = ori.shift;

% Set options
opts.vox_cor = vox_cor; opts.fft_cor = fft_cor;

% Calculate coordinates for CSI volume.
csi.ori = CSI_coordinates_calculate(...
    ori.res, ori.offcenter, ori.dim, def_shift, opts);
% ---------------------------------------------------- %

% DAT-file Flip corrections % ------------------------------------------ %
if ori.flipCorrection == 1

    % Flip over k-space
    ind_to_flip = csi_findDimLabel(csi.data.labels,{'ky','kz','y','z'});
    ind_to_flip(isnan(ind_to_flip)) = [];
    for kk = 1:size(ind_to_flip,2)
        csi.data.raw = flip(csi.data.raw,ind_to_flip(kk));
    end
    % 
    csi.ori.flipCorrection = -1;

    CSI_Log({'Applied CSI-volume orientation-correction for dat-file: '},...
        { sprintf('%s | ', csi.data.labels{ind_to_flip}) });
else
   csi.ori.flipCorrections = -1;
end

% DAT-file Voxel Shift corrections % ----------------------------------- %
if ori.voxShiftCorrection == 1
    
    % This should shift according to the readout-dimension. Which could be
    % either X or Y... 

    ind_to_shift = csi_findDimLabel(csi.data.labels,{'kx', 'x'});
    ind_to_shift(isnan(ind_to_shift)) = [];
    csi.data.raw = circshift(csi.data.raw, -1, ind_to_shift);
    csi.ori.voxShiftCorrection = -1;

    CSI_Log({'Applied voxel shift correction for dat-file: '},...
        { sprintf('%s | ', csi.data.labels{ind_to_shift}) });
else
   csi.ori.voxShiftCorrection = 0;
end
% ---------------------------------------------------------------------- %

% Save to appdata
setappdata(gui.CSIgui_main,'csi',csi);   

% --- Executes to calculate CSI coordinates from offcenter and voxelsize
function ori = CSI_coordinates_calculate(res, offcenter, dim, shft, opts)
% Calculate MRSI/CSI coordinates and additional spatial information.
%
% [col row slice] == [RL AP FH] == [X Y Z]
%
% Input:
% res           Voxel resolution in mm [RL AP FH];
% dim           CSI array size [RL AP FH];
% offcenter     Offcenter in milimeters [RL AP FH];
% shft          Voxel shift in each direction [RL AP FH]; 
%               Default: [0 0 0];
% opts          Voxel and fft corrections applied yes(1) or no(0).
%               Default: vox(1) and fft(0).
% 
% Output: 
% ori-struct with all volumetric, coordinates, and spatial information of
% the CSI grid (MRS-data) with the requested shift.


% Handle nargin
if nargin < 5, opts.vox_cor = 1; opts.fft_cor = 0; end
if nargin < 4, shft = [0 0 0]; end

% Set resolution
ori.res = res; ori.offcenter = offcenter; ori.dim = dim;

% Correction factors saved.
if isfield(opts, 'vox_cor'), ori.vox_cor = opts.vox_cor; end
if isfield(opts, 'fft_cor'), ori.fft_cor = opts.fft_cor; end
    

% Apply any given additional shift
for kk = 1:3
    ori.offcenter(kk) = ori.offcenter(kk) + (res(kk).*shft(kk)); 
end
ori.shift = shft;

% Coordinates of CSI data. % ------------------------ %
% Fields: 
%   Coordinate vector (vector) 
%   Coordinate limits (limit) 
%   Volume limit      (limit_vol)
ori = csi_coordinates(ori, 'center', ori.vox_cor, ori.fft_cor); 

% Rewrite FOV; redundant
ori.fov = ori.lim_vol(:,2)' - ori.lim_vol(:,1)' ;

% Calculate volume grid % ----------- %
% Meshgrid
[ori.mesh.x, ori.mesh.y, ori.mesh.z] = ...
    meshgrid(ori.vector{1} , ori.vector{2}, ori.vector{3});


% Clean up % ------------------------ %
% Update LOG
CSI_Log({ '',...
    'CSI-parameters ----------------------------------------',...
    'Direction:', 'Dimensions:', 'Resolution:',...
    'Offcenter:', 'FOV:', ...
    'Voxel shift due odd/even: ','Voxel shift due FFT-method:',...
    'Voxel limit (Min):', 'Voxel limit (Max):', ...
    'Volume limit (Min):', 'Volume limit (Max)','Voxel Shift:' '',''},...
   { '','','[AP/LR/FH]',...
     ori.dim,             ori.res,...
     ori.offcenter(1,:),  ori.fov,...
     ori.vox_cor,         ori.fft_cor,...
     ori.lim(:,1)',       ori.lim(:,2)',...
     ori.lim_vol(:,1)',   ori.lim_vol(:,2)',...
     shft,...
     '--------------------------------------------------------------',''});




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
elseif strcmp(mri.ext, 'ima')
    MRI_coordinates_IMA(gui, mri);
end

% Get updated appdata.
mri = getappdata(gui.CSIgui_main,'mri');

% Show data to user.
if ~isfield(mri,'ori'), return; end
mid_slice = ceil(size(mri.ori.offcenter,1)/2);

% Show details to user: Res, dims, FOV, Limits.
CSI_Log({'MRI-parameters -------------------------------', ...
    'Direction: [AP/LR/FH]', ...
    'Dimensions', 'Resolution: ', ...
    'Full FOV: ',...
    'Voxel limit (min): ','Voxel Limit (max): ',...
    'Offcentre (Middle): ', '', ''},...
    {'','', mri.ori.dim, mri.ori.res, ...
    sprintf(' %3.2f',mri.ori.fov),...
    sprintf(' %3.2f',mri.ori.lim(:,1)),...
    sprintf(' %3.2f',mri.ori.lim(:,2)),...
    sprintf(' %3.2f',mri.ori.offcenter(mid_slice,:)),...
    '---------------------------------------------------',' '});

% --- Executes on button press in button_MRI_PlotIMG_Grid.
function button_MRI_PlotIMG_Grid_Callback(~, ~, gui)
% Plot converted images with the grid superposed on the images. No spectra
% are shown.

% Plot.
MRI_plotImage_tabbed(gui);

% --- Executes on button press in button_MRI_PlotIMG_inCSI.
function button_MRI_PlotIMG_inCSI_Callback(hobj, ~, ~)
% Show image in current displayed MRS array.
MRI_plotImage_current_CSI(hobj);

% --- Executes on button press in button_MRI_Convert.
function button_MRI_Convert_Callback(~, ~, gui)
%
% Convert MRI images to CSI data and create appdata CONV.
%
% Requires:
%       image data and the coordinates meshgrid, 
%       csi limit, resolution and dimensions.

MRI_to_CSIspace(gui);


% MRI Functions % ------------------------------------------------------ %

% --- Executes on buttons press in button_IMGcoordinates *if IMA
function MRI_coordinates_IMA(gui, mri)
% Get all parameters from the header file and set up to calculate the
% image coordinates.

% Image orientation structure
imgori = struct;

if ~isfield(mri, 'par')
    CSI_Log({'No parameter field present in MRI-data.'},{'Aborted.'});
    return;
end

% Image Positient Patient: the coordinates of the top-left voxel of the
% images. The row and column ipp are identical for every slice. The slice 
% ipp changes. 
imgori.ipp = extractField(mri.par,'ImagePositionPatient');
tmp = cell2mat(imgori.ipp);

imgori.tlv = unique(tmp(1:2,:))';
imgori.tlv([2 1]) = imgori.tlv([1 2]); % Flip AP/LR (row, col)
imgori.tlv(3) = tmp(3,1);

% This is the coordinate vector for the slices.
% imgori.tlv_slice = tmp(3,:);

% Image resolution: the resolution of a voxel in  the images.
% First value is the row, second row the column as per DICOM-description
tmp = extractField(mri.par,'PixelSpacing');
tmp = cell2mat(tmp);
imgori.res = unique(tmp', 'rows');
imgori.res(3) = mri.par{1}.SliceThickness;

imgori.row = mri.par{1}.Rows;
imgori.col = mri.par{1}.Columns;

% Image dimensions: number of pixels in row and col.
imgtype = fields(mri.data);
imgori.dim = size(mri.data.(imgtype{1}));

% Image Field of View
imgori.fov = imgori.dim.*imgori.res;

% Slice gap-size
% Dont use this as is the spacing between the center of each adjacent slice
% and not the gap between two slices after adding the slicethickness. This
% also includes any overlap of slices during acquisition.
% imgori.gap = mri.par{1}.SpacingBetweenSlices;


% Process remaining ranges % -------------------------------------------- %
% Create range vectors of the image-volume (mm);
        
% For the Row and Column
vox_correction = 1;
for kk = 1:2
    % Select resolution and offcenter of current axis (row/col).
    N = imgori.dim(kk); res = imgori.res(kk);
    offc  = imgori.tlv(1,kk);
                
    % Top Left Voxel (TLV) Coordinate Definition:    
    % TLV is relative to top left corner wrt the axis origin and
    % increases towards the positibe side of the image. This requires
    % the addition of half a voxel to describe the centre of a voxel.
    if vox_correction == 0
        % NO center-of-voxel correction
        tmp = offc;                
    else
        % YES center-of-voxel correction
        tmp  = offc + (0.5*res);    
    end
    
    % Volume limits of current direction
    Vbegin = tmp; 
    Vend   = tmp + (res*N) - res;

    % Create grid vector and limits field.
    imgori.vec{kk}     = Vbegin:res:Vend;
    imgori.lim(kk,1:2) = imgori.vec{kk}([1 end]);
end

% Slice coordinate-vector.
imgori.vec{3} = cell2mat(cellfun(@(x) x(3),imgori.ipp, 'uniform', 0));
imgori.lim(3,1:2) = imgori.vec{3}([1 end]);

% Get offcenter values from coordinate vector.
imgori.offcenter = cellfun(@median, imgori.vec); 

% MESH
imgori.mesh.nb = 'Meshgrid: row and column ARE swapped.';
[imgori.mesh.x, imgori.mesh.y,  imgori.mesh.z] = ...
    meshgrid(imgori.vec{2}, imgori.vec{1}, imgori.vec{3});

% Save data
mri.ori = imgori;
setappdata(gui.CSIgui_main, 'mri', mri);

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

% --- Executes on button press in button_IMGcoordinates *if DCM
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
    img = NaN; 
    CSI_Log({'MRI Match Slices: No converted image data present.'},...
            {''});
    return;
end

if ~isappdata(gui.CSIgui_main,'csi')
    CSI_Log({'MRI Match Slices: No CSI data present.'},...
            {''});
    img = NaN; return; 
end

csi = getappdata(gui.CSIgui_main,'csi');
if ~isfield(csi,'ori')
    CSI_Log({'MRI Match Slices: No geometry info for MRS data present.'},...
            {''});
    img = NaN; return; 
end


        % ----- % Create MRS matching images % ----- %

nSlices = size(csi.data.raw, 4);
szData = size(conv.data);
img = NaN(szData(1),szData(2),nSlices);
for sli = 1:nSlices              % For every CSI slice

    % Slice coordinates for CSI and CONV
    cz = unique(csi.ori.mesh.z); mz = unique(conv.mesh.z);

    % Find CONV slices matching to CSI slices.
    img_range = CSI2MRI(cz(sli), mz, csi.ori.res(3), conv.res(3));
    if sli == 1, img_all_slice_range = NaN(nSlices, size(img_range,2)); end
    img_all_slice_range(sli,:) = img_range;
    
    % Minimum and maximum index
    indMn = img_range(1); indMx = img_range(2);
    
    % 08/01/24 QH
    if indMx > size(conv.data,3), indMx = size(conv.data,3); end

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

% Let user pick image type if possible.
if size(imtmp,2) > 1
    
    % Get imagetype of interest; imtoi;
    uans = getUserInput_Popup({'Image type to convert to MRSI space:'},...
                           {imtmp});
    if isempty(uans), CSI_Log({'Skipped image convert.'},{''}); return; end
    % Image type of interest
    imtoi = uans{1}; 

else
    imtoi = imtmp{1};
end
   

             % ---------- % Convert IMG % ---------- %

% Converted resolution equals original MR image resolution. (3D)
% However, to fit correctly in csi grid, the resolution is changed (below).
conv.res = mri.ori.res; % Initial... May change!!

% Calculate a resolution fitting the CSI space such that there is a integer
% amount of image pixels fitted in each CSI direction of space.
res_fit = csi.ori.res ./ conv.res;      % #MRpix / CSIpix
% res_rem = res_fit - floor(res_fit);     % Pixel change
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
[x,y,z] = meshgrid(conv.vec{2} ,conv.vec{1}, conv.vec{3});
conv.mesh.x = x; conv.mesh.y = y; conv.mesh.z = z; 

% Interp values @ CSI space % --------------------------------- %
% Interp3 is used to convert images to csi-grid in conv-struct.

% Check for stack availability - these images need seperate inclusion!
if isfield(mri.ori,'stack_of_interest')
    stack_ind = mri.ori.stack_of_interest;
else
    stack_ind = 1; 
end

if strcmp(mri.ext ,'par') || strcmp(mri.ext ,'ima')   
    image_convert_data = mri.data.(imtoi)(:,:,:);
else
    image_convert_data = mri.data.(imtoi)(:,:,stack_ind,:);
end

% Interp
conv.data = interp3(mri.ori.mesh.x, mri.ori.mesh.y, mri.ori.mesh.z,...      % Original MRI coordinates
                    squeeze(image_convert_data),...                         % Original MRI values
                    conv.mesh.x, conv.mesh.y, conv.mesh.z,'Linear',0);      % Requested coordinates
conv.dim  = size(conv.data);

% Export any previous contrast info for image display.
if isfield(mri,'contrast'), conv.contrast = mri.contrast; end


% Post processing MRSI % ------------------------------------ %
% Apply required rotations or flips.

% % Apply additional rotation if it is list/data file
if strcmp(csi.ext,'.list') || strcmp(csi.ext,'.data')

end

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
    {'Which images to save: '},{{'All','Converted','MRS Matching'}}, [],...
    'Save MRI');
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
% clr_tx = gui.colors.text_main; 
clr_tb = gui.colors.text_title;

% PLOT_PAR: STRUCT WITH ALLL 2D PLOT OPTIONS
plot_par.colors = gui.colors;
                   % -------- % Figure: create window % -------- %

% Create figure
fh = figure('Tag', tag ,'Name', 'CSIgui - Tabs',...
            'Color', clr_bg, 'Toolbar', 'None', 'MenuBar', 'None',...
            'NumberTitle', 'Off');                   

% Create tab group
tabg = uitabgroup(fh); % tabh = cell(1,nSlices);

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

% --- Executes on button press in button_MRI_Delete_IMG_Data.
function button_MRI_Delete_IMG_Data_Callback(~, ~, gui)

% Delete MRI and CONV appdata 
if isappdata(gui.CSIgui_main, 'mri')
    rmappdata(gui.CSIgui_main, 'mri'); 
    gui.txt_fnIMG.String = '';
end
if isappdata(gui.CSIgui_main, 'conv')
    rmappdata(gui.CSIgui_main, 'conv');
end





% PLOT 2D CSI % -------------------------------------------------------- %


% --- Executes on button press in button_plotCSI.
function button_plotCSI_Callback(hobj, evt, gui)
% Show CSI button; to display CSI and possible converted MRI data.
% Launches:
%       CSIpanel_2D_DataToDisplay for navigation through all dimensions
%       CSI_plot2D to open up the CSIgui 2D figure and plot MRS data

if isappdata(gui.CSIgui_main, 'csi')
    panel_2D_DataSliders(hobj, evt, gui);
    CSI_2D_initiate2D(hobj, evt, gui);
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

% ------- % Plot data using options from plot_par
CSI_2D_plotVoxels(plot_par,gui);

% ------ % Set the figure ratio to voxel size
% CSI_2D_setFigure_ratio(gui);

% ------ % Let Scroll-panel follow 2D CSIgui plot figure
if strcmp(gui.menubar.MRSI.snapWindow.main.Checked,'on')
    panel_2D_followPlot2D_initiate(plot_par.fh);
end

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

% SNAP 2 PLOT % ------------------------------------ % DEV Experimental
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


% Get FULL data cell indexing % --------- %
% Used to get all data including all slices.

% Create cell for indexing outside the spatial-dimensions 
% if numel(plot_par.dim) > 3, nDimC = num2cell(plot_par.dim(4:end));
% else, nDimC = {1};
% end
% % To linear vector per cell.
% nDimC = cellfun(@(x) 1:x, nDimC,'UniformOutput',0);

% Store all other indexes..
if numel(plot_par.dim) > 2, nDimCtmp = num2cell(plot_par.dim(3:end));
else, nDimCtmp = {1};
end
nDimCtmp = cellfun(@(x) 1:x, nDimCtmp,'UniformOutput',0);
plot_par.select_all_dim = nDimCtmp;


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
if strcmp(gui.menubar.MRSI.ColorScale.ScalebyWindow.Checked,'on')
    plot_par.scale_by_window_color = 1;    
else
    plot_par.scale_by_window_color = 0;
    plot_par.scale_range_color = [1 size(data_volume,1)];
end


% Color scaling by SLICE, VOLUME or STATIC.
vol_data = CSI_getUnit(data_volume, plot_par.data_unit);
full_scale_range =  ...
    plot_par.scale_range_color(1):plot_par.scale_range_color(2);

% Correct volume data used to calculate specific bins for different colors
% scaled by minimum and maximum values: In x-axis window or from full
% spectral field of view.
if plot_par.scale_by_window_color
    % #N, X, Y, Slice (Z), non-spatial dimensions.
    vol_data = vol_data(full_scale_range,:,:,plot_par.select_all_dim{:}); 
end

scaleby = CSI_2D_Scaling_Color_Get(gui);
switch scaleby
    case 'vol' % Scale by volume
        % Plot color scaled for all voxels in the *volume*
        max_per_voxel = max(vol_data,[],1);
        data_ylimits_color = [min(max_per_voxel(:)), max(max_per_voxel(:))]; 
        plot_par.scale_color = 0;
    case 'sli' % Scale by slice
        % Calc limits in slice        
        max_per_voxel = max(vol_data(:,:,:,plot_par.plotindex{:}),[],1);
        data_ylimits_color = [min(max_per_voxel(:)) max(max_per_voxel(:))];           
        plot_par.scale_color = 1;
    case 'sta' % Static color
        % Color scaling limited to 1 color: static line color
        data_ylimits_color = NaN;
        plot_par.scale_color = 2;
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
    % tmp = vol_data(:,:,:,plot_par.plotindex{:}); 
    plot_par.clrs = gui.colors.lines1; 
    plot_par.clrs_data_range = max(vol_data(:)); 
end

% --- % Executed by CSI_plot2D_initiate: get plot2D settings
function plot_par = CSI_2D_getPlotSettings(plot_par, gui, data_volume)
% Add to structure plot_par the following plot-settings and plot-data
% fields: 
% 
% Scaling plot color, axis-y and x limit (by volume/static/voxel) and
% more.
% Output fields to plot_par:
%               scale_by_window_color   Color scale data by window or full 
%                                       spectra
%               scale_color             Scale type of colors: 
%                                       vol(0/sli(1/sta(2
%               clrs                    Plot colors for color scaling 
%                                       spectra
%               xlimit                  Visual x-limits of spectra plots
%               xaxisdata               Xaxis plot data (ppm/unitless)
%               scale_by_window_axis    Scale y-axis by window or full
%                                       spectrum
%               scale_range             Index scale range for previous
%               scale                   Scale type: vox(0/sli(1/vol(2
%               axScale_ylimit          Y-limits for scaling vol/sli
%               voxel_grid              Enable or disable grid in each vox


% Axis: scale by window and X-limits% ---------- %

% Get visual and index range of x-axis data to use for scale by window of
% the y-axis.

% X-axis visual limits
plot_par.xlimit = plot_par.xaxis.xlimit;

% Get the frequenty axis
% unitless(none) or frequency(ppm) and get unitless range of data to 
% scale y-axis; to full spectrum or to visible part of spectrum.
if isfield(plot_par.xaxis, 'ppm')
    plot_par.xaxisdata = plot_par.xaxis.ppm;
    
    [~,scale_range(1)] =  ...
        min(abs(plot_par.xaxis.ppm - plot_par.xaxis.xlimit(1)));
    [~,scale_range(2)] =  ...
        min(abs(plot_par.xaxis.ppm - plot_par.xaxis.xlimit(2)));
   
else         
    
    plot_par.xaxisdata = plot_par.xaxis.none; 
    [~,scale_range(1)] =  ...
        min(abs(plot_par.xaxis.none - plot_par.xaxis.xlimit(1)));
    [~,scale_range(2)] =  ...
        min(abs(plot_par.xaxis.none - plot_par.xaxis.xlimit(2)));
end

% If only a single value... 
if diff(scale_range) > size(data_volume,1)
    scale_range = [1 size(data_volume,1)];
end

% Save scale range - for color window of spectra.
plot_par.scale_range_color = scale_range;

% If user request Y-axis scaling by full spectrum - set full index as new
% scaling range.
if strcmp(gui.menubar.MRSI.AxisScale.ScalebyWindow.Checked,'off')
    plot_par.scale_range_axis = [1 size(data_volume,1)];
    plot_par.scale_by_window_axis = 0;
else
    plot_par.scale_range_axis = scale_range;
    plot_par.scale_by_window_axis = 1;
end


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
scale_range = plot_par.scale_range_axis;
switch scaleby 
    case 'vox', plot_par.scale = 0; 
    case 'sli', plot_par.scale = 1; 
        % Unit
        tmp_data = CSI_getUnit(plot_par.data2D, plot_par.data_unit);
        % Correct x-axis window
        tmp_data = tmp_data(scale_range(1):scale_range(2),:,:);
        plot_par.axScale_ylimit = [min(tmp_data(:)) max(tmp_data(:))]; 
    case 'vol', plot_par.scale = 2; 
        % Unit
        tmp_dataVol = CSI_getUnit(data_volume, plot_par.data_unit);
        % Correct x-axis window
        tmp_dataVol = tmp_dataVol(scale_range(1):scale_range(2),...
            :,:,plot_par.select_all_dim{:});
        plot_par.axScale_ylimit = [min(tmp_dataVol(:)) max(tmp_dataVol(:))];
end

 
% Voxel grid on or off % ---------------- %
if strcmp(gui.menubar.MRSI.AxisScale.Grid.Checked, 'on')
    plot_par.voxel_grid = 1;
else
    plot_par.voxel_grid = 0;
end


% Scaling: Plot Color % ------------- %
% Create - plot_par.clrs, plot_par.clrs_data_range
plot_par = CSI_2D_getPlotSettings_ColorScaling(plot_par, data_volume);



% --- % Executed by CSI_plot2D_initiate: plot images
function plot_par = CSI_2D_plotImages(plot_par, csiguiObj)
% GUI of 2D figure

            % --------------% IMAGES TO PLOT % --------------%
                  % %%% MRI MATCHED TO CSI PLOT %%% %

% Get MRI image data and set to CSI space to plot.
if isappdata(csiguiObj, 'conv')
    conv = getappdata(csiguiObj, 'conv');
    try 
        
        if isfield(plot_par, 'img_all')
            plot_par.img = plot_par.img_all(:,:,plot_par.plotindex{1});
        else % No image-data available in plot-parameters
            
            % Get image @ this slice: uses CSI2MRI
            img = MRI_matchSlices(csiguiObj);
                         
            if (length(img) == 1) && isnan(img)
                CSI_Log({'Converting image data to MRS-space failed.'},...
                        {'Canceled image plot.'});
                return;
            else                
                plot_par.img_all = img;              
            end

            % Set current image
            plot_par.img = plot_par.img_all(:,:,plot_par.plotindex{1});
        end

        % Plot IMG at CSI-coordinates! %%%%
        if ~isfield(plot_par, 'imax')
            plot_par.fh.Units = 'Pixels';
            ax_sz = plot_par.fh.Position(3:4);
            plot_par.imax = ...
                axes('parent',plot_par.fh,...
                'Unit','pixels','InnerPosition',[1 1 ax_sz], 'Color', 'blue',...
                'XTick', [],'YTick', []);
            
        end
        plot_par.imax.Units = 'normalized';
        
        % hold(plot_par.imax, 'on'); 
        if (sum(plot_par.img(:)) == 0) % Image is only zeroes.
            imagesc(plot_par.imax,plot_par.img); 
            colormap(plot_par.imax,gray(255)); 
            set(plot_par.imax,'Color', 'Black'); alpha(plot_par.imax,0); 
        else
            % Image plotting
            % Imscale plots over the entire figure and does not
            % imply any border issues as with imshow-function.
           imagesc(plot_par.imax,plot_par.img); 
                      
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
                % Calculate contrast
                if min(plot_par.img(:)) >= max(plot_par.img(:)) 
                    caxis(plot_par.imax,...
                        [min(plot_par.img(:)) min(plot_par.img(:))+1]);
                else
                    caxis(plot_par.imax,...
                        [min(plot_par.img(:))-1 0.5*max(plot_par.img(:))]);
                end
            end
            colormap(plot_par.imax, gray(255));
        end
        
        % Remove axis-decorations
        grid(plot_par.imax,'off');
        plot_par.imax.Box = 'off';
        plot_par.imax.XTick = [];
        plot_par.imax.YTick = [];
        
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

% Scale range window of x-axis
full_scale_range =  ...
    plot_par.scale_range_axis(1):plot_par.scale_range_axis(2);

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
        if createFig % Only if figure did not exist before...
            plot_par.ax{ri,ci} = axes('parent',plot_par.fh,'position',pos);
        end

        
        % VOXEL DATA % ----------------- %
        plot_data_vox = plot_par.data2D(:,ci,ri);
        plot_data_vox = CSI_getUnit(plot_data_vox, plot_par.data_unit);
        
        
        % COLOR SCALING % -------------- %
        
        % This is NOT the Ylimit eventually used for axis-display.
        % See below plot() @ YAXIS SCALING
        ylimit = [min( plot_data_vox( : ) ) max( plot_data_vox( : ) )];
        
        % Relative to maximum Y-data.
        [~, clr_ind] = min(abs(plot_par.clrs_data_range-ylimit(2)));
        plot_color = plot_par.clrs(clr_ind,:); % See before ri/ci for-loops.
  
       
        %%%% PLOT CSI data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(plot_par.ax{ri,ci}, plot_par.xaxisdata, plot_data_vox,...
                                'color', plot_color, 'LineWidth', 1.25);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % YAXIS SCALING % -------------- %
        
        % Scale Y axis by limits of voxel, slice or volume.
        if plot_par.scale == 0                                              % Voxel 
            ylimit = ...
                [min( plot_data_vox( full_scale_range ) )...
                 max( plot_data_vox( full_scale_range ) )];
        elseif plot_par.scale == 1 || plot_par.scale == 2                   % Slice or Volume
            ylimit = plot_par.axScale_ylimit;                               % Use pre-calculated limit.
        end
            
        % Largest absolute limit: min/max.
        if     abs(ylimit(2)) > abs(ylimit(1)), ylimfac = abs(ylimit(2));
        elseif abs(ylimit(1)) > abs(ylimit(2)), ylimfac = abs(ylimit(1));
        else, ylimfac = 1;
        end
       
        % Create visualised Ylimit. Centered around zero.
        ylimit = [-1.05*ylimfac  1.05.*ylimfac];
        if ylimit(2) <= ylimit(1), ylimit(2) = ylimit(1)+1; end % Safety 
        ylim(plot_par.ax{ri,ci}, ylimit); % Set Ylimit

        
        % X LIMIT % -------------------- %
        
        % Use limits of ppm axis(1) and (end) plus reverse the x-axis
        % direction to correctly display MRSI data.
        xlim(plot_par.ax{ri,ci}, sort(plot_par.xlimit(1:2))); 
        set(plot_par.ax{ri,ci}, 'xdir', 'reverse');
        
        
        % VOXEL GRID % ------------------ %
        
        if plot_par.voxel_grid
            plot_par.ax{ri,ci} = ...
                CSI_2D_plotVoxels_grid(plot_par.ax{ri,ci},...
                                       plot_par.colors.grid2);
        end
        

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
% Why: Due to a axis linewidth error in Matlab, the anchor point of the
% position of the axis is not the bottom left point but the center of 
% the axis linewidth. This create descripancies in the plot layout due to
% pixel/normalize/inches conversions to on the screen display - especially 
% when saving figures. Using seperate grids, this error is not visible.
if createFig
    CSI_2D_grid(plot_par.fh,  plot_par.fh.Position(3:4), ...
                plot_par.dim, plot_par.range, gui.colors.grid);
end

% Bring figure to front.
figure(plot_par.fh);

% --- Calculate and add a grid in each voxel in 2D-plots
function ax = CSI_2D_plotVoxels_grid(ax, clr)
% Add a grid to a voxel in the 2D CSI plot. Used in csi_2D_plotVoxels and
% requires the axes handle to the voxel and a grid color;

% Prep options                                  ----- !Can Be Edited!
% Number of grid lines, grid linewidth and grid alpha opacity; latter
% overwrites userinput.
gridN = 5; gridW = 0.5; 
gridAlpha = 0.75; clr(4) = gridAlpha;


% Get current axis limits
yl = ax.YLim; xl = ax.XLim;

% Calulate (mesh) grid-cooridnates for voxel
xstep = diff(xl)/gridN; ystep = diff(yl)/gridN;
xv = xl(1):xstep:xl(2); yv = yl(1):ystep:yl(2);
[xgh, ygh] = meshgrid(xv,yv); xgv = xgh'; ygv = ygh';

% Plot grid
hold(ax,'on');
plot(ax,xgh(:,2:end-1),ygh(:,2:end-1),...
    'r-','Linewidth',gridW,'Color',clr);
plot(ax,xgv(:,2:end-1),ygv(:,2:end-1),...
    'r-','Linewidth',gridW,'Color',clr);    
hold(ax,'off');

% Reverse plot-layering e.g. plotted data on top.
ax.Children = flip(ax.Children,1);



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
wv = 1; wh = wv; wv = wv./w; wh = wh./h; 
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
    private_instance = 1 ; 
    private_tag = string(datetime('now','Format','HHmmss'));
elseif mouse_button == 3 % Right
    private_instance = 0 ; 
    private_tag = string(datetime('now','Format','HHmmss'));
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

% Because of JavaFrame, surpress any warnings about it...
warning('off', 'all');

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
end

warning('on', 'all'); 

function panel_2D_followPlot2D_stop(hobj)
% Disable the slider window to follow the plot2D window by removing the
% callback from the plot2D window. 

% Because of JavaFrame, surpress any warnings about it...
warning('off', 'all');

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
    set(jbh,'ComponentMovedCallback','');
catch err
    CSI_Log({'Error in followPlot2D_stop'},{err.Message});
end

warning('on', 'all');

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


% --- % Executed if user click snap Window in menubar
function panel_2D_followPlot2D_menubar(~, ~)
% Set checkmark on or off in menubar when clicking menubar > snap Window
% Enable or disable the corresponding function that snaps the slider 
% window to the plot2D figure.

% CSIgui data
hobj = findobj('Tag', 'CSIgui_main', 'Type', 'Figure');
if isempty(hobj), return; end
gui = guidata(hobj);

% Get state of button
state = gui.menubar.MRSI.snapWindow.main.Checked;
switch state
    case 'on',  gui.menubar.MRSI.snapWindow.main.Checked = 'off';
        state = 'off';
    case 'off', gui.menubar.MRSI.snapWindow.main.Checked = 'on';
        state = 'on';
    otherwise,  gui.menubar.MRSI.snapWindow.main.Checked = 'on';
        state = 'on';
end

fobj = findobj('Tag', 'CSIgui_plot2D', 'Type', 'Figure');
if ~isempty(fobj) && strcmp(state,'on')
    panel_2D_followPlot2D_initiate(fobj);
elseif  ~isempty(fobj) && strcmp(state,'off')
    panel_2D_followPlot2D_stop(fobj);
end



% ---------------------------------------------------- %
% PLOT 2D SCALING % -------------------------------------------------- %
% ---------------------------------------------------- %

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

% xaxis data.
plot_par.xaxis = csi.xaxis;

% Axis size to fit figure
plot_par.dim      = size(csi.data.raw);   % Data dimensions
plot_par.dim(1)   = [];                   % Remove time index e.g. 1
plot_par.data_dim = numel(plot_par.dim);  % 3D/2D/1D volume. 

% Get all necessary data.
plot_par = CSI_2D_getData(plot_par, gui, csi.data.raw);
plot_par = CSI_2D_getPlotSettings(plot_par, gui, csi.data.raw);

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
    
% --- Executed by user in menubar: CSI > Scaling Color > Scale by window
function CSI_2D_Scaling_Color_ScaleByWindow(hObject, ~)
% Enabel or disbale scale by Color.
gui = guidata(hObject);

if strcmp(get(gui.menubar.MRSI.ColorScale.ScalebyWindow,'Checked'), 'off')
    set(gui.menubar.MRSI.ColorScale.ScalebyWindow,'Checked', 'on');
else
    set(gui.menubar.MRSI.ColorScale.ScalebyWindow,'Checked', 'off');
end

% 3. Update GUI data.
guidata(hObject,gui);

function CSI_2D_Scaling_Axis_ScaleByWindow(hObject, ~)
% Enabel or disbale scale by Color.
gui = guidata(hObject);

if strcmp(get(gui.menubar.MRSI.AxisScale.ScalebyWindow,'Checked'), 'off')
    set(gui.menubar.MRSI.AxisScale.ScalebyWindow,'Checked', 'on');
else
    set(gui.menubar.MRSI.AxisScale.ScalebyWindow,'Checked', 'off');
end

% 3. Update GUI data.
guidata(hObject,gui);

function CSI_2D_Scaling_GridVoxels(hObject, ~)
% Enabel or disbale scale by Color.
gui = guidata(hObject);

if strcmp(get(gui.menubar.MRSI.AxisScale.Grid,'Checked'), 'off')
    set(gui.menubar.MRSI.AxisScale.Grid,'Checked', 'on');
else
    set(gui.menubar.MRSI.AxisScale.Grid,'Checked', 'off');
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
% Creates struct xaxis which contains all frequency information and
% corresponding axis for visualization.
%
% If auto equals 1, no user answer is required; 
%   Works if specific frequency info is available. (Userinput or spar)
%   Without this data, only arbitrary units are returned.

% If no automatic option is given, turn if off.
if nargin <= 2, auto = 0; end

% Get gui-data
gui = guidata(hObject); 

% Require csi: data, samples, csi
csi = getappdata(gui.CSIgui_main,'csi'); 
if isempty(csi), return; end

% x-axis structure
xaxis = struct;

% Store previous static data
if isfield(csi, 'xaxis') % If already present, copy some data.
    foi = {'BW','nucleus','trans','gyro','tesla', 'shift'};
    for kk = 1:size(foi,2)
        tmp = foi{kk};
        if isfield(csi.xaxis, tmp)
            xaxis.(tmp) = csi.xaxis.(tmp);
        end
    end
end

% READ SPAR-header (Philips) % ------------------------------- %
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

% READ DAT-Header (Siemens/TWIX) -------------------------------- %
if isfield(csi, 'twix') && ~isfield(csi, 'xaxis')
    % Nucleus
    xaxis.nucleus = csi.twix.Config.Nucleus;
    
    % Magnet strength
    if isfield(csi.twix, 'Meas')
        if isfield(csi.twix.Meas,'flNominalB0')
            xaxis.tesla = csi.twix.Meas.flNominalB0;
        end
    elseif isfield(csi.twix, 'Dicom')
        if isfield(csi.twix.Dicom,'flMagneticFieldStrength')
            xaxis.tesla = csi.twix.Dicom.flMagneticFieldStrength;
        end
    end
    
    % Bandwidth (Receiver)
    if isfield(csi.twix.Config, 'DwellTime')
        dwelltime = csi.twix.Config.DwellTime * 1e-9;
    else       
        vals = getFieldValues(csi.twix, 'dwelltime');
        ind = ~cellfun(@ischar, vals); vals = cell2mat(vals(ind));
        vals = unique(vals);
        if numel(vals) > 1, vals = mode(vals); end
        dwelltime = vals * 1e-9;
    end
    xaxis.BW = 1/dwelltime;
    
    % Oversampling Correction    
    if isfield(csi.twix.Config, 'ReadOSFactor')
        OSfactor = csi.twix.Config.ReadOSFactor;
    elseif isfield(csi.twix.Config,'ReadoutOSFactor')
        OSfactor = csi.twix.Config.ReadoutOSFactor;
    else
        OSfactor = 1;
    end
    OSvectorsz = csi.twix.Config.VectorSize*OSfactor;           
    if OSvectorsz ~= size(csi.data.raw,1)
        xaxis.BW = xaxis.BW ./ OSfactor;
    end
    
    switch xaxis.nucleus 
        case '1H',  xaxis.gyro = 42.57747892*10^6;
        case '2H',  xaxis.gyro = 6.536*10^6;
        case '31P', xaxis.gyro = 17.235*10^6;
        case '23Na',xaxis.gyro = 11.262*10^6;
        case '19F', xaxis.gryo = 40.052*10^6;
        otherwise,  xaxis.gyro = 42.57747892*10^6;
    end 

    % Get imaging frequency
    vals = getFieldValues(csi.twix, 'Frequency', 0);
    ind = ~cellfun(@ischar, vals); vals = cell2mat(vals(ind));
    vals = unique(vals);
    if numel(vals) > 1, vals = mode(vals); end
    if vals ~= 0
        xaxis.trans = vals;
    else
        xaxis.trans = xaxis.tesla * xaxis.gyro;    
    end
end

% USER INPUT % ------------------------------------------ %
% If auto == 1, this is skipped

if auto == 0
    % Request required values
    qry = { 'Nucleus (1H): ', 'Magnet strength(T): ', ...
            'Bandwidth (Hz): ', 'PPM Shift (Num): '} ;

    % Create possible answers if exist: previous input or from header file.
    if isfield(xaxis,'nucleus'),an{1} = xaxis.nucleus; 
    else, an{1}= '31P';  
    end    
    if isfield(xaxis,'tesla'),  an{2} = xaxis.tesla;   
    else, an{2}= '7';   
    end
    if isfield(xaxis,'BW'),     an{3} = xaxis.BW;      
    else, an{3}= '8000'; 
    end
    if isfield(xaxis,'shift'),  an{4} = xaxis.shift;   
    else, an{4}= '0';    
    end

    % Display UI to get user input
    inp = getUserInput(qry, an);
    if isempty(inp)
        CSI_Log({'Skipped setting parameters.'},{''}); return; 
    end

    % Set user input.
    xaxis.nucleus = inp{1};
    xaxis.tesla   = str2double(inp{2});
    xaxis.BW      = str2double(inp{3});
    xaxis.shift   = str2double(inp{4});

    switch inp{1} % Gyromagnetic constants in Hz
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

% -------------------------------------------------------- %
% CALC AXES % -------------------------------------------- %
% -------------------------------------------------------- %

%%% Calculate axis-parameters

% Find dimension with FID/Spectrum
fid_dim_label = {'sec','fid', 'spec', 't'}; 
ind = csi_findDimLabel(csi.data.labels,fid_dim_label);
ind(isnan(ind)) = [];

% Set length of spectro-data
if ~isempty(ind), xaxis.N = csi.data.dim(ind);
else % If not found - use largest dimension            
    [~, tmp] = max(csi.data.dim);
    xaxis.N = csi.data.dim(tmp);
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
    if ~isfield(xaxis,'xlimit') || (xaxis.xlimit(2) ~= csi.data.dim(1)-1)
        xaxis.xlimit = [xaxis.none(1) xaxis.none(end)];
    end
end

% CLEAN UP % -------------------------------------------- %

% Add CSI-output.
csi.xaxis = xaxis; % Contains frequency details of MRS data

% Save in appdata
setappdata(gui.CSIgui_main,'csi',csi);

% ---------------------------------------------------- %
% (Additional window for editting 2D plot options.) ------------------ %
% ---------------------------------------------------- %

% --- Executes on button press in button_CSI_DisplayOptions.
function button_CSI_DisplayOptions_Callback(~, ~, gui)
CSI_2D_Scaling_Options(gui); % Set all options for 2D plot using a GUI.

% --- See button_CSI_DisplayOptions_Callback
function CSI_2D_Scaling_Options(gui)


% Get CSI app data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

% CURRENT 2D-PLOT SETTINGS % ------------------------------------------- %

% Color scaling % --- %
% Create option per setting and sort with current setting at top.
current_clr = CSI_2D_Scaling_Color_Get(gui);
switch current_clr
    case 'vol', current_clr = 'Volume';
    case 'sli', current_clr = 'Slice';
    case 'sta', current_clr = 'Static';
end
% Create option per setting and sort with current setting at top.
inp_clr = {'Voxel', 'Volume', 'Static'};
ind = contains(inp_clr, current_clr); inp_clr = {inp_clr{ind} inp_clr{~ind}};

% Scale by x-axis window: on/off
cl_scaleWindow = gui.menubar.MRSI.ColorScale.ScalebyWindow.Checked;
if strcmp(cl_scaleWindow,'on'), cl_scaleWindow = 1; 
else,                           cl_scaleWindow = 0; 
end


% Y-axis scaling % --- %
% Create option per setting and sort with current setting at top.
current_axs = CSI_2D_Scaling_Axis_Get(gui);
switch current_axs
    case 'vox', current_axs = 'Voxel';
    case 'sli', current_axs = 'Slice';
    case 'vol', current_axs = 'Volume';
end
inp_axs = {'Voxel','Slice','Volume'};
ind = contains(inp_axs, current_axs); inp_axs = {inp_axs{ind} inp_axs{~ind}};

% Scale by x-axis window: on/off
ax_scaleWindow = gui.menubar.MRSI.AxisScale.ScalebyWindow.Checked;
if strcmp(ax_scaleWindow,'on'), ax_scaleWindow = 1; 
else,                           ax_scaleWindow = 0; 
end

% Voxel Grid % --- %
ax_grid = strcmp(gui.menubar.MRSI.AxisScale.Grid.Checked, 'on');

% X-Limit % --- %
xlimit = num2str(csi.xaxis.xlimit);     


% ---------------------------------------------------------------------- %


% 2D PLOT OPTIONS GUI % --------------------- %
% Create figure and the UI elements.

% Figure
tg = 'CSI_2D_OptionsGUI';
nm = 'CSIgui - 2D-Plot Display Options';
fh = figure('Name',nm,'Tag',tg,'Menubar','none','Toolbar','none',...
    'Color', gui.colors.main,'Unit','Pixels','NumberTitle', 'Off');
gdat = guidata(fh); axis off;

% Figure position
fig_wh = [400 250];
sz = gui.CSIgui_main.Position; fig_pos = sz(1:2)+((sz(3:4)-fig_wh)./2);
fh.Position = [fig_pos fig_wh];



% --- % Title static text % --- % 

gdat.title.displayOptions = ...
    uicontrol(fh, 'Style', 'Text','Unit', 'Pixels',...
    'Position',[10 230 240 15],...
    'String','Plot 2D - Display Options',...
    'Fontsize', 10, 'FontWeight','Bold',...
    'Foregroundcolor', gui.colors.hilight1,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');


% --- % Axis scaling % --- % 

%Popup menu title
gdat.title.scaleAxis = ...
    uicontrol(fh, 'Style', 'Text','Unit', 'Pixels',...
    'Position',[50 202 100 15],...
    'String','Y-Axis scaling','Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_title,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');

% Dropdown menu title    
gdat.popup.scaleAxis = ...
    uicontrol(fh, 'Style', 'popupmenu','Unit', 'Pixels',...
    'Position',[50 180 150 20],...
    'String',inp_axs,'Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_main,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');

% Radio button scale by x-axis window.
gdat.radio.ax_scaleWindow = ...
    uicontrol(fh, 'Style','radiobutton','Unit', 'Pixels',...
            'position', [50 155 125 20],'value',ax_scaleWindow,...
            'String','Scale by Window',...
            'Foregroundcolor', gui.colors.text_main,...
            'BackgroundColor', gui.colors.main);

% Radio button grid.
gdat.radio.grid = ...
    uicontrol(fh, 'Style','radiobutton','Unit', 'Pixels',...
            'position', [50 135 125 20],'value',ax_grid,...
            'String','Grid per voxel',...
            'Foregroundcolor', gui.colors.text_main,...
            'BackgroundColor', gui.colors.main);

% --- % Color Scaling % --- % 

%Popup menu title
gdat.title.scaleColor = ...
    uicontrol(fh, 'Style', 'Text','Unit', 'Pixels',...
    'Position',[50 110 100 15],...
    'String','Color scaling','Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_title,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');

% Dropdown menu title    
gdat.popup.scaleColor = ...
    uicontrol(fh, 'Style', 'popupmenu','Unit', 'Pixels',...
    'Position',[50 90 150 20],...
    'String',inp_clr,'Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_main,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');
   
% Radio button scale by x-axis window.
gdat.radio.cl_scaleWindow = ...
    uicontrol(fh, 'Style','radiobutton','Unit', 'Pixels',...
            'position', [50 65 125 20],'value',cl_scaleWindow,...
            'String','Scale by Window',...
            'Foregroundcolor', gui.colors.text_main,...
            'BackgroundColor', gui.colors.main);
             

  
% --- % X-Limits % --- %

% X-Limit edits title
gdat.title.xlimits = ...
    uicontrol(fh, 'Style', 'Text','Unit', 'Pixels',...
    'Position',[250 202 100 15],...
    'String','X-Axis limits','Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_title,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Left');

gdat.edit.xlimit = ...
    uicontrol(fh, 'Style', 'edit','Unit', 'Pixels',...
    'Position',[250 180 100 22],...
    'String',xlimit,'Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_main,...
    'BackgroundColor', gui.colors.main,...
    'HorizontalAlignment', 'Center');


% --- % Buttons % --- % 

% Apply button
gdat.button.Apply = ...
    uicontrol(fh, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[275 30 75 20],'String','Apply','Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_main,...
    'BackgroundColor', gui.colors.main,...
    'Callback', @CSI_2D_OptionsGUI_Apply);

% Update button
gdat.button.Update = ...
    uicontrol(fh, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[200 30 75 20],'String','Update','Fontsize', 8,...
    'Foregroundcolor', gui.colors.text_main,...
    'BackgroundColor', gui.colors.main,...
    'Callback', @CSI_2D_OptionsGUI_Update);

% Store gui-data.
guidata(fh, gdat);
% --- See CSI_2D_Scaling_Options GUI
function CSI_2D_OptionsGUI_Apply(hobj, ~)
% Process data from Plot 2D - Display GUI.

% Get GUIdata
gdat = guidata(hobj);

csiobj = findobj('Name','CSIgui','Tag', 'CSIgui_main', 'Type', 'Figure');
if isempty(csiobj), return; end
gui = guidata(csiobj); csi = getappdata(gui.CSIgui_main,'csi');

% PROCESS INPUT % -------------------------------------- %

% Axis scaling
ax_ans = gdat.popup.scaleAxis.String{gdat.popup.scaleAxis.Value};
evt.Source.Text = ax_ans; evt.Source.Label = ax_ans;
CSI_2D_Scaling_Axis_Set(gui.CSIgui_main, evt);

% Color scaling
cl_ans = gdat.popup.scaleColor.String{gdat.popup.scaleColor.Value};
evt.Source.Text = cl_ans; evt.Source.Label = cl_ans;
CSI_2D_Scaling_Color_Set(gui.CSIgui_main, evt);

% Grid per voxel
if gdat.radio.grid.Value
    gui.menubar.MRSI.AxisScale.Grid.Checked = 'on';
else
    gui.menubar.MRSI.AxisScale.Grid.Checked = 'off';
end

% Process x-axis range answer
xl = gdat.edit.xlimit.String;
csi.xaxis.xlimit = sort(str2double(strsplit(xl,' ')));

% Scale by window for Y-axis and Color scaling
if gdat.radio.cl_scaleWindow.Value
    gui.menubar.MRSI.ColorScale.ScalebyWindow.Checked = 'on';
else
    gui.menubar.MRSI.ColorScale.ScalebyWindow.Checked = 'off';
end
if gdat.radio.ax_scaleWindow.Value
    gui.menubar.MRSI.AxisScale.ScalebyWindow.Checked = 'on';
else
    gui.menubar.MRSI.AxisScale.ScalebyWindow.Checked = 'off';
end

% Store
guidata(gui.CSIgui_main, gui);
setappdata(gui.CSIgui_main,'csi',csi);
% --- See CSI_2D_Scaling_Options GUI
function CSI_2D_OptionsGUI_Update(hobj,~)

% Get GUIdata
gdat = guidata(hobj);

csiobj = findobj('Name','CSIgui','Tag', 'CSIgui_main', 'Type', 'Figure');
if isempty(csiobj), return; end
gui = guidata(csiobj); csi = getappdata(gui.CSIgui_main,'csi');

% Current CSI 2D Plot Options % ----------------------------------- %

% Color scaling % --- %
% Create option per setting and sort with current setting at top.
current_clr = CSI_2D_Scaling_Color_Get(gui);
switch current_clr
    case 'vol', current_clr = 'Volume';
    case 'sli', current_clr = 'Slice';
    case 'sta', current_clr = 'Static';
end
% Create option per setting and sort with current setting at top.
inp_clr = {'Voxel', 'Volume', 'Static'};
ind = contains(inp_clr, current_clr); inp_clr = {inp_clr{ind} inp_clr{~ind}};

% Scale by x-axis window: on/off
cl_scaleWindow = gui.menubar.MRSI.ColorScale.ScalebyWindow.Checked;
if strcmp(cl_scaleWindow,'on'), cl_scaleWindow = 1; 
else,                           cl_scaleWindow = 0; 
end

% Y-axis scaling % --- %
% Create option per setting and sort with current setting at top.
current_axs = CSI_2D_Scaling_Axis_Get(gui);
switch current_axs
    case 'vox', current_axs = 'Voxel';
    case 'sli', current_axs = 'Slice';
    case 'vol', current_axs = 'Volume';
end
inp_axs = {'Voxel','Slice','Volume'};
ind = contains(inp_axs, current_axs); inp_axs = {inp_axs{ind} inp_axs{~ind}};

% Scale by x-axis window: on/off
ax_scaleWindow = gui.menubar.MRSI.AxisScale.ScalebyWindow.Checked;
if strcmp(ax_scaleWindow,'on'), ax_scaleWindow = 1; 
else,                           ax_scaleWindow = 0; 
end

% Voxel Grid % --- %
ax_grid = strcmp(gui.menubar.MRSI.AxisScale.Grid.Checked, 'on');

% X-Limit % --- %
xlimit = num2str(csi.xaxis.xlimit);  

% Set Options % ----------------------------------- %

% Axis scaling
gdat.popup.scaleAxis.String = inp_axs;
gdat.radio.ax_scaleWindow.Value = ax_scaleWindow;
% Voxel grid
gdat.radio.grid.Value = ax_grid;
% Color scaling
gdat.popup.scaleColor.String = inp_clr;
gdat.radio.cl_scaleWindow.Value = cl_scaleWindow;
% X-axis limits.
gdat.edit.xlimit.String = xlimit;

% Save handles.
guidata(hobj,gdat);
  
% -------------------------------------------------------------------- %

  


% ---------------------------------------------------- %
% PLOT 2D PANEL % ---------------------------------------------------- %
% ---------------------------------------------------- %

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

% ---------------------------------------------------- %
% DATA RETRIEVE FUNCTIONS % ------------------=----------------------- %
% ---------------------------------------------------- %

% --- Get CSI domain; frequency or time
function domain = CSI_getDomain(gui)

% Get both domain states
tstate = gui.menubar.MRSI.domain.time.Checked;
fstate = gui.menubar.MRSI.domain.frequency.Checked;

% Set output argument domain
if     strcmp(tstate,'on'), domain = 'time';
elseif strcmp(fstate,'on'), domain = 'freq';
else,                       domain = 'none';
end

% --- Set CSI domain; frequency or time
function CSI_setDomain(hObj, evt, varargin)
% Set the data domain of MRSI data to time or frequency e.g. spectrum or
% FID. hobj refers to CSIgui object.
% Input: 
%       evt.Source.Text (>MATLAB2017) = 'Frequency' or 'Time' 
%       evt.Source.Label(<MATLAB2016) = 'Frequency' or 'Time' 
% If varargin: 'Frequency'/'freq' or 'Time'/'time' variable evt will 
% not be used.

    
% Get gui appdata.
gui = guidata(hObj);
    
if nargin > 2
    evt_click = varargin{1}; 
else
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
end

% Set clicked property to on in menubar
evt_click  = lower(evt_click(1:4));
switch evt_click
    case 'freq'
        gui.menubar.MRSI.domain.frequency.Checked = 'on';
        gui.menubar.MRSI.domain.time.Checked = 'off';
        gui.popup_domain.Value = 2;
    case 'time'
        gui.menubar.MRSI.domain.time.Checked = 'on';
        gui.menubar.MRSI.domain.frequency.Checked = 'off';
        gui.popup_domain.Value = 3;
    otherwise
        gui.menubar.MRSI.domain.time.Checked = 'off';
        gui.menubar.MRSI.domain.frequency.Checked = 'off';
        gui.popup_domain.Value = 1;
end


% Get appdata
if isappdata(gui.CSIgui_main, 'csi')
   csi = getappdata(gui.CSIgui_main, 'csi');
   csi.data.domain = evt_click;
   setappdata(gui.CSIgui_main, 'csi', csi);
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
% Input 2: xaxis can be empty and an arbitrary range will be used otherwise 
%          a ppm-axis for the spec data is expected.
% Input 3: range can be empty and the user will be asked to enter a range
%          otherwise the index range is requested (low to high - not ppm)

% Get peak of interest from user if absent
if nargin <= 2 || isempty(range)
    % Get range from user.
    range = CSI_getPeakOfInterest(xaxis); 
    
    % If no range is found
    if isempty(range)
    CSI_Log({'CSI_getDataAtPeak: No data in given range.'},{'Returning.'});
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
% Input: 
% #1    xaxis struct with unitless index of data and/or corresponding ppm
%       or other unity.
%
% #2    OPTIONAL - if not given; no tag will be displayed in the questions
%                  dialog screen.
%       Tag for questions dialog
%
% Output: unitless index range, ! not the ppm range !
% The output can be used as an index in xaxis.unitless or ppm;
%
% Will be changed so the storage of POI is taken care of automagically
%
% USES CSIpar(ameters) appdata.


% Process input
if nargin == 1, poi_tag = ''; end

% Check for stored POI
CSIguiObj = findobj('Name','CSIgui'); poi_def = NaN;
if ~isempty(CSIguiObj) && ishandle(CSIguiObj)
    if isappdata(CSIguiObj,'CSIpar')
        CSIpar = getappdata(CSIguiObj,'CSIpar');
        poi_def = NaN;
        if isfield(CSIpar,'poi') && isfield(xaxis,'ppm')           
            if CSIpar.poi <= numel(xaxis.ppm), poi_def = CSIpar.poi; end
        end
    end       
end

% Checksum for PPM-field or unitless!
if isfield(xaxis, 'ppm')
    ppm = 1; unit_str = '(ppm):'; 
    if isnan(poi_def)
        unit_ans = [xaxis.ppm(1) xaxis.ppm(end)];
    else
        unit_ans = xaxis.ppm(poi_def);
    end
else
    ppm = 0; unit_str = '(Unitless):'; 
    if isnan(poi_def)
        unit_ans = [1 xaxis.none(end)+1];
    else
        unit_ans = poi_def;
    end
end
% For UserInput dialog NFO.
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

% --- CSIparameters --- %

% Store userinput POI
% 1. GET CSIgui-settings if not available 2. SET poi to CSIpar.POI
if ~exist('CSIpar', 'var')
    CSIguiObj = findobj('Name','CSIgui'); 
    if ~isempty(CSIguiObj) && ishandle(CSIguiObj)
        CSIpar = getappdata(CSIguiObj,'CSIpar'); % Empty or struct;
    else
        CSIpar = NaN; % No CSIguiObj available - set to NaN; dont store.
    end
end
if isempty(CSIpar) || isstruct(CSIpar) % If empty or not NaN;
    CSIpar.poi = poi;
    setappdata(CSIguiObj,'CSIpar',CSIpar);
end

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



% CSIgui-1D: DISPLAY & PANEL % --------------------------------------- %


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

if instance.private
    instance.tag = string(datetime('now','Format','HHmmssSSS'));
else                
    instance.tag = '';
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
if isfield(appdata1D, 'ylimit')
    appdata1D.axes.YLim =  appdata1D.ylimit;
else
    limy = [min(plot_y) max(plot_y)].*1.05;
    if limy(1) >= limy(2), limy(1) = limy(1)-1; end
    ylim(appdata1D.axes, limy);
end

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
w = 180; h = 420;
% Set panel figure size: linked to the 1D plot!
set(fh_1D, 'position', [ceil(scrsz(3).*22/32) ceil(scrsz(4).*2/4) w h]);

% Save panel instance data
setappdata(fh_1D,'instance',instance);

% ADD BUTTONS % ------------------------------------ %

% Add buttons
bw = w-20; bh = 20; % Normalised w and h of button

% Buttons, handles and their info description.
bName = {'Phasing','FFT','iFFT','Apodization', 'Zero Fill', 'SNR',...
         'Linewidth','Baseline','Shift', 'Smoothing','Data Display',...
         'Replace Voxel', 'Export'};
bCall = {@panel_1D_PhasingMethod,...
         @panel_1D_FFT, ...
         @panel_1D_iFFT, ...
         @panel_1D_Apodization, ...
         @panel_1D_ZeroFill, ...
         @panel_1D_SNR,  ...
         @panel_1D_FWHM, ...
         @panel_1D_Baseline, ...
         @panel_1D_Shift,...
         @panel_1D_Smooth,...
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
         'Apply a time/phase shift to the FID or spectrum.',...
         'Denoise data by applying a Sgolay smoothing filter.',...
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


% CSIgui-1D: FUNCTIONS % --------------------------------------------- %

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
%
% This function acts like a phase-corrections manager.


% Manual VS. Automatic % -------------------------------------------- %

% Manual or automatic phasing?
qry = {'Automatic or manual phasing:',...
       'Apply phase single spectrum or parts of MRS-volume:'};
opt = {{'Automatic', 'Manual'}, {'Phase Spectrum', 'Phase Parts'}};
uans = getUserInput_Popup(qry, opt, [], 'Phase Corrections (1D)');
if isempty(uans), CSI_Log({'CSIgui-1D: Aborted phase corrections.'},{''});
    return;
end

% Process userinput
switch uans{1}
    case 'Manual'
        panel_1D_PhaseCorrection_Manual(hObj); 
        do_auto = 0;
    case 'Automatic'
        do_auto = 1;
end

switch uans{2}
    case 'Phase Spectrum', do_parts = 0;
    case 'Phase Parts', do_parts = 1;
end

% Automatic phase-corrections % ------------------------------------ %

% Only continue this part if automatic is chosen
if do_auto    
    % User input
    qry = {'Automatic zeroth- or first- order phase corrections:'};
    opt = {{'Both', 'Zeroth', 'First'}};
    uans = getUserInput_Popup(qry, opt, [], 'Phase Corrections (1D)');
    if isempty(uans), CSI_Log({'CSIgui-1D: Aborted phase corrections.'},{''});
        return;
    end
    
    switch uans{1}
        case 'Zeroth'        
            panel_1D_PhaseCorrection_Zero_Auto(hObj);
        case 'First'
            panel_1D_PhaseCorrection_Auto_First(hObj);
        case 'Both'
            panel_1D_PhaseCorrection_Zero_Auto(hObj);
            panel_1D_PhaseCorrection_Auto_First(hObj);
    end
end

% Phase Parts % ---------------------------------------------------- %

if do_parts
    panel_1D_PhaseCorrection_ApplyToAll(hObj);
end

% --- Executes when user in "Phasing" -> "Automatic" & "First"
function panel_1D_PhaseCorrection_Auto_First(hObj, ~, ~)
% Apply automatic first-order phase corrections

% PREP DATA % ------------------------------------ %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);


% USER INPUT % ----------------------------------- %

% Ask userinput
qry = {'Apply additional zeroth-order corrections afterwards: ',...
       'Accuracy: '};
def = {{'Yes', 'No'}, {'Default','Custom'}}; 
uans = getUserInput_Tick(qry, def, ...
        'First-order phase corrections', [1, 1]);    
if isempty(uans)
    CSI_Log({'Aborted first-order phase corrections.'},{''});
    return; 
end

% Process userinput
if uans(2), acc = 4e-2; 
else
    qry = {'First-order phase correction accuracy:'};
    def = {'3e-2'};
    uans_acc = getUserInput( qry,def, [], 'First-order phase corrections');
    if isempty(uans_acc)
        CSI_Log({'CSIgui-1D: Aborted first-order phase corrections.'},{''});
        return;
    end
    acc = str2double(uans_acc{1}); 
end
do_zero = uans(1);

% APPLY CORRECTION % ----------------------------- %
[data_phased, dphase, cphase] = csi_autoFirstPhase(data, do_zero, acc, 0);

% PLOT & SAVE % ---------------------------------- %

% Add 1D data structure to the 1D figure.
appdata1D.voxel.processed = data_phased;
appdata1D.phasing.zero    = cphase;
appdata1D.phasing.first   = dphase;

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D:'},{'Applied automatic first-order phase correction.'});
if do_zero, uans = 'Yes'; else, uans = 'No'; end
CSI_Log({'CSIgui-1D: Additional zeroth-order phasing: '},{uans});
CSI_Log({'CSIgui-1D: Calculation accuracy:'},{acc});

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

% zeroloc = find(first_ph == 0);
% linone = first_ph(1):abs(first_ph(1))/(zeroloc):0;
% lintwo = 0:abs(first_ph(1))/(zeroloc):first_ph(end);


CSI_Log({'Phase change 0th order: '},{rad2deg(zero_ph)});
CSI_Log({'Phase change 1th order: '},{[first_ph(1) first_ph(end)]});

% CLEAN UP % ----------------------------------------- %

% Save
setappdata(CSI_1D_obj,'data1D',appdata1D);

% Plot the data
CSI_1D_displayData(CSI_1D_obj);

% Update info to user.
CSI_Log({'CSIgui-1D:'},{'Applied manual phase correction.'});

% --- Executes if user chooses in "Phasing" -> "Phase All"
function panel_1D_PhaseCorrection_ApplyToAll(hObj,~)
% See if there is any phase correction data. Else call the manual 
% phasing function first and then apply it to all voxel in the data 
% set or a slice.

% INITIATE % ------------------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D figure object, data_1D appdata and data array
[CSI_1D_obj, appdata1D] = CSI_1D_getData(obj1D);

% 2D/3D appdata
CSImain_obj = findobj('Tag','CSIgui_main');
if isempty(CSImain_obj)
    fprintf('Error: CSIgui appears to be closed. Returning.\n');
    return; % CSIgui is closed!
end
csi = getappdata(CSImain_obj, 'csi');

% Get user input for volume or slice phasing
uans = getUserInput_Popup({'Apply phasing to: '},...
                         {{'Full volume', 'Current slice'}},...
                         [], 'Phase Corrections (1D)');
if uans{1} == 0, return; end
if strcmpi(uans{1}, 'Full volume'), vol_or_sli = 0;
else, vol_or_sli = 1;
end
% vol_or_sli = 0 for volume and vol_or_sli = 1 for slice.

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

% Backup of CSI data.
gui = guidata(CSImain_obj); % This is CSIgui-gui data
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui,'Before 1D-phase-corrections to data.'); end
csi = getappdata(CSImain_obj, 'csi');

if vol_or_sli == 0

    % Add zero order phase correction
    pha_new = pha + appdata1D.phasing.zero;
    % Add first order phase correction
    pha_new = pha_new + appdata1D.phasing.first;
    
    log_msg = ...
        'Applied the phase correction to all voxel in the MRSI dataset.';
elseif vol_or_sli == 1    
    
    % Get current visualized slice
    if sum(appdata1D.voxel.index) == 0, return; end 
    ind = num2cell([size(pha,1) appdata1D.voxel.index]); % C/R/D4/D5 etc
    
    % fid, x/c, y/r, slice and the rest
    dims = num2cell(size(pha));    
    dims_range =cellfun(@(x) 1:x, dims,'uniform',0);
    for kk = 4:size(dims,2)
        dims_range{kk} = ind{kk};
    end
    
    % Add zero order phase correction
    pha_new_cut = pha(dims_range{:}) + appdata1D.phasing.zero;
    % Add first order phase correction
    pha_new_cut = pha_new_cut + appdata1D.phasing.first;
    
    pha_new = pha;
    pha_new(dims_range{:}) = pha_new_cut;

    log_msg = ...
        'Applied phase correction to current slice.';
end

% Create complex data
csi.data.raw = complex(magn.*cos(pha_new), magn.*sin(pha_new));

% Save the appdata
setappdata(CSImain_obj, 'csi', csi);

% Replot the 2D plot!
CSI_2D_initiate2D();

% Update info to user.
CSI_Log({'CSIgui-1D:'},{log_msg});

% --- Executes if user chooses in "Phasing" -> "Automatic" & "Zeroth"
function panel_1D_PhaseCorrection_Zero_Auto(hObj, ~)
% Apply auto zeroth order phase correction.

% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);

% USER INPUT % ---------------------------------- %

% Get peak of interest.
range = CSI_getPeakOfInterest(appdata1D.axis, ...
                              'CSIgui-1D: Automatic phase correction');

qry = {'Phase corrections method:'};
def = {'Maximize the real-component.','Match real-compontent to absolute.'}; 
uans = getUserInput_Popup(qry,{def}, [], 'Phase Corrections (1D)');
if isempty(uans), return; end
switch uans{1}
    case 'Maximize the real-component.', ph_meth = 1;
    case 'Match real-compontent to absolute.', ph_meth = 2;
end

% APPLY CORRECTION % ----------------------------- %

% POI from user.
if length(range) > 1, poi = range(1):range(2); end

% Apply autophasing function
[data_phased, phase_change] = csi_autoZeroPhase(data, poi, ph_meth, 0);

% PLOT & SAVE % ---------------------------------- %

% Add 1D data structure to the 1D figure.
appdata1D.voxel.processed = data_phased;
appdata1D.phasing.zero = phase_change;
appdata1D.phasing.first = zeros(1,size(data_phased,1));

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
uans = getUserInput_Popup({'Append direction:'},{{'Post','Pre','Both'}},...
                          [], 'Zerofill');
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

function [peak_names, peak_pos, peak_on] = csi_loadPeakInfo()
% Loads text-file peaks.txt to get peak-info. Returns the peak names, peak
% position and on/off for queries.
%
% Expected line-structure in peaks.txt:
% peak-name peak-lowerbound peak-upperbound peak-on/off
% First line of text-file will be ignored.
%
% peak on/off is used for tick-box user input, allowing enabling and
% disabling peaks of interest whilst most peaks of a species are denoted in
% the text-file.

% Get file-identifier
fid = fopen('Scripts\Peaks.txt','r'); 


% Get #lines in file
n = 0; while ~feof(fid), fgetl(fid); n=n+1; end
frewind(fid); % Return to start of file

% Create variables and read line-data
peak_names = cell(1,n-1); peak_pos = NaN(n-1, 2); peak_on = NaN(1,n-1);
n = 0;
while ~feof(fid)
    n = n+1;
    line = fgetl(fid);
    if n > 1
        line = strsplit(line);
        peak_names{n-1} = line{1};
        peak_pos(n-1,1:2) = str2double([line(2) line(3)]);
        peak_on(n-1) = str2double(line(4));
    end
end
fclose(fid); % Close file

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
method = getUserInput_Popup({'Signal unit: '},{{'Real','Magnitude'}},...
                            [], 'SNR (1D)');
if isempty(method), return; end
switch method{1}, case 'Real', method = 1; case 'Magnitude', method = 0; end


% Mask size
masksz = getUserInput({'Mask size: '},{round(size(data,1)/12)});
if isempty(masksz), return; end
masksz = str2double(masksz);

% ----- % Calculate SNR  

% Offset in y-axis used for plotting at peak positons
yst = diff(appdata1D.axes.YLim)/50; 

% Automatic or for single peak
if strcmp(uans,'Automatic') % All peaks
    
    % Automatic with prior knowledge list or via local-maximum peak search
    qry = {'Which automatic method? (See peaks.txt):'};
    def = {'Use peaks info text-file', 'Use local maximum'};
    uans = getUserInput_Popup(qry,{def});
    if isempty(uans), return; end

    % Switch method for automatic peak localization
    switch uans{1}
        % -------------------------------------------------------------- %
        case def{1} % --- Use Peaks.txt nfo ---------------------------- %

            % Get peak-names, bounds and peak-on/off for userinput
            [peak_name, peak_boundary, peak_on] = csi_loadPeakInfo();
        
            % Get userinput which peaks to include.
            space_string = repmat({repmat(char(32),1,24)}, size(peak_name));
            qry = cellfun(@sprintf, repmat({'%s%s'}, size(peak_name)), ...
                space_string, peak_name, 'UniformOutput', false);
            opt = {'Yes', 'No'}; opt = repmat({opt}, 1, size(qry,2));
            def = peak_on;
            uans = getUserInput_Tick(qry, opt, 'Peaks to include:', def);
            if isempty(uans), return; end
        
            % Remove unselected peaks.
            peak_boundary(uans == 0, :) = [];
            peak_name(uans == 0) = [];
                        
            if isfield(appdata1D.xaxis,'ppm')
                ppm = appdata1D.axis.ppm; 
            else
                % Converts ppm-boundaries wrt center spectrum length.
                % This will not be perfect!
                
                % Peak furthest from center
                max_peak = max(abs(peak_boundary(:,2)));
                
                % use that to create a ppm-axis
                half_samples = round(size(data,1)./2);
                axis_max = ((half_samples*1.05)./max_peak);
                ppm = linspace(-axis_max,axis_max,size(data,1))';
            end
    
            % Get peak-maximum within boundary using ppm-axis
            boundary_unit = nan(size(peak_boundary));
            peak_pos_unit = NaN(1,size(peak_boundary,1));
            for kk = 1:size(peak_boundary,1) % ppm 2 unitless
                % Convert peak-boundary in ppm to nearest value in ppm-axis
                [~, boundary_unit(kk,1)] = ...
                    min(abs(ppm - peak_boundary(kk,1)));
                [~, boundary_unit(kk,2)] = ...
                    min(abs(ppm - peak_boundary(kk,2)));
        
                % Peak-boundaries unitless
                boundary_unit(kk,:) = sort(boundary_unit(kk,:));
                
                % Peak of interest - maximum within boundaries
                [~, nn] = max( ...
                    real( data( boundary_unit(kk,1):boundary_unit(kk,2))));
                peak_pos_unit(kk) = boundary_unit(kk,1) + nn-1;
            end
        
            % Peak position for plotting, i.e. max value within peak boundaries
            if isfield(appdata1D.axis,'ppm')
                peak_pos_plot = ppm(peak_pos_unit); 
            else                       
                peak_pos_plot = peak_pos_unit;
            end
    
        % -------------------------------------------------------------- %
        case def{2} % --- Use local maximum method --------------------- %
    
            % ------ % Find peaks
            def = 1:15; def(1) = 7; def(2:7) = 1:6; def(8:end) = 8:15;
            uans = getUserInput_Popup({'Number of peaks:'},{{def}});
            if isempty(uans), return; end
        
            nPeaks = str2double(uans{1});
        
            % Positions of peaks
            peak_pos_unit = csi_findPeaks(data, nPeaks)';
        
            % Convert peak position to ppm
            if isfield(appdata1D.axis,'ppm')
                peak_pos_plot = appdata1D.axis.ppm(peak_pos_unit); 
            else                       
                peak_pos_plot = peak_pos_unit;
            end
            
            peak_name = repmat({''}, size(peak_pos_unit));
    end

    % Plot peak locations
    hold(appdata1D.axes,'on'); 
    plot(appdata1D.axes, peak_pos_plot, real(data(peak_pos_unit)), 'or');

    % ------ % Calculate SNR
    
    % Calculate SNR for found peaks
    np = size(peak_pos_unit,2); snr = NaN(1,np);
    
    % Calculate SNR per peak
    for kk = 1:np
        range = peak_pos_unit(kk) + [-5 5]; % Range of the peak
        snr(kk) = csi_SNR(data, masksz, method, range);
        
        % Plot SNR as text at peak-max location
        txt_plot_y = real(data(peak_pos_unit(kk)))+yst;
        txt_string = sprintf('%3.1f',snr(kk));
        text(appdata1D.axes, peak_pos_plot(kk), txt_plot_y, txt_string,...
            'FontSize',8,'FontWeight','Bold');
    end


else % Individual peak
    
    % ----- % Calculate SNR
    
    % Get range from user (easy to use the dataAtPeak function)
    range = CSI_getPeakOfInterest(appdata1D.axis,'CSIgui-1D: SNR');
    if isempty(range), return; end
    
    % Calculate SNR
    snr = csi_SNR(data, masksz, method, range); % Output SNR
    
    % ----- % Display SNR 
    
    % At max value in range
    tmp = range(1):range(2); [val, ind] = max(real(data(tmp,:))); 
    ind = tmp(ind);
    
    % Convert position to ppm values
    if isfield(appdata1D.axis,'ppm'), ind = appdata1D.axis.ppm(ind); end
    
    % Plot marker at peak location
    hold(appdata1D.axes,'on'); plot(appdata1D.axes, ind, real(val), 'or');
    
    % Plot SNR as text at peak location    
    text(appdata1D.axes, ind, val+yst, ...
       sprintf('%3.1f',snr),'FontSize',8,'FontWeight','Bold');
    peak_name = {'N.A.'};
end
hold(appdata1D.axes,'off'); % Turn off hold

% Show SNR to user.
CSI_Log({'CSIgui-1D: Voxel -','Name:', 'SNR:'},...
    {appdata1D.voxel.index, strjoin(peak_name,'   '), snr});

% --- Executes when user presses "Linewidth" in panel_1D
function panel_1D_FWHM(hObj, ~)
% Calculate linewidth at FWHM at a certain peak of for all peaks found
% automatically.

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, appdata1D, data] = CSI_1D_getData(obj1D);

% Automatic or manual calculation?
qry = {'FWHM estimation method:',...
       'Calculate FWHM for one individual peak or all automagically:'};
def = {{'Intersection', 'Local minima'},...
       {'Individual', 'Automatic'}};
uans = getUserInput_Popup(qry, def, [], 'FWHM (1D)');
if isempty(uans), return; end

if strcmp(uans{1}, 'Intersection')
    panel_1D_FWHM_method_intersect(hObj);
    return;
else
    uans = uans(2);
end

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
    range = CSI_getPeakOfInterest(appdata1D.axis,'CSIgui-1D: FWHM');
    if isempty(range), return; end
    
    % Maximum in this range
    [~,mi] = max(real(data(range(1):range(2))));
    % Set max as peak centre.
    peak_pos = mi+range(1);
    
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
        csi_lineWidth(data, ax, range(1):range(2), 1);
    
    % Plot data
    hold(appdata1D.axes,'on'); % Hold 1D plot figure
    plot(appdata1D.axes, lwp(pp,:), real(lwv(pp,:)),'or');    % Plot marker
    text(appdata1D.axes, lwp(pp,1)-xst, lwv(pp,1)+yst, ...    % Plot text
       sprintf('%2.3f',lw(pp)),'FontSize',8,'FontWeight','Bold');
    
end
hold(appdata1D.axes,'off');

% Update in UI
CSI_Log({'Lowest position to highest, FWHM: '},{lw'});

% --- Executes by choice of FWHM-method
function panel_1D_FWHM_method_intersect(hObj,~,~)

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end
% Do checks and get CSI_1D object, data_1D appdata and data array
[~, appdata1D, data] = CSI_1D_getData(obj1D);


% Get range from user (Unitless range given)
range = CSI_getPeakOfInterest(appdata1D.axis,'FWHM');
if isempty(range), return; end

% Prep xaxis % -------- %

% Use unitless or ppm
if isfield(appdata1D.axis, 'ppm')
    ax = appdata1D.axis.ppm;
else                         
    ax = appdata1D.axis.none;
end


% Baseline | Samples
N = size(data,1); 
bv = mean(real([data(1:round(N/20)); data(N-round(N/20):N)])); 
bvline = cat(2,ax',repmat(bv,N,1));

% MAXIMUM
data_ranged = real(data(range(1):range(2)));
[mv, mi] = max(data_ranged);
% Set max as peak centre.
peak_pos = mi+range(1)-1; % Data index
xv = ax(peak_pos);        % In unit (ppm/unitless etc.)
ax_ranged = ax(range(1):range(2));
% FWHM ESTIMATE
fwhm = mv - ((mv-bv)/2);
fwhmline = cat(2,ax',repmat(fwhm,N,1));

% PREPARE FIGURE % -------- %

% Create new figure for user input
fh = figure(); 
% Plot the spectrum (reverse xaxis)
spObj = plot(ax, real(data)); spObj.Parent.XDir = 'reverse'; hold on; 
% Axes
axObj = spObj.Parent;

% PLOT MAXIMUM
plot(axObj, xv,  mv,'or');
    % datapoint
cursorMode = datacursormode(fh);hdtip = cursorMode.createDatatip(spObj);
hdtip.Position = [xv mv 0]; updateDataCursors(cursorMode);


% PLOT BASELINE
baseObj = plot(axObj, bvline(:,1),bvline(:,2),'--k');
    % datapoint
cursorMode = datacursormode(fh); hdtip = cursorMode.createDatatip(baseObj);
hdtip.Position = [bvline(1,1) bv 0]; updateDataCursors(cursorMode);


% PLOT FWHM ESTIMATE
fwhmObj = plot(axObj, fwhmline(:,1),fwhmline(:,2),'--c');
    % datapoint
cursorMode = datacursormode(fh); hdtip = cursorMode.createDatatip(fwhmObj);
hdtip.Position = [fwhmline(1,1) fwhm 0]; updateDataCursors(cursorMode);


% Intersect FWHM line and DATA(range) (Estimate)
[xest,yest] = polyxpoly(ax_ranged',data_ranged, fwhmline(:,1)', fwhmline(:,2)');

% Plot results of found intersects
plot(axObj, xest, yest, '-om');

% The answer in PPM or unit-less: ALWAYS use the first two - if more
% intersects are found it is shown in the plot. 
if size(xest,1) == 1 % The peak is not within the range probabaly
    fwhm_outp = NaN; 
else
    xest = [xest(1) xest(end)]; yest = [yest(1) yest(end)];
    fwhm_outp = diff(xest); 
end


% Plot data in CSIgui 1D
yst = diff(appdata1D.axes.YLim)/50; xst = diff(appdata1D.axes.XLim)/100;
hold(appdata1D.axes,'on'); % Hold 1D plot figure
plot(appdata1D.axes, xest, yest,'or');    % Plot marker
text(appdata1D.axes, xest(1)-xst, yest(1)+yst, ...    % Plot text
   sprintf('%2.3f',fwhm_outp),'FontSize',8,'FontWeight','Bold');

% Update in UI
if isnan(fwhm_outp)
    fwhm_outp = 'Single intersect with FWHM found - redefine peak range';
end
CSI_Log({sprintf('FWHM @ [%2.2f ; %2.2f]: ',ax_ranged(1), ax_ranged(2))},...
        {fwhm_outp});

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


% USER DATA % ---------------------------------- % 

% Get user input
qry =  {'Correct for N term?',...
            'Correct for onesided fft in spectroscopy?',...
            'Shift before and after FFT? (Echo)'};
def = {{'No','Yes'}, {'Yes','No'}, {'No','Yes'}};
uans = getUserInput_Popup(qry, def, [], 'FFT (1D)');
if isempty(uans), CSI_Log({'Skipped FFT.'},{''}); return; end

if strcmpi(uans{1},'yes'), correct_N = 1; else, correct_N = 0; end
if strcmpi(uans{2},'yes'), onesided = 1;  else, onesided = 0; end
if strcmpi(uans{3},'yes'), dbl_shift = 1; else, dbl_shift = 0; end


% FFT % ---------------------------------------- %

% Create a N x 1 vector!
if size(data2fft,2) > size(data2fft,1), data2fft = data2fft'; end

% Apply 1D fourier transform
appdata1D.voxel.processed = csi_fft(data2fft,correct_N,onesided,dbl_shift);

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

% USER DATA % -------------------------- %

uans = getUserInput_Popup(...
           {'Correct for onesided fft in spectroscopy?'},...
           {{'Yes','No'}}, [], 'iFFT (1D)');
if isempty(uans), CSI_Log({'Skipped FFT.'},{''}); return; end
if strcmpi(uans{1},'yes'), onesided = 1;  else, onesided = 0; end

% iFFT % --------------------------------------- %

% Create a N x 1 vector.
if size(data2ifft,2) > size(data2ifft,1), data2ifft = data2ifft'; end

% Apply 1D fourier transform
appdata1D.voxel.processed = csi_ifft(data2ifft,onesided);

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
qry = {'Filter Method', 'Flip window: '};
def = {{'Gaussian', 'Hamming','Hann', 'Exponential', 'Blackman', ...
        'Flattop'}, {'No', 'Yes'}};
uanstype = getUserInput_Popup(qry, def, [], 'Apodization (1D)');
if isempty(uanstype), return; end

% Additional options for specific filters
switch uanstype{1}
    case 'Gaussian'
        % ((1/(1/BW x nSamples) * (ans/nSamples)) / pi) == (bw/ans)/pi 
        % Hz = bw / (N*pi) && N = (Hz * Pi)/BW
        if isfield(appdata1D.axis,'BW')
            uopts = getUserInput({'Apodization factor: (Hz)'},{20});    
            if isempty(uopts)
                CSI_Log({'Skipped apodization.'},{''}); return; 
            end
            hz = str2double(uopts{1});
            uansopts{1} = num2str(appdata1D.axis.BW ./ (hz * pi));
        else
            uansopts = getUserInput(...
            {'Standard deviation e.g. half length at FWHM: (Samples)'},...
            {(appdata1D.xaxis.N/4)});
        end
    case 'Exponential'
        uansopts = getUserInput(...
            {'Exponential decay strength (Samples): '},{0.5});
end

if strcmp(uanstype{2}, 'Yes'), opts(3) = 1; end

% Set options if applicable
if exist('uansopts', 'var'), opts(1) = str2double(strsplit(uansopts{1},' '));
else, uansopts = {''}; opts(1) = 0;
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
[CSI_1D_obj, appdata1D, data] = CSI_1D_getData(obj1D);


% USER INPUT % -------------------------------- %

% X-axis ------------ %
uans = getUserInput_Popup({'X-axis type: ','Plot unit: '},...
                          {{'PPM', 'Frequency', 'Time','None'},...
                           {'Real', 'Magnitude', 'Phase', 'Imaginary'}},...
                           [], 'Data Display (1D)');
if isempty(uans), return; end

% User input handling
axis_unit = lower(uans{1}); plot_unit = uans{2};

% Edit some x axis labels...
if strcmp(axis_unit,'frequency'), axis_unit = 'freq'; end 

% Get current x-axis settings
xlimit = appdata1D.axes.XLim;                     

% Ask user x-axis limits
xans = getUserInput({'X-axis limits: (ppm or unitless)'},{xlimit}); 
if isempty(xans),return;end
appdata1D.axis.xlimit= str2double(strsplit(xans{1}));


% Y-axis ------------ %
data_new_unit = CSI_getUnit(data, plot_unit);
uans = getUserInput({'Y-axis limits:'},...
    {[min(data_new_unit) max(data_new_unit)]});

% APPLY SETTINGS % ---------------------------- %
appdata1D.voxel.unit = plot_unit;
appdata1D.axis.unit  = axis_unit;
appdata1D.ylimit     = str2double(strsplit(uans{1}));

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

% --- Executes when user presses "Smooth" in panel_1D.
function panel_1D_Smooth(hObj, ~)
% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end

% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, dataOri] = CSI_1D_getData(obj1D);

% Shift % ---------------------------------------- %
data  = dataOri;

datasmo = smoothdata(real(data),'sgolay', 'degree', 2);

appdata1D.voxel.processed = datasmo;
disp ERROR!!!!!QH20202
% Save the complex result to the 1D-plot figure's appdata.
setappdata(CSI_1D_obj, 'data1D', appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj)

% Show info to user
CSI_Log({'CSIgui-1D:'},{'Denoised data by smoothing.'});

% --- Executes when user presses "Shift" in panel_1D.
function panel_1D_Shift(hObj, ~)
% Shift the FID to its first sin-amplitude or first maximum 
%                         OR...
% Use the TE to shift.
%
% BETA - Experimental function

CSI_Log({'WARNING:'},...
    {'Automatic "CSIgui-1D Shift" is an experimental function'});


% PREP DATA % ---------------------------------- %

% Get CSIgui 1D figure object
obj1D  = panel_1D_getInstanceData(hObj); if ~isobject(obj1D),return; end

% Do checks and get CSI_1D object, data_1D appdata and data array
[CSI_1D_obj, appdata1D, dataOri] = CSI_1D_getData(obj1D);

% Shift % ---------------------------------------- %
data  = dataOri;

qry = {'Shift method: '};
def = {{'Automatic','Use echotime (TE)'}};
uans = getUserInput_Popup(qry, def, [] , 'Phase Shift (1D)');
if isempty(uans), CSI_Log({'CSigui-1D: Skipped 1D Shift'},{''}); return; 
end  

switch uans{1}
    case 'Automatic'
        % Cheap approach - move data to position of maximum FID signal
        sz = size(data,1); 
        if isfield(appdata1D.xaxis, 'time')
            ax = appdata1D.xaxis.time;
        else
            ax = appdata1D.xaxis.none;
        end

        % Smooth data
        winsz = round(size(data,1)./60);
        datasmo = smoothdata(real(data),'gaussian', winsz);

        % Plot original data.
        figure(); plot(ax, datasmo, '-b' ); 
        hold on; plot(ax,real(data),'--y'); 

        % Find peaks A
        [pks, locs, w] = findpeaks(datasmo,'NPeaks', 5); 
        pks_smooth = cat(2,locs,pks);

        % Maximum width-peak
        [~, max_w_loc] = max(w);
        shiftval = locs(max_w_loc);
        mm = pks(max_w_loc);
    
        % Save new data.
        appdata1D.voxel.processed = zeros(sz,1);
        appdata1D.voxel.processed(1:(sz-shiftval+1)) = dataOri(shiftval:end);

        % Plot all nicely.
        figure(); 
        plot(ax,real(data),'-b','linewidth',1.8);  hold on; 
        plot(ax,datasmo,'-y','linewidth',1.8); 
        plot(ax(pks_smooth(:,1)),pks_smooth(:,2),'og','linewidth',1.8); 
        plot(ax([1 shiftval]),[mm mm], '-om','linewidth',2);
        pl = plot(ax,real(appdata1D.voxel.processed),'-r');
        legend({ 'Original', 'Smooth', 'Used maxima', 'Shift-line', 'New data'});


        dt = mean(diff(ax)); str = sprintf('Shifted: %1.6f | %3.0f', dt.*shiftval,...
            shiftval);
        text(ax(round(shiftval./2)), mm.*1.05, str, 'Color', [0.8 0.8 0.8]);

        pl.Parent.Color = [0 0 0]; 


    case 'Use echotime (TE)'
        uans = getUserInput({'Echotime, TE (ms):'},{0.35});
        if isempty(uans), return; end
        
        TE = str2double(uans{1});
        shiftval = round(TE./(mean(diff(appdata1D.xaxis.time)).*1000));
        sz = size(data,1);
        
        % Output
        CSI_Log({'CSIgui-1D: FID shifted by '},{shiftval});

        % Save new data.
        appdata1D.voxel.processed = zeros(sz,1);
        appdata1D.voxel.processed(1:(sz-shiftval+1)) = dataOri(shiftval:end);
        
end


% PLOT AND SAVE % ------------------------------ %

% Save the complex result to the 1D-plot figure's appdata.
setappdata(CSI_1D_obj, 'data1D', appdata1D);

% Plot the 1D data
CSI_1D_displayData(CSI_1D_obj)

% Show info to user
CSI_Log({'CSIgui-1D:'},{'Applied shift.'});



% -------------------------------------------------------------------- %

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
if nargin >= 1
    infolabel = varargin{1}; 
    if nargin < 2, infoval = ''; else, infoval  = varargin{2}; end

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


% -------------------------------------------------------------------- %
% Noise Data Manager and Processing % -------------------------------- %
% -------------------------------------------------------------------- %

% --- View noise component executed in menu

% --- Noise data manager for separate noise-data
function CSI_Noise_ViewManager(hObject, ~, ~)
% Gui-data struct.
gui = guidata(hObject);


% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

if isfield(csi.data, 'noise')
    
    if isfield(csi.data, 'backup_toview_noise')         % REVERT TO DATA
        
        % Store noise
        csi.data.noise.raw = csi.data.raw;
        csi.data.noise.labels = csi.data.labels;

        % Return to original data
        csi.data.raw = csi.data.backup_toview_noise.raw;
        csi.data.labels = csi.data.backup_toview_noise.labels;

        % Dimensions update
        csi.data.dim = size(csi.data.raw);

        % Remove the backup field.
        csi.data = rmfield(csi.data,'backup_toview_noise');
        
        % Message
        msg = 'Stored noise data and reverted to backup of original data.';
        % Set correct menu-label
        gui.menubar.MRSI.noise.Label = 'Load Noise...';
    else                                                % SET NOISE ACTIVE
        % Create backup raw and labels
        csi.data.backup_toview_noise.raw = csi.data.raw;
        csi.data.backup_toview_noise.labels = csi.data.labels;

        % Replace raw with noise.
        csi.data.raw = csi.data.noise.raw;
        % Dimensions update
        csi.data.dim = size(csi.data.raw);
        
        % Labels
        if isfield(csi.data.noise,'labels')
            csi.data.labels = csi.data.noise.labels;
        else
            new_label = cell(1,numel(size(csi.data.raw)));
            new_label{1} = 'fid';
            empty_lab = cellfun(@isempty, new_label) == 1;

            % String characters to double
            labelAsNum = cellfun(@double, new_label, 'UniformOutput', false);
            % Number of characters
            nChar = cellfun(@size, new_label, repmat({2},1,size(new_label,2)));
            % Index of single charactes
            strt = max(cell2mat(labelAsNum(nChar == 1)))+1;
            % If no single character label found
            if isempty(strt), strt = 65; end
            
            % Create new label names for empty indexes.
            ind = find(empty_lab == 1);
            for kk = 1:size(ind,2)
                new_label(ind(kk)) = {char(strt+kk-1)}; 
            end
            
            csi.data.labels = new_label;
        end

        % Info for user
        msg = 'Noise data inserted. Backup of original created.';
        % Set correct menu-label
        gui.menubar.MRSI.noise.Label = 'Revert to Data...';
    end
    
    % Save appdata.
    setappdata(gui.CSIgui_main, 'csi', csi);
    
    CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

    %Show user.
    CSI_Log({'% ------------------------------------------------ %'},{''}); 
    CSI_Log({msg},{''}); 
    CSI_Log({'% ------------------------------------------------ %'},{''});
else
    CSI_Log({'No noise data present.'},{''}); return;
end
   
% --- Prepare noise-data for use in calculations
function [csi, hobj, gui] = CSI_Noise_Prepare(hobj, gui)
% Process noise from noise-measurement. This is stored in csi.data.noise
% and loaded separately or extracted from a list-data-file itself.
%
% This method uses the CSIgui built-in switch noise/data view option to
% process the data automatically.

% Apodization | FFT | Delete Channel
qry = {'Process noise data:', 'Delete a channel from noise-data:', ...
       'Average noise', 'Remove OS', 'Apodize noise:', 'FFT noise:', };
opt = repmat({{'Yes', 'No'}}, size(qry));
dans = ones(1,size(qry,2)); dans(3) = 0; dans(4) = 0;
uans = getUserInput_Tick(qry, opt, 'Noise Preparation', dans);
if isempty(uans)
    CSI_Log({'Noise Preperation:'},{'Aborted.'}); csi = NaN; return; 
end

% UANS 1: process noise data? or no options enabled
if uans(1) == 0 || sum(uans(2:end)) == 0
    csi = getappdata(gui.CSIgui_main,'csi'); return; 
end

% Process noise by setting it as active data set    
CSI_Noise_ViewManager(hobj, [], []);

% Noise - Delete channel
if uans(2)
    button_CSI_Delete_Callback([],[],gui, 0);
    gui = guidata(hobj);
end

% Noise - Average
if uans(3)
    % Apodize
    button_CSI_Average_Callback([],[], gui, 0);
    gui = guidata(hobj);
end

% Noise - Remove oversampling
if uans(4)
    button_CSI_RemoveOS_Callback([],[],gui,0);
    gui = guidata(hobj);
end

% Noise - Apodization
if uans(5)
    % Apodize
    button_CSI_Apodization_FID_Callback([],[],gui,0);
    gui = guidata(hobj);
end

% Noise - FFT
if uans(6)
    button_CSI_FFT_Callback([], [], gui, 0);
    gui = guidata(hobj);
end
            
% Revert to data
CSI_Noise_ViewManager(hobj, [], []);

% Reload csi-app-data
csi = getappdata(gui.CSIgui_main,'csi');


% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %


% --- Executed via menubar to save data to file
function CSI_saveData(hObject, ~, ~)
% Gui-data struct.
gui = guidata(hObject);

% End if no CSI data present
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Find default path to open
if isfield(csi, 'filepath'), fp = csi.filepath;
else, fp = [];
end
           
% UI for file from user
[fn, fp] = uiputfile({'*.mat','MATLAB file';'*.sdat','Sdat file';...
    '*.txt','Text file';},'Save MRSI data to file...',fp); 
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
        if isfield(csi,'filepath'), fp = csi.filepath;
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
    
    % Reduce file size? Save Backup?
    uans = getUserInput_Popup(...
        {'Reduce file size by converting to single?',...
         'Include backup-data?'},...
        {{'Yes','No'},{'No', 'Yes'}}, [], 'Save Data');
    
    if strcmpi(uans{1},'yes')
        csi.data.raw = single(csi.data.raw);
    end
    if strcmpi(uans{2},'No') && isfield(csi,'backup')
        csi = rmfield(csi,'backup');
    end
    
    % Add log
    log = get(gui.listbox_CSIinfo,'String');
    csi.log = char(log);


    % % Get MRSI data
    % csigui = csi.data;
    % if strcmpi(reduce_size{1},'yes')
    %     csigui.raw = single(csigui.raw);
    % end
    % 
    % if isfield(csi,'twix')
    %     csigui.twix = csi.twix;
    % end
    % if isfield(csi,'list')
    %     csigui.list = csi.list;
    % end
    % if isfield(csi,'voxelmask')
    %     csigui.voxelmask = csi.voxelmask;
    % end
    % 
    % % Add log
    % log = get(gui.listbox_CSIinfo,'String');
    % csigui.log = char(log);
    % 
    % % Add xaxis structure
    % csigui.xaxis = csi.xaxis;
    % 
    % Add Ori
    % if isfield(csi,'ori'), csigui.ori = csi.ori; end
end

                  % ------- % Save CONV data % ------- %

% Check for converted data
if isappdata(gui.CSIgui_main,'conv')
    % csigui.conv = getappdata(gui.CSIgui_main,'conv');
    csi.conv = getappdata(gui.CSIgui_main,'conv');
end

                  % ------- % Save MRI data % ------- %

% Check for mri data
if isappdata(gui.CSIgui_main,'mri')    
    % csigui.mri = getappdata(gui.CSIgui_main,'mri');
    csi.mri = getappdata(gui.CSIgui_main,'mri');
end

                   % ------- % Save file  % ------- %

% Save it.
save([fp fn ext], 'csi','-v7.3');

% Done.
CSI_Log({'Saved MRSI data to MAT-file:'},{[fp fn ext]});


% --- Save CSIgui 2D plot to file       
function CSI_saveFig(hObject, ~, ~)
% Save the 2D CSIgui plot to file supporting multiple formats;
% jpg, png, eps, tiff, fig.



% Get guidata
gui = guidata(hObject);

% Check for CSI data
if ~isappdata(gui.CSIgui_main,'csi'), return; end
csi = getappdata(gui.CSIgui_main,'csi'); 

if isfield(csi, 'filepath'), fp = csi.filepath; else, fp = cd; end
% Get filepath destination
% Get file path and extension from user
[fn, fp, fi] = uiputfile({'*.png','Portable network graphic';...
                          '*.jpg', 'JPEG image';'*.eps','EPS file';...
                          '*.tiff','TIFF image';'*.bmp','Bitmap file';...
                          '*.fig','MATLAB Figure'},...
                          'Save CSIgui 2D figure...', fp);
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

uans = getUserInput_Popup({'Specify figure(s) to save:',...
                           'Resolution (DPI):',...
                           'Transparency:',...
                           'Show index in figure:',...
                           'Fast mode: '},...
                         {{'Current figure', 'All', 'All of current slice'},...
                           cat(2,num2cell(0:200:1200),'Custom'),...
                           {'No','Yes'},{'Yes','No'},{'No','Yes'}}, ...
                         [], 'Save Figure');
if isempty(uans), CSI_Log({'Skipped exporting figures.'},{''}); return; 
end  

% Process what to save
switch uans{1}
    case 'Current figure', save_all = 0;
    case 'All',            save_all = 1;
    case 'All of current slice', save_all =  2;
end

% Process resolution
if strcmp(uans{2},'Custom')
        uans{2} = getUserInput({'Custom DPI scaling:'},{''});
end
dpi = str2double(uans{2});

% Process transparency
if strcmp(uans{3},'Yes'), transparency = 1; else, transparency = 0; end

% Process plotting index in figure
if strcmp(uans{4},'Yes'), showIndex = 1; else, showIndex = 0; end

% Process fast saving
if strcmp(uans{5}, 'Yes'), fast = 1; else, fast = 0; end

% FIGURE OBJECT % ------------------------------------------------------ %
% Get 2D-plot figure its object, apply chosen options and save as image.

% Check if 2D Plot is active (open).
fig_obj =  findobj('type','figure','tag','CSIgui_plot2D');
% If 2D Plot not active, open it.
if isempty(fig_obj)
    panel_2D_DataSliders([],[],gui); CSI_2D_initiate2D(); 
end

% Get 2D-Plot and 2D-Panel figure
fig_obj = findobj('type','figure','tag','CSIgui_plot2D');
pan_obj = findobj('type','figure','tag','CSIpanel_2D_DataToDisplay');
if ~isempty(pan_obj)
    pan_gui = guidata(pan_obj);
else
    % Open slider panel
    panel_2D_DataSliders([],[],gui);
    % Get the object
    pan_obj = findobj('type','figure','tag','CSIpanel_2D_DataToDisplay');
    pan_gui = guidata(pan_obj);
    % Get plot index for 2D plot figure gui data.
    gui2D = guidata(fig_obj); gui2D.plotindex
    
    % Set sliders
    for sli = 1:size(pan_gui.sliders,2)
        pan_gui.sliders{sli}.Value = gui2D.plotindex{sli};
        pan_gui.texts{sli}.String = sprintf('%i/%i', gui2D.plotindex{sli},...
            gui2D.dim(sli+2));
    end
end



% FIGURE DATA % ------------------------------------------- %

if save_all == 1
    % Data slice dimensions
    dim = csi.data.dim(4:end); 
    % All figure indexes to cover
    indArray = slice2rowIndex(num2cell(dim));
elseif save_all == 0
    % Data slice dimensions
    dim = csi.data.dim(4:end);
    % Get index of current figure - Loop over each available slider
    indArray = NaN(1,numel(dim));
    for sli = 1:size(pan_gui.sliders,2)
        indArray(sli) = pan_gui.sliders{sli}.Value;
    end
    
elseif save_all == 2
    % Data slice dimensions
    dim = csi.data.dim(4:end);
    dimtmp = dim;
    
    % Get current slice    
    sloi = pan_gui.sliders{1}.Value;
    
    % Create indexes to loop over keeping "slice" static
    dimtmp(1) = [];
    if ~isempty(dimtmp)
        indArray = slice2rowIndex(num2cell(dimtmp));
        indArray_Slice = repmat(sloi,size(indArray,1)',1);
        indArray = cat(2,indArray_Slice, indArray);
    end
    
end


% LOOP & SAVE % -------------------------------------------- %
% Loop each index & Save the figure
for indi = 1:size(indArray,1)
    % 1. Set index at the slider
    for sli = 1:size(pan_gui.sliders,2)
        pan_gui.sliders{sli}.Value = indArray(indi,sli);
        pan_gui.texts{sli}.String = ...
            sprintf('%i/%i', indArray(indi,sli),dim(sli));
    end
    
    % 2. Replot 2D CSI plot.
    CSI_2D_initiate2D();
    
    % 3. Get new 2D CSI plot figure object
    fig_obj = findobj('type','figure','tag','CSIgui_plot2D');
    
    % 4. Add a index-numbers to the figure.
    if showIndex        
        nr = sprintf('%02i|', indArray(indi,:)); nr(end) = [];
        ax = axes(fig_obj, 'position', [0 0 0.10 0.05], 'Color','None');
        ax.Visible = 'off'; ax.YLim = [0 1]; ax.XColor = 'none'; 
        ax.YColor ='none'; ax.XTick =[]; ax.YTick = []; 
        text(ax,0,0.5,nr, 'Color', gui.colors.text_main,'FontSize',12);
    end
    
    number_title = sprintf('%03i', indArray(indi,:));

    % 4. Save the figure to image file    
    if strcmp(ext,'.fig')        
        savefig(fig_obj,[fp fn '_' number_title ext]);
    else
        
       % TRANSPARENCY % --------------------------------------------------- % 
       if transparency
               export_fig([fp fn '_' number_title ext], ...
                '-transparent',['-' ftype],'-nocrop', '-m2', fig_obj);
       else
           scrsz = get(0,'ScreenSize'); scrsz = scrsz(3:4);
           if fast
               pause(0.3); % Make sure display is updated.

               % Solution is screenshot!
               figpos = fig_obj.Position; % [l b w h]
               figpos(2) = abs(figpos(2) + figpos(4) - scrsz(2));
               img = screenshot(figpos); % [l t w h]

               % Save to file
               imwrite(img, [fp fn '_' number_title ext]);
           else
               export_fig([fp fn '_' number_title ext], ...
                           ['-' ftype],'-nocrop', '-m2', fig_obj);
           end
       end
    end
    % Remove axes used to show index-numbers
    if showIndex, delete(ax); end
    pause(0.005)
    
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


% CSI Backup System % ----------------------------------------------- %
% ---------------------------------------------------------------------- %

% --- Executes on button press in button_backupSet.
function button_backupSet_Callback(~, ~, gui)
% Store the csi.data field into csi.backup
% csi.backup.(bu time).data = csi.data;
% csi.backup.(bu time).tag  = info_str;
CSI_backupSet(gui, 'User backup.');

% --- Executes on button press in button_backupGet.
function button_backupGet_Callback(~, ~, gui)
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
bup_time = sprintf('%s', string(datetime('now','Format','HH:mm:ss')));
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
loc = contains(bup_list, 'lastBackup'); bup_list(loc) = [];

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
    uans = getUserInput_Popup({'Available backups: '}, {bup_disp},...
                              [], 'Backup Manager');
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

% --- Executes on button press in button_backupDel.
function button_backupDel_Callback(~, ~, gui)
% Delete specific backups made to clear up memory.
% Check for CSI data
CSI_backupDel(gui, 1)

% --- Delete chosen backup or oldest one.
function CSI_backupDel(gui, askUser)

% Without askUser input, automatically delete oldest backup.

if nargin == 1, askUser = 0; end

if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main,'csi');

% Check for backup data
if ~isfield(csi,'backup')
    CSI_Log({'No MRSI backup available.'},{''}); return;
end

% AVAILABLE BACKUP CHOICE % -------------------- %

% Get list of available backup fields
bup_list = fieldnames(csi.backup);
if numel(bup_list) < 1, return; end

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
    uans = getUserInput_Popup(...
        {'Available backups to delete:'},{bup_disp},[], 'Backup Manager');
    if isempty(uans), return; end

    % Process user choice: backup of interest
    boi_ind = strcmp(uans, bup_disp); boi = bup_list{boi_ind};
    tagoi = tag{boi_ind};
    
else % NO user input required.
    % Delete oldest backup
    boi = bup_list{1}; tagoi = csi.backup.(boi).tag;        
end

% DELETE BACKUP % --------------------------- %
csi.backup = rmfield(csi.backup, boi);

% SAVE AND CLEAN UP % --------------------------- %
% Store appdata.
setappdata(gui.CSIgui_main, 'csi', csi);
% Show user.
CSI_Log(...
    {['Deleted backup ' strrep(boi(8:end),'_',':') ]},{tagoi});

% --- Revert to previous backup: ctrl+z or menubar
function CSIgui_menubar_backupUndo(hObj,~,~)
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

% --- Executes on button press in checkbox_backup.
function checkbox_backup_Callback(~, ~, ~)
% This checkbox is used to turn auto-backup on or off.
% Use the following code to enable this in a function
%
% Create backup
% backup = gui.checkbox_backup.Value; 
% if backup, CSI_backupSet(gui, 'Before "this function"'); end



% IMG BUTTONS % ----------------------------------------------------- %
% ------------------------------------------------------------------- %

% ------------------------------------------------------------------- %
% Some IMG coodinate functions are located at COORDINATES paragraph which
% also containes the CSI coordinates function(s).
% ------------------------------------------------------------------- %

% --- Executes on button press in button_MRI_PlotIMG.
function button_MRI_PlotIMG_Callback(hObj, ~, ~)

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
function button_MRI_setContrast_Callback(~, ~, gui)
% Set the contrast for plotting the image behind the CSI data

% Get data.
if ~isappdata(gui.CSIgui_main, 'mri'), return; end
mri = getappdata(gui.CSIgui_main, 'mri');
if ~isappdata(gui.CSIgui_main, 'conv'), return; end
conv = getappdata(gui.CSIgui_main, 'conv');

% Get current contrast or calculate.
if isfield(conv, 'contrast')
    disp_cont =  conv.contrast; 
else
%     fn = fieldnames(conv);
%     disp_cont = [min(conv.data.(fn{1})(:)) max(conv.data.(fn{1})(:))];
    disp_cont = [min(conv.data(:)) max(conv.data(:))];
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


% Close CSIgui % ---------------------------------------------------- %


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


% GUI COLOR THEME % ------------------------------------------------- %

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
    colors.grid2      = [0.4 0.4 0.4];       % Grid color 2D plot
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
    colors.grid2      = [0.4 0 0];           % Grid color 2D plot
elseif strcmp(custom, 'on')
   
   % Read colors from file
   colors = setGUIcolor_custom_read;

end

% Special object to apply highlight colors
obj_hilight1 = {'txt_infoCSI',  'txt_infoIMG', 'button_ws'};   % Text hili
obj_hilight2 = {'button_plotCSI',  'button_MRI_plotIMG'};          % BG hili

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
            
        % If its a uicontrol object: GUI control
        elseif strcmp((gui.(fn{fi}).Type),'uicontrol') 
            
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
        elseif strcmp((gui.(fn{fi}).Type),'uipanel')
            gui.(fn{fi}).BackgroundColor = colors.main;
            gui.(fn{fi}).ForegroundColor = colors.text_title;
            gui.(fn{fi}).HighlightColor = colors.text_title;

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



%%% MERGE VOXELS % -------------------------------------------------- %

% --- Executed if user presses merge voxels button
function MergeVoxels_Initiate(hObj, gui)
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
MergeVoxels_plotVoxels(hObj, gui, plot_par);

% --- Executes on button press in button_CSI_MergeVoxels.
function button_CSI_MergeVoxels_Callback(hObj, ~, gui)
MergeVoxels_Initiate(hObj, gui);

% --- Executes by CSI_MergeVoxels_Initiate
function MergeVoxels_plotVoxels(hObj, gui, plot_par)
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
              'Callback',@MergeVoxels_Button_Merge);                  
    % Save selected voxel data or figure
    uicontrol('Parent',save_data.tabh{sli},'Style','pushbutton','String', 'Save',...
              'BackgroundColor',plot_par.colors.main,...
              'ForegroundColor',gui.colors.text_main,...
              'Position',[5 20 50 15],...
              'Callback',@MergeVoxels_Button_Save);

end % ------ % End of Slice Loop % ------ %

save_data.highlight = gui.colors.hilight1;
% Store data to figure.
guidata(plot_par.fh, save_data);

loadBar(NaN);

% --- Executed by functions related to CSI_MergeVoxels
function selected = MergeVoxels_getSelected(hObj)
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
function MergeVoxels_Button_Merge(hObj,~)

% --- % Get the selected voxels data
sel = MergeVoxels_getSelected(hObj);

% --- % Userinput

quest = {'SNR Filtering','Aligning','Weighting'};
uans = getUserInput_Radio(quest,ones(1,size(quest,2)));
if isempty(uans), CSI_Log({'Skipped merging voxels.'},{''}); return; 
end  


% --- % Apply options
if sum(uans) > 0 
    
    % --- % Get POI
    [sel.poi, sel.pex] = MergeVoxels_POI(sel.xaxis);
    
    % --- % SNR Filtering
    if uans(1), sel = MergeVoxels_SNR_Filter(sel);  end

    % --- % Aligning
    if uans(2), sel = MergeVoxels_Align(sel); end

    % --- % Weighting
    if uans(3), sel = MergeVoxels_Weighted(sel); end

else
    % No Options
end

% --- % Average & display
MergeVoxels_Average(sel);

% --- % Executes on button press: Save @ CSI_MergeVoxel figure
function MergeVoxels_Button_Save(hObj,~)

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
    MergeVoxels_SaveData(hObj, fp, fn); 
elseif idx ==2 % Save Figure
    MergeVoxels_SaveFig(hObj,  fp, fn); 
elseif idx ==3 % Save both
    MergeVoxels_SaveData(hObj, fp, fn);
    MergeVoxels_SaveFig(hObj,  fp, fn); 
end

% --- % Executed by CSI_MergeVoxels_Button_Save
function MergeVoxels_SaveData(hObj, fp, fn)
% Get the selected voxel data and save to file.

% --- % Get selected voxels.
selected = MergeVoxels_getSelected(hObj);

% --- % Store data
save([fp '\' fn '.mat'], 'selected');

% --- % LOG
CSI_Log({'Saved selected voxel data to:'},{[fp '\' fn '.mat']})

% --- % Executed by CSI_MergeVoxels_Button_Save
function MergeVoxels_SaveFig(hObj, fp, fn)

% ---- % Figure object
fh = hObj.Parent; 

% ---- % Save figure

matver = version('-release'); matyr = str2double(matver(1:end-1));
% Compact for MATLAB 2015 and later versions only.
if (matyr <= 2014), savefig(fh, [fp '\' fn '.fig']);
else,               savefig(fh.Parent.Parent, [fp '\' fn '.fig'], 'compact');
end

% LOG
CSI_Log({'Saved selected voxel figure to:'},{[fp '\' fn '.fig']});

% --- % Executes is user requests merge with options none
function MergeVoxels_Average(sel)
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
function [poi, pex] = MergeVoxels_POI(xaxis)
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
function sel = MergeVoxels_SNR_Filter(sel)
% Include and/or exclude voxels using SNR.

voxels = sel.voxels; pex_on = ~isnan(sel.pex); index = sel.index;

% ---- % Userinput
% SNR Limits

% SNR Limits to include or exclude
tags = {'Alignment: if > SNR Limit',...
        'Exclusion: exclude if < SNR Limit'}; 
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

% SNR limit for POI (minimum)
poi_snr_lim = str2double(uans{1});                  
% SNR limit for PEX
if pex_on, pex_snr_lim = str2double(uans{2}); else, pex_snr_lim = 0; end

% ---- % Phasing and SNR

% SNR of POI 
poi_snr = MergeVoxels_SNR(sel, 'alignment');         

% SNR of PEX
if pex_on
    seltmp = sel; seltmp.poi = sel.pex;
    pex_snr = MergeVoxels_SNR(seltmp, 'inclusion'); 
end

% ---- % Select Voxels by SNR

% For POI: peak to align to.
poi_incl = (poi_snr >= poi_snr_lim); 
poi_excl = (poi_snr <  poi_snr_lim);

% For PEX: peak to use for exclusion of voxel
if pex_on, pex_excl = (pex_snr <= pex_snr_lim); else, pex_excl = 0; end

% Merge SNR incl/excl from POI and PEX
vox_incl = poi_incl; vox_excl = poi_excl;
if pex_on
    vox_incl(pex_excl == 1) = 0; 
    vox_excl(pex_excl == 0) = 1;
end

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
         '#Included by SNR of POI/PEX:','-----------------------'},...
         {'',size(voxels,2),...
         [num2str(poi_snr_lim) ' | ' num2str(pex_snr_lim)],...
         sum(vox_incl), ...
         [num2str(sum(poi_incl)) ' | ' num2str(sum(pex_excl))],''});

% ---- % Set Output
sel.voxels = voi; sel.snr = voi_snr; sel.index = voi_ind;

% --- % Executed to calculate SNR and apply phasing
function snr = MergeVoxels_SNR(sel, poi_tag)
% Calculate SNR and apply phasing prior to calculations if necessary.
%
% poi_tag = 'Name for peak'

if nargin == 1, poi_tag = ''; end

voxels = sel.voxels; % Get voxel data.
poi = sel.poi;       % Peak of interest index

% ---- % Userinput
% Phasing

uans = getUserInput_Popup({['Apply phasing to ' poi_tag ' peak?']},...
                          {{'No','Yes'}}, [], 'Merge Voxels');
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
% Calc SNR for POI   
snr = csi_SNR(poi_vox_4snr, round(size(poi_vox_4snr,1)*0.1), 1, sel.poi); 

% --- % Executed to weight voxels to SNR
function sel = MergeVoxels_Weighted(sel)
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
function sel = MergeVoxels_Align2(sel)
% Frequency alignement of spectra from  multiple voxel.
%
% Sel-structure requires field voxels and peak of interest poi.
% This versions used the old MergeVoxels algorithm based on maximum value
% only.


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

% --- % Executed by CSI_MergeVoxels_Button_Align- OFF
function sel = MergeVoxels_Align(sel)
% Frequency alignement of spectra from  multiple voxel.
%
% Sel-structure requires field voxels and poi, in a 1x2 vector representing
% the integer unit-less index range of the peak of interest.
%
% Aligning is based on the spectrum with the largest signal at the POI
% (ref) and other spectra aligned to this spectrum using the correlation
% (2D) between ref and voxel.


plot_on = 1;
% Undress the input structure:
voxels = sel.voxels; poi = sel.poi; 

                   % ---- % Aligning voxels % ---- %
                    
% Spectrum with maximum signal at peak of interest (POI)
[voi_poi_max , voi_poi_max_pos] = ...
    max( real( voxels(poi(1):poi(2),:)), [], 1);
[~,template_ind] = max(voi_poi_max);


% Correlation % ------------- %
% Ref chosen by max value of spectrum/peak of interest

% If signal processing toolbox available; use xcorr from Matlab.
% else use xcorrq; simplified version which does the job.
tb = isToolbox('Signal Processing Toolbox');

template = voxels(:,template_ind); voi_shifts = size(voxels,2);
for si = 1:size(voxels,2)
    if tb
        [corval, lags] = xcorr(template, voxels(:,si));
    else
        [corval, lags] = xcorrq(template, voxels(:,si));
    end
    % Find most significant correlation for lag-distance
    [~, ind] = max(corval); voi_shifts(si) = lags(ind);
end

% Shift each voxel by its voi_shift amount
voi_cell = mat2cell(voxels, size(voxels,1), ones(1,size(voxels,2)));
voi_aligned_cell = ...
    cellfun(@circshift, voi_cell, num2cell(voi_shifts), 'uniform',0);
voi_aligned = cell2mat(voi_aligned_cell);

% ---- % Set Output
sel.voxels = voi_aligned;


% ---- % Additional plots

if plot_on == 1
    nsli = unique(sel.index(:,3)); nfig = numel(nsli);
    % Plot for each slice seperately.
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
            plA = plot(xax, real(voi_aligned(:,vox_ind)),'-r'); hold on;
            % Plot voxel
            plV = plot(xax, real(voxels(:,vox_ind)),'--b'); 

            % Axis cosmetics
            axh = plV.Parent; % 
            plV.Parent.XDir = 'reverse';

        end

        % Display data
        str = ['Slice ' num2str(nsli(kk)) ...
               ' | Avg shift: ' num2str(mean(voi_shifts(vox_ind_all)))];      
        uicontrol(fh, 'Style','text','Units','Normalized',...
            'Position',[0.35 0.925 0.3 0.05],...
            'string', str);    
    end
end

% ------------------------------------------------------------------- %
% --------------------------------------------------------------------- %

                            % MISCELLANEOUS %

                                           
% --- Executes on selection change in listbox_CSIinfo.
function listbox_CSIinfo_Callback(~, ~, ~)
% --- Executes during object creation, after setting all properties.
function listbox_CSIinfo_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in listbox_MRIinfo.
function listbox_MRIinfo_Callback(~, ~, ~)
% --- Executes during object creation, after setting all properties.
function listbox_MRIinfo_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popup_plotUnit.
function popup_plotUnit_Callback(~, ~, ~)
% --- Executes during object creation, after setting all properties.
function popup_plotUnit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in popup_plotIMG.
function popup_plotIMG_Callback(~, ~, ~)
% --- Executes during object creation, after setting all properties.
function popup_plotIMG_CreateFcn(hObject, ~, ~)
% hObject    handle to popup_plotIMG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ------------------------------------------------------------------- %
            %   Everything below here is new OR beta!  %
                      % Possibly not finished %
% ------------------------------------------------------------------- %

% --- Executes on button press in button_TestSomething.
function button_TestSomething_Callback(hObj, evt, gui)

warndlg('Watch out! Developer testing button. Panic!');

if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

csi.data.raw = single(csi.data.raw);
setappdata(gui.CSIgui_main, 'csi', csi);


% --- Simple sinc function
function an = sincQ(x)
an = sin(x)./x;

% --- Executes on button press in button_CSI_VoxelsShift.
function button_CSI_VoxelsShift_Callback(~, ~, gui)
% Subtract or add voxels to the CSI coordinates grid or by manipulating the
% data indexing.

uans = ...
    getUserInput_Buttons('Shift MRS data by manipulating:',...
                        {'Coordinate grid', 'Data index'});
if isempty(uans), return; end

switch uans
    case 'Data index'
        CSI_voxelShift_Indexing(gui);
    case 'Coordinate grid'        
        CSI_voxelShift_Coordinates(gui);
end
             
% --- Shift coordinate grid of CSI volume to apply a voxel shift
function CSI_voxelShift_Coordinates(gui)
% Shift N-voxels by manipulating the cooridnate grid of the CSI data. 

% Get csi data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');
if ~isfield(csi, 'ori'), return; end

% Get orientation data
ori = csi.ori;

% Get requested shift
uans = getUserInput({'Voxel shift for each direction?'},{'0.5 0.5 0.5'});
if isempty(uans), return; end
shft = str2double(strsplit(uans{1}));

% Get previous voxel-coordinate calculation shift options.
opts.fft_cor = ori.fft_cor; opts.vox_cor = ori.vox_cor;

% Calculate coordinates.
csi.ori = CSI_coordinates_calculate...
    (ori.res, ori.offcenter, ori.dim, shft, opts);

% Save to appdata
setappdata(gui.CSIgui_main,'csi',csi);   
gui = guidata(gui.CSIgui_main);

% Recalculate all conversion data.
MRI_to_CSIspace(gui);

% --- Executes by selecting this method in CSI_VoxelsShift_Callback
function CSI_voxelShift_Indexing(gui)

% Get csi data
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before voxel shift.'); end

% Dimension to shift
spat_dim = csi_findDimLabel(csi.data.labels,{'kx','ky','kz','x','y','z'});
spat_dim(isnan(spat_dim)) = [];
if isempty(spat_dim), spat_dim = [2 3 4]; end

% Ask user #shift and which dimension
uans = getUserInput(...
    {   {'Spatial index in MRSI data: (kx ky kz) '},...
        {'Number of voxels to shift data: (x, y, z | col, row, slice'}},...
        {spat_dim, '1 0 0'} );
if isempty(uans), return; end    
% Proces input
kspace = str2double(strsplit(uans{1}));
shifts = str2double(strsplit(uans{2}));  

% Apply shift
for kk = 1:size(kspace,2)
    % Kk +1 to correct for the time dimension at index 1.
    if shifts(kk) ~= 0
    csi.data.raw = circshift(csi.data.raw,shifts(kk),kspace(kk));
        CSI_Log({'Circular shift #voxels | N-dimension:'},...
                {[num2str(shifts(kk)) ' | ' num2str(kspace(kk))]})
    end
end


% Save to appdata
setappdata(gui.CSIgui_main,'csi',csi); 

% --- Executes on selection change in popup_domain.
function popup_domain_Callback(hobj, ~, ~)
str = hobj.String{hobj.Value}; CSI_setDomain(hobj, [],  str);

% --- Executes during object creation, after setting all properties.
function popup_domain_CreateFcn(hobj, ~, ~)
if ispc && isequal(get(hobj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hobj,'BackgroundColor','white');
end

% --- Executes on button press in button_CSI_FrequencyAlignment.
function button_CSI_FrequencyAlignment_Callback(~, ~, gui)
% Align the seperate voxels to a specific peak e.g. frequency alignment or
% peak alignment.

disp('This function is under construction')

domain = CSI_getDomain(gui);
if strcmp(domain, 'time')
    CSI_Log({'MRS data is in time domain; '},...
            {'change it to frequency domain to apply peak alignment.'});
    return;
end 

% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before alignment.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Userinput % -------------------------------------------- %

% Get data of peak of interest
[doi, ~, ~] = CSI_getDataAtPeak(csi.data.raw, csi.xaxis);

% Which voxel to use as reference? Or highest SNR?
qry = {'Reference voxel or voxels for shift: '};
def = {{'Single voxel','Along a dimension'}};
uans = getUserInput_Popup(qry, def, [], 'Frequency Alignment');
if isempty(uans), CSI_Log({'Skipped alignment.'},{''}); return; end

switch uans{1}
    case 'Single voxel'
        
        % Get reference voxel
        uans = getUserInput(...
           {'Voxel index: '},{ones(1,numel(size(csi.data.raw))-1)});
        vox_ref = num2cell(str2double(uans{1}));
        
        % Get peak indices
        [~, ind_peaks] = max(real(doi),[],1);

        % Get reference index
        ind_ref = ind_peaks(:,vox_ref{:});
        
        % Shift values
        shiftval = ind_peaks - ind_ref;
               
        % Create cell of raw data
        sz = size(csi.data.raw);
        cell_layout = ...
        arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
        % Creat cell of data.
        data = mat2cell(csi.data.raw, sz(1), cell_layout{:}); 
        
        % circular shift
        csi.data.raw = cell2mat(cellfun(@circshift, data, num2cell(-1.*shiftval),...
            repmat({1},size(data)),'uniform', 0));

        
     case 'Along a dimension'
         
        % Get reference dimension
        uans = getUserInput_Popup({'Reference dimension:'},...
                                  {csi.data.labels(2:end)},...
                                  [], 'Frequency Alignment');
        if isempty(uans), return; end
        % Dimension of interest for voxel-reference
        ind_dim = find(strcmp(csi.data.labels, uans{1})==1);
         
        % Get reference voxel
        uans = getUserInput_Popup({'Reference voxel in this dimension:'},...
                                  {{1:size(csi.data.raw,ind_dim)}},...
                                  [], 'Frequency Alignment');
        if isempty(uans), return; end        
        % Reference voxel for each voxel along ind_dim to shift to
        ind_vox = str2double(uans{1});
        
          
             
        % Get peak indices
        [~, ind_peaks] = max(real(doi),[],1);
        
        % Get reference voxel over ind_dim
        sz = size(csi.data.raw);
        tmp_ind = arrayfun(@(x) 1:x,sz,'uniform', 0);
        tmp_ind{ind_dim} = ind_vox; tmp_ind{1} = 1;
        % Reference to shift to  
        ind_ref = ind_peaks(tmp_ind{:});
        
        % Shift values
        shiftval = ind_peaks - ind_ref;
        
        % Convert data to cell
          % Create cell of raw data
        sz = size(csi.data.raw);
        cell_layout = ...
        arrayfun(@ones, ones(1,size(sz(2:end),2)),sz(2:end),'Uniform',0);
        % Creat cell of data.
        data = mat2cell(csi.data.raw, sz(1), cell_layout{:}); 
        
        % Apply shift
        csi.data.raw = cell2mat(cellfun(@circshift, ...
            data, num2cell(-1.*shiftval), repmat({1},size(data)), ...
            'Uniform', 0));
end

        
% Store data
setappdata(gui.CSIgui_main,'csi', csi);

% Show nfo
CSI_Log({'Alignment unique shift values: '}, {unique(-1.*shiftval)'});

% --- Executes on button press in button_CSI_AutoProcessing_openScript.
function button_CSI_AutoProcessing_openScript_Callback(hObj, evt, gui)
% Load and execute autoscripting... 

% Load text file
rt = mfilename('fullpath');
fp = fileparts([rt '.m']);

if exist([fp '\Scripts'], 'dir') == 7, fp = [fp '\Scripts']; end

[fn,fp,ind] = uigetfile('*.txt', 'defname', fp);
if ind == 0, return; end

fid =  fopen([fp fn]);
inp = fscanf(fid,'%s\n');
fclose(fid);
inp = strsplit(inp,';'); inp(end)=[];

butfunc = fieldnames(gui); butfunc(~contains(butfunc,'button')) = [];

% Process text file
for kk = 1:size(inp,2)
    tmp = inp{kk}; tmp = strsplit(tmp,'='); %1=func 2=inp
    
    inpfnc = tmp{1}; 
    if length(tmp)>1, inpvar = tmp{2}; else, inpvar=[]; end


    if strcmp('plotCSI', inpfnc) || strcmp(inpfnc(1:3), 'MRI')
        qrystr = ['button_' inpfnc];
    else
        qrystr = ['button_CSI_' inpfnc];
    end
    qry = strcmp(butfunc, qrystr);
    
    if sum(qry)==0
        switch inpfnc % Put all exceptions here
            case 'setDomain'
                CSI_setDomain(hObj, evt, inpvar); 
            case 'saveData'
                CSI_saveData(hObj,[],[]);
            case 'backupDel'
                CSI_backupDel(gui);
            case 'setBackup'
                gui.checkbox_backup.Value = str2double(inpvar);
        end
    else
       btn = butfunc(qry);
       if ~isempty(btn) 
           fnc2run = [btn{:} '_Callback'];
           eval([ fnc2run '(hObj,evt,gui)']);
       else
           CSI_Log({'Could not parse scripted input:'},...
                   {[inpfnc ' | ' btn{:}]});
       end
        
    end

end

% --- Executes on button press in button_CSI_Delete.
function button_CSI_Delete_Callback(~, ~, gui, backup)
% Delete a part of the loaded MRS data matrix.

if nargin < 4, backup = gui.checkbox_backup.Value; end

% Create backup
if backup, CSI_backupSet(gui, 'Before deleting data.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get index to delete
qry = {'Data index/dimension to delete from:'};
def = {csi.data.labels(2:numel(size(csi.data.raw)))};
uans = getUserInput_Popup(qry, def, [], 'Delete MRS data');
if isempty(uans), CSI_Log({'Aborted deleting data.'},{''}); return; end
ind = find(strcmp(uans{1},csi.data.labels)); lab = uans{1};

qry = {'Specify data index to delete: (M:N or mutliple allowed)'};
def = {'1'};
uans = getUserInput(qry, def, [], 'Delete MRS data');
if isempty(uans), CSI_Log({'Aborted deleting data.'},{''}); return; end
tmp = strfind(uans{1},':');
if isempty(tmp)
    tbdeleted = str2double(strsplit(uans{1}));
else
    tmp = strsplit(uans{1},':');
    tbdeleted = str2double(tmp{1}):str2double(tmp{2}); 
end

% Delete!
sz = size(csi.data.raw); dimind = arrayfun(@(x) 1:x, sz, 'uniform', 0);
tmp = dimind{ind}; tmp(tbdeleted) = []; dimind{ind} = tmp;
% Update data and dim
csi.data.raw = csi.data.raw(dimind{:});
csi.data.dim = size(csi.data.raw);

% Store data
setappdata(gui.CSIgui_main,'csi', csi);

% Update X-axis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show nfo
CSI_Log({'Deleted data (dimension|index): '}, {[lab ' | ' int2str(tbdeleted)]});

% --- Executes on button press in button_CSI_Concatenate.
function button_CSI_Concatenate_Callback(~, ~, gui)
% Concatenate specific dimensions of the MRS data matrix.
% Create backup
backup = gui.checkbox_backup.Value; 
if backup, CSI_backupSet(gui, 'Before concatenating data.'); end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get indexes to merge
qry = {'Dimension to merge #1:','Dimension to merge #2:'};
def =  { csi.data.labels(2:numel(size(csi.data.raw))),...
                            csi.data.labels(2:numel(size(csi.data.raw)))};
uans = getUserInput_Popup(qry, def, [], 'Concatenate CSI');
if isempty(uans), CSI_Log({'Aborted concatening data.'},{''}); return; end

% Process userinput
ind1 = find(strcmp(uans{1},csi.data.labels)); lab1 = uans{1};
ind2 = find(strcmp(uans{2},csi.data.labels)); lab2 = uans{2};
if isequal(ind1,ind2)
    CSI_Log({'Cant merge the same indexes!'},{''}); return; 
end

% Merge!
sz = size(csi.data.raw); 
sznew = sz; sznew(ind1) = sz(ind2).*sz(ind1); sznew(ind2) = 1;
nDimC = cellfun(@(x) 1:x, num2cell(sz),'UniformOutput',0); tmp = [];
for aa = 1:sz(ind2)
    tmpind = nDimC; tmpind{ind2} = aa;
    if aa == 1
        tmp = csi.data.raw(tmpind{:});
    else
        tmp = cat(ind1,tmp,csi.data.raw(tmpind{:}));
    end    
end
csi.data.dim = sznew; csi.data.raw = tmp;

% Store data
setappdata(gui.CSIgui_main,'csi', csi);


% Show nfo
CSI_Log({'Merged data (dimension|index): '}, ...
    {[lab1 ' & ' lab2 ' | ' int2str(ind1) ' & ' int2str(ind2)]});

% --- Executes on button press in button_CSI_FAdynamic.
function button_CSI_FAdynamic_Callback(~, ~, gui)
% Process and start fit of FA dynamic data.
%
% Required userinput: peak of interest, sine or MR-signal equation, plot 
% and save images, dimension of FA-dynamic and add S(0) = zero.

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'),return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Data
data = csi.data.raw;

% USERINPUT % -------------------------------------------------------- %
% Settings and data-variables

% Get peak of interest - (data of interest; doi)
[doi, doi_axis] = CSI_getDataAtPeak(data, csi.xaxis);
if isnan(doi), CSI_Log({'Aborted FA dynamic.'},{''}); return; end

% ----------------------------- %
% // User Input: settings ----- %
% ----------------------------- %

qry = {'Flip angle dynamics index:',...
       'Zeroth-order phase corrections:', ...
       'Fit equation:',...
       'Add data point S(FA=0) = 0:',...
       'Display graphs:', 'Save figures' };
def = {csi.data.labels(2:numel(size(csi.data.raw))),...
      {'No','Full','After zero-crossing'},...
      {'S(FA): Simple Sine',...
       'S(FA,TR,T1): MRS Signal equation'},...
      {'Yes','No'},{'Yes','No'}, {'Yes','No'}};
uans = getUserInput_Popup(qry, def, [], 'FA Dynamic');
if isempty(uans), CSI_Log({'Aborted FA dynamic.'},{''}); return; end

% ---------------------- %
% // Process Input ----- % 
% ---------------------- %

% USER: Index of FA-dynamic
ind = find(strcmp(uans{1},csi.data.labels)); lab = uans{1};

% USER: Zeroth-order phasing
switch uans{2}
    case 'Full',phase_data = 1; 
    case 'No', phase_data = 0;
    case 'After zero-crossing', phase_data = 2;
end

% USER: Fit Equation
switch uans{3}
    case 'S(FA): Simple Sine',               fit_method = 0;
    case 'S(FA,TR,T1): MRS Signal equation', fit_method = 1;
end

% USER: Add data zero @ time zero
switch uans{4}, case 'Yes', add_zero = 1; case 'No', add_zero = 0; end

% USER: Plot
switch uans{5}, case 'Yes', plot_data = 1; case 'No', plot_data = 0; end

% USER: Plot
switch uans{6}, case 'Yes', save_plot = 1; case 'No', save_plot = 0; end

% -------------------------------- %
% // User Input: Flip angles ----- %
% -------------------------------- %

% Load previous FA-input if available or calculate
if isfield(csi.xaxis,'FA')
    FAdef = int2str(csi.xaxis.FA);
elseif isfield(csi,'nfo') && isfield(csi.nfo, 'startingFA_deg')
    FAdef = csi.nfo.startingFA_deg + ...
        ( (1:csi.nfo.nrFAs)-1 ).*(csi.nfo.FASpacing_deg);
else 
    FAdef = '';
end

% Get FA-values from user
uans = getUserInput(...
    {'Start FA:','Stepsize:','Custom FA:'}, {'50','50', FAdef});
if isempty(uans), CSI_Log({'Aborted FA dynamic.'},{''}); return; end

% ----------------------------------- %
% // Process Input: Flip angles ----- % 
% ----------------------------------- %

% Calculate flip angles for FA-series
FAstart = str2double(uans{1}); FAstep = str2double(uans{2});
if ~isempty(uans{3})
    if contains(uans{3}, ':')
        xdat_parts = strsplit(uans{3}); perPart = cell(size(xdat_parts));
        for kk = 1:size(xdat_parts,2)
            perPart{kk} = str2double(strsplit(xdat_parts{kk}, ':'));        
            if numel(perPart{kk}) == 3
                perPart{kk} = perPart{kk}(1):perPart{kk}(2):perPart{kk}(3);
            end
        end
        xdat = cell2mat(perPart);
    else
        xdat = str2double(strsplit(uans{3}));
    end 
else
    FAstop = ((FAstep.*(size(doi,ind)-1) ) + FAstart);
    xdat  = (FAstart:FAstep:FAstop);
end
csi.xaxis.FA = xdat; % Store xaxis data.
setappdata(gui.CSIgui_main,'csi', csi);   

% ------------------------------ %
% // User Input: TR and T1 ----- %
% ------------------------------ %

if fit_method == 1
    
    % --- % Check if protocol nfo is available
    if isfield(csi,'nfo')
        if isfield(csi.nfo, 'TR')
            if strcmpi(csi.nfo.TR, 'user defined')
                fntmp = fieldnames(csi.nfo); 
                TRind = find(strcmp(fntmp,'TR') == 1);
                TR = (csi.nfo.(fntmp{TRind+1}));
            else
                TR = (csi.nfo.TR);
            end
        end
    else
        TR = '50';
    end
    
    uans = getUserInput({'TR(ms):','T1(ms):'}, {TR,'6000'});
    if isempty(uans), CSI_Log({'Aborted FA dynamic.'},{''}); return; end
    TR = str2double(uans{1}); T1 = str2double(uans{2});
else
    TR = []; T1 = [];
end

% PREPARE DATA % ------------------------------------------------------ %
                     
% MRSI data to cell format
% {#samples x #spectra} with  #spectra == #FA;

% 1. Set FA-index to second dimension.
sz = size(doi); permv = 1:numel(sz); permv(2) = ind; permv(ind) = 2;
doi = permute(doi, permv);

% 2. Create cell-arrays
sz = size(doi); cell_layout = arrayfun(@ones,...
    ones(1,size(sz(3:end),2)),sz(3:end),'UniformOutput',0);
doi_cell = mat2cell(doi, sz(1),sz(2), cell_layout{:});
doi_cell = squeeze(doi_cell);

% Data counter
iterNumber = num2cell(find(ones(size(doi_cell)) == 1));
iterNumber = reshape(iterNumber, size(doi_cell));

% Directory for output: data and figures
daystr = string(datetime('now','Format','yyMMddHHmm'));
dirstr = strjoin(['FAdynamic_' daystr],'');
% Create output directory
if exist(dirstr,'dir') ~= 7, mkdir(dirstr); end

% FIT FUNC % -------------------------------------------------------- %

% Call the fit-function CSI_FAdynamic
outp = cellfun( @csi_fit_FASeries, ...
    repmat({xdat}, size(doi_cell)),...        % xdata (FA)
    doi_cell,...                              % data    
    repmat({fit_method},size(doi_cell)),...   % fit_method
    repmat({phase_data},size(doi_cell)),...   % phase_data boolean
    repmat({add_zero},size(doi_cell)),...     % Add S(FA=0)=0;
    repmat({plot_data},size(doi_cell)),...    % plot_data boolean
    repmat({save_plot},size(doi_cell)),...    % save_plot boolean
    iterNumber,...                            % Iterator number
    repmat({dirstr}, size(doi_cell)),...      % Directory for output
    repmat({TR},size(doi_cell)),...           % TR
    repmat({T1},size(doi_cell)),...           % T1
    repmat({doi_axis},size(doi_cell)),...     % ppm
    'UniformOutput', 0);      

% PROCESS OUTPUT % ---------------------------------------------------- %

% Save output-data to file
save(strjoin([dirstr '\FAdynamic_OutputAll.mat'],''), 'outp');

% Display output in CSI-log
outp = squeeze(outp); outp = reshape(outp,[],1);
fns = fieldnames(outp{1});
for kk = 1:size(outp,1) % Outp loop
    CSI_Log({''},{''});
    CSI_Log({'% ------------------------------- %'},{''});
    for li = 1:size(fns,1) % Field loop
        tmp = outp{kk};
        CSI_Log({[fns{li} ': ' ] },{tmp.(fns{li})});
    end
    CSI_Log({sprintf('%i. ----------------------------- %%',kk)}, {''});                          
end
CSI_Log({''},{''});

% Show statistics
zc_ext = cell2mat(extractField(outp, 'zero_crossing_fit'));
if numel(zc_ext) > 1
stats = csi_statistics_of_volume(zc_ext);
CSI_Log({'Zero-crossing statistics ---------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});
end

% --- Executes on button press in button_CSI_RemoveOS.
function button_CSI_RemoveOS_Callback(~, ~, gui, backup)
% A very crude way of removing oversamling from the data:
% Given a new value for dimension of interest doi, cut the data from start
% to the new value and store it.
%
% Possible upgrade: allow removing/cutting of a range.
% Example: column of size 11 requires OS removal of first data point and
% last two data points to end up with column size 8.

if nargin < 4, backup = gui.checkbox_backup.Value; end

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Check data domain
domain = CSI_getDomain(gui); doFFT = 0;
if ~strcmp(domain, 'time')
    CSI_Log({'MRS data is possibly in frequency domain;'},...
            {'Removing oversampling requires time domain data.'});
    qry = {'MRS-data is possibly in frequency domain:'};
    def = {{'Convert', 'Ignore', 'Abort'}};
    uans = getUserInput_Popup(qry,def,[],'Remove OS');
    if isempty(uans) || strcmp(uans, 'Abort'); return; end
    if strcmp(uans{1}, 'Convert'), doFFT = 1; end
end 

% Create backup
if backup, CSI_backupSet(gui, 'Before removing oversampling.'); end

% ------------------------------------------------- %

% Get dimensions to remove OS: ind
uans = getUserInput_Popup({'Index:'}, {csi.data.labels}, [], 'Remove OS');
if isempty(uans), CSI_Log({'Aborted OS removal.'},{''}); return; end
ind = find(strcmp(uans{1},csi.data.labels)); ind_lab = uans{1};

% Get method of removal
qry = {'Split data by:'};
def = {{'Split','Odd values','Even values'}};
uans = getUserInput_Popup(qry, def, [], 'Remove OS');
if isempty(uans), CSI_Log({'Aborted OS removal.'},{''}); return; end
switch uans{1}
    case 'Split',       method = 'Split';
    case 'Odd values',  method = 'Odd';
    case 'Even values', method = 'Even';
end

% Get value of OS to remove: osval
uans = getUserInput({'New size of dimension:'},...
                    {round(csi.data.dim(ind)./2)});
if isempty(uans),CSI_Log({'Aborted OS removal.'},{''}); return; end
osval = str2double(uans{1});

% Apply on given dim: csi.data.raw @ ind = 1:osval
sz = csi.data.dim;
if ~strcmp(method, 'Split')
    % Bandwidth changes if you take odd or even values from data.
    csi.xaxis.BW = csi.xaxis.BW / (sz(1) / osval);
end
% Index vector for selected method
switch method
    case 'Split'        
        sz(ind) = osval;
        rng = cellfun(@(x) 1:x, num2cell(sz), 'uniform', 0);
    case 'Odd'               
        rng = cellfun(@(x) 1:x, num2cell(sz), 'uniform', 0);
        rng{1} = 1:2:osval*2;
    case 'Even'        
        rng = cellfun(@(x) 1:x, num2cell(sz), 'uniform', 0);
        rng{1} = 2:2:osval*2;
end

% Possible iFFT
if doFFT, csi.data.raw = csi_ifft(csi.data.raw); end

% Cut the data with the new value
csi.data.raw = csi.data.raw(rng{:});

% Update csi.data.dim size 
csi.data.dim = size(csi.data.raw);

% Possible FFT
if doFFT, csi.data.raw = csi_fft(csi.data.raw); end

% ------------------------------------------------- %

% If a peak was selected - this needs to be removed as its point changed.
if isappdata(gui.CSIgui_main,'CSIpar')
    rmappdata(gui.CSIgui_main,'CSIpar');
end

% Store data %
setappdata(gui.CSIgui_main,'csi', csi);

% Update X-axis data
CSI_2D_Scaling_calc_xaxis(gui.CSIgui_main,[],1);

% Show nfo
CSI_Log({'Removed oversampling (dimension | size | Method): '}, ...
        {[ind_lab, ' | ', num2str(osval) ' | ' method]});

% --- Executes on button press in button_CSI_MapB1.
function button_CSI_MapB1_Callback(~, ~, gui)
% Choose a method to calculate B1.
% Double CSI - use FA = invcos(M(2a) / 2M(a));
% Multi CSI - use signal equation and fit FA, T1 and M0.

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
% csi = getappdata(gui.CSIgui_main, 'csi');


% User Input % --------------------------------------------------------- %
qry = 'Choose a B1-method:';
opt = {'Double-CSI', 'Multi-CSI'};
uans = getUserInput_Popup({qry},{opt}, [], 'Calculate B1-Map');
if isempty(uans), return; end

switch uans{1}
    case 'Multi-CSI'
        CSI_MapB1_MultiCSI(gui);
    case 'Double-CSI'
        CSI_MapB1_DoubleCSI(gui);
end

% --- Executes by button_CSI_MapB1
function CSI_MapB1_MultiCSI(gui)
% Calculate B1 maps using multiple-CSI protocol and fitting the CSI-signal
% equation for FA, T1 and M0.
%
% Requires TR. 
% T1 range of upper lower bound
% FA range of upper lower bound
%
%
% Input: gui-handle of CSIgui

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% User Input % ------------------------------------------------- %

% Read TR from file, Get stored TR, if available.
if isfield(csi.xaxis,'TR'), TR = csi.xaxis.TR;
elseif strcmp(csi.ext, 'dat') || isfield(csi,'twix')
    fn = strtrim(strsplit(csi.filename, '|'));
    TR = NaN(1,numel(fn));
    for kk = 1:numel(fn)
        twix = mapVBVD([csi.filepath fn{kk}]);
        TR(kk) = twix.hdr.Config.TR./(1e6);        
    end
    clear('twix'); 
else, TR = [0.3 0.5 0.65 0.8 1 1.2 1.5 2 4 6 8 10 12 14]; 
end

% Get TR, T1, Zero, Parallel - Add FA?
ele = {'popup', 'edit', 'edit', 'edit', 'edit', 'edit', 'popup',...
       'popup', 'popup', 'popup', 'popup', 'popup'};
qry = {'TR Index:', 'Give a TR range [s]:',...
       'Initial T1 [s]:', 'T1 boundary [s]:'...
       'Initial FA [deg]:','FA Boundary [deg]:', ...
       'Add S(0) = 0:', 'Extend TR by copying S(end): ',...
       'Double fit cycle', 'Save data:', 'Display data:'};
opt = {flip(csi.data.labels), num2str(TR, '%2.2f '), ...
      '8.9', '8.8999 8.9001', '20', '1 90', {'Yes', 'No'},...
      {'No', 'Once','Twice', 'Thrice'},...
      {'Yes', 'No'}, {'Yes', 'No'}, {'Yes', 'No'}};
uans = getInput(ele, qry, opt, 'Fit FA-maps');
if isempty(uans), return; end

% UANS - TR-index
ind_TR = csi_findDimLabel(csi.data.labels,uans(1));        

% UANS - TR
TR = str2double(strsplit(uans{2}));
csi.xaxis.TR = TR; setappdata(gui.CSIgui_main, 'csi', csi);

% UANS - T1
T1init = str2double(uans{3});

% UANS - T1 Boundary
T1bound = str2double(strsplit(uans{4}));

% UANS - FA
FAinit = str2double(uans{5});

% UANS - FA Boundary
FAbound = str2double(strsplit(uans{6}));

% UANS - Add Zero
if strcmp(uans{7},'Yes'), doZero = 1; else, doZero = 0; end

% UANS - Extend TR(end) by S(end)
if strcmp(uans{8},'No'), doExtend = 0; doExtendVal = 0; 
else, doExtend = 1;  doExtendVal = numel(uans{8}) - 3; 
end

% UANS - Double Fit Cycle
if strcmp(uans{9},'Yes'), doDouble = 1; else, doDouble = 0; end

% UANS - Save data to file
if strcmp(uans{10},'Yes'), doSave = 1; else, doSave = 0; end

% UANS - Display data
if strcmp(uans{11},'Yes'), doDisp = 1; else, doDisp = 0; end

% Prepare FIT loop % ------------------------------------------- %

% Add zero to TR if enabled
if doZero, TR = [0 TR]; end

% Add TR if doExtend S(end) for TR(end)*N
if doExtend
    TR = [TR repmat(TR(end), 1, doExtendVal)]; 
    for kk = 0:doExtendVal-1, TR(end-kk) = TR(end-kk) * (2.5 - (0.5*kk));
    end
end

% Get k-space index
sz = size(csi.data.raw); % For Kx,Ky,Kz
lab_spat = {'kx','ky','kz', 'x', 'y', 'z'};
ind_spat = csi_findDimLabel(csi.data.labels,lab_spat);
ind_spat(isnan(ind_spat)) = [];        

% Permute: #S x Spatial Dimensions x TR-dimensions x Rest
permv = [1 ind_spat ind_TR];
add_ind = find(~ismember(1:numel(csi.data.dim), permv) == 1);
permv = [permv add_ind];
csi.data.raw = permute(csi.data.raw, permv);

% Create voxel index for loop
vec = arrayfun(@(x) 1:x, sz(ind_spat), 'UniformOutput', false);
voi = allCombinations(vec);
ind_full = arrayfun(@(x) 1:x, size(csi.data.raw), 'UniformOutput', false);
ind_full(ind_spat) = {0};

% Double iteration value
if doDouble, double_iter = 2; else, double_iter = 1; end
fit_par_di = cell(1,double_iter);

% CSI-fitting model: PAR(1) = FA, PAR(2) = T1 and PAR(3) = M0.
csi_model = @(PAR, TR) PAR(3) * (1-exp(-TR./PAR(2))).*sind(PAR(1))./ ...
                                (1-(exp(-TR./PAR(2)).*cosd(PAR(1))));

% Fit twice - using first calculated T1 average and FA for initial fit
% parameters.

for di = 1:double_iter
        
    if di == 1 % Initiate output-var
        ind_spat_sz_c = num2cell(csi.data.dim(ind_spat));
        fit_par = NaN(ind_spat_sz_c{:}, 3);         
        rsq = NaN(1,ind_spat_sz_c{:}); 
        confint = NaN(ind_spat_sz_c{:}, 3, 2);               
        M0bound = [NaN NaN]; opts = cell(ind_spat_sz_c{:}, double_iter);
        output = cell(ind_spat_sz_c{:}, double_iter);  
    end
    
    % @Second iteration: Set new boundaries, get fitted FA/T1.
    if doDouble && di == 2 
        FA = fit_par(:,:,:,1); T1 = fit_par(:,:,:,2); 
        T1init = mean(T1(:)); T1std = std(T1(:)); fac = 0.001;
        T1bound = T1init + [-fac*T1std +fac*T1std];
    end
    
    % Display fit-variable limits   
    CSI_Log({''},{''});
    CSI_Log({'% ----------------------------------- %'},{''});
    CSI_Log({'M0:'},{'Voxel and T1/FA dependent.'});
    CSI_Log({'TR:', 'FA-Initial:','FA-boundary:',...
                    'T1-Initial:', 'T1-boundary:', },...
            {TR, FAinit, FAbound, T1init, T1bound});    
    CSI_Log({sprintf('Fit iteration #%i ------------------- %', di)},{''});
    
    yval = NaN(1, numel(TR));
    for kk = 1:size(voi,2) % Loop voxels

        % Continue-second fit dialog box.
        if di == 2 && kk == 1
            uans = getUserInput_Buttons({'Start second iteration?'},...
                {'Continue', 'Stop'}, [], 'B1-map fitting');
            if isempty(uans), loadBar(NaN); return; end
            switch uans
                case 'Continue' % Continue script...
                case 'Stop', loadBar(NaN); return;
            end
        end
        
        % Progress bar
        loadBar(kk./size(voi,2), 'Fitting voxels...');

        % Use index-lookup to get voxel-data
        voitmp = num2cell(voi(:,kk)); ind_full(ind_spat) = voitmp;

        % Convert data to double
        spoi = squeeze(double(csi.data.raw(ind_full{:})));
        
        %  Peak Value
        ydat = max(real(spoi),[],1);
        
        % lsqcurvefit properties
        opt = optimoptions(@lsqcurvefit,'Display', 'off', ...
            'TolFun', 1e-8, 'TolX', 1e-8, ...
            'MaxFunEvals', 600, 'MaxIter', 600);

        % Add zero to the data
        if doZero, yval(1:end-doExtendVal) = [0 ydat]; 
        else, yval(1:(end-doExtendVal)) = ydat; end

        % Extend by S(end)
        if doExtend 
            yval(end-doExtendVal+1:end) = ...
                repmat(yval(end-doExtendVal), 1, doExtendVal); 
            % opt.Weights = ...
            % [opt.Weights repmat(opt.Weights(end), 1, doExtendVal)];
        end
        
        % Get fitted FA from first iteration if second iteration
        if doDouble && di == 2, FAinit = FA(voitmp{:}); end
        
        % Initial M0
        M0init = max(yval) .* ...
            (1- (exp(-TR(end-doExtendVal)./T1init).*cosd(FAinit))) ./ ...
            (sind(FAinit).*( 1-exp(-TR(end-doExtendVal)./T1init)));
        M0bound(1) = M0init * 0.75; M0bound(2) = M0init * 1.25;        
    
        % Upper and lower bound of fitting parameters.
        % M0 needs to be adjusted according T1/FA init
        lb = [FAbound(1)  T1bound(1)  M0bound(1)];
        ub = [FAbound(2)  T1bound(2)  M0bound(2)];
        starting_values =  [FAinit T1init M0init];
                
        % Fit using nonlinear least square curve fitting fcn        
        [beta, ~, R, ~, outp, ~, J] = ...
        lsqcurvefit(csi_model, starting_values, TR, yval, lb, ub, opt);    
        
        % Store output
        outp.startpoint = starting_values;
        outp.lower = lb; outp.upper = ub;  
        outp.TR = TR; outp.yval = yval;
        outp.beta = beta; outp.model = csi_model;

        % Store fit-options and output
        output{voitmp{:}, di} = outp; opts{voitmp{:}, di} = opt;

        % Store fit-parameters per voxel
        fit_par(voitmp{:}, :) = beta;    
         
        % R-Squared
        yfit  = csi_model(beta,TR);
        rsq(:, voitmp{:}) = ...
            1 - (sum((yval(:)-yfit(:)).^2)/sum((yval-mean(yval)).^2));
                
        % Confidence Interval  
        confint(voitmp{:}, :,:) = nlparci(beta, R, 'jacobian', J);                
    end
    loadBar(NaN);

    % Store fit-par of volume per double-iter
    fit_par_di{di} = fit_par;

    % Save data to file? % -------------------------------------- %
    if doSave
        if isfield(csi.data,'filepath'), fp = csi.data.filepath;
        else, fp = [];
        end
        tstr = sprintf('Save data iteration #%i', di);

        [fn_out, fp_out, idx] = ...
            uiputfile({'*.mat','MATLAB File';}, tstr, fp);
        if idx ~= 0
            fpn = [fp_out fn_out]; 
            fit_par_label = {'FA', 'T1', 'M0'};
            save(fpn, 'fit_par', 'TR', 'rsq', 'confint', 'opts', ...
                'csi_model', 'fit_par_label', 'output');
        end
    end
    
    % Display data? % ------------------------------------------- %
    if doDisp               
        % Fit parameters
        FAt = permute(fit_par(:,:,:,1), [4 1 2 3]);
        CSI_dataAs_Initiate(FAt,  ['FA - ' int2str(di)], gui,...
            csi.data.labels);        
        T1t = permute(fit_par(:,:,:,2), [4 1 2 3]);
        CSI_dataAs_Initiate(T1t,  ['T1 - ' int2str(di)], gui,...
            csi.data.labels);
        M0t = permute(fit_par(:,:,:,3), [4 1 2 3]);
        CSI_dataAs_Initiate(M0t,  ['M0 - ' int2str(di)], gui,...
            csi.data.labels);
        % Goodness of fit
        CSI_dataAs_Initiate(rsq,  ['R2 - ' int2str(di)], gui,...
            csi.data.labels);        
    end

end % End of double fit-iterations

% --- Executes by button_CSI_MapB1
function CSI_MapB1_DoubleCSI(gui)
% Given two CSI measurements with M1(a) and M2(2a), calculate the flipangle
% using: invcos(m2/2*m1) = FA.
%
% Input: gui-handle of CSIgui

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% ------------------------------------------------- %

% -- USERINPUT -- %

% \\ Get peak of interest
range = CSI_getPeakOfInterest(csi.xaxis, 'Calculate B1-maps');
if isempty(range), return; end
doi = CSI_getDataAtPeak(csi.data.raw, csi.xaxis, range);

% \\ Index of interest and data type
uans = getUserInput_Popup({'Data type (calculations): ',...
                           'Data type (results): '},...
                         {{'Real', 'Magnitude','Complex'},...
                          {'Real', 'Magnitude', 'Imaginary',...
                           'Real raw','Imaginary raw', 'Complex'}},...
                          [], 'B1-mapping');
if isempty(uans), CSI_Log({'Skipped B1-mapping.'},{''}) ; return; end

% \\ B1 index
B_ind = find(strcmpi(csi.data.labels,'b1'));
if isempty(B_ind)
    B_ind = getUserInput_Popup({'Data B1 index'},{csi.data.labels},...
                               [], 'B1-mapping');
    if isempty(B_ind)
        CSI_Log({'Skipped B1-map calculations.'},{''}) ; return; 
    end 
    B_ind  = find(strcmp(csi.data.labels,B_ind)==1);
end

% \\ SNR Filtering
filter_snr = getUserInput(...
    {'Minimum value for SNR Filtering: [0 = off]', 'SNR window:'},...
    {'0',round(csi.data.dim(1)./10)});
if isempty(filter_snr)
    CSI_Log({'Skipped B1-map calculations.'},{''}) ; return; 
end

% -- PROCESS USERINPUT -- %

% -- Datatype: DATA CONVERTED TO DATATYPE
switch uans{1}
    case 'Real',      doi = real(doi);
    case 'Magnitude', doi = abs(doi);
end

% -- CALCULATIONS  -- % B1

% \\ Get maximum value of peak of interest
CSI_Log({'Calculating B1 per voxel. Using data-type: '}, uans(1));

% Max value
mx_val = max(doi,[],1);

% \\ Split the data by B_ind
doi_ranges = cellfun(@(x) 1:x, num2cell(csi.data.dim),'uniform', 0);
doi_ranges{1} = 1;
Mone_range = doi_ranges; Mone_range{B_ind} = 1;
Mtwo_range = doi_ranges; Mtwo_range{B_ind} = 2;

Mone = mx_val(Mone_range{:});
Mtwo = mx_val(Mtwo_range{:});

% \\ Calculate using acos(M2/2*M1)
inner = Mtwo ./ (Mone.*2);
FA = acosd(inner);


% -- CALCULATIONS  -- % SNR
if str2double(filter_snr{1}) ~= 0
    
    snr_limit = str2double(filter_snr{1});
    snr_mask = str2double(filter_snr{2});
    
    % Display Info %
    CSI_Log({['Calculating SNR per voxel, '... 
              'and filter data using minimum SNR limit: ']},...
            {snr_limit});

    % Noise mask
    mask_size = snr_mask;

    % Using high-FA data for SNR
    snr_dims = doi_ranges; snr_dims(1) = [];
    snr_dims = [{1:csi.data.dim(1)}, snr_dims]; 
    
    % Use highest flip angle data to calculate SNR
    snr_dims{B_ind} = 2;

    % Calculate using noise mask
    SNR_all = csi_SNR(csi.data.raw(snr_dims{:}), mask_size, 1, range);
    
    % Convert NaNs to zero
    SNR_all(isnan(SNR_all)) = 0; 
    
    % Filter boolean 
    SNR_bool = SNR_all < snr_limit;
    
    % Filter data
    FA(SNR_bool) = NaN;
end


% \\ Process output FA-maps
% FA_main = FA; % Store FA
imag_FA = imag(FA); % Imaginary FA
real_FA = real(FA); % Real FA
abso_FA = abs(FA); % Absolute FA

% When the ratio of M(2a) / 2xM(a) is larger than 1, the resulting inverse
% cosine will  result in imaginary results, which is fals. We filter that
% here:

% Where does the result contain imaginary-output
index_imag_output = (imag_FA ~= 0);
index_real_output = (index_imag_output == 0);

% Real without any voxel containing imaginary results
real_FA_noImag = real_FA;
real_FA_noImag(index_imag_output) = NaN;

% Imaginary without any voxel containing correct results
% So the imaginary part of the data, with correctly calculated values
% removed such that you only see the incorrect results.
imag_FA_noReal = imag_FA;
imag_FA_noReal(index_real_output == 1) = NaN;

% Process result to data-type 
switch uans{2}
    case 'Real',          FA_out = real_FA_noImag;
    case 'Magnitude',     FA_out = abso_FA;
    case 'Imaginary',     FA_out = imag_FA_noReal;
    case 'Real raw',      FA_out = real_FA;
    case 'Imaginary raw', FA_out = imag_FA;
    case 'Complex',       FA_out = FA;
end

% Show statistics nfo
stats = csi_statistics_of_volume(FA_out);
CSI_Log({'B1-map statistics ------------------------------- %',...
         'Mean: ', 'Mode: ', 'Median: ', 'Min | Max: '},...
     {'', sprintf('%.2f +/- %.2f',stats.mean, stats.std), ...
          sprintf('%.2f | freq. %3.0f || ',cat(1,stats.mode, stats.freq)),...
          sprintf('%.2f', stats.median), ...
          sprintf('%.2f | %.2f', stats.min, stats.max)});
% -------------------------------------------- Display
% \\ Display Data
CSI_dataAs_Initiate(FA_out, 'B1-Map', gui, csi.data.labels);

% --- Executes on button press in button_CSI_VoxelMask.
function button_CSI_VoxelMask_Callback(~, ~, gui)
% Create a voxel mask to exclude voxels from calculations. 
%
% A window to select voxels to exclude, with option to copy selection
% to other slices, export indexes of excluded voxels and redo or remove
% mask.
%
% Require: csi-dimensions, image-data.
% Create: figure with slice-dimensions only, ability to select.


% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Get image-data
plot_img = 0;
if isappdata(gui.CSIgui_main,'conv')
    conv = getappdata(gui.CSIgui_main,'conv'); plot_img = 1;
    img = MRI_matchSlices(gui.CSIgui_main); 
    if isfield(conv,'contrast')
        contrast = conv.contrast;    
    else
        contrast = [min(img(:)) max(img(:))];
        if contrast(2) <= contrast(1), contrast(2) = contrast(1)+1; end
    end
end

% Ask user to load or use editor
uans = getUserInput_Popup({'Create voxel-mask:'}, {{'Editor', 'Import'}});
if isempty(uans{1}), return; end

switch uans{1}
    case 'Editor'

    % Set plot parameters % ------------------- %
    % Plot_par requires data-dimensions and colors. Other parameters are
    % calculated by VoxelMask_Editor.
    plot_par = struct; 
    
    plot_par.colors = gui.colors; % GUI Color-nfo: main, text_title, hilight1
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
        CSI_Log({'Mask conversion error.'},{'Mask set to empty.'});        
    end
end

% Clean up
csi.voxelmask = mask;
setappdata(gui.CSIgui_main, 'csi', csi);

% Message user
CSI_Log({'Exclusion voxel mask created.'},...
        {'Applicable to data visualization'});  

% --- Executes on button press in button_CSI_Statistics.
function button_CSI_Statistics_Callback(~, ~, gui)
% Calculate statistics over the volume or current slice. Includes location
% of max and min signal voxel.

% Get appdata
if ~isappdata(gui.CSIgui_main, 'csi'), return; end
csi = getappdata(gui.CSIgui_main, 'csi');

% Ask user data-range of interest
qry = {'Calculate statistics over slice(s):', ...
       'Apply exclusion voxel-mask:', 'Select peak/frequency range:'};
def = {{'All', 'Current','Displayed'},{'Yes','No'},{'No','Yes'}};
uans = getUserInput_Popup(qry, def);
if isempty(uans), return; end

% Analyse user input
if strcmp(uans{2}, 'Yes'), doMask = 1; else, doMask = 0; end
if strcmp(uans{3}, 'Yes'), doPeak = 1; else, doPeak = 0; end

% Set data
doi = csi.data.raw;

% Select peak % ------------------------------------------------------- %
if doPeak
    poi = CSI_getPeakOfInterest(csi.xaxis,'Calculate statistics');
    doi = CSI_getDataAtPeak(doi, csi.xaxis, poi);
    if isempty(poi), return; end
end

% Get data of interest % ----------------------------------------------- %
dim = size(doi);
doi_cell_layout = arrayfun(@(x) 1:x, dim, 'UniformOutput', false);
switch uans{1}
    case 'All', msg = 'Full volume.';
    case 'Current'
        % Get the panel-slider object and 2d-plot object
        [~, gui2D] = CSI_2D_getDataSliders(gui);
    
        % Plotted slice index of interest
        sloi = gui2D.plotindex(1);
    
        % Get data of interest.            
        doi_cell_layout{4} = sloi{:};       

        % Log message        
        msg = sprintf(...
        'Current slice (%i) including all non-spatial dimensions.',sloi{:});

    case 'Displayed'
        % Get the panel-slider object and 2d-plot object
        [~, gui2D] = CSI_2D_getDataSliders(gui);    
        % Plotted slice index of interest
        sloi = gui2D.plotindex;
            
        % Get data of interest.    
        doi_cell_layout(4:4+numel(sloi)-1) = sloi; 

        % Log message
        addStr = repmat(' %i |', 1, numel(sloi));
        msg = sprintf(['Displayed slice only: ' addStr], sloi{:});

end
% Data of interest
doi = doi(doi_cell_layout{:});

% Voxel mask % --------------------------------------------------------- %
if doMask
    % Apply voxel-mask
    if ~isfield(csi, 'voxelmask')
        button_CSI_VoxelMask_Callback([], [], gui)
    end
    csi = getappdata(gui.CSIgui_main, 'csi');
    % User can quit/cancel voxel-mask creation - this catches that error
    % and continues without creating masked-data
    if isfield(csi, 'voxelmask')
        mask = csi.voxelmask;  
    
        % Data-size - is expected to have spatial dimensions on ind(2:4);
        dsz = size(doi); msz = size(mask);     
    
        % Retrieve mask for current size: if non-spatial dimensions are 
        % used to calculate a parameter, that index is not present in the 
        % data-volume.
        if numel(msz) ~= numel(dsz)
            cell_ind = arrayfun(@(x) 1:x, dsz, 'UniformOutput', false);
            % First dimension of mask is only single value
            cell_ind{1} = 1;
            % Set specific slice of interest (see % get DOI % --- %)
            if exist('sloi','var'), cell_ind(4:4+numel(sloi)-1) = sloi; end
            % Cut mask
            mask = mask(cell_ind{:});
        end
    
        % Apply mask to data
        doi(mask) = NaN;
    else
        CSI_Log({'Voxel-Mask was not loaded into memory.'}, ...
                {'Continuing with full data-volume without masking.'})
    end
end



% Calculate statistics % ----------------------------------------------- %
stats = csi_statistics_of_volume(doi);
Statistics_Viewer(stats);

% Log to user
CSI_Log({'Statistics calculated over:',msg}, {' ',' '});

