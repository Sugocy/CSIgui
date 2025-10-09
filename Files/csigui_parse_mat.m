function [csi, mri, conv, log, gui, success] = csigui_parse_mat(fp, fn, gui)
% Read mat-file and apply checksums for CSIgui data-handling. Returns the
% loaded CSI structure. Returns all appdata structures related to CSIgui
% present in the mat-file.

% --- Check mat-file integrety for CSIgui  
mat_cont = whos('-file',[fp '\' fn '.mat']); format_old = 0; format_new = 0;
if strcmp(mat_cont(1).name,'csigui'),  format_old = 1;
elseif strcmp(mat_cont(1).name,'csi'), format_new = 1;
end

success = 1; csi = nan; mri = csi; conv = csi; log = nan;
if (format_old + format_new) == 0
    success = 0; 
    CSI_Log(...
    {'Incorrect mat-file. Use by CSIgui generated mat-file.',...
     'Expected fields:'},...
    {'Required structure: csigui (old) or csi (new).',...
     'raw, dim, filepath and name, noise, split, conv, mri,'}); 
    return;
end

% Integrety verified, read matfile
if format_old
    inp = load([fp '\' fn '.mat'], 'csigui');
    csigui = inp.csigui; clear inp;
elseif format_new
    inp = load([fp '\' fn '.mat'], 'csi');
    csi = inp.csi; clear inp;
end

             
if format_old            
    
    % Create csi-structure.
    csi = struct;

    % Find MRI struct and store   
    if isfield(csigui,'mri')
        mri = csigui.mri; csigui = rmfield(csigui,'mri');
    end

    % Find conv struct and store
    if isfield(csigui,'conv')
        conv = csigui.conv; csigui = rmfield(csigui,'conv');
    end

    % Twix (Siemens) header raw data
    if isfield(csigui,'twix')
        csi.twix = csigui.twix; csigui = rmfield(csigui,'twix');
    end

    % List (Philips) header raw data 
    if isfield(csigui,'list')
        csi.list = csigui.list; csigui = rmfield(csigui,'list');
    end
    
    % Exclusion voxel-mask.
    if isfield(csigui,'voxelmask')
        csi.data.voxelmask = csigui.voxelmask;
    end   

    % Set remaining CSI input to structure;
    csi.data = csigui; csi.ext = 'mat';    
    
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
                      
    % Get log data and parse
    if isfield(csigui,'log')
        log = cellstr(csigui.log);        
        log = log(3:end); % Exclude first two lines          
        csi.data = rmfield(csi.data,'log'); % Remove log-field
    end

    % Clear memory
    clear csigui
    
elseif format_new
        
    if isfield(csi, 'conv')
        conv =  csi.conv; csi = rmfield(csi, 'conv');
    end
   
    if isfield(csi, 'mri')
        mri = csi.mri;  csi = rmfield(csi, 'mri');
    end

    if isfield(csi, 'log')
        log = cellstr(csi.log);                
        log = log(3:end); % Exclude first two lines        
        csi = rmfield(csi,'log'); % Remove log-field
    end

    if isfield(csi, 'voxelmask')
        csi.data.voxelmask = csi.voxelmask; 
        csi = rmfield(csi, 'voxelmask');
    end

    % Set extensions
    csi.ext = 'mat';    
end