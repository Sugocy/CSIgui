function [img, nfo] = dicomreadSiemens(varargin)
% Read dicom files from Siemens platforms, ima-files.
%
% Input: none - opens a select-file(s) user interface.
%        path - a full path to file
%        filepath, filename - full path to directory of file(s) and
%                             corresponding filename(s). Filename must be a 
%                             cell arry with 1xN filenames. The image size
%                             for each filename is expected equal.
%                             
%
% Creator: Quincy van Houtum, PhD. 2023
% quincyvanhoutum@gmail.com


% Parse input-arguments.
if nargin == 0
    [fn, fp, id] = uigetfile({'*.ima', 'Siemens DICOM files (*.ima)'},...
                    'Select one or multiple dicom files.',...
                    'MultiSelect', 'on');
    if id == 0, return; end    
    fn = {fn};
elseif nargin == 1 % Only file-path to file
    [fp, fn, ext] = fileparts(varargin{1}); 
    fn = {[fn ext]};    
elseif nargin == 2 % Directory path and filename.
    fp = varargin{1}; fn = varargin{2};
    if ~iscell(fn), fn = {fn}; end        
end

% Correct filepath
if ~strcmp(fp(end),'\'), fp = [fp '\']; end

% Set dicom-dictionary
setDicomDict(fp, fn{1});

% Check mr-data type
[type, sop] = checkDicomModality(fp, fn{1});
if strcmpi(type,'mrs') || strcmpi(type,'nan')
    warning('dicomreadSiemens | Incompatible dicom-type (%s) file: %s',...
        sop, fn{1});
    img = nan; nfo = {nan}; return;
end

% Image size and storage
matrix_size = getMatrixSize(fp,fn{1});
if isnan(matrix_size)
    warning('dicomreadSiemens | No col/row nfo in dicom header: \n%s',...
    fn{1});
    img = nan; nfo = {nan}; return; 
end
img = NaN(matrix_size(2),matrix_size(1),size(fn,2));

% Loop files; we assume all are MR-image dicoms.
nfo = cell(1,size(fn,2));
for fni = 1:size(fn,2)
    fnt = fn{fni}; % Temp filename of iteration.
    nfo{fni} = dicominfo([fp fnt],'UseDictionaryVR', true);

    % Read dicom file
    img(:,:,fni) = single(dicomread([fp fnt]));
    
end % End of file(s) for-loop.


end

function [type, sopclass] = checkDicomModality(fp, fn)
% Returns MR data type: mri, mri enhanced dicom, mrs or nan for unknown
% SOPclassUID.
    nfo = dicominfo([fp fn],'UseDictionaryVR', true); % Load dicom-nfo
    sopclass = nfo.MediaStorageSOPClassUID; 
    
    if strcmpi(sopclass,'1.2.840.10008.5.1.4.1.1.4')   % MRI        
        type = 'mri';
    elseif strcmpi(sopclass,'1.2.840.10008.5.1.4.1.1.4.1') % MRI enhanced
        type = 'mri_enhanced';
    elseif strcmpi(sopclass,'1.2.840.10008.5.1.4.1.1.4.2') % MRS   
        type = 'mri';
    elseif strcmpi(sopclass,'1.3.12.2.1107.5.9.1 ')
        type = 'CSA Non-Image';
    else
        type = 'nan';
    end
    
end

function setDicomDict(fp, fn)
    % Set correct dicom-dictionary        
    % dicomdict('factory')
    nfo_tmp = dicominfo([fp fn],'UseDictionaryVR', true); % Load dicom-nfo
    [~,~,ext] = fileparts([fp fn]);
    if strcmpi(nfo_tmp.Manufacturer,'Siemens') || strcmpi(ext,'.ima')
        dicomdict('set', 'dicom-dict-siemens.txt');
    elseif strcmpi(nfo_tmp.Manufacturer,'Philips')
        dicomdict('set', 'dicom-dict-philips.txt');
    end
    clear('nfo_tmp');    
end

function matrix_size = getMatrixSize(fp,fn)
    % dicomdict('factory'); % Revert to factory settings
    nfo_tmp = dicominfo([fp fn],'UseDictionaryVR', true); % Load dicom-nfo
    if isfield(nfo_tmp, 'Columns') && isfield(nfo_tmp, 'Rows')
        matrix_size = double([nfo_tmp.Columns nfo_tmp.Rows]);
    else
        matrix_size = NaN;
    end
    clear('nfo_tmp'); 
end