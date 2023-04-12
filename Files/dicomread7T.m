function [images, sliceinfo, examinfo] = dicomread7T(varargin)
%%%% Description:                     Read Philips DICOM files from MRI
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2015-09
%%% --------------------------------------------------------------------
%%%
%%% [img, slice info, examinfo] = dicomread7T({....});
%%% Open a DICOM MRI file and load the images and image information. All
%%% images are stored in an array. The information of all N slice is stored 
%%% as a struct per slice in a cell-array of size 1xN while general exam
%%% information is returned as a single structure.
%%%
%%% Default output mode is set to floating point values and scaled using 
%%% the MR rescale interscept and slope parameters in the dicom-header.*
%%%
%%% Input: 
%%%         None     (Open via GUI)
%%%         Filepath (string)
%%%         Filepath (string) + image display (1 or 0)        
%%%         Image display (1 or 0) + output mode (1-4)
%%%
%%% *Source: Quantitative Image Analysis due to Platform-Dependent
%%% Image Scaling. - Chenevert 2014 
%%%
%%%                                             % WRITTEN FOR 7T
%%%                                             % Compatible with 3T
fprintf('\n');

switch size(varargin{:},2) 
    case 0, [filename, filepath, fi] = ...
            uigetfile('*.DCM', 'Open DICOM-file'); if fi == 0,return; end; 
            filepath = [filepath '\' filename]; imgdisp = 0; outp = 1;
    case 1, filepath = varargin{1}{1}; imgdisp = 0; outp = 1;
    case 2
        if ischar(varargin{1}{1})                                           % If first input argument is filepath string 
            filepath = varargin{1}{1}; imgdisp = varargin{1}{2}; outp = 1;
        else                                                                % Else open via UI
            imgdisp = varargin{1}{1};  outp = varargin{1}{2};
            [filename, filepath, fi] = ...
            uigetfile('*.DCM', 'Open DICOM-file');    
            if fi == 0, return; end; filepath = [filepath '\' filename];
        end
    case 3, filepath = varargin{1}{1}; imgdisp = varargin{1}{2};
            outp = varargin{1}{3};
    otherwise
        disp('To many input arguments...'); return; 
end


[fp, fn, ext] = fileparts(filepath);
if isempty(ext), ext ='.dcm'; end; filepath = [fp '\' fn ext];


%% Read DICOM-header & get image(s) information
dcmdict = 'dicom-dict-philips.txt';
if exist(dcmdict,'file'), dicomdict('set', dcmdict);
else, fprintf('No dicom-dictionary found, default by Matlab used.\n');
end
% Load dicom header/info
if exist(filepath, 'file')
    dcminf = dicominfo(filepath);
else
    images = NaN; sliceinfo = NaN; examinfo = NaN; return; 
end


% Check if MR image(s)
if ~strcmp(dcminf.Modality, 'MR')
    fprintf('Error: Image modality not equal to MRI\n'); return; 
end

% Extract nr of image(s) if not possible, abort.
if isfield(dcminf, 'NumberOfFrames'), nrimg = dcminf.NumberOfFrames; 
else, disp('Error: No images found!'); images = 0; sliceinfo = 0; examinfo=0; 
    return;
end

% Load all files
imgdcm = dicomread(filepath); 

%% Frame specific floating point calculation/frame
dcmsz  = size(imgdcm); slinfo = cell(1, nrimg); 
if imgdisp == 1, figh = figure(); end

switch outp
    case 4,    imgout = NaN(dcmsz(1), dcmsz(2), size(imgdcm,4),3);
    case 5,    imgout = NaN(dcmsz(1), dcmsz(2), size(imgdcm,4),4);
    otherwise, imgout = NaN(dcmsz(1), dcmsz(2), size(imgdcm,4));      
end

for qq = 1:nrimg
    % Rescale intercep 2005,100d - Sclope 2005, 100e
    % Get normalization factors: RI and RS
    
    % Scale intercept - RI
    RI = double(dcminf.PerFrameFunctionalGroupsSequence...
        .(['Item_',num2str(qq)])...
        .PrivatePerFrameSq.Item_1.RescaleIntercept);
    % Scale slope - RS
    RS = double(dcminf.PerFrameFunctionalGroupsSequence...
        .(['Item_',num2str(qq)])...
        .PrivatePerFrameSq.Item_1.RescaleSlope);
    % Scale Type - RT
    RT = (dcminf.PerFrameFunctionalGroupsSequence...
        .(['Item_',num2str(qq)])...
        .PrivatePerFrameSq.Item_1.RescaleType);
    
    % Get normalization factors: SI and SS
    % Other fields with scaling factors are not equal!
    SS = double(dcminf.PerFrameFunctionalGroupsSequence...
        .(['Item_',num2str(qq)])...
        .PrivatePerFrameSq.Item_1.MRScaleSlope);
    SI = double(dcminf.PerFrameFunctionalGroupsSequence...
        .(['Item_',num2str(qq)])...
        .PrivatePerFrameSq.Item_1.MRScaleIntercept);
    
    % Take image and covert to double
    PV = double(imgdcm(:,:,1,qq));      % Raw data from DICOM.     
    
    % =================================================================== %
    % BE WARNED: Vendor specific scaling is applicable to certain
    % calculated maps suchs as ADC/B1/T1 etc.
    % Example: ADC maps on philips MR systems require use of the DV-values
    % to obtain real world values. Evenso, the values are still off by a
    % factor of 1000.
    % =================================================================== %
    
    DV  = (PV.*RS+RI);                   % Display values
    FP  = (PV.*RS+RI)./(RS*SS);          % Floating point values
    FP2 = (DV .* SS + SI);               % FP values: both formulas hold!
    
    % Save image(s) according to outp type.
    switch outp
        case 1 % PV to FP - scaling 1 - Rescale to floating point values
            imgout(:,:,qq) = FP;
       
        case 2 % PV to DV - Scaling 2 - Rescale to display values
            imgout(:,:,qq) = DV;
            
        case 3 % RAW data - No scaling applied.
            imgout(:,:,qq) = PV;
           
        case 4 % All PV-scaling data requested - see above.
            imgout(:,:,qq,1) = FP;  % 1  = fp
            imgout(:,:,qq,2) = DV;  % 2  = dv            
            imgout(:,:,qq,3) = PV;  % 3  = pv, rawdata
                            
        case 5
            % THIS OUTP IS FOR TESTING PURPOSES ONLY
            % ====================================== %
            imgout(:,:,qq,1) = PV;  
            imgout(:,:,qq,2) = DV;            
            imgout(:,:,qq,3) = FP;  
            imgout(:,:,qq,4) = FP2;
            % 3 and 4 should be equal outputs!
             
    end            
    

    
    % Get slice specific information
    sl = dcminf.PerFrameFunctionalGroupsSequence.(['Item_',num2str(qq)]); 
    
    % Define TR per slice.
    sl.TR = dcminf.MRSeriesRepetitionTime;
    
    % Add scaling parameters.
    sl.Scaling.rs_ri = [RS ; RI]; sl.Scaling.ss_si = [SS;  SI];
    sl.Scaling.rtype = RT;
    
    
    % Determine orientation from dicom header and store per frame
    % Stack, volume and orientation information!
    if isfield(dcminf, {'Stack'})
        % Save stack structure
        sl.Stack = dcminf.Stack;
        % Get orientation from stack.
        if isfield(dcminf.Stack, 'Item_1')
              sl.MRStackViewAxis = dcminf.Stack.Item_1.MRStackViewAxis;
              ViewAxis = dcminf.Stack.Item_1.MRStackViewAxis;               % View axis e.g. slice-direction
        else,  ViewAxis = NaN; sl.MRStackViewAxis = 'NaN';
        end
        
        % Set orientation according to slice-dimension
        switch ViewAxis
            case 'FH', sl.PlaneOrientation = 3; % transversal
                sl.PlaneOrientationStr = 'TRA';    
            case 'RL', sl.PlaneOrientation = 2; % sagittal
                sl.PlaneOrientationStr = 'SAG';       
            case 'AP', sl.PlaneOrientation = 1; % coronal
                sl.PlaneOrientationStr = 'COR'; 
            otherwise, sl.PlaneOrientation = 0; % unknown 
                sl.PlaneOrientationStr = NaN;
        end
    else
        fprintf('Plane orientation not found! Set to NaN.'); 
        sl.MRStackViewAxis  = NaN;
    end


    
    % Add MRFOVGeo field
    sl.PrivatePerFrameSq.Item_1.PatientPosition  = dcminf.PatientPosition;
    
    % Save slice information!
    slinfo{qq} = sl;
    
    % Display info
    msg = sprintf('%d/%d',qq,nrimg); msgsz = size(msg);
    if qq == 1,    fprintf([ '\n' msg] );
    else           fprintf([repmat('\b',1,msgsz(2)) msg]);
    end
    
    % Display image if activated
    if imgdisp == 1, imshow(PV(:,:), [min(PV(:)) max(PV(:))]); 
        drawnow expose; pause(0.1);
    end
end
if imgdisp == 1, close(figh); end

%% Set image-output
images = imgout; sliceinfo = slinfo;

%% Generate and save general data/info
examinfo = struct;

% General exam info
examinfo.patientname      = dcminf.PatientID;
examinfo.examname         = dcminf.PatientName.FamilyName;
examinfo.protname         = dcminf.ProtocolName;
examinfo.examdatetime     = dcminf.AcquisitionDateTime;
examinfo.scanmode         = dcminf.PulseSequenceName;
examinfo.scandur          = dcminf.MRSeriesScanDuration;
examinfo.patientposition  = dcminf.PatientPosition;

%% Anatomical plane orientation to general exam data.
if sum(isfield(dcminf, {'SeriesVolume', 'Stack'})) > 0
    
    % Get plane orientation from SeriesVolume or Stack struct.
    if isfield(dcminf, 'SeriesVolume') && ...
            isfield(dcminf.SeriesVolume, 'Item_1')
        ViewAxis = dcminf.SeriesVolume.Item_1.MRVolumeViewAxis;             % For now, only 1 stack operations.
    elseif isfield(dcminf, 'Stack') && isfield(dcminf.Stack, 'Item_1')
        ViewAxis = dcminf.Stack.Item_1.MRStackViewAxis;                     % View axis e.g. slice-direction
    else
        ViewAxis = NaN;
    end
    
    % Set orientation according to slice-dimension
    switch ViewAxis
        case 'FH', examinfo.PlaneOrientation = 3; % transversal
                   examinfo.PlaneOrientationStr = 'TRA';           
        case 'RL', examinfo.PlaneOrientation = 2; % sagittal
                   examinfo.PlaneOrientationStr = 'SAG';          
        case 'AP', examinfo.PlaneOrientation = 1; % coronal
                   examinfo.PlaneOrientationStr = 'COR'; 
        otherwise, examinfo.PlaneOrientation = 0; % unknown 
                   examinfo.PlaneOrientationStr = NaN;
    end 
else
    fprintf('Plane orientation not found! Set to NaN.\n'); 
    examinfo.PlaneOrientationStr = NaN;
end

%% Miscellaneous general parameters
examinfo.acqnr       = dcminf.AcquisitionNumber;
examinfo.maxnrslice  = dcminf.MRSeriesNrOfSlices;
examinfo.maxnrdyn    = dcminf.MRSeriesNrOfDynamicScans;
examinfo.maxnrechoes = dcminf.MRSeriesNrOfEchoes;
examinfo.dimensions  = [dcminf.Columns dcminf.Rows];
examinfo.flipangle   = dcminf.MRSeriesFlipAngle;
examinfo.NSA         = dcminf.MRSeriesNumberOfAverages;
examinfo.TR          = dcminf.MRSeriesRepetitionTime;
examinfo.pxBW        = dcminf.PixelBandwidth;
examinfo.WFshift     = dcminf.MRSeriesWaterFatShift;
examinfo.Stacks      = dcminf.MRSeriesNrOfStacks;
examinfo.StacksInfo  = dcminf.Stack;                                        % All stacks: item_1, item_2 etc.
examinfo.DCMdump     = dcminf;                                              % The entire DICOM head is dumped in here - just in case thingy;)                              
                                                                            % Missing max nr of cardiac cycles.

fprintf('\n')