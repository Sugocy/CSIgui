function examinfo = dicomreadSiemens_sortParameters(nfo)
% Given a cell array with the nfo per slice of an image array, retuns a
% sorted struct 'examinfo' with specific important parameters.
%
% The resulting examinfo structure is comparable to dicomread7T for Philips
% dicom files.
%
% Fields which are Philips dicom dependent will be set to NaN.


% Go...
examinfo = struct;

% General information
examinfo.patientname      = nfo{1}.PatientID;
examinfo.examname         = nfo{1}.PatientName.FamilyName;
examinfo.protname         = nfo{1}.ProtocolName;
examinfo.examdatetime     = nfo{1}.AcquisitionDate;
examinfo.scanmode         = nfo{1}.SequenceName;
examinfo.scandur          = nfo{1}.SeriesTime;
examinfo.patientposition  = nfo{1}.PatientPosition;


% Check for image-orientation: sag/cor/tra
if isfield(nfo{1}, 'ImageOrientationPatient')
    % Describes orientation of image to patient.
    ori = extractField(nfo,'ImageOrientationPatient');
    ori = cellfun(@iop2plane, ori, 'uniform', false);
    % Compare for all slices.
    ind = cellfun(@strcmpi,ori,repmat(ori(1),size(ori)));
    if sum(ind) == size(ori,2)
        ori = ori(1);
    else
        ori = join(ori,'\');
    end
end
examinfo.planeorientationstr = ori{:};


switch examinfo.planeorientationstr
        case 'tra', examinfo.planeorientation = 3; % transversal        
        case 'sag', examinfo.planeorientation = 2; % sagittal
        case 'cor', examinfo.planeorientation = 1; % coronal
        otherwise,  examinfo.planeorientation = 0; % unknown or multiple. 
end


examinfo.acqnr       = nfo{1}.AcquisitionNumber;
examinfo.maxnrslice  = size(nfo,2);
examinfo.maxnrdyn    = NaN;
examinfo.maxnrechoes = nfo{1}.EchoTrainLength;
examinfo.dimensions  = double([nfo{1}.Columns nfo{1}.Rows]);
examinfo.flipangle   = nfo{1}.FlipAngle;
examinfo.NSA         = nfo{1}.NumberOfAverages;
examinfo.TR          = nfo{1}.RepetitionTime;
examinfo.TE          = nfo{1}.EchoTime;
examinfo.pxBW        = nfo{1}.PixelBandwidth;
examinfo.WFshift     = NaN;
examinfo.Stacks      = NaN;
examinfo.StacksInfo  = NaN;
examinfo.DCMdump     = nfo; 

end


function view = iop2plane(iop)
% Using the crossproduct of the iop, determine the plane of view of data.
% ['1', '0', '0', '0', '0', '-1'] Cor
% ['0', '1', '0', '0', '0', '-1'] Sag
% ['1', '0', '0', '0', '1',  '0'] Tra    
    iopr = round(iop);
    plane = abs(cross(iopr(1:3),iopr(4:6)));
    ind = find(plane(1:3) == 1);
    if     ind == 1, view = 'sag';
    elseif ind == 2, view = 'cor';
    elseif ind == 3, view = 'tra';
    end

end

