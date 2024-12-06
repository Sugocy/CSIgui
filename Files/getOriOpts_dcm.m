function opts = getOriOpts_dcm(im, nf)
% Input: {Image array , info cell-structs}
% Create opts. Struct with all geometric information of the DICOM file.
% 
% NB. Image array indexing requires order: height, width, slice.
%     Data will be squeezed.
% 
% 1. Image volume details 
% 2. Orientation 
% 3. Limits and Ranges
% 4. Transformation matrices

% Output structure and the images
opts = struct; opts.img = squeeze(im);

opts = getImageDetails(opts,nf);   % Dim, resolution, position, orientation
opts = getStackDetails(opts,nf);   % Stack- FOV, offcenter. 
opts = getCoordinateGrid(opts,nf); % Image and Ori pos, limits, meshgrid.
opts = getCRP2NWV(opts);           % Matlab correction row=y, col=x;
opts = getMagnetToPatient(opts);   % Tpom - See philips documentation
opts = getScaling(opts);           % Transformation matrix scaling
opts = getTranslation(opts);       % Transformation matrix translation
opts = getRefObject(opts);         % RefObject of transformation matrices
opts = getAngulation(opts);        % Angulation transformation matrix

% Get Tsom matrix: Slice Orientation Matrix
Tsom = getTsom(opts.slice_plane); opts.array.Tsom = Tsom;


function opts = getImageDetails(opts, nf)
% Image volume details: height, width, slice, dimensionality, voxel size.
%                       plane orienation, 
opts.row    = size(opts.img,1);   % Height/Width/St/Sl/Ech/Dyn/Misc.
opts.col    = size(opts.img,2);   % Height/Width/St/Sl/Ech/Dyn/Misc.
opts.slices = size(opts.img,3);   % Nr of slices.
opts.imSz   = size(opts.img)  ;

% 3D or 2D 
opts.dim    = nf{1}.PrivatePerFrameSq.Item_1.MRAcquisitionType;

% Patient orienation: FF/HF and Supine/Prone e.g. HFP, FFS, etc.
opts.pos    = nf{1}.PrivatePerFrameSq.Item_1.PatientPosition;

% Image plane orientation
opts.slice_dir = nf{1}.Stack.Item_1.MRStackViewAxis;
switch opts.slice_dir
    case 'RL', opts.slice_plane = 'SAG';  
    case 'FH', opts.slice_plane = 'TRA';
    case 'AP', opts.slice_plane = 'COR';
end

% Voxel sizes
opts.vox(1:2) = nf{1}.PixelMeasuresSequence.Item_1.PixelSpacing;
opts.vox(3)   = nf{1}.PrivatePerFrameSq.Item_1.SliceThickness;
opts.gap      = nf{1}.PrivatePerFrameSq.Item_1.SpacingBetweenSlices;

function opts = getStackDetails(opts,nf)
try
% FOV from stack parameter.
opts.stack.fov_dim  = {'RL', 'AP','FH'};
opts.stack.fov(1)   = double(nf{1}.Stack.Item_1.MRStackFovRL);
opts.stack.fov(2)   = double(nf{1}.Stack.Item_1.MRStackFovAP);
opts.stack.fov(3)   = double(nf{1}.Stack.Item_1.MRStackFovFH);

% Stack offcenter - Select knowing LR/AP/FH! 
opts.stack.offcenter_dim = {'RL', 'AP','FH'};
offc_tmp(1,:) = cell2mat(extractField(nf,...
                        'Stack.Item_1.MRStackOffcentreRL'))';
offc_tmp(2,:) = cell2mat(extractField(nf,...
                        'Stack.Item_1.MRStackOffcentreAP'))';
offc_tmp(3,:) = cell2mat(extractField(nf,...
                        'Stack.Item_1.MRStackOffcentreFH'))';
opts.stack.offcenter = double(unique(offc_tmp','rows')); clear('offc_tmp');
catch err
    err.message
    fprintf('Did not succeed in getting opts.stack.offcenter.')
    return;
end

function opts = getCoordinateGrid(opts,nf)
% Create meshgrid for voxels in image volume by using coordinate limits,
% range and voxelsizes.

% 1. Image orientation Patient and Position.

% Image Orientation Patient
opts.iop = cell2mat(extractField(nf,...
    'PlaneOrientationSequence.Item_1.ImageOrientationPatient'))';

% Image Patient Position
opts.ipp = cell2mat(extractField(nf,...
    'PlanePositionSequence.Item_1.ImagePositionPatient'))';

%%% NB %%%
% If you have TRA recored protocols the following is true:
% 1. IPP(:,3)   == z-coordinate of the slices top left pixels.
% 2. IPP(:,1:2) == equal for all slices (X/Y doesnt change)


%% Inital limits and span

% Get limits of row/col/slice coordinates
opts.lim.row    = [opts.ipp(1,1) ...
                   opts.ipp(1,1) + (opts.row-1).*opts.vox(1) ];   % row
opts.lim.col    = [opts.ipp(1,2) ...
                   opts.ipp(1,2) + (opts.col-1).*opts.vox(2) ];   % col
opts.lim.slice  = [opts.ipp(1,3) opts.ipp(end,3)];

               
% NOW things get interesting as vendors do their own thing in a 
% STANDARDIZED file-type.
% 1. https://itk.org/pipermail/insight-users/2008-November/027903.html
% 2. https://groups.google.com/forum/#!msg/comp.protocols.dicom...
%          /gyvsgAj4y6o/KOjJywTWFlAJ
%
% To conclude:
% IOP and IPP are safer for use of DICOM slice reconstruction than the
% Spacing between Slices or SliceThickness fields in the dcm-header.
%
% Use IPP(:,3)                      : coordinates/span of the slices
% Use SpacingBetweenSlices (Gap)    : measure for center 2 center of slices
% 
               
% Get span vectors: no angulation taken into account.
opts.range.row   = opts.lim.row(1)  :opts.vox(1): opts.lim.row(2);
opts.range.col   = opts.lim.col(1)  :opts.vox(2): opts.lim.col(2);

% This is dangerous!  BE AWARE NOT TO USE THIS FORMULA/METHOD
% opts.span.slice = opts.lim.slice(1):opts.gap : opts.lim.slice(2)+opts.gap;
%
% So why is this? Floating point accuracy, gap/spacingbetweenslices and
% slicethickness which do not cooperate all the time. Problems with span
% going from lim(1) to lim(2) but excludes the lim(2) due floating point
% errors. Correcting by +1-gap works only for 3D. 
%
% THUS solution:            Use IPP(:,3) as slice span.
opts.range.slice = opts.ipp(:,3)';



% Create meshgrid

opts.mesh.nb = 'Meshgrid: row and column are NOT swapped.';
[opts.mesh.x, opts.mesh.y,  opts.mesh.z] = ...
    meshgrid(opts.range.col, opts.range.row, opts.range.slice);

function opts = getCRP2NWV(opts)
% Transformation matrices
% Tmat  :   Matlab to MRI transformation
% Action:   Transform Col Row Page (CRP) orientation to North West
%           View (NWV). First two index dimensions are switched.
%           Already in a 3D 4x4 format.
opts.array.Tmat_nb = 'Tmat: Transform CRP to NWV orientation.';
opts.array.Tmat    = [ 0 -1  0  0 ; 
                      -1  0  0  0 ; 
                       0  0  1  0 ; 
                       0  0  0  1 ]; 

function opts = getMagnetToPatient(opts)
% Input:       opts.pos e.g. string HFS/FFS/HFP/FFP
% Description: Create tform for geometric operations: Tpp and Tpo       
% Output:      Tpom in opts. ... array. 4x4 matrix.


% How a point is given in magnet coordinates is given by PatientOrientation
% and position: HF/FF - S/P.
% Tpo = orientation invariant axis permutation
% Tpp = orientation changing axis permutation going from right to left
%       handed system.
switch opts.pos
    case 'HFS', Tpp = [0 1 0 0 ; -1 0 0 0; 0 0 -1 0; 0 0 0 1];
                Tpo = [1 0 0 0 ;  0 1 0 0; 0 0  1 0; 0 0 0 1]; 
    case 'FFS', Tpp = [0 -1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1];
                Tpo = [1 0 0 0 ;  0 1 0 0; 0 0  1 0; 0 0 0 1];
    case 'HFP', Tpp = [0 1 0 0 ; -1 0 0 0; 0 0 -1 0; 0 0 0 1];
                Tpo = [-1 0 0 0; 0 -1 0 0; 0 0 1 0 ; 0 0 0 1];
    case 'FFP', Tpp = [0 -1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1];
                Tpo = [-1 0 0 0; 0 -1 0 0; 0 0 1 0 ; 0 0 0 1];
    otherwise 
        fprintf('Orientation patient unknown: HFS/FFS/HFP/FFP?\n');
end

try
    opts.array.Tpom_nb = ...
        'Tpom: Transform magnet point to patient coordinates.';
    opts.array.Tpom = Tpo * Tpp; 
catch 
    opts.array.Tpom_nb = ...
        'Tpom: Transform magnet point to patient coordinates.';
    fprintf('Tpom not found!');
end

function opts = getScaling(opts)
%% Scaling
% Input:       opts.vox(1:3) e.g. resolution.
% Description: Create tform for geometric operations: SCALING        
% NB:          Switch of vox(1) and vox(2) -> e.g. Col/Row.
% Output:      scaling in opts. ... tform&array. 4x4 matrices.

opts.array.scaling = [opts.vox(2) 0 0 0; 0 opts.vox(1) 0 0;...
                      0 0 opts.vox(3) 0; 0 0 0 1];
opts.tform.scaling = affine3d(opts.array.scaling); 

function opts = getTranslation(opts)
%% Translation
% Input:       opts.ipp(1:3) e.g. image position patient/slice.
% Description: Create tform for geometric operations: TRANSLATION       
% NB:          Switch of ipp(1) and ipp(2) -> e.g. Col/Row.
% Output:      translation in opts. ... array&tform. 4x4 matrices.
%              Both as cell (1xNslices)

for kk = 1:size(opts.ipp,1)
    ipp = opts.ipp(kk,:);
    opts.array.translation{kk} = [1 0 0 ipp(2); 0 1 0 ipp(1);...
                                  0 0 1 ipp(3); 0 0 0 	1];
    opts.tform.translation{kk} = affine3d(opts.array.translation{kk}');                      
end
                        
function opts = getRefObject(opts)
%% Create reference objects:             ROW AND COL ARE SWAPPED
% Swap of image dimensions NOT image size as proposed by Matlab.

% An imref3d object encapsulates the relationship between the intrinsic 
% coordinates anchored to the columns, rows, and planes of a 3-D image and
% the spatial location of the same column, row, and plane locations in a 
% world coordinate system.

% The intrinsic coordinate values (x,y,z) of the center point of any pixel 
% are identical to the values of the column, row, and plane subscripts for
% that pixel. For example, the center point of the pixel in
% row 5, column 3, plane 4 has intrinsic coordinates x = 3.0, y = 5.0, 
% z = 4.0. Be aware, however, that the order of the coordinate 
% specification (3.0,5.0,4.0) is reversed in intrinsic coordinates 
% relative to pixel subscripts (5,3,4).

opts.refobj.nb = ...
'Imref3D: Swap of voxel sizes NOT the image size as proposed by Matlab.';
opts.refobj.vox = imref3d([opts.row opts.col opts.slices],...
                           opts.vox(2),    opts.vox(1),    opts.vox(3));
          
opts.refobj.size=  imref3d([opts.row opts.col opts.slices]);
try 
    opts.refobj.lim = imref3d([opts.row opts.col opts.slices],...          
              opts.lim.col, opts.lim.row, opts.lim.slice);
catch
    opts.refobj.lim = imref2d([opts.row opts.col opts.slices],...          
              opts.lim.col, opts.lim.row);
end
    
function opts = getAngulation(opts)

%% Angulation

% Source: http://dicomiseasy.blogspot.nl/2013/06/...
%                           getting-oriented-using-image-plane.html          
% The First vector is (0.5,0,-0.8660254) is the direction of the image 
% rows in the patient coordinate system. 0.5 is cos(600) and -0.8660254 
% is cos(1500). So the image rows are rotated 60  from the patient's X 
% direction (right-to-left) and 150 from the 
% patient's Z direction (feet-to-head) and the patient's Y direction is
% perpendicular to the image X axis.

% Source: http://cmic.cs.ucl.ac.uk/fileadmin/cmic/Documents/
%                                       DavidAtkinson/DICOM_6up.pdf


% Angulation in degrees.
tmp_ang          = acosd(opts.iop(1,:))'; 

% Add to options.
opts.ang.info = 'No swap of row&col.'; 
opts.ang.x       = tmp_ang(1:3); opts.ang.y = tmp_ang(4:6); 
opts.ang.dircos  = opts.iop(1,:)'; 

% Direction Cosine                    
% Vector n: cross of F columns: iop.                    Row/Col swapped.
% Source:
% http://nipy.org/nibabel/dicom/dicom_orientation.html#dicom-orientation

opts.ang.F_dim   = 'F-direction cosine matrix - row&col are swapped.';
opts.ang.F       = [opts.ang.dircos(4) opts.ang.dircos(1);...
                    opts.ang.dircos(5) opts.ang.dircos(2);...
                    opts.ang.dircos(6) opts.ang.dircos(3)];
% Unit vector of dir cosines
opts.ang.Fcross         = cross(opts.ang.F(:,1), opts.ang.F(:,2));

% Transformation matrix Trot, rotation/angulation
opts.array.ang  = [opts.ang.F(1,1) opts.ang.F(1,2) opts.ang.Fcross(1) 0;...
                   opts.ang.F(2,1) opts.ang.F(2,2) opts.ang.Fcross(2) 0;...
                   opts.ang.F(3,1) opts.ang.F(3,2) opts.ang.Fcross(3) 0;...
                       0                  0                  0          1];  
try                         
   % Affine3d
   opts.tform.ang = affine3d(opts.array.ang);
catch 
   fprintf('Affine3D of angulation arrays not executed.\n');
end



