function data = PARcoordinates(data, origin, varargin)
%
% Data:         Offcenter, dimensions, resolution.
%
% Origin:       Offcenter origin set to topleft (default) 
%               or center of image volume.
% 
%
% Offcenter: 1=LR; 2=AP; 3=FH;
% Assume it describes the center of first top left voxel transmitted (equal
% definition to dicom header def. of the offcenter tag).
%
% For CSI data be aware of the extra dimension in data.res and data.dim for
% the samples itself. Remove them;)
%
% PARREC and DICOM share equal offcenter definition!
%
%
%

if nargin == 1, origin = 'topleft'; end
if nargin == 3, vox_correction = varargin{1}; else vox_correction = 0; end

switch origin
    case 'topleft'                  % Offcenter == center of top left voxel 
        % Get fields from data struct;
        offc  = data.offcenter; N = data.dim; res = data.res;
        
        % For the Row and Column
        for kk = 1:2

            if vox_correction == 0
                % NO center-of-voxel correction
                tmp = offc(1,kk);                
            else
                % YES center-of-voxel correction
                tmp  = offc(kk) + (0.5*res(kk));    
                
                % Is it minus or plus half a voxel? 
                % Depends on origin of axis-scheme....
                % ASSUME: LTV is relative to top left corner being 
                % axis origin and increasing towards right side of image
                % thus ADD HALF VOXEL. 
            end
            
            Vbegin = tmp; 
            Vend   = tmp + (res(kk)*N(kk)) - res(kk);

            % Create grid vector and limits field.
            data.vec{kk}     = Vbegin:res(kk):Vend;
            data.lim(kk,1:2) = data.vec{kk}([1 end]);
        end
        % SLICE
        data.vec{3} = data.offcenter(:,3);
        data.lim(3,1:2) = data.vec{3}([1 end]);
        
    case 'center'
        for kk = 1:2                       % Offcenter == center of volume.
                       
            N     = data.dim(kk); res = data.res(kk);
            offc  = data.offcenter(1,kk);
            
            if vox_correction == 1
                odd = mod(N,2);
                if odd % There is an actual middle voxel
                    Vbegin = offc-(res*(N-1) *0.5);
                    Vend   = offc+(res*(N-1) *0.5);
                else   % There is no middle voxel
                    Vbegin = (offc-(0.5*res))-(res*(N-1)*0.5);
                    Vend   = (offc+(0.5*res))+(res*(N-1)*0.5);
                end
                vec = Vbegin:res:Vend;
            else
                % Around zero.
                vec = (-1.*N/2 ): 1 : (N-1)/2;
                vec = vec .* data.res(kk) + offc;
            end

            data.vec{kk}     = vec; 
            data.lim(kk,1:2) = data.vec{kk}([1 end]);
        end
        
        % SLICE
        data.vec{3} = data.offcenter(:,3);
        data.lim(3,1:2) = data.vec{3}([1 end]);
        
    otherwise
        disp('Wrong orientation given.');
end
