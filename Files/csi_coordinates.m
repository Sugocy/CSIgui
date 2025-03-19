function data = csi_coordinates(data, origin, vox_cor, fft_cor)
% COMPATIBLE WITH CSI DATA LOADED FROM LIST/DATA FILES.
%
% Data:         Offcenter, dimensions, resolution.
%
% Origin:       Offcenter origin set to topleft (default) 
%               or center of image volume.
%
% Optional - voxel correction: subtract 0.5 a voxel to offcenter
% Optional - fourier transform correction: add 0.5 voxel
%
% Offcenter: 1=LR; 2=AP; 3=FH;
% Assume it describes the center of first top left voxel transmitted (equal
% definition to dicom header def. of the offcenter tag). NOT TRUE..
%
% NB Be aware of the extra dimension in data.res and data.dim for
% the sample dimension itself. Clear this index. 
%
% contact:
% qhoutum2@umcutrecht.nl | quincyvanhoutum@gmail.com


if     nargin == 1, origin = 'topleft'; vox_cor = 0; fft_cor = 0;
elseif nargin == 2, vox_cor = 0; fft_cor = 0;
elseif nargin == 3, fft_cor = 0;
end

% Field of View
% data.fov = data.res .* data.dim;

switch origin
    case 'topleft'      % Offcenter defined at center of top left voxel 
        
        % Get fields from data struct;
        offc  = data.offcenter; N = data.dim; res = data.res;
        
        % For the Row and Column
        for kk = 1:3

            if vox_cor == 0
                % NO center-of-voxel correction
                Vbegin = offc(kk);                
            else
                % YES center-of-voxel correction
                Vbegin = offc(kk) + (0.5*res(kk));    
            end

            Vend  = Vbegin + (res(kk)*N(kk)) - res(kk);

            % Create grid vector and limits field.
            data.vec{kk}     = Vbegin:res(kk):Vend;
            data.lim(kk,1:2) = data.vec{kk}([1 end]);
        end

        
    case 'center' % Offcenter defined at center of volume.
       
        negfac = [1 1 1];
        for kk = 1:size(data.dim,2) 
                
            % Parameters for direction N
            N = data.dim(kk); res = data.res(kk); 
            offc = data.offcenter(kk);
            
            % Grid vector
            % For odd; the vector is always moved "down" w.r.t axis.
            vec_form = @(N) ceil((N/2)*-1 : 1 : (N-1)/2);
            vec = vec_form(N) ; 
            
            % Calculate coordinates per point and add offcenter.
            vec = (vec .* res) + offc;               
            
            % Correct half a voxel shift due FFT method
            if fft_cor, vec = vec + negfac(kk).* (0.5 .* res); end


            % Store the vector
            data.vector{kk} = vec; 
            % Limits of voxel i.e. the center
            data.lim(kk,1:2) = data.vector{kk}([1 end]);
            % Limits of volume i.e. the grid
            data.lim_vol(kk,1:2) = data.lim(kk,1:2) + (0.5*[-res res]);
        end

    otherwise
        disp('Wrong orientation given.');
end
