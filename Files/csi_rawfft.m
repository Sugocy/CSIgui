function fft_data = csi_rawfft(csi, dim, shift_method, loop)
%%%% Description:                            Forward Fourier of CSI volume
%%% Creator: Dr. Q. van Houtum       Version: 1.5          Date: 2023-08
%%% --------------------------------------------------------------------
%%%
%%% Applies fourier over the Kx, Ky and Kz dimensions of the given
%%% volume. 
%%%
%%% Input: csi              - data representing raw-kspace-csi data.
%%%        dim              - index dimension of kx, ky, kz in csi. (1xN)
%%%        shift_method     - switch between circshift(0) or fftshift(1) or
%%%                           automatic(2) based on odd/even matrix-size.
%%%        loop             - process data per dimension (1) or step by
%%%                           step (0, default). Latter is slightly faster.
%%%
%%% NFO: FFT-shift for odd sized matrices, the circular shift is used 
%%% for even sized matrices.
%%%                           
%%% Contact: quincyvanhoutum@gmail.com


% Process input
if nargin < 3 , shift_method = 2; loop = 0; end
if nargin < 4 , loop = 0; end

% Prepare variables
fft_data = csi; sz = size(fft_data);

% odd even bool
odd_bool = mod(sz(dim),2);    
% Number of positions to shift for circshift-fnc
shift_val = ceil( ( sz(dim) ./ 2 ) + 1 );
% shift_val = ceil( ( sz(dim) ./ 2 ) + 1);



% NFO output
fprintf('K-space Matrix: %i x %i x %i.\n', sz(dim));
if ~loop 
    fprintf('Shifting all k-space dimensions before iFFT.\n');
else 
    fprintf('Processing all k-space dimensions consecutively.\n');
end

% CIRC SHIFT
if shift_method == 0
    if ~loop
        
        % Shift every dimension circular
        for kk = 1:size(dim,2)
            fft_data = circshift(fft_data,shift_val(kk), dim(kk)); 
        end
    
        % Apply iFFT
        for kk = 1:size(dim,2), fft_data =  ifft( fft_data,[], dim(kk) ); end
    
        % Shift back every dimension
        for kk = size(dim,2):-1:1
            fft_data = circshift(fft_data,-1.*shift_val(kk), dim(kk)); 
        end

    else
        % Process every dimension consecutively.
        for kk = 1:size(dim,2)
            fft_data = circshift( ifft( circshift( ...
                fft_data, shift_val(kk), dim(kk))...
                    , [], dim(kk))...
                        , -1.*shift_val(kk), dim(kk));
        end
    end

% FFT SHIFT
elseif shift_method == 1
        
    if ~loop 
        
        % Shift everything first.
        for kk = 1:size(dim,2)
            fft_data =  ifftshift(fft_data, dim(kk));
        end

        % FFT all dimensions in one loop.
        for kk = 1:size(dim,2)
            fft_data =  ifft(fft_data,[], dim(kk));
        end

        % Shift everything back
        for kk = 1:size(dim,2)
            fft_data =  fftshift(fft_data, dim(kk));
        end
    
    else

        % FFT and FFTshift in one loop.
        for kk = 1:size(dim,2)
            fft_data =  ...
            fftshift(ifft(ifftshift(fft_data, dim(kk)), [], dim(kk)),dim(kk) );
        end
    end
    
% CIRC/FFTSHIFT based on odd/even    
elseif shift_method == 2
    
    if ~loop
        
        % SHIFT 1/2
        for kk = 1:size(dim,2)
            if odd_bool(kk) % Odd
                fft_data =  ifftshift(fft_data, dim(kk));
            else % Even
                fft_data = circshift(...
                    fft_data, shift_val(kk), dim(kk)); 
            end
        end

        % FFT
        for kk = 1:size(dim,2)
            fft_data = ifft(fft_data,[], dim(kk));
        end

        % SHIFT 2/2
        for kk = 1:size(dim,2)
            if odd_bool(kk) % Odd
                fft_data =  fftshift(fft_data, dim(kk));
            else % Even
                fft_data = circshift(...
                    fft_data,-1 .* shift_val(kk), dim(kk)); 
            end
        end

    else
        
        % Loop each dimension and process fully.
        for kk = 1:size(dim,2)
            % If odd use fft-shift
            if odd_bool(kk)
                fft_data = fftshift( ifft( ifftshift( ...
                    fft_data, dim(kk) ), [], dim(kk)), dim(kk));
            else
            % If even use circ-shift
                fft_data = ...
                    circshift(fft_data,ceil(sz(dim(kk))/2+1), dim(kk)); 
                fft_data = ifft(fft_data,[],dim(kk));
                fft_data = ...
                    circshift(fft_data,-1.*ceil(sz(dim(kk))/2+1), dim(kk));
            end
        end
    end
    
end

