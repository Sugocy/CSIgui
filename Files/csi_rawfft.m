function fft_data = csi_rawfft(csi, dim, shift_method)
%%%% Description:                            Forward Fourier of CSI volume
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%%
%%% Applies fourier over the Kx, Ky and Kz dimensions of the given
%%% volume. 
%%%
%%% Input: csi              - data representing raw-kspace-csi data.
%%%        dim              - index dimension of kx, ky, kz in csi. (1xN)
%%%        shift_method     - switch between circshift(0) or fftshift(1)
%%%                           
%%%        NB:  Circshift shifts all dimensions first before FFT.
%%%             FFTshift shifts the dimension, applies FFT and shifts back
%%%             before moving to the next dimension to FFT over.
%%%                           
%%%
%%% Contact: qhoutum2@umcutrecht.nl

% Process input
if nargin ~= 3, shift_method = 0; end
fft_data = csi; sz = size(fft_data);



% CIRC SHIFT
if shift_method == 0
    
    % Shift every dimension circular
    for kk = 1:size(dim,2)
        fft_data = circshift(fft_data,ceil(sz(dim(kk))/2+1), dim(kk)); 
    end

    % Apply iFFT
    for kk = 1:size(dim,2), fft_data =  ifft( fft_data,[], dim(kk) ); end

    % Shift back every dimension
    for kk = size(dim,2):-1:1
        fft_data = circshift(fft_data,-1.*ceil(sz(dim(kk))/2+1), dim(kk)); 
    end

% FFT SHIFT
elseif shift_method == 1
    
    
    % FFT and FFTshift in one loop.
    for kk = 1:size(dim,2)
        fft_data =  ...
        fftshift(ifft(ifftshift(fft_data, dim(kk)), [], dim(kk)),dim(kk) );
    end
    
end

