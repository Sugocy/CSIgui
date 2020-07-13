function w = window1D(L, type, state, ac)
%%%% Description:                    Create filter window of specific type
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-06
%%% --------------------------------------------------------------------
%%%
%%% Input: window size  - L, integer and non-zero size of the window
%%%        window type  - Hamming, hann, blackman or flattop.
%%%        window state - Symmetric (default) or periodic.
%%%        window coeff - Actual coefficients (ac) for general window 
%%%                       formula. If this input is given, type is set 
%%%                       to custom and these ac-values are used.
%%%
%%% Output: window      - 1D window if length L.
%%%
%%% Generalised window equation:
%%% wineq = @(ac, x) ac(1) - ac(2)*cos(2*pi*x) + ac(3)*cos(4*pi*x) - ...
%%%                          ac(4)*cos(6*pi*x) + ac(5)*cos(8*pi*x);

% If window size equals 1.
if L == 1, w = 1; return; end

% Handle input
if     nargin == 1, state = 'symmetric'; type = 'hamming'; 
elseif nargin <= 2, state = 'symmetric'; 
elseif nargin == 4, state = 'custom';
end

% fix uppercase-issues
type = lower(type); state = lower(state);

%% Set window type

switch type
    case 'hamming',  ac = [ 25/46 0.46 0 0 0 ];
    case 'hann'   ,  ac = [ 0.50 0.50 0 0 0 ];
    case 'blackman', ac = [ 0.42 0.50 0.08 0 0 ];
    case 'flattop',  ac = [ 0.21557895  0.41663158  ...
                            0.277263158 0.083578947 0.006947368 ];
    case 'custom',   disp('Using custom coefficients for window.');
    otherwise,       disp('Incorrect window type selected. Abort.'); 
       return;
end


%% Set window state e.g. sflag.
switch state
    case 'symmetric', x = (0:L-1)'./(L-1);
    case 'periodic',  x = (0:L)'./(L);
end

%% Create window

% Generalised window equation.
wineq = @(ac, x) ac(1) - ac(2)*cos(2*pi*x) + ac(3)*cos(4*pi*x) - ...
    ac(4)*cos(6*pi*x) + ac(5)*cos(8*pi*x);

% Done. Define output window w.
w = wineq(ac, x);
w = w(1:L); w(w > -1e-15 & w <-1e-15) = 0;


