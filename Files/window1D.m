function win = window1D(L, type, state, ac)
% win = window1D(L, type, state, ac)
%
% Return a 1D filter window of size L in symmetric or periodic state. 
% Multiple filter types are available.
%
% Input: window size  - L, integer and non-zero size of the window
%        window type  - "hamming_ideal" (default), "hamming" classic, 
%                       "hann", "blackman" or "flattop". See below for
%                       details on the ideal hamming window.
%        window state - Symmetric (default) or periodic.
%        window coeff - Actual coefficients (ac) for general window 
%                       formula. "Type" will be ignored (set to custom) if 
%                       the 5 ac are given.
%
% Output: window      - 1D window if length L.
%
% Generalised window equation:
% wineq = @(ac, x) ac(1) - ac(2)*cos(2*pi*x) + ac(3)*cos(4*pi*x) - ...
%                          ac(4)*cos(6*pi*x) + ac(5)*cos(8*pi*x);
%
% Hamming_ideal:
% In the equiripple sense, the optimal values for the coefficients are 
% a0 = 0.53836 and a1 = 0.46164. Equiripple: the ripple in the frequency 
% response is spread out evenly across, avoiding large peaks. The zero
% crossing at 5pi/(N-1) will cancel the first sidelobe present in the Hann
% window, which is the main goal/aim of the hamming windows. 
% NB. This is NOT equal to equiripple Parks–McClellan methods and it has
% only as small (neglible) <0.03dB change in sidelobe level.
%
% Quincy van Houtum. v06.2017
% quincyvanhoutum@gmail.com


% // --- Handle input

% If window size equals 1.
if L == 1, win = 1; return; end

% Check input-arguments
if     nargin == 1, state = 'symmetric'; type = 'hamming_ideal'; 
elseif nargin <= 2, state = 'symmetric'; 
elseif nargin == 4, state = 'custom';
end
% If no type input - empty - set to default.
if isempty(type), type = 'hamming_ideal'; end

% fix uppercase-issues
type = lower(type); state = lower(state);

% // --- Set window type

switch type
    case 'hamming_ideal',  ac = [ 0.53836 0.46164 0 0 0 ];
    case 'hamming',        ac = [ 0.54 0.46 0 0 0 ];
    case 'hann'   ,        ac = [ 0.50 0.50 0 0 0 ];                         % ac = [1.78 1.73 0 0 0 ]; % https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.1110       
    case 'blackman',       ac = [ 0.42 0.50 0.08 0 0 ];
    case 'flattop',        ac = [ 0.21557895  0.41663158  ...
                                  0.277263158 0.083578947 0.006947368 ];
    case 'custom',   fprintf('Using custom coefficients for window.\n');
    otherwise,       fprintf('Incorrect window type selected. Abort.\n'); 
                     win = NaN; return;
end


% Set window state e.g. sflag.
switch state
    case 'symmetric', x = (0:L-1)'./(L-1);
    case 'periodic',  x = (0:L)'./(L);
end

% // --- Create window

% Generalised window equation.
wineq = @(ac, x) ac(1) - ac(2)*cos(2*pi*x) + ac(3)*cos(4*pi*x) - ...
    ac(4)*cos(6*pi*x) + ac(5)*cos(8*pi*x);

% Done. Define output window w.
win = wineq(ac, x);
win = win(1:L); 


