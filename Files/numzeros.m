function [ndecimals, acc] = numzeros(val, acc)
% Find number of zero-decimals using a addition/subtraction method
%
%
% Usefull for automated printing of values, example:
% random variable val = 0.00000001424145; 
% numzero = numzeros(val,20); 
% fprintf('%3.3fe-%i\n', val*10.^(numzero+1), numzero+1)
% Number of zeros is 7, in engineering format - 1e-8 to have one integer
% before the decimal-separator.
if nargin < 2, acc = 32; end


val = abs(val); sign = -1; 

val_bool = 1; dt = 10.^-(acc); n = 0;
while val_bool        
    val = val + (sign * dt * (10.^n));
    if sign, val_bool = val > 0; else, val_bool = val < 0; end
    n = n+1;
end

ndecimals = acc - n + 1;