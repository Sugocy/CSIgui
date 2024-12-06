function [window, formula] = window1D_exp_type1(sample_size, target)
% Given the sample size of the data and a target between 0 and 1, a
% exponential decay is returned with end value "target".
%
%
% Source: http://zone.ni.com/reference/en-XX/help/370051V-01/cvi/libref/
% analysisconcepts/characteristics_of_different_smoothing_windows/

if nargin == 2

% 1. Define exponential decay function
expw = @(f,s) f.^( (1:s)./ (s-1) );

% 2. Define exponential decay function - Equal to the above.
% expw = @(f,s) exp(1).^( (1:s .* log(f)) ./ (s - 1));

% 3. Calculate window
window = expw(target,sample_size);

% 4. Set formula as output
if nargout > 1, formula = expw; end

else
    
    
    
    
    
    
end