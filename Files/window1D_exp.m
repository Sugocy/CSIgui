function window = window1D_exp(L, varargin)
% Enable two different exponential decay functions:
%
% Type 1: L = window size;   var(1) = final target value
% --------------------------------------------------------------------
% Given the sample size of the data and a target between 0 and 1, a
% exponential decay is returned with end value "target". If target set to
% >1 decay becomes growth ;)
%
% Type 2: L = window size; var(1) = Exp Length; Var(2) = Stretch
% --------------------------------------------------------------------
% Create expontential apodization window given window size L and
% exponential length (EL) and strength (EL) or stretch.


% Type 1
if nargin == 2
    window = window1D_exp_type1(L, varargin{1});
% Type 2
else
    window = window1D_exp_type2(L, varargin{1},varargin{2});
end

