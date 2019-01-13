function out = gaussmf(x,opts)
% Opts equals vector [sig c], the std and center.
func = @(x,sig,c) exp( (-1.*(x-c).^2) ./ (2.*sig.^2) );
out = func(x,opts(1),opts(2));

% sqrt(a/pi) .* exp(-a.*x.^2)
% exp( -((ii-1)^2) / (2*apo_f_gauss^2))