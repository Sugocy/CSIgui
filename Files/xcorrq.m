function [corval, lags, pval] = xcorrq(x,y)
%%%% Description:                    Simplified cross correlation xcorr()
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2019-06
%%% --------------------------------------------------------------------
% [] = xcorrq(x,y) 
%                       returns the cross correlation of x and y.
%
% Input
%       X = Nx1 vector;       Y = Nx1 vector
%       Default convert to Nx1 implemented. Therefor only 1 column
%       possible.
%
% MATLAB DESCRIPTION:
% ..returns the cross-correlation of two discrete-time sequences, x and y.
% Cross-correlation measures the similarity between x and shifted (lagged) 
% copies of y as a function of the lag.
%
% xcorr() is part of the signal processing toolbox. This is a simplified
% version to workaround version incompability. Advised is using the default
% xcorr() from the signal processing toolbox whenever possible.
%
% Usage:
% crval = xcorrq(x,y); [mx, ind] = max(crval); shift = lags(ind)
% Best probable shift/lag  = lags(ind); 
% Can be refined using p-value!
%
% GNU GENERAL PUBLIC LICENSE - Version 3, 29 June 2007
%
% Contact: quincyvanhoutum@gmail.com


% Prepare ---- %

% Confirm Nx1 vector
szx = size(x); szy = size(y);
if szx(1) < szx(2), x = x'; szx = fliplr(szx); end
if szy(1) < szy(2), y = y'; szy = fliplr(szy); end

% Equal N
if     szx(1) > szy(1), y( end+1 : szx(1)-szy(1) ) = 0; 
elseif szy(1) > szx(1), x( end+1 : szy(1)-szx(1) ) = 0; 
end

% Correlate ---- %

% Define shifts e.g. lags
n = size(x,1)-1; lags = -n:n; m = size(y,2);
corval = zeros(m,1); pval = zeros(m,1);
% Loop each shift
for kk = 1:size(lags,2)
	[corval(kk), pval(kk)] = corr(x, circshift(y,lags(kk),1));
end


