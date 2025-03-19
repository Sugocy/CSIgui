function N = findDecimals(val)
% Returns the number of decimal places actually has a value other than 0.
% e.g. 0.0000054321 == 5.4321e-06 thus N = 6;
%
% Be aware, if integer values are in "val", this is subtracted and some 
% accuracy issues mighta rise (eps).


if round(val) > 0, val = abs(val - floor(val)); end

lim = 32; rnd = NaN; aim = 0; n = -1;
while rnd ~= aim
    n = n + 1;
    rnd = round(val,lim-n);
    if n > lim, break; end
end


if n > lim
    val_str = num2str(val);
    ind = strfind(val_str,'.');
    N = numel(val_str) - ind;
else
    % #decimals
    N = lim - n + 1;
    
    % Correction for 0.5-cutoff
    if (val*10^N) <= 1, N = N + 1; end
end

