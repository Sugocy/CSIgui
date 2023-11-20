function [spec, dphase, cphase] = csi_autoFirstPhase(spec, method, acc, disp)
% Quick and dirty but functional method to apply first order phase
% corrections. 
%
% Input:    spectrum, method
% Output:   spectrum, first order phase changes, constant phase-change.
%           step-accuracy for calculating UAC (default: 1e-3)
%
% Method 0 (default)
% Maximizes the area under the curve of the real part of the spectrum by 
% adding a linear phase offset to the spectrum.
%
% Method 1
% Uses method 0 but also applies a zeroth order phase correction after the
% first order corrections again.
%
% quincyvanhoutum@gmail.com

if nargin < 2, method = 0; end
if nargin < 3, acc = 2e-2; end % 4e-2 would also work
if nargin < 4, disp = 0; end

% Split spectrum
m = abs(spec);
p = angle(spec);

% Original UAC
uac_ori = trapz(real(spec),1);

% Preperation
ns = size(p,1);
nv = (-ns/2:1:(ns/2)-1)';

% Spectrum function
spec_fnc = @(ph) m .* cos(ph) + m .* sin(ph) * 1i;

% Real-part function
real_fnc = @(ph) m .* cos(ph);

% Linear function: the slope is a multiplication of 1-rad
lin_fnc = @(sl) sl .* (pi./180 .* nv);

% Calculate AUC for few slopes
sl_inp = -(pi/180)*360:acc:(pi/180)*360; % Slope input
sl_inp_uac = trapz(real_fnc(p + lin_fnc(sl_inp)),1); % Slope input UAC

% Find maximum UAC
[uac_max, loc] = max(sl_inp_uac);

% Set first order phase change
dphase = lin_fnc(sl_inp(loc));

% If method == 1, apply additional zero-order phase corrections (constant)
if method
    cp_inp = (-180:1:180) .* (pi/180);
    cp_inp_uac = trapz(real_fnc(p + cp_inp + dphase),1);
    [uac_max_zero, loc_zero] = max(cp_inp_uac);
    cphase = cp_inp(loc_zero);
else
    cphase = 0; uac_max_zero = uac_max;
end

% Reconstruct spectrum
spec = spec_fnc(p + dphase + cphase);

% Display change
if disp
    npl = 3;
    f = figure('NumberTitle','Off',...
        'Name', 'CSI First-order Phase Corrections.');

    subplot(npl,1,1); plot(real_fnc(p)); 
    legend(['Original: ' sprintf('%.2f',uac_ori)]);
    
    subplot(npl,1,2); plot(real_fnc(p + dphase)); 
    legend(['First order corrections: ' sprintf('%.2f',uac_max)]);
   
    subplot(npl,1,3); plot(real(spec)); 
    legend(['First order corrections: ' sprintf('%.2f',uac_max_zero)]);
    
    for kk = 1:size(f.Children)
        if strcmp(f.Children(kk).Type, 'axes')
            f.Children(kk).XDir = 'reverse';
        end
    end

    % fprintf('UAC original: %.2f\n', uac_ori);
    % fprintf('UAC first-order corrected: %.2f\n', uac_max);
    % fprintf('UAC first- and zeroth-order corrected: %.2f\n', uac_max_zero);
end
