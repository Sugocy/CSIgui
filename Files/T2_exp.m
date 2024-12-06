function [T2, R2, CI, Ampl] = T2_exp(echoes, TE, W, show_plot) 
% Echo data, echotimes or echo time difference, show plot and weight.
% If no weighting required - set W to [];
% Echoes (nFits x nEchoes);
% TE     (1 x nEchoes) || (1 x 1);
% W      (1 x nEchoes);
%
% fit function = @(param,x) param(2)*exp(-x/param(1));       

switch nargin
    case 2, W = 0; show_plot = 0;
    case 3,        show_plot = 0;
end

% Calculate TE
if size(TE,2) == 1, TE = TE:TE:((size(echoes,2))*TE); end

% Fit function: param(1) = T2 - param(2) = A: A*e^(-x/T2)
fitfunc = @(param,x) param(2).*exp(-x./param(1));                            

% Get number of echoes for all requested fit and # of fits.
nEcho = size(echoes,2);
nFits = size(echoes,1);

T2 = NaN(nFits,1); R2 = NaN(nFits,2); Ampl = NaN(nFits,1);
CI = NaN(nFits,2);
for fi = 1:nFits
    vals = echoes(fi,:);    
    
    % Nonlinear fit estimates for parameters
    T2est   = -1*(TE(end)-TE(1))/(log(vals(end))-log(vals(1)));             % Math: e^rx = e^-x/T2 => r = 1/T2;   r =  ln(Y2)-ln(Y1)/X2-X1 
    Aest    = vals(1)/exp( -1*TE(1)/T2est );                                % Math: y = A*e^rx = A*e^-x/T2;       A = y/e^rx;
   
    % Safety: NO imaginary estimates due bad experimental data.
    if ~isreal(T2est), T2est = real(T2est); end
    if ~isreal(Aest),  Aest  = real(Aest);  end
    
    
    init_param = [T2est Aest];
    
    %     Aest2      = vals(2)/exp( -1*TE(2)/T2est );
    %     initparam2 = [T2est Aest2]; Try again... Later.
    if W == 0
        [fit_par, R, ~, CovB, ~, ~] = ...
            nlinfit(TE, vals, fitfunc, init_param); 
    else
        [fit_par, R, ~, CovB, ~, ~] = ...
            nlinfit(TE, vals, fitfunc, init_param, 'Weights', W);
    end
    
    % Store parameters
    T2(fi,1) = fit_par(1); Ampl(fi,1) = fit_par(2);
    
    % Values using fit-parameters
    vals_fit = fitfunc(fit_par, TE);
        
    % Method I
    SSres   = sum(R.^2);                                                    % Sum of squares residual
    SStot   = sum((vals - mean(vals)).^2);                                  % Sum of squares total
    R2(fi, 1)   = 1 - (SSres/SStot);                                        % Method 1: R-square calculation 
    
    % Method II
    C = corrcoef(vals,vals_fit);                                            % Correlation coefficients of y/yfit
    R2(fi, 2) = C(1,2)^2;                                                   % Method 2: R-square calculation 
    
    % Confidence interval
    CItmp    = nlparci(fit_par, R,'Covar',CovB);
    CI(fi,:) = CItmp(1,:);
    
    % Show plot if requested.
    if show_plot == 1
        figh = figure(8484); clf(figh, 'reset');
        plot(TE, vals, '-b'); hold on; plot(TE, vals_fit, '-r');
        title([num2str(fi) ' - R^2: ' num2str(R2(fi,1))]);
        drawnow; pause(1); if fi == nFits, close(figh); end
    end
            
end