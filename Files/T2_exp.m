function [T2, R2, CI, Ampl, Nconst, fitfunc] = T2_exp(echoes, TE, W, show_plot) 
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

% Fit function: param(1) = T2 - param(2) = A: A*e^(-x/T2) + C
fitfunc = @(param,x) param(2).*exp(-x./param(1));                            

% Find a value for C - noise component
% Using std between last echoes
noise = std(echoes(:,end), [], 'all');

% Get number of echoes for all requested fit and # of fits.
nEcho = size(echoes,2);
nFits = size(echoes,1);


% Containers
T2 = NaN(nFits,1); R2 = NaN(nFits,2); Nconst = NaN(nFits,1);
Ampl = NaN(nFits,1); CI = NaN(nFits,2);
R2a = NaN(nFits,1); R2b = NaN(nFits,1);


if isempty(gcp('nocreate'))
    c = parcluster('Processes'); 
    nCores = ceil(feature('numcores') .* 0.9);
    c.NumWorkers = nCores; saveProfile(c);
    parpool(c);
end
 
TE_end = TE(end); TE_start = TE(1); TE_three = TE(3);
% Using 3rd echo --> more signal than last - thus noise in init-param
parfor fi = 1:nFits
    vals = echoes(fi,:);    
    
    % Nonlinear fit estimates for parameters
    T2est   = -1*(TE_three-TE_start)/(log(vals(3))-log(vals(1)));             % Math: e^rx = e^-x/T2 => r = 1/T2;   r =  ln(Y2)-ln(Y1)/X2-X1 
    Aest    = vals(1)/exp( -1*TE_start(1)/T2est );                            % Math: y = A*e^rx = A*e^-x/T2;       A = y/e^rx;
    % Nest    = noise;                                                        % Math: need to add C to the equations above

    % Safety: NO imaginary estimates due bad experimental data.
    if ~isreal(T2est), T2est = real(T2est); end
    if ~isreal(Aest),  Aest  = real(Aest);  end
    
    % Initial parameters
    init_param = [T2est Aest];

    warning('off', 'all')
    % if W == 0
    %     try
    %     [fit_par, R, ~, CovB, ~, ~] = ...
    %         nlinfit(TE, vals, fitfunc, init_param); 
    %     catch
    %         fit_par = [NaN, NaN, NaN]; R = NaN; CovB = NaN;
    %     end
    % else
    try
    [fit_par, R, ~, CovB, ~, ~] = ...
        nlinfit(TE, vals, fitfunc, init_param, 'Weights', W); 
    catch
        fit_par = [NaN, NaN, NaN]; R = NaN; CovB = NaN;
    end
    % end
    

    % Store parameters
    T2(fi,1) = fit_par(1); Ampl(fi,1) = fit_par(2);
    % Nconst(fi,1) = fit_par(3);
    
    % Values using fit-parameters
    vals_fit = fitfunc(fit_par, TE);
        
    % Method I
    SSres   = sum(R.^2);                                                    % Sum of squares residual
    SStot   = sum((vals - mean(vals)).^2);                                  % Sum of squares total
    % R2(fi, 1)   = 1 - (SSres/SStot);                                      % Method 1: R-square calculation 
    R2a(fi) = 1 - (SSres/SStot);                                        

    % Method II
    C = corrcoef(vals,vals_fit);                                            % Correlation coefficients of y/yfit
    % R2(fi, :) = C(1,2)^2;                                                 % Method 2: R-square calculation 
    R2b(fi) = C(1,2).^2;
    
    % Confidence interval
    CItmp    = nlparci(fit_par, R,'Covar',CovB);
    CI(fi,:) = CItmp(1,:);
    
    % Show plot if requested.
    if show_plot == 1
        figh = figure(8484); clf(figh, 'reset');
        plot(TE, vals, '-b'); hold on; plot(TE, vals_fit, '-r');
        title([num2str(fi) ' - R^2: ' num2str(R2(fi,1))]);
        % This doesnt work with parfor - using correct R2b
        % drawnow; pause(1); if fi == nFits, close(figh); end
    end
            

end

% Catenate R2 calc for output - needed bc of parfor loop.
R2 = cat(3, R2a, R2b);