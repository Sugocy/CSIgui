function outp = csi_fit_FASeries(xdata, ydata, fit_method, phase_data, ...
    add_zero, plot_data, save_plot, iter, dirstr, varargin)
% Given a flip angle dynamic data set with increasing flip angle for each
% spectrum at a specific index, this function...
%
% ...will return a graph with the peak signal over FA 
% including a fit and the corresponding zero-crossing which translates 
% to the 180deg flip angle. The fit itself will be weighted using the SNR
% of the given spectral-data.
%
% Input:
%
%       1. FA           = flip angles // voltage // etc (x-axis)
%       2. data         = spectrum for each FA
%       3. phase_data   = 0, 1 or 2; 
%                    |1: phase all data prior to processing
%                    |2: phase spectra only after initial zero-crossing 
%                    |0: do not phase data at all.
%       4. fit_method   = 0, 1;
%                    |(1) -REQUIRES varargin(1)=TR and varargin(2)=T1;
%                    |Sine (0) or MR signal equation (1)
%       5. add_zero     = Adds S(FA=0) = 0 to data before fitting.
%       6. plot_data    = Boolean for plotting fit and result
%       7. save_plot    = Boolean for saving plot
%       8. counter      = Number for data-iteration
%       9. dirstr       = directory to save data
%      10. variables    = TR and T1, required for fit-method 1
%                       = ppm-axis
%
% NFO Phasing:
% If the zero-crossing is present in the data, it needs to be phased
% correctly for a correct fit. This means the MR signal should be negative
% after the zero-crossing. This is important for both the sine and
% MR-equation method.
%
% NFO Method:
% If TR is many orders larger than the T1, T1 weighting will be
% neglible to almost zero; use the sine function.
% If there is no zero-crossing in the data, use the sine-function.
% If the T1 is unknown, it is better to use the sine-function.
% If the flip-angles (xdata) is not in degrees, use the sine-function.
%
% This function uses: csi_SNR, csi_findZeroCross, csi_autoZeroPhase
%
%
% vThesis: https://doi.org/10.1002/nbm.4178 (2020)
% V3: updated 12/2023
% Q. van Houtum, PhD; quincyvanhoutum@gmail.com


% INPUT/SETUP % -------------------------------------------------------- %

% Process variable input arguments.
if (nargin >= 10)
    if numel(varargin) == 2
        TR = varargin{1}; T1 = varargin{2}; 
    elseif numel(varargin) == 1
        ppm = varargin{1};
    elseif numel(varargin) == 3
        TR = varargin{1}; T1 = varargin{2}; ppm = varargin{3};
    end
end



% Create output variable
outp = struct; 

% Ease of coding: data of interest
doi = double(ydata);

% Filename for use of saving data and figures
strfn = [sprintf('%02i', iter), '_FAdynamic'];

% PHASING % ------------------------------------------------------------ %


% Phase-correct all spectra
doi_main = doi;
if phase_data == 1, doi_phased = phaseCorrection(doi);
else,               doi_phased = doi;    
end

% ------------------------------------------ %
% Expected: All spectra are in absorption mode.

% Phase-correct after zero-crossing only
if phase_data >= 1
    [doi_main, ~] = phaseCorrection_zeroCrossing(doi_phased);
end


% FIT PREP % ----------------------------------------------------------- %
% If zero-crossing is present in the data, the maximum value of the spectra
% before the zero-crossing and minimum values after the zero-crossing need
% to be calculated accordingly.

% Find initial zero-crossing in data: using abs-magnitude of min/max of a
% spectrum.
mx = max(real(doi_main)); mn = min(real(doi_main));
mxmn = (abs(mx) <= abs(mn));

% Get max/min according to initial zero-crossing
vals = NaN(1,size(mx,2));
vals(mxmn==0) = mx(mxmn == 0); vals(mxmn) = mn(mxmn);

% Add S(FA = 0) = 0 to data
if add_zero, vals = [0 vals]; xdata = [0 xdata]; end

% X values at high res for zero-crossing and plot.
xdat_HR = xdata(1):1:xdata(end);

% Weighting for fit: SNR
snrw = csi_SNR(doi_main, round(size(doi_main,1)./8), 0, [1 size(doi_main,1)]);
if add_zero, weights = [max(snrw) snrw]; else, weights = snrw; end

% FIT % ---------------------------------------------------------------- %
% According to fit-method, fit equation through data-points.

if fit_method == 0
    % --- % FIT SINE

    % // Initial parameters

    % Where values increase and decrease again
    dy = diff(vals); grad_decr = find(dy < 0); 
    if numel(grad_decr) > 1
        [init_amp, init_amp_ind] = max(vals(1:grad_decr(1)));
    else
        [init_amp, init_amp_ind] = max(vals); 
    end
    
    % [init_amp, init_amp_ind] = max(vals); 
    init_fre = 1./(4*xdata(init_amp_ind)); init_shift = 0;
    init_par = [init_amp, init_fre, init_shift];

    % Sine function: radians with frequency
    func = @(par,x) (par(1) .* sin( 2*pi*par(2).*x + par(3) )) ;
    
    % FIT    
    beta = nlinfit(xdata, vals, func, init_par, 'Weights', weights);

    % Save data to file.
    % fn = string(datetime('now','Format','yyMddHHmmssSSS'));
    save(strjoin([dirstr '\' strfn '_Output.mat'],''), 'beta', 'func',...
        'xdata', 'vals', 'xdat_HR', 'weights', 'init_par');


elseif fit_method == 1
    % --- % Fit MRS equation 
    
    % Ratio TR/T1
    rat = (TR.*10^-3)./(T1.*10^-3);

    % Fit function for MR signal: par(1:2) = [Amplitude Period RatioTR/T1]
    func = @(par, FA) par(1) .*(sind(FA + par(2)) .* ...
    ((1-exp(-1.*par(3))) ./ (1 - exp(-1.*par(3)).*cosd(FA + par(2)))));

    % Initial parameters
    [init_amp, init_per] = max(vals); 
    init_per = xdata(init_per) * 4; init_rat = rat; 

    % FIT
    init_par = [ init_amp init_per init_rat];
    [beta, R, J, ~,  MSE] = ...
        nlinfit(xdata, vals, func, init_par, 'Weights',weights);
    
    % Statistics
    ratio_CI = nlparci(beta,R,'jacobian',J);
    T1_fit = (TR./beta(3)); T1_CI = sort(TR./ratio_CI(3,:));
    
    % T1 and TR to output
    outp.T1_fit = T1_fit;
    outp.ratio_CI = ratio_CI;
    
    % Save data to file.   
    save(strjoin([dirstr '\' strfn '_Output.mat'],''), ...
        'func', 'xdata', 'vals', 'xdat_HR', 'beta', 'R', 'J','MSE',...
        'weights', 'T1_fit', 'T1_CI');

end


% --- CALCULATE FIT-VALUES
% Input resolution fit-values
fitval = func(beta, xdata);
% High resolution fit-values
fitval_HR = func(beta, xdat_HR);

% --- STASTISTICS

% R-squared: 1 - (SSR/SST)
RSQ = 1 - (sum((vals - fitval).^2) ./ sum((vals - mean(fitval)).^2));
outp.r_square = RSQ;


% --- % ZERO CROSSING FROM FIT-DATA-range
zc = findZeroCross(fitval_HR); zc_data_ind = zc;
zc_data = xdat_HR(zc) + xdat_HR(1);
    
% --- % ZERO CORSSING FROM EXTENDED FIT
% Period/Frequency from fit-data
if     fit_method == 0, period = 1/beta(2);
elseif fit_method == 1, period = beta(2);
end

% Exteneded data
xdat_HR_ext = xdata(1):1:abs(period);
fitval_HR_ext = func(beta, xdat_HR_ext);

% Zero-crossing calcualtion using HR extended fit-data.
if fit_method == 1

    % Closest to zero-method
    zc_ext = findZeroCross(fitval_HR_ext); zc_ext_ind = zc_ext;
    zc_ext = xdat_HR_ext(zc_ext) + xdat_HR_ext(1);

elseif fit_method == 0    
    
    % Calculating zero-crossing via sine-function   
    zc_ext = -1; N = 1;
    while zc_ext < 0
        zc_ext = (-beta(3) ./ (2*pi*beta(2))) + abs( (period/2)*N );
        fprintf('Zero-crossing #%i at x = %7.2f: %3.3g\n', N, zc_ext,...
            func(beta, zc_ext));    
        zc_ext_ind = round(zc_ext)+1;
        N = N +1; 
        
        % Safety data-indexing
        if (zc_ext_ind > size(xdat_HR_ext,2))
            zc_ext = (-beta(3) ./ (2*pi*beta(2)));
            zc_ext_ind = round(zc_ext);
        end

        % Safety while-loop
        if N > 10, break; end        
    end
    fprintf('---------------------------------------\n');
end

% --- % Output structure
outp.zero_crossing_data = zc_data; outp.zero_crossing_fit = zc_ext;
outpstr = {'T1_fit','T1_CI', 'MSE'};
for kk = 1:size(outpstr,2)
    if fit_method == 1, outp.(outpstr{kk}) = eval(outpstr{kk}); end
end


% --- % Plot % --------------------------------------------------------- %
if plot_data

    try
    fh = figure(); 
    
    % Graph limits.
    ylim_val(1) = min(real(cat(1, doi(:), doi_main(:))));
    ylim_val(2) = max(real(cat(1, doi(:), doi_main(:))))*1.05;
    
    % ------------------------------------------ %
    % Bottom Row Plot 

    % The color orange predefined.
    orange = [0.9290 0.6940 0.1250];

    % Plot fit-values, fit-result and zero-crossings
    ax_bottom = subplot(2,1,2); 
    % Plot original fit-values
    plot(xdata, vals,'--ob','LineWidth', 1.25); hold on;     
    % Plot extended fit
    hold on; plot(xdat_HR_ext, fitval_HR_ext, '-r', 'LineWidth', 1.5);  
    ylabel('Amplitude [a.u.]'); xlabel('Dynamic Range');

    % Plot calculated zero-crossing.
    zc_list = [zc_data zc_ext]; zc_ind = [zc_data_ind zc_ext_ind];
    for kk = 2:2 % ZC from data is disabled.
        if kk == 1, tmpy = fitval_HR(zc_ind(kk));
        else,       tmpy = fitval_HR_ext(zc_ind(kk));
        end
        zcoi = zc_list(kk);
        plot(zcoi, tmpy,'om','MarkerSize', 6,'LineWidth',2);        
        text(zcoi, (max(fitval_HR_ext) - tmpy)./10,...
            sprintf('%3.2f\n', zcoi), ...
            'Fontweight', 'bold','FontSize', 8);        
    end
    % Zero-line
    plot([xdata(1) xdat_HR_ext(end)],[0 0], '--k','Linewidth', 1.5);

    % Legend
    lg = legend('Data','Fit', 'Zero-crossing'); 
    lg.FontSize =  6; lg.Location = 'southwest';        

    % Fontsizes
    ax_bottom.FontSize = 7;
    ax_bottom.YLabel.FontSize = 9; ax_bottom.XLabel.FontSize = 9;

    % ------------------------------------------ %
    % Top Row Plot 

    % Plot the top-spectra / phased and raw
    % Manual subplots of top-row
    nsubs = size(doi_main,2); sbplot_n = nsubs; sbplot_m = 2;

    % Create x-position for each axis its normalized location using the 
    % positon info from the bottom plot to keep it aligned.
    ax_pos_bottom = ax_bottom.Position;
    ax_dx = ax_bottom.Position(3) ./ nsubs; 
    ax_x_end = (ax_bottom.Position(3) + ax_bottom.Position(1)) - ax_dx; 
    ax_x_vec = linspace(ax_pos_bottom(1), ax_x_end, nsubs);
    
    % Y-location of subplots
    ax_tmp = subplot(sbplot_m, sbplot_n, 1); ax_tmp_pos = ax_tmp.Position;
    delete(ax_tmp); ax_y = ax_tmp_pos(2);
    
    % Size of sub-plots
    ax_sz = [ax_dx ax_pos_bottom(4)];

    % Create ppm-axis if it does not exist.
    if exist('ppm', 'var') ~= 1, ppm = linspace(-4,4,size(ydata,1)); end
    xlim_val = [min(ppm) max(ppm)];

    % Prepare x-tick values
    center = ppm(round(size(ppm,2)./2));
    rb = round( ppm( round(size(ppm,2).*(5/24)) ),1);
    lb = round( ppm( round(size(ppm,2).*(19/24))+1),1);
    xtcks = sort([lb center rb]);
    

    ax = gobjects(1,nsubs);
    for kk = 1:nsubs
        ax(kk) = axes('Position', [ax_x_vec(kk) ax_y ax_sz]);
        plot(ax(kk), ppm, real(doi_main(:,kk)), 'color', orange);
        hold on;
        plot(ax(kk), ppm, real(doi(:,kk)), '--r')

        cf = 0; if add_zero, cf = 1; end
        xlabel(sprintf('%i', xdata(kk+cf))); 
        ylim(ylim_val); xlim(xlim_val);
        

        % Remove x-limit xticks values.
        ax(kk).XTick = xtcks; xtickangle(0);
        ax(kk).FontSize = 7; ax(kk).XDir = 'Reverse';
        ax(kk).XLabel.FontSize = 8;

        % Cosmetics
        if kk == 1
            ax(kk).YLabel.String = 'Amplitude [a.u.]'; 
            ax(kk).YLabel.FontSize = 9;
        end
        if kk ~= 1, ax(kk).YTick = []; end
        if kk == nsubs
            lg = legend('Phased', 'Original');
            lg.FontSize = 6; lg.Location = 'northeast';
            lg.Position = [ax_x_vec(1) (ax_y + ax_sz(2)) - lg.Position(4) ...
                lg.Position(3) lg.Position(4)];
        end
    end

    % Write to file.
    if save_plot
        figname = strjoin([dirstr '\' strfn '.png'],'');
        export_fig(figname{1}, '-transparent','-nocrop','-m1', fh);
    end

    catch err
        warning(err.message);
    end
% ---------------------------------------------------------------------- %

end

function doi_phased = phaseCorrection(doi)
% Phase-correct data of interest (doi) into absorption mode

% MRSI data to cell format
sz = size(doi); N = sz(1);
cell_layout = arrayfun(@ones,...
    ones(1,size(sz(2:end),2)),sz(2:end),'UniformOutput',0);
doi_cell = mat2cell(doi, sz(1), cell_layout{:});

% Apply auto zerophase to each cell spectrum
doi_phased = cellfun(@csi_autoZeroPhase, ...
            doi_cell, ...                              % data
            repmat({1:N}, size(doi_cell)),...          % range
            repmat({2},size(doi_cell)),...             % method
            repmat({0},size(doi_cell)),...             % plot
            'UniformOutput', 0);      
doi_phased = cell2mat(doi_phased); 

function [doi_phased, indmin ]= phaseCorrection_zeroCrossing(doi)    
% Phase-correct all spectra after zero-crossing, reversing the full
% absoprtion mode by 180degrees.

% Find zero-crossing in the data
[~, mxmn] = findZeroCross(max(real(doi)));
indmin = find(mxmn==1); indmax = find(mxmn==0);

% Which spectrum has the lowest maximum signal if all are in absorption?
mnph = max(real(doi)); [~, lowest_max_signal_ind] = min(mnph);
indminPh = lowest_max_signal_ind+1:size(mnph,2);
indmaxPh = 1:size(mnph,2); indmaxPh(indminPh) = [];
if numel(indmin) > numel(indminPh) || isempty(indmin)
    indmin = indminPh; indmax = indmaxPh;
end   

% I) Apply 180 degree phasechange
%    ...to all spectra after the zerocrossing
doi_min_mag = abs(doi(:,indmin)); doi_min_pha = angle(doi(:,indmin));
% Add 180-deg phase change
doi_min_pha_corr = doi_min_pha + pi; 
doi_min_phased = complex(doi_min_mag.*cos(doi_min_pha_corr),...
                         doi_min_mag.*sin(doi_min_pha_corr));

% Phased data                     
doi_phased(:,indmax) = doi(:,indmax); 
doi_phased(:,indmin) = doi_min_phased;

