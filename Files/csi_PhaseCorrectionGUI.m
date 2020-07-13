function varargout = csi_PhaseCorrectionGUI(varargin)
%%%% Description:   GUI to apply first and zero order phasing corrections.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-12
%%% --------------------------------------------------------------------
%%% 
%%% [spec, zero_ph, first_ph] = csi_PhaseCorrectionGUI(...);
%%%
%%% INPUT:
%%% 1. Spectral data (1xN).
%%% 2. Spectral data (1xN) & PPM axis (1xN).
%%% 3. Spectral data (1xN) + Struct with fields:
%%%                          bw, TE, freq, offset (optional).
%%%
%%% bw        = bandwidth
%%% TE        = echotime
%%% freq_tx   = Transmit frequency (Proton/Phosphorus etc.)
%%% ppm-shift = Offset of the ppm axis (in ppm!)
%%%
%%% OUTPUT:
%%% 1. Spec   = phase corrected spectrum (1xN)
%%% 2. Zero
%%% 3. First
%%%
%%% If the structure is present with above mentioned fields, an initial
%%% zero-order phasing will be calculated and applied at start-up.

% TBD
% 1. Process input - Check (19/12/17)
% 2. Show degrees of phase change in GUI somewhere.
% 3. Show pivot peak ppm location in GUI somewhere.
% 4. Get phase corrections as seperate output.

% INPUT analysis ======================================================== %             
% ======================================================================= %

mockup = 0; % This boolean is used to create a mockup x-axis if not given.
switch nargin
    case 0
        fprintf('No input arguments csi_PhaseCorrection, returning.\n'); 
        return;
    case 1 % Spectrum only.
        spec = varargin{1}; 
        mockup = 1;
    case 2 % PPM-axis OR structure.
        spec = varargin{1}; 
        
        % STRUCT: BW/TE/IMFREQ
        if isstruct(varargin{2})
            x = varargin{2};
            
            % Boolean for required fields.
            xbool = isfield(x, {'bw','TE','freq','offset'});
            if sum(xbool(1:3)) ~= 3
                % missing required fields.
                fprintf('Missing required fields: bw, TE or freq. \n');
                mockup = 1;
            end
            
            % Optional offset of ppm-axis.
            if xbool(4) == 0, x.offset = 0; end
        else
            x.ppm = varargin{2}; x.N = size(spec,1);
            mockup = 2; % Skip mockup parts.
        end
end

% INPUT Processing ====================================================== %
% ======================================================================= %

if mockup == 0
    % X-Axis (PPM/Hz/Sec)
    x.N = size(spec,1);
    x.dwell = 1/x.bw; x.acqT  = x.dwell * x.N; x.freq_tx = x.freq * 10^-6;
    % Calculations
    x.time = (0 : (x.N)-1) .* x.dwell;           % N -> Time
    x.hz   = ( -x.N/2 : (x.N/2) - 1) ./ x.acqT;  % N -> HZ
    x.ppm  = (x.hz / x.freq_tx) - x.offset;      % HZ -> PPM
    fprintf('Proceeding with calculated ppm-axis.\n');
elseif mockup == 1 % Create mockup x-axis for data.
    x.N = size(spec,1); x.ppm = 1:x.N;
    fprintf('Proceeding with unitless x-axis.\n');
elseif mockup == 2
    % Nothing!
    % x.N & x.ppm already exist;
    fprintf('Proceeding with given x-axis.\n');
end


% Initiate GUI ========================================================== %
% ======================================================================= %

% Create figure
W = 680; H = 520;
figh = figure('Toolbar', 'None', 'Menubar', 'None','NumberTitle','Off',...
       'Color', 'Black', 'Name', 'csi_PhaseCorrectionGUI',...
       'Tag','csi_PhaseCorrectionGUI',...
       'WindowScrollWheelFcn', @scrollAdjust,...
       'WindowButtonDownFcn',  @test, 'CloseRequestFcn',@closeReqFcn);
% Test doesnt do anything ATM: click on axis to set pivot peak. In
% development.
       
% Figure size and position
pc_res = get(0,'screensize'); fig_pos = [pc_res(3)/3 pc_res(4)/3 W H];
figh.Position = fig_pos;

% GUI-data struct from figure
gdat = guidata(figh); gdat.figh = figh;

% Apply menu to figure.
menu.File = uimenu(figh,'Label','File');
menu.Save = uimenu(menu.File, 'Label', 'Save as..','Accelerator', 'S',...
                              'Callback', @saveSpec);                          
menu.Close = uimenu(menu.File,'Label','Close','Callback', @closeReqFcn,...
                              'Separator', 'on', 'Accelerator', 'Q');                          
% Store menu in guidata.
gdat.menu  = menu;

   
% Add axes
gdat.ax = axes(gdat.figh,'Color','Black',...
                         'YColor',[0.60 0 0.25],'XColor',[0.60 0 0.25],...
                         'Xdir','Reverse'); hold on;
gdat.ax.Position = [0.125 0.2 0.75 0.70]; % [X,Y,W,H]

% Add slider
% Initialized for zero order phasing
gdat.sl = uicontrol('Style','Slider','BackGroundColor','Black',...
                    'Units','Normalized','Callback',@sliderChange,...
                    'Value', 0, 'min', -180, 'max', 180,...
                    'SliderStep',[1/(360) 1/(36)]); 
gdat.sl.Position = [0.4 0.075 0.48 0.04]; % [X,Y,W,H]

% Add radio buttons
% To select either first or zero order correction with srollbar.
% or to select/change the pivot peak.
rwh           = [0.08 0.025 0.085]; % Radio width and height and yposition.
rx1           = 0.125;              % X-Position of first radio button
radio_clr_txt = [0.35 0.35 0.35]; 
gdat.radio_Zero  = ...
    uicontrol(figh,'Style','radiobutton','Tag','radioZero',...
             'Unit','Normalized','Position',[rx1 rwh(3) rwh(1) rwh(2)],...
             'String','Zero','BackGroundColor', 'Black',...
             'ForeGroundColor', radio_clr_txt,'FontWeight','Bold',...
             'FontSize', 10, 'Value', 1,'Callback', @radioSwitch);
gdat.radio_First  = ...
    uicontrol(figh,'Style','radiobutton','Tag','radioFirst',...
             'Unit','Normalized','Position',[rx1+rwh(1) rwh(3)  rwh(1) rwh(2)],...
             'String','First','BackGroundColor', 'Black',...
             'ForeGroundColor', radio_clr_txt,'FontWeight','Bold',...
             'FontSize', 10, 'Value', 0,'Callback', @radioSwitch);
gdat.radio_Pivot  = ...
    uicontrol(figh,'Style','radiobutton','Tag','radioPivot',...
             'Unit','Normalized','Position',[rx1+(rwh(1)*2) rwh(3) rwh(1) rwh(2)],...
             'String','Pivot','BackGroundColor', 'Black',...
             'ForeGroundColor', radio_clr_txt,'FontWeight','Bold',...
             'FontSize', 10, 'Value', 0, 'Callback', @radioSwitch);
             
         
% Store input in GUI ==================================================== %        
% ======================================================================= %

% SAVE original X-data struct;
gdat.x = x; gdat.spec = spec;
% SAVE original;
gdat.ori.spec = gdat.spec;


% INITIAL PARAMETERS ==================================================== %
% ======================================================================= %
% Zero, first and pivot peak values.

% Zero order
% Formula(BW/ndimf) * TE * 2Pi; only if parameters are given.
if mockup == 0
    gdat.phase.zero  = x.bw./x.N * (x.TE/1000) * (2*pi);
    gdat.sl.Value    = gdat.phase.zero;
    % Show user initial phase change
    fprintf('Applied initial 0th order phase correction: %03.3f rad.\n',...
            gdat.phase.zero);
else
    gdat.phase.zero = 0; gdat.sl.Value = gdat.phase.zero;
end
    
% First order
gdat.phase.first = zeros(x.N, 1); gdat.phase.first_slidervalue = 0;

% Pivot peak
gdat.pivot.position   = (x.N/2)+1;
gdat.pivot.slider_val = gdat.pivot.position;

% Save data and plot spectrum
guidata(gdat.figh, gdat);  plotSpectrum(gdat);


% OUTPUT + CLOSEDOWN ==================================================== %
% ======================================================================= %
% Wait for user to close down GUI: closeReqFcn.
uiwait(gdat.figh); % The closeReqFcn resumes this function.

% Set output(s)
if nargout
    % Reload gui-data and apply phase changes
    gdat = guidata(gdat.figh); 
    
    % Apply phase change to spectrum, output this spectrum.
    spec = applyPhase(gdat); 
    
    % Set as output.
    switch nargout 
        case 1, varargout{1} = spec;
        case 2, varargout{1} = spec; 
                varargout{2} = gdat.phase.zero; 
        case 3, varargout{1} = spec;
                varargout{2} = gdat.phase.zero;
                varargout{3} = gdat.phase.first;
    end
    
    % Close down: Delete figure.
    delete(gdat.figh);
end

% Close down: Delete figure.
delete(gdat.figh); % Close always if UI resumes. (See closeReq)

function closeReqFcn(hObj, ~)

try % FAST
    % Resume UI, see initialisation, to create output and close figure.
    gdat = guidata(hObj); uiresume(gdat.figh);
catch % SLOW
   fobj = findobj('Type','Figure','Tag','csi_PhaseCorrectionGUI');
   if ~isempty(fobj), uiresume(fobj); delete(fobj); end
end





function plotSpectrum(gdat)
% 1. Apply phase changes to spectrum
% 2. Calc x and y limits
% 3. Plot zero-signal line
% 4. Plot spectrum and edit axis
% 5. Calc pivot peak
% 6. Plot pivot peak
% 7. Flush java.

% Colors to plot
plotclr = 'm'; pivotclr = 'y';

% 1. Apply phasing.
% 2. gdat.spec will contain the spectrum with the corrected phases.
applyPhase(gdat); gdat = guidata(gdat.figh);


% 2. Y LIMIT
ymax = (max(real(gdat.spec))); 
ymin = (min(real(gdat.spec)));
if      abs(ymin) > abs(ymax)
    gdat.ylim = [-1*abs(ymin)*1.05 abs(ymin)*1.05];
elseif  abs(ymin) < abs(ymax)
    gdat.ylim = [-1*abs(ymax)*1.05 abs(ymax)*1.05];
end

% 3. X LIMIT
if ~isfield(gdat,'xlim')
    xmax = max(gdat.x.ppm); xmin = min(gdat.x.ppm);
    gdat.xlim = [xmin*1.01 xmax*1.01];
end

% 3. PLOT AND AXIS COSMETICS
% a. Showing the zero-axis to user.
hold(gdat.ax,'off'); % By setting this to on ...
ydat_zero_ax = zeros(size(gdat.x.ppm,2),1);
gdat.plot.xzero = plot(gdat.ax, gdat.x.ppm, ydat_zero_ax,...
    '-','Color', [0.35 0.35 0.35]);
hold(gdat.ax,'on');  % ... the axis are cleared.

% b. Spectrum
gdat.plot.spec = plot(gdat.ax, gdat.x.ppm, real(gdat.spec), plotclr);
set(gdat.ax,'Color','Black','YColor',[0.60 0 0.25],...
                            'XColor',[0.60 0 0.25],...
                            'Xdir','Reverse','Box', 'Off');
% Apply limits.

gdat.ax.YLim = gdat.ylim; gdat.ax.XLim = gdat.xlim;

% c. Pivot calculation
guidata(gdat.figh, gdat); pivot_peak_calc(gdat.figh); 
gdat = guidata(gdat.figh);

% d. Pivot plot
gdat.plot.pivot = ...
    plot(gdat.ax,gdat.pivot.line(:,1), gdat.pivot.line(:,2), pivotclr);

% e. Java flush
vr = version('-release'); vr(end) = []; vr = str2double(vr);
if vr >= 2016, drawnow limitrate;
else,           drawnow expose;
end
   
% // END update.
guidata(gdat.figh, gdat);

function spec = applyPhase(gdat)
% Sets the spectrum corrected according phasing options set by user 
% in gdat.spec;

% Get zero and first order phase corrections.
phcor0 = gdat.phase.zero; phcor1 = gdat.phase.first;

% DEV
% figure(3); plot(gdat.x.ppm, phcor1); ylim([-2*pi 2*pi]);
% set(gca,'Xdir', 'Reverse');

% Get original data.
spec = gdat.ori.spec;

% Get phase and magnitude.
phase = angle(spec); magn = abs(spec);

% Apply phase corrections.
new_phase = phase + phcor0 + phcor1;

% %%%%%%%%%% Calculate new complex spectrum. %%%%%%%%%%%
spec = complex(magn.*cos(new_phase), magn.*sin(new_phase));
gdat.spec = spec;

% Update guidata.
guidata(gdat.figh,gdat);

function zero_order(hObj)

% 1. Get appdata struct with spectrum/fid data
gdat       = guidata(hObj);    % get GUI data
% 2. Get slider value --> degree of change
slider_val = gdat.sl.Value;    % get slider value

% 3. Phase correction using slider value
phase_cor_zero = (slider_val * (pi/180)); % Convert to radians.

% 4. Display output
% fprintf('Applied 0th Phase correction factor: %03.3f degrees, %03.3f radians.\n',...
% phase_cor_zero*(180/pi), phase_cor_zero);

% 5. Store spectrum.
gdat.phase.zero = phase_cor_zero;

% 6. Update GUI data and PLOT
guidata(hObj, gdat); plotSpectrum(gdat)

function first_order(hObj)

% 1. Get appdata struct with spectrum/fid data
gdat       = guidata(hObj);    % get GUI data
slider_val = gdat.sl.Value;    % get slider value


% 2. Set up the linear phase change over the full spectrum
% i. Get x-axis data
x = gdat.x;

% ii. Add -dPhase to +dPhase linearly. dPhase is set by user between 0 and
% 2pi.
slope = ( 1 * slider_val * ((pi/180))) /(x.N) ;                             % Lin. Slope dPhase/(N)

% Container for range over spectrum.
range = zeros(1,x.N); 

% Get pivot peak dependent x-axis over the entire spectrum range
% This is the x-axis value in ppm.
piv_pos = gdat.pivot.position; 

% Create range accordingly.
if piv_pos ~= 0
    range(piv_pos:x.N) = 0:(x.N - piv_pos);
    range(1:piv_pos-1) = (-1*(x.N - (x.N - piv_pos))) +1 : 1 :-1;
else
    range = 0:x.N-1;
end
   

% Calculate the phase change over the spectrum range.
phase_cor_range = range.*slope; % Range over spectrum.
                                        

% 3. Save phase correction range
gdat.phase.first = phase_cor_range';
gdat.phase.first_slidervalue = slider_val;

% 4. Display output
% fprintf('Applied 1st phase correction factor: %03.2f degrees, %2.2f radians.\n',...
% phase_cor_range(end)*(180/pi), phase_cor_range(end));
% fprintf('Applied 1st phase correction slope:  %03.2f degrees, %2.2f radians.\n',...
% ((slider_val)/(x.N/2)), ((slider_val * ((pi/180)))/(x.N/2)));

% 5. Update GUI data and PLOT
guidata(hObj, gdat); plotSpectrum(gdat);

function pivot_peak_set(hObj)
% Sets the position of the pivot peak line.

% 1. Get appdata struct with spectrum/fid data
gdat       = guidata(hObj);    % get GUI data
slider_val = floor(gdat.sl.Value);

% 2. Set slider value as new pivot position.
gdat.pivot.position   = abs(slider_val - gdat.x.N);
gdat.pivot.slider_val = gdat.sl.Value;
% Here the values of the slider are reversed as the plotted x-axis is
% reversed. E.g. moving the slider to left means moving pivot peak line to 
% the left and vice versa.

% 3. Update GUI data and PLOT: this calls the pivot_peak_calculation and
% draws the newly set pivot peak :)
guidata(hObj, gdat); plotSpectrum(gdat);

function pivot_peak_calc(hObj)
% Called by plotSpectrum to calculate a vertical line representing the
% pivot peak line set by the user.

gdat = guidata(hObj); % Appdata
res  = 2;             % Resolution of pivot line.

% Pivot position on the axis thus 1:size(x-axis);
pos = gdat.pivot.position; % This is 0*slope in first order phasing
pos = round(pos); if pos == 0, pos = 1; end
xvalue = gdat.x.ppm(pos);

% Limits of the Y-axis and its range.
lim = [gdat.ax.YLim]; ydat = linspace(lim(1), lim(2), res)';
% Create vertical line data for x.
xdat = repmat(xvalue,[res 1]); 

% Set and save
gdat.pivot.line = cat(2,xdat, ydat); guidata(hObj, gdat);

function sliderChange(hObj, ~)
% Change of slider position will execute correct function according set
% radio-boxes.

gdat = guidata(hObj);

if gdat.radio_Zero.Value == 1
    zero_order(hObj);
elseif gdat.radio_First.Value == 1
    first_order(hObj);
elseif gdat.radio_Pivot.Value == 1
    pivot_peak_set(hObj);
end

function sliderSwitch(hObj)
% If radio boxes are changed, this function is called to set correct
% slider-options for its set operation.
gdat = guidata(hObj);

if gdat.radio_Zero.Value == 1
    % Settings ZERO ORDER PHASE SLIDER
    gdat.sl.Min = -360; gdat.sl.Max = 360;
    gdat.sl.SliderStep = [1/(720) 1/(72)];
    
    gdat.sl.Value = gdat.phase.zero .* (180/pi); 
elseif gdat.radio_First.Value == 1
    % Settings FIRST ORDER PHASE SLIDER
    gdat.sl.Min = -7200; gdat.sl.Max = 7200;
    gdat.sl.SliderStep = [1/(14400*2) 1/(14400*2)];
    
    gdat.sl.Value = gdat.phase.first_slidervalue;
elseif gdat.radio_Pivot.Value == 1
    % Settings PIVOT SLIDER
    gdat.sl.Min = 1; gdat.sl.Max = gdat.x.N;
    gdat.sl.SliderStep = [1/(gdat.x.N) 1/(gdat.x.N/10)];
    
    gdat.sl.Value = gdat.pivot.slider_val;
    %NB: the pivot peak values are reversed when calculated such that
    %shifting the slider to the left makes the pivot peak shift to the left
    %and vice versa.
end

% Update.
guidata(hObj, gdat);

function radioSwitch(hObj,evnt)
gdat = guidata(hObj);
tag = evnt.Source.Tag;
switch tag
    case 'radioZero'
        if gdat.radio_Zero.Value == 1 
            gdat.radio_First.Value = 0;
            gdat.radio_Pivot.Value = 0;
        elseif gdat.radio_Zero.Value == 0
            gdat.radio_First.Value = 1;
            gdat.radio_Pivot.Value = 0;
        end
    case 'radioFirst'
        if gdat.radio_First.Value == 1
            gdat.radio_Zero.Value = 0;
            gdat.radio_Pivot.Value = 0;
        elseif gdat.radio_First.Value == 0
            gdat.radio_Zero.Value = 1;
            gdat.radio_Pivot.Value = 0;           
        end
    case 'radioPivot'
        if gdat.radio_Pivot.Value == 1
            gdat.radio_Zero.Value = 0; 
            gdat.radio_First.Value = 0;
            
        elseif gdat.radio_Pivot.Value == 0 
            gdat.radio_Zero.Value = 0;
            gdat.radio_First.Value = 1;
        end
end
guidata(hObj,gdat);

% Update the slider accordingly!
sliderSwitch(hObj)

function scrollAdjust(hObj, evnt)
% Applies correct action upon mouse-scroll-wheel change.

gdat = guidata(hObj);
scrInp = evnt.VerticalScrollCount;

% Get current functionality of GUI
currFunc = [gdat.radio_Zero.Value, gdat.radio_First.Value,...
            gdat.radio_Pivot.Value];
% Get current slider value
currSli = gdat.sl.Value;
if currFunc(2) == 1 % first order
    newVal  = currSli + (scrInp*20);
else
    newVal  = currSli + scrInp;
end

% Dont surpass min and max values of slider.
if  (newVal > gdat.sl.Max) ||  (newVal < gdat.sl.Min), return; end

% Save the new value.
gdat.sl.Value = newVal; guidata(hObj,gdat);

if     currFunc(1) == 1 % Zero order
    zero_order(hObj);
elseif currFunc(2) == 1 % First order
    first_order(hObj);
elseif currFunc(3) == 1 % Pivot peak.
    pivot_peak_set(hObj);
end

function test(hObj, evnt)
% [x,y] = ginput(1);
% 
% disp('w');
% Something...Scaling if figure is resized so using click pos is dangerous.
% uicontrol('Style','Text','String', 'TEST', 'BackGroundColor', 'red',...   
%          'Position', [click_pos(1) click_pos(2) 10 10] );
%
%
%
%click_pos = evnt.Source.CurrentPoint

function saveSpec(hObj, ~)

gdat = guidata(hObj);
applyPhase(gdat); gdat = guidata(hObj);

% Get user to designate path and name.
[fn,fp] = uiputfile({'*.sdat','SDAT file';'*.txt','Text file'},'Save as');
if fn == 0, return; end
[~, ~, ext] = fileparts(fn);

switch ext
    case '.txt'
        % Write real and imaginary to txt file.
        repart = real(gdat.spec);
        impart = imag(gdat.spec);
        tosave = cat(2,repart, impart);
        
        % Create file
        txtfid = fopen([fp fn], 'w+');
        for kk = 1:size(tosave,1)
            fprintf(txtfid,'%f %f\n', tosave(kk,:));
        end
        
    case '.sdat'
        % Write as the sdat file
        csi_writeSDAT([fp fn], gdat.spec, 'spec');
end
