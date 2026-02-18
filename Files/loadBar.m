function barObj = loadBar(perc, info, barObj)
% Opens a load bar if not existant and sets the load-bar percentage to
% perc. If perc equals NaN, the bar is closed.
%
% Input info, a string, is optional and can be used to set text in the
% window of the loadbar. 

% If no info string is given
if nargin == 1, info = 'Busy...'; end

% Normalize percentage
if ~isnan(perc) && perc > 1, perc = perc/100; end

% Check for barObject if necessary
if nargin < 3    
    % Check for bar-figure
    barObj = findobj('Type','Figure','Tag','loadBar');
elseif ~ishandle(barObj)
    barObj = findobj('Type','Figure','Tag','loadBar');
end

                % ------- % Create loadBar % ------- %
if isempty(barObj)
    % Create figure
    figh = figure('Tag', 'loadBar', 'Name', info,...
           'Color', 'Black','Toolbar', 'None', 'MenuBar', 'None',...
           'NumberTitle', 'Off','Resize','off');  
      
    % Set figure position
    w = 480; h = 20;
    scrsz = get(0,'Screensize'); 
    figpos = round(scrsz(3:4)./2) - ([w h]./2);
    set(figh,'Position', [figpos w h]);
    
    % GUI data
    bgui = guidata(figh); bgui.fig = figh;
    
    % Set a (text) bar
    bgui.bar = uicontrol('style','text', ...
                'units','normalized','position',[0 0 0 1],...
                'BackgroundColor', [0.9 0.0 0.0]);
            
    % Save guidata
    guidata(bgui.fig, bgui);
else
    % Get loadbar guidata
    bgui = guidata(barObj);
end

                % ------- % UPDATING loadBAR % ------- %

% Close if NaN
if isnan(perc), delete(bgui.fig); return; end

% Set bar progress
bgui.bar.Position = [0 0 perc 1];

% Set info string
bgui.fig.Name = ['CSIgui: ' info];

% Flush
drawnow;                            
           