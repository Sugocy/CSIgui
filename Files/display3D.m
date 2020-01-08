%%%% Description:                            Quickly view ND image arrays.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2018-01
%%% --------------------------------------------------------------------
%%% 
%%% display3D(data, varargin)
%%% 
%%% data        ND data array with widht and height of "image" data on
%%%             first and second index.
%%%
%%% Variable input arguments: 'label', label input;
%%% Accepted labels: 
%%%    'colormap'       Standard matlab colormaps or matrix
%%%    'limit'          Axis limits of colormap to show i.e. contrast
%%%    'tag'            Tag to show @ bottom of image display
%%%    'pos'            Position of display3D on the screen.
%%%                     If (1x2) positioned, if (1x4) sizes and positions
%%%                     expecting [x y w h];
%%%
%%% Control options:
%%%    'Middle mouse'   When hold, mouse movements changes contrast
%%%    'Right mouse'    Revert to default contrast of the current array

% ------- Initiate display 3D
function display3D(data,varargin)
% Display 3D data quickly in a 2D plot with option to navigate through the
% data.

if nargin == 0
    data = rand(25,25).*10;
end

% --- Process input

% Squeeze and reshape data if necessary
data = squeeze(data); sz = size(data);
data = reshape(data, sz(1), sz(2), []); 


% Process any variable options
opt = struct;
if nargin > 1
   inp_arg = {'colormap','limit','tag','pos'};
   for ii = 1:2:size(varargin,2)
       opt_type = lower(varargin{ii});
       if sum(contains(inp_arg,opt_type))
           opt.(opt_type) = varargin{ii+1};
       end
   end
end

% --- Create figure
fh = figure('Name','Display 3D in 2D','Tag', 'Display3D',...
            'NumberTitle','Off','Unit','Pixels',...
            'Menubar','None','Toolbar','None','Color','Black',...
            'WindowScrollWheelFcn', @scrollWheel,... % Scrolling arrays
            'WindowButtonMotionFcn',@hoverMouse,...  % Hovering over axis
            'WindowButtonDownFcn',  @clickDown,...   % Click initiate
            'WindowButtonUpFcn',    @clickUp);       % Click release

% --- Calculate position

% 30 percent of screensize
scrsz = get(0,'screensize');

% Figure size
fig_sz = scrsz(3:4)./3; 

% Check size ratio of data
sz = size(data); 
if sz(1) >= sz(2)
    ratio = sz(2)/sz(1); % R=Y=H/C=X=W
    fig_sz(1) = fig_sz(2)*ratio;
elseif sz(2) > sz(1)
    ratio = sz(1)/sz(2); % C=X=W/R=Y=H
    fig_sz(2) = fig_sz(1)*ratio;
end

% At screen center
fig_ps = (scrsz(3:4)./2) - ([fig_sz(1)  fig_sz(2)]./2); 
fh.Position = [fig_ps fig_sz];

% User Input: Position
if isfield(opt,'pos')
    if size(opt.pos,2) == 4
        fh.Position = opt.pos;
    else
        fh.Position(1:2) = opt.pos;
    end
end


% --- Add UI elements

% 1. AXIS
ax = axes(fh,'Position',[0 0 1 1]);
ax.XTick = []; ax.YTick = []; ax.Box = 'off';
ax.Color = 'Black'; ax.YColor = 'Black'; ax.XColor = 'Black';


% 2. SCROLLBAR
sb = uicontrol(fh,'Style','Slider','Units', 'Normalized',...
                  'Position',[0.05 ax.Position(2) 0.03 ax.Position(4)],...
                  'ForegroundColor','Black','BackgroundColor','Black');              
% Set min and max value of the slider + step size 
sb.Value = 1; sb.Max = size(data,3); sb.Min = 1; 
% Safet for array size 1 slider steps
step = ([1 2]./(sb.Max-1)); if ~isfinite(step), step = [1 1]; end
sb.SliderStep  = step;

% Function to call when using slider
sb.Callback = @SB_control;
sb.Visible = 'off';

% 3. INFO TEXT
txt = uicontrol(fh,'Style','Text','Units','Pixels',...
                'Position',[0 0 60 10],...
                'ForegroundColor','Yellow','BackgroundColor','Black',...
                'HorizontalAlignment','Center','FontSize',7);              
txt.String = sprintf('%i/%i',1,size(data,3));
txt.FontWeight = 'bold';

% 4. TAG TEXT
if isfield(opt,'tag')
    tag_w = 80;
    tag = uicontrol(fh,'Style','Text','Units','Pixels',...
                   'Position',[txt.Position(3) 0 tag_w 10],...
                   'ForegroundColor','Yellow','BackgroundColor','Black',...
                   'HorizontalAlignment','Center','FontSize',7);    
    tag.String = opt.tag;
    tag.FontWeight = 'bold';
else
    tag = [];
end


% EDIT CONTAST
edit_w = 80;
edit = uicontrol(fh,'Style','edit','Units','Pixels',...
               'Position',[0 fh.Position(4)-9 100 10],...
               'ForegroundColor','Yellow','BackgroundColor','Black',...
               'HorizontalAlignment','Center','FontSize',7,'String','TEST');    
edit.String = 'contrast'; edit.FontWeight = 'bold';



           
% Normalize UI elements
fh.Units = 'Normalized'; sb.Units  = 'Pixels';

% Save GUI components and data
gui = guidata(fh); gui.fh = fh; gui.data = data; gui.ax = ax; gui.sb = sb;
gui.txt = txt; gui.opt = opt; gui.tag = tag; gui.edit = edit;
guidata(fh,gui);

% Display data
AX_display(gui);

% ------- Control data display
function AX_display(gui)


% ---- Display Data
sl = round(gui.sb.Value); % Get array# to display from slider
imshow(gui.data(:,:,sl),'Parent',gui.ax);

% ---- Txt
gui.txt.String = sprintf('%i/%i',sl,size(gui.data,3));

% ---- Set any options if given

% Contrast limits
if ~isfield(gui.opt,'limit')
    gui.opt.limit = getContrast(gui.data);
end
if gui.opt.limit(1) >= gui.opt.limit(2)
    gui.opt.limit(2) = gui.opt.limit(1)+1;
end
gui.ax.CLim =  gui.opt.limit;
gui.edit.String = num2str(gui.opt.limit);


% Colormap type
if isfield(gui.opt,'colormap')
    colormap(gui.ax, gui.opt.colormap);
end
    
% ------- Control slider
function SB_control(hObj, ~)
gui = guidata(hObj); AX_display(gui);

% ------- Control scrolling
function scrollWheel(hObj, evt)
% Act on scroll wheel action, moving the slider and thus showing a
% different slice in the data.

% GUI handle
gui = guidata(hObj); 

% If scrolled up, increase scrollbar.
if evt.VerticalScrollCount < 0 
    if gui.sb.Value+1 <= gui.sb.Max
            gui.sb.Value = gui.sb.Value+1; 
    end 
% If scrolled up, decrease scrollbar.
elseif  evt.VerticalScrollCount > 0 
    if gui.sb.Value-1 >= gui.sb.Min
            gui.sb.Value = gui.sb.Value-1; 
    end
end

% Update app-structure
guidata(hObj,gui);

% Display using new slider value
AX_display(gui);

% ------- Mouse Down
function clickDown(hObj, ~)

% Get clicked mouse button.
mouse_button = get(hObj,'selectiontype');

gui  = guidata(hObj);
switch lower(mouse_button)
    case 'normal'
    case 'alt'
        % Default contrast set
        gui.opt.limit = getContrast(gui.data(:,:,gui.sb.Value));
        caxis(gui.ax, gui.opt.limit);
        
        
    case 'extend'
        % Add additional callback to Mouse Motion Fcn.
        gui.fh.WindowButtonMotionFcn = ... 
            @(hObj,evt)( cellfun( @(x) feval(x,hObj,evt),...
           {@(hObj,evt)hoverMouse(hObj,evt),...
            @(hObj,evt)hoverMouse_changeContrast(hObj,evt)}));
end
guidata(hObj,gui);

% ------- Mouse Up
function clickUp(hObj, ~)

% Get clicked mouse button.
mouse_button = get(hObj,'selectiontype');

gui  = guidata(hObj);
switch lower(mouse_button)
    case 'normal'
    case 'alt'
    case 'extend'
        % Set @hoverMouse as the only callback from mouse movements.
        gui.fh.WindowButtonMotionFcn = @hoverMouse;
        % Remove any mouse-position data
        if isfield(gui,'mousePos'), gui = rmfield(gui, 'mousePos'); end
end
% Update GUI data
guidata(hObj,gui);

% ------- Upon hovering if middle mouse is down
function hoverMouse_changeContrast(hObj, ~)

% Get GUI data
gui = guidata(hObj);  
% Get current mouse position
curPos = gui.ax.CurrentPoint; curPos = round(curPos(1,1:2));

% Calculate mouse-position change and edit contrast
if isfield(gui, 'mousePos')
    
    % Calculate change and flip the x-axis change
    % to right/up   --> positive 
    % to left/down  --> negative
    change = (gui.mousePos(1,1:2) - curPos(1,1:2)) .* [-1 1]; 
    
    % Mouse change over axis in percentages i.e. "voxel change"
    sz = size(gui.data); changeSpatial = (abs(change) ./ sz(1:2)).*100;
    % No change if only 10% spatial movement of mouse
    change(changeSpatial <= 10) = 0;
    
    % Converty to a binary change
    change(change > 0) = 1; change(change < 0) = -1;
    

    % Calculate new contrast
    currSlice    = gui.data(:,:,gui.sb.Value);
    contrastStep = (abs( max(currSlice(:))) - abs(min(currSlice(:)))) / 256;
    contrastNew  = gui.ax.CLim + (contrastStep .* change);
    
    % Axis limit safety check
    if contrastNew(1) >= contrastNew(2), contrastNew = gui.ax.CLim; end

    % Set contrast and save it in the options
    gui.ax.CLim = contrastNew; gui.opt.limit = contrastNew;
    gui.edit.String = num2str(gui.opt.limit);
    
else
    gui.mousePos = curPos;
end

% Update GUI data
guidata(hObj,gui);

% ------- Control mouse hover
function hoverMouse(hObj, ~)
% Display the value of the pixel which the mouse hovers on in titlebar of
% the figure.

% Get data handle
gui = guidata(hObj);  
C = gui.ax.CurrentPoint;

% Only if mouse over axis only == positive values
if ( round(C(1,1))> 0 ) && ( round(C(1,2)) >0 )

    % X and Y coordinate
    xc = (round(C(1,1))); yc = (round(C(1,2)));
    
    % Get image pixel data
    % Reverse X and Y coordinate for Row/Column discrepancy in Matlab!
    if   (size(gui.data,2) >= xc) && (size(gui.data,1) >= yc)
        px = gui.data(yc, xc, gui.sb.Value);
    else
        gui.fh.Name = 'Display 3D in 2D';
        return;
    end

    text_str = sprintf('Row: %i Col: %i - Value: %5.5f', ...
                        round(C(1,2)), round(C(1,1)), px);

    % Set name of figure title
    gui.fh.Name = text_str;
else
    gui.fh.Name = 'Display 3D in 2D';
end


guidata(hObj,gui);

% ----- % Calculate contrast
function contrast = getContrast(data)

% Calculate the contrast once.
contrast = double( [min(data(:)) max(data(:))] ); 

% Safety check contrast
if contrast(1) >= contrast(2), contrast(2) = contrast(1)+1; end
