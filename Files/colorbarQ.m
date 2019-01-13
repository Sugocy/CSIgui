function colorbarQ(colors, values)

% Create figure
figh = figure('Tag', 'ColorBar','Name', '',...
       'Color', 'Black','Toolbar', 'None', 'MenuBar', 'None',...
       'NumberTitle', 'Off','Resize','on');  

% Set figure position
w = 320; h = 25;
scrsz = get(0,'Screensize'); 
figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(figh,'Position', [figpos w h]);   

% GUI data
bgui = guidata(figh); bgui.fig = figh;

nText = size(colors,1);
nText_w = w/nText; nText_h = h.*1.25;
nText_x = nText_w; nText_w = w/nText.*1.25;

bgui.ax = axes(bgui.fig, 'units','Pixels','Position',[0 0 w h]);
bgui.ax.XLim = [0 nText];

for kk =1:nText
    
 % Set a (text) bar
 bgui.bar{kk} = uicontrol('style','text', ...
            'units','pixels',...
            'position',[nText_x*(kk-1) 0 nText_w nText_h],...
            'BackgroundColor', colors(kk,:));
  bgui.bar{kk}.Units = 'normalized';
end   
    
set(figh,'WindowButtonMotionFcn',@mouseHover);

if size(values,1) == 1, values = values'; end
bgui.values = values;

guidata(figh, bgui)


function mouseHover(hObj,~)

% Get data handle
gui = guidata(hObj); C = gui.ax.CurrentPoint;

% Get value of mouse position
n = ceil(C(1));
if n <= size(gui.values,1) && n > 0
    val = gui.values(n) ;
    hObj.Name = (sprintf('%4.4E', val));
end



