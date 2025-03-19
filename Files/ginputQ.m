function pos = ginputQ(N, fh)
% ginput made by Quincy. Beta.

% Get current fig if no obj given.
if nargin ==1, fh = gcf; end

cursor_on = 1;

% Get gui-structure, add number of positions
data = guidata(fh); data.N = N;
data.fh = fh; 
% Find axis-obj in figure
kids = fh.Children;
for kk = 1:size(kids,1)
    if strcmp(kids(kk).Type,'axes'), ind = kk; break; end
end
% Get data object
data.axh = fh.Children(ind);
% Set title to show N-clicked/N;
old.title = data.axh.Title; title(data.axh,sprintf('%i/%i',0,data.N));

% --- Set Callback: motion and mouse down

% Store old fcn
old.fcn_button = data.axh.ButtonDownFcn; 
old.fcn_motion = data.fh.WindowButtonMotionFcn;

% Apply motion
if cursor_on, data.fh.WindowButtonMotionFcn = @mouse_hover; end

% Apply mouse-down
data.axh.ButtonDownFcn = @mouse_down; 
% Turn hittest for all line children off 
% This enables click-through! (in axis)
for kk = 1:size(data.axh.Children,1),data.axh.Children(kk).HitTest = 'off';
end


% Axis unit
old.unit = data.axh.Units;data.axh.Units = 'Pixels'; 

% Update GUI-data
guidata(fh,data); pos = NaN; 

                    % --- % Pause until user is done % --- %
% Wait here             
uiwait(data.fh);
                    % --- % Pause until user is done % --- %

% Get updated gui-data
if ~ishandle(data.fh),return; end % If user closed figure before finishing
data = guidata(data.fh);

% Set created output
pos = data.pos;

% Set default values
data.axh.ButtonDownFcn = old.fcn_button; 
data.fh.WindowButtonMotionFcn = old.fcn_motion;
data.axh.Title = old.title;

% Turn hittest of axis children on
for kk = 1:size(data.axh.Children,1),data.axh.Children(kk).HitTest = 'on';
end


function mouse_down(hObj, evt)

% Get gui structure data
data = guidata(hObj);

% Save position of mouse down
if isfield(data,'pos'), data.pos(end+1,:) = evt.IntersectionPoint;
else,                   data.pos = evt.IntersectionPoint;
end

% Break if requested #mouse_down is reached.
if size(data.pos,1) >= data.N, uiresume(data.fh);  end

% Update title
title(data.axh,sprintf('%i/%i',size(data.pos,1),data.N));

% Update gui structure data
guidata(hObj,data);

function mouse_hover(hObj, evt)

% Mouse position
mouse_pos = evt.Source.CurrentPoint;

% Guidata
data = guidata(hObj);

% Axis location

ax_pos = data.axh.Position;
ax_pos_inner = data.axh.InnerPosition;

% Inside axis-box
if mouse_pos(1) >= ax_pos(1) &&   mouse_pos(1) <= (ax_pos_inner(3)+ ax_pos(1))
 if mouse_pos(2) >= ax_pos(2) &&  mouse_pos(2) <= (ax_pos_inner(4)+ ax_pos(2))
    
    hold(data.axh,'on'); 
    if isfield(data, 'mouse_cursor')
        for kk = 1:size(data.mouse_cursor,2)
            delete(data.mouse_cursor{kk});
        end
    end
    
    curr_final = 12; curr_sz = (curr_final/2);
    curr_offset = 6;
    
% Vertical

    data.mouse_cursor{1} = uicontrol('Style','Text',...
        'Position',...
        [mouse_pos(1) mouse_pos(2)-curr_sz-curr_offset 2 curr_sz],...
        'HandleVisibility', 'off', 'HitTest', 'off',...
        'BackGroundColor','Blue');

    data.mouse_cursor{2} = uicontrol('Style','Text',...
        'Position',[mouse_pos(1) mouse_pos(2)+curr_offset 2 curr_sz],...
        'HandleVisibility', 'off', 'HitTest', 'off',...
        'BackGroundColor','Blue');

    % Horizontal

    data.mouse_cursor{3} = uicontrol('Style','Text',...
        'Position',...
        [mouse_pos(1)-curr_sz-curr_offset mouse_pos(2)-1.5 curr_sz 2],...
        'HandleVisibility', 'off', 'HitTest', 'off',...
        'BackGroundColor','Blue');
    
    data.mouse_cursor{4} = uicontrol('Style','Text',...
        'Position',...
        [mouse_pos(1)+curr_offset mouse_pos(2)-1.5 curr_sz 2 ],...
        'HandleVisibility', 'off', 'HitTest', 'off',...
        'BackGroundColor','Blue');
   
 end
end

guidata(hObj,data);






