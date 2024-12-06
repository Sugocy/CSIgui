function userInput = getInput(e_name, e_query, e_input, fig_title)
% Get userinput using some input elements. Available input-elements are: 
% edit or popup.
%
% e_name is the element name, edit or popup.
% e_query is the query or string for the element.
% e_input is the default input seperate per element as a cell-array.
%           popup e_input: {{'1','2','3'},{'Banana','Apple','Lime'}};
%            edit e_input: {{'5'}, {'Lime'}}
% fig_title is the figure its name shown in the window-bar.
%
% userInput is empty when user skips or closes the window.
%
%
% Quincy van Houtum, PhD. 2023/10
% quincyvanhoutum@gmail.com

% Input processing
if nargin < 4, fig_title = 'Input'; end

% Colors of the GUI
clrs = [0.000 0.000 0.000; 0.941 0.941 0.941];

% Create sub-figure
fig_userinp = figure('ToolBar','None','Menubar','None','Unit','Pixels',...
                     'NumberTitle', 'Off', 'Color', [0 0 0],...
                     'Name', fig_title,'Tag', 'getInput'); 
subdat = guidata(fig_userinp); axis off;

% Nr of dropdown menu's
Ninput = size(e_query,2);

% Total size of figure.
dw = 10;  dy = 10; w = 320+2*dw; 
% Define values for EDIT and TEXT sizes.
sz_popup  = [w-2*dw 25]; sz_text = [w-2*dw 20]; 
sz_button = [(w-2*dw)/2 20];

% Define height
h = Ninput * (sz_popup(2)*2) + (sz_button(2)*2) + 2*dy;

% Size tune parameters and set figure size and window-position
scrsz = get(0,'Screensize'); 
figpos = round(scrsz(3:4)./2) - ([w h]./2);
set(fig_userinp,'Position', [figpos w h]);

% --------------------- % Set UI control elements

% Set each title text and popup menu
for kk = 1:Ninput
    %Popup menu title
    subdat.h_title{kk} = ...
        uicontrol(fig_userinp, 'Style', 'Text','Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_popup(2)*kk + sz_text(2)*(kk-1)) ...
                    sz_text(1) sz_text(2)],...
        'Tag', 'h_popup','String',e_query{kk},'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left',...
        'TooltipString', e_query{kk});
    
    % Dropdown menu title    
    switch e_name{kk}
        case 'popup', elem = 'popupmenu'; 
        case 'edit', elem = 'edit'; 
    end

    subdat.h_elem{kk} = ...
        uicontrol(fig_userinp, 'Style', elem,'Unit', 'Pixels',...
        'Position',[dw h-dy-(sz_popup(2)*kk + sz_text(2)*(kk))...
                    sz_popup(1) sz_popup(2)],...
        'Tag', ['h_' elem], 'String', e_input{kk}, 'Fontsize', 8,...
        'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
        'HorizontalAlignment', 'Left',...
        'KeyPressFcn', @do_KeyPress);   
end


% Set button CONTINUE
subdat.h_buttonContinue = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw        dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_CONTINUE','String','Continue','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @setOutput);

% Set button SKIP
subdat.h_buttonSkip = uicontrol(fig_userinp, 'Style', 'Pushbutton','Unit', 'Pixels',...
    'Position',[dw+sz_button(1)  dy/2  sz_button(1) sz_button(2) ],...
    'Tag', 'button_SKIP','String','Skip','Fontsize', 12,...
    'Foregroundcolor', clrs(2,:),'BackgroundColor', clrs(1,:),...
    'Callback', @setOutput);

% --------------------- % Wait for user input

% Update gui-data of this function' figure.
guidata(fig_userinp, subdat); 
% Wait for user t2bdone
uiwait(fig_userinp); pause(0.1);

% --------------------- % Process user input

% Get updated gui-data
if ishandle(fig_userinp)
    subdat = guidata(fig_userinp);
    if isfield(subdat, 'ans')
        userInput = subdat.ans;
    else  
        userInput = [];
    end
    % Close figure
    close(fig_userinp);
else
    userInput = [];
end


% -. Executed by button in getUserInput_Popup to save answer
function setOutput(hObject,~,~)
% When user closes, skips or continues the getUserInput dlg.

userInputString = get(hObject,'String'); % What did get us here? 

% Get guidata from the getUserInput window
subobj = findobj('type','figure','tag','getInput');
if iscell(subobj), subobj = subobj{1}; end
subdat = guidata(subobj);

if strcmp(userInputString, 'Skip')
    val = {};
else
    % Get strings input from popup menu
    val = cell(1,size(subdat.h_elem,2));
    for kk = 1:size(subdat.h_elem,2)
        pval = get(subdat.h_elem{kk},'Value');
        pstr = get(subdat.h_elem{kk},'String');
        if pval == 0
            val{kk} = pstr;
        else
            val{kk} = pstr{pval};
        end
    end
end
% Set userinput value and Update storage
subdat.ans = val; guidata(subobj,subdat);

% Close figure - getUserInput script resumes
uiresume(subobj);



% -. Executed by button press in edits
function do_KeyPress(H, E) 
drawnow;
if strcmp(E.Key, 'return')
    setOutput(H, E);    
else
    return;
end
