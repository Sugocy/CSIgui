function auto_click

% fprintf('Running!\n')

N = 0; boi = {'CSIgui_UI_popup', 'CSIgui_UI_edits', 'CSIgui_UI_tick', 'CSIgui_UI_button'}; 
while N ~= 1
    
    for kk = 1:numel(boi)
        obj = findobj('type', 'Figure', 'Tag', boi{kk});        
        if ~isempty(obj), break; end
    end  
        
    if ~isempty(obj)
        for kk = 1:numel(obj.Children)
            if strcmp(obj.Children(kk).Type, 'uicontrol')
                if strcmp(obj.Children(kk).Tag, 'button_CONTINUE') ...
                    || strcmp(obj.Children(kk).Tag, 'buttonA')
                    obj.Units = 'Pixels'; obj.Children(kk).Units = 'Pixels';
                    
                    fsz = obj.Position;
                    loc = obj.Children(kk).Position; % x  y w h
                    click_pos = fsz(1:2) + loc(1:2) + (loc(3:4)./2);
                    moveMouse(click_pos);                
                    
                end
            end

        end
   
    end
    N = N + 1;
end




