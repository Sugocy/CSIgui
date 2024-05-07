function Statistics_Viewer(stats)
% Plot the fields and values of struct "stats" in a listbox gui, with
% option to copy using ctrl + c.
%
% Dr. Q. van Houtum, quincyvanhoutum@gmail.com
% Version 1.0 - 04/2024

% Figure
fh = uifigure('MenuBar','none', 'ToolBar','none',...
              'Name', 'Statistics Viewer',...
              'Color','k', 'NumberTitle','off',...
              'WindowKeyPressFcn', @keyPressFcn);
fh.Position(3:4) = round(fh.Position(3:4) ./ [2 1.5]);

dxy = 5;
lb_sz = fh.Position(3:4)-(dxy*2); lb_sz(1) = lb_sz(1)./2.75;
lbv_sz = lb_sz; lbv_sz(1) = fh.Position(3) - lb_sz(1) - (dxy*2);

% Listbox fieldname
gui.lb = uilistbox(fh, 'Items', {}, 'Value', {});
gui.lb.Position = [dxy dxy lb_sz];
gui.lb.BackgroundColor = [0 0 0]; gui.lb.FontColor = [0.88 0.88 0.88];
gui.lb.Tag = 'fields';
gui.lb.ValueChangedFcn = @deselect;
gui.lb.FontWeight = 'bold';

% Listbox value
gui.lbv = uilistbox(fh, 'Items', {}, 'Value', {});
gui.lbv.Position = [lb_sz(1)+dxy dxy lbv_sz];
gui.lbv.BackgroundColor = [0 0 0]; gui.lbv.FontColor = [0.94 0.94 0.94];
gui.lbv.Tag = 'values'; 

flds = fieldnames(stats);
vals = cell(1,numel(flds)); 
for kk = 1:numel(flds)
    val = stats.(flds{kk});
    if iscell(val), val = cell2mat(val); end
    
    if numel(val) > 1
        nvals = numel(val);
        str = repmat('%g | ', 1, nvals);
        vals{kk} = sprintf( str, val(:));
    else
        vals{kk} = sprintf( '%f', val);    
    end
    flds{kk} = sprintf('%s:', flds{kk});
end

gui.lb.Items = flds; gui.lbv.Items = vals;
guidata(fh, gui);

end

function keyPressFcn(hobj, evt)
% To enable copy-data via ctrl+c

modf = evt.Modifier{:}; key = evt.Key;
switch modf
    case 'control'
        switch key
            case 'c'
                for kk = 1:numel(hobj.Children)
                    if strcmp(hobj.Children(kk).Tag, 'values')
                        data2exp = hobj.Children.Value;
                    end
                end
                clipboard("copy", data2exp)
        end
end
end

function deselect(hobj, ~)
% Deselect any line in fields-listbox
    hobj.Value = {};
end