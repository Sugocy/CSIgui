function Statistics_Viewer(stats)
% Plot the fields and values of struct "stats" in a listbox gui, with
% option to copy using ctrl + c.
%
% Dr. Q. van Houtum, quincyvanhoutum@gmail.com
% Version 2.0 - 06/2024

% Matlab year/version
ver = version('-release'); ver = str2double(ver(1:end-1));


% Figure
if ver > 2018
    fh = uifigure('MenuBar','none', 'ToolBar','none',...
                  'Name', 'Statistics Viewer',...
                  'Color','k', 'NumberTitle','off',...
                  'WindowKeyPressFcn', @keyPressFcn);
else
    fh = uifigure('MenuBar','none', 'ToolBar','none',...
                  'Name', 'Statistics Viewer',...
                  'Color','k', 'NumberTitle','off');
end
fh.Position(3:4) = round(fh.Position(3:4) ./ [1.9 1.25]);
scrsz = get(0,'screensize');
if scrsz(3) < 1280, scrsz = [0 0 1280 720]; end
fh.Position(1:2) = (scrsz(3:4) .* [0.15 0.85]) - [0 fh.Position(4)];



dxy = 5;
lb_sz = fh.Position(3:4)-(dxy*2); lb_sz(1) = lb_sz(1)./2.5;
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
    
    if ischar(val)
        vals{kk} = val;
    else

        if numel(val) > 1
            nvals = numel(val);
            str = repmat('[%g] ', 1, nvals);
            vals{kk} = sprintf( str, val(:));
        else
            
            [nZeros, acc] = numzeros(val);
            if nZeros <= 0 || nZeros == acc
                vals{kk} = strip(sprintf('%6.2f', val));    
            else
                prefix = '0'; if nZeros >= 9, prefix = ''; end
                vals{kk} = sprintf('%3.3fe-%s%i', ...
                                   val*10.^(nZeros+1), prefix, nZeros+1);
            end

        end
        
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