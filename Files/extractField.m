function fieldvals = extractField(istruct, ifield)
%%%% Description:                     Extract multiple struct fields 
%%% Creator: Ir. Q. van Houtum       Version: 2.1          Date: 2017-04
%%% --------------------------------------------------------------------
%%%
%%% If istruct is a cell-struct (e.g. istruct{:,:,:,:}), extractField will
%%% get all values from the ifield from each cell-index.
%%%
%%% See also: isfieldfull();
%%%
%%% Supported for up to 5 index-dimensions.
%%% Supported for graphics object in Matlab 2016a. Uses fieldnames and
%%% isfield combined to figure out existance of subfields in struct-cell.

% Get field-depth field of interest in structure!
% 1. Get all subfields in requested ifield.
dotindex = strfind(ifield, '.'); subfields = cell(1,size(dotindex,2));
for qq = 1:size(dotindex,2)+1
    if     qq == 1
        if ~isempty(dotindex)
            subfields{qq} = ifield(1:dotindex(qq)-1);
        else
            subfields{qq} = ifield;
        end
    elseif qq == size(dotindex,2)+1
        subfields{qq} = ifield(dotindex(end)+1:end);
    else
        subfields{qq} = ifield(dotindex(qq-1)+1:dotindex(qq)-1);  
    end
end

% Check if all fields exist in istruct.
tmpstruct = istruct{1,1,1,1,1};
for qq = 1:size(subfields,2)
    % Get fieldnames if possible
    fnames = fieldnames(tmpstruct); 
    fnindex = find(strcmp(fnames,subfields{qq}),1);
    if isfield(tmpstruct, subfields{qq}) || ~isempty(fnindex) 
        tmpstruct = tmpstruct.(subfields{qq});
    else warning('RealWorldSys:NoExtractField',...
        'Extracted field(s) not in structure - values not extracted'); 
        fieldvals = []; return;
    end
end

% Get size of field-values ---> use tmpstruct
valsz = size(tmpstruct);
% Storage variable
fieldvals = cell(size(istruct)); % 1 x N values - stslecdy

% Loop through expected 4 dimensions
for st = 1:size(istruct,1)
    for sl = 1:size(istruct,2)
        for ech = 1:size(istruct,3)
            for dyn = 1:size(istruct,4)
                for dum = 1:size(istruct,5) % Dummy dim.
                    tstruct = istruct{st,sl,ech,dyn,dum};
                    for qq = 1:size(subfields,2)
                        tstruct = tstruct.(subfields{qq});
                    end
                    fieldvals{st,sl,ech,dyn,dum} = tstruct; 
                end
            end
        end
    end
end


end