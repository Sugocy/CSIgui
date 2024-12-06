function spar = csi_readSPAR(fpn)
%%%% Description:                                       Read SPAR-file
%%% Creator: Ir. Q. van Houtum       Version: 1.0           Date: 2017-07
%%% --------------------------------------------------------------------
%%% Read SPAR-file into structure. Will load all parameters found.
%%% Input: full filepath.
%%%
%%% Contact: qhoutum2@umcutrecht.nl


% Handle input
[fp, fn, ~] = fileparts(fpn); fpn = [fp '\' fn '.SPAR'];


% Check if file exists
if ~exist(fpn,'file'), warning('SPAR-file does not exist.'); 
    spar = NaN; return; 
end

% Open file
fid = fopen(fpn,'r'); 
if fid == -1, warning('SPAR-file not found missing!'); return; end
% Loop each line in SPAR-file
n = 1;  while ~feof(fid), spar_lines{n,1} = fgetl(fid); n=n+1; end
% Close file ID
fclose(fid);

% Create output structure
spar = struct;

% Loop each line and save only parameter lines.
for li = 1:size(spar_lines,1)
    tmp = spar_lines{li};
    if ~isempty(tmp) && sum(strcmp(tmp(1),{'!','',' ',})) == 0
        attr = strrep(genvarname(tmp(1:strfind(tmp, ':')-2)),'0x2E','_');
        val  = tmp(strfind(tmp, ':')+1:end);
        
        % If conversion to double results in NaN e.g. its a char.
        if isnan(str2double(val)), ...
                spar.(attr) = strtrim(val);                                        % Save as string.
        else                      ...
                spar.(attr) = str2double(val);                              % Save as double.
        end
        
    end
end


