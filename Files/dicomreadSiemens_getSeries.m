function files = dicomreadSiemens_getSeries(varargin)
% Given a single ima file. It will look into its directory to find
% ima-files that are part of the scanning protocol. This script looks at
% the filenames only and does not load any dicom-related information.
%
% input: full filepath to file
%               Returns all files that match to the input file.
%        filepath and filename seperatly
%               Returns all files taht match to the input file.
%        filepath only ending with backwards dash '\'
%               Returns a cell array with groups of matching ima-files.
%
% Returns NaN if no group is found.

if nargin == 1
    % Either a directory or full file path
    if isfolder(varargin{1})
        fp = varargin{1}; fn = [];
    elseif isfile(varargin{1})
        [fp, fn, ext] = fileparts(varargin{1});
        fn = [fn ext];
    end
elseif nargin == 2
    fp = varargin{1}; fn = varargin{2};
end

% Correct filepath
if ~strcmp(fp(end),'\'), fp = [fp '\']; end

if isempty(fn)
    files = getProtocols_directory(fp);    
else
    % get all groups
    files = getProtocols_directory(fp);
    
    % Find filename in groups
    inCell = ...
    cellfun(@contains, files, repmat({fn},size(files)), 'uniform', false);
    % Select this group for output.
    bool = cellfun(@sum, inCell);
    
    if sum(bool) ~= 0    
        % Output
        files = files{logical(bool)};       
    else
        files = NaN;
    end
end


function files = getProtocols_directory(fp)
% This reads the full given directory and groups ima-file according to
% their file-names.

% Directory content
content = dir(fp);

% Create dicom-image groups for all files in content.
imgcont = {};
protcont = struct;
for kk = 1:size(content,1)
    if ~content(kk).isdir % Skip directories
        [~,~,ext] = fileparts(content(kk).name); % Get extension
        if strcmpi(ext,'.ima')
            snip = strsplit(content(kk).name,'.'); % Find prot/img-nr
            nr_iter = 0;
            for sni = 1:size(snip,2)
                if ~isnan(str2double(snip{sni}))
                    if nr_iter == 0
                        protnr = snip{sni}; % Protocol-nr
                    elseif nr_iter == 1
                        imgnr = snip{sni}; % Image-nr in protocol
                        break; 
                    end
                    nr_iter = nr_iter + 1;
                end
            end
            
            % Add to group
            protnrf = ['p' protnr];
            if isfield(protcont,protnrf) % Existing group
                protcont.(protnrf).nfiles = protcont.(protnrf).nfiles + 1;               
                imgcont{protcont.(protnrf).ind} = ...
                    [imgcont{protcont.(protnrf).ind} kk];             
            else % New group
                protcont.(protnrf).nfiles = 1;
                protcont.(protnrf).ind = size(imgcont,2)+1;
                imgcont{size(imgcont,2)+1} = kk;             
            end            
        end
    end    
end

% convert groups to output
files = cell(1,size(imgcont,2));

for grpi = 1:size(imgcont,2)   
    for pri = 1:size(imgcont{grpi},2)
        if pri == 1
            files{grpi} = {content(imgcont{grpi}(pri)).name};
        else
            files{grpi} = [files{grpi}{:} {content(imgcont{grpi}(pri)).name}];
        end
    end % End of group-dcm for-loop.    
end


