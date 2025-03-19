function [data, nfo] = csi_readText_jmrui(fpn)
%%%% Description:                    Read wJMRUI exported data text-file.
%%% Creator: Dr. Ir. Q. van Houtum     Version: 1.0          Date: 2020-10
%%% --------------------------------------------------------------------
%%% Read MRS data from text-files exported with JMRUI. Returns both the
%%% time and frequency domain data as complex values.
%%%
%%% Returned NFO-fields: N (#samples), nfids (#fids), dwelltime, 
%%%                      BW (bandwidth), filename, 
%%%                      PhasingZero (zero order phasing),
%%%                      PhasingFirst (first order phasing).
%%%                 
%%% Contact: quincyvanhoutum@gmail.com    



%% File processing

if nargin == 0
    % UI for file selection
    [fn, fp, fi] = uigetfile({'*.txt', 'JMRUI text file'},...
        'Select JMRUI exported text file...'); if fi == 0, return; end
    % Create one variable for full path.
    fpn = [fp fn];     
end

% Read text-file ID
fid = fopen(fpn); if fid==-1, data = NaN; return; end
if ismac || isunix % No timepenalty compared to CMD method
    % Nlines in files (loop through file ID lines and rewind file ID)
    nlines=0; while ~feof(fid), fgetl(fid); nlines=nlines+1; end
    frewind(fid); % Beginning of file
else
    % Nlines in file (uses Windows - CMD)
    [~, cmdout] = system(['find /c /v "" ', fpn]); % #lines in txtfile
    nlines = strsplit(cmdout,' '); nlines = str2double(nlines{end});
end
% Read all lines in text file.
kk = 1; flines = cell(1,nlines);
while ~feof(fid), tmp = fgetl(fid); flines{kk} = tmp; kk=kk+1; end
fclose(fid);


%% Parse Header % --- %
% Process header information of txt file with prior knowledge of JMRUI
% text-file header. 

% Fields of interest in the JMRUI header are renamed. Renaming depends on
% strings foi (jmrui-labels) and foi_new (CSIgui-labels).
foi = {'PointsInDataset', 'DatasetsInFile', 'SamplingInterval',...
       'MagneticField','SignalNames','ZeroOrderPhase','BeginTime'};
foi_new = {'N', 'nfids', 'dwelltime', 'tesla', 'filename', ...
    'PhasingZero','PhasingFirst'};
nfo = struct;
for kk = 1:size(foi,2)
    ind = find(contains(flines,foi{kk}));
    if ~isempty(ind)
        tmp = strsplit(flines{ind},':');
        if size(tmp,2) > 1
            if ~isnan(str2double(tmp{end}))
                nfo.(foi_new{kk}) = str2double(tmp{end});
            else
                nfo.(foi_new{kk}) = strtrim(tmp{end});
            end 
        end
    end
end

% Post parsing
if isfield(nfo,'dwelltime')
    nfo.dwelltime = nfo.dwelltime./1000; % In seconds
    nfo.BW = 1./nfo.dwelltime;           % Bandwidth in Hz
end

% ------ Output Parse Header
% NFO-structure created containing fields:
% N (samples), nfids (#fids), dwelltime, BW (bandwidth), filename, zero
% order phasing (PhasingZero) and PhasingFirst.


%% Parse Data % --- %
% Load txt-data-lines indexed into matrix

% Find data cell index in flines
ind_data = find(contains(flines,sprintf('out of %i in file',nfo.nfids)));
ind_data = ind_data + 1;

% Container
data = NaN(nfo.N, 2 ,nfo.nfids);
% Loop each FID (ease of programming)
for kk = 1:nfo.nfids
    tmp = flines(1,ind_data(kk):ind_data(kk)+nfo.N-1); % Lines/Fid
    tmp = cellfun(@strsplit, tmp,'uniform',0);         % Split string
    tmp = cellfun(@str2double,tmp,'uniform',0);        % Convert to dbl
    tmp = cat(1,tmp{:});                               % Concatenate
    data(:,:,kk) = cat(2,complex(tmp(:,1),tmp(:,2)),...% Store FID
                         complex(tmp(:,3),tmp(:,4)));  % Store Spectrum
end


% ----- Output Parse Data
% Data-variable created with all FIDS/Spectra loaded as complex data into
% memory. Both FID and Spectra data stored in the JMRUI text file is
% present and indexed to dimension 2;



