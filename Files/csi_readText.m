function data = csi_readText(fp, data_size)
%%%% Description:                    Load MRSI data text files.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2018-04
%%% --------------------------------------------------------------------
%%%
%%% Input:
%%%         fp:         Full filepath to MRSI data text file.
%%%         data_size:  1. Array    - with data size, number of elements 
%%%                                   after reshape loaded must not change.
%%%                     2. Filepath - to data_size text file.
%%%                     3. None     - mrsi_loadText will check for a text
%%%                                   file with path and name fp appended
%%%                                   with _ArraySize.txt.
%%%
%%% If no input is given, a UI for file selection is opened. If wrong data
%%% size or no data size file is present, the output is not reshaped.
%%% 
%%%
%%% Contact: qhoutum2@umcutrecht.nl

%% Process input

% If no data size is given.
if nargin == 1
    
    % Data size file:
    if exist([fp(1:end-4) '_ArraySize.txt'],'file')
        % Check if ArraySize file is present
        data_size = [fp(1:end-4) '_ArraySize.txt'];
    else
        % Else dont reshape the data
        data_size = 0;
    end
    
elseif nargin == 0
    
    % UI for file selection
    [fn, fp, fi] = uigetfile({'*.txt', 'MRSI data text file'},...
        'Select MRSI data text file...'); if fi == 0, return; end
    % Create one variable for full path.
    fp = [fp fn]; 
    
    % Data size file:
    if exist([fp(1:end-4) '_ArraySize.txt'],'file')
        % Check if ArraySize file is present
        data_size = [fp(1:end-4) '_ArraySize.txt'];
    else
        % Else dont reshape the data
        data_size = 0;
    end
end

%% Read data file
fid = fopen(fp);                            % Open file ID 
if fid == -1, data = NaN; warning('Text-file does not exist.'); return; end
data_raw = fscanf(fid,'%f %f \n',[2 inf])'; % Read file
if isempty(data_raw), data = NaN; return; end
fclose(fid);                                % Close file ID

%% Read data size file
if ischar(data_size)
    fid = fopen(data_size);                 % Get file ID
    data_size = fscanf(fid,'%f')';          % Read file
    fclose(fid);                            % Close file
end

%% Create output and reshape data
data = complex(data_raw(:,1), data_raw(:,2));
try % Using try, output is guaranteed if wrong reshape size is given.
    if data_size ~= 0, data = reshape(data, data_size); end
catch err
    fprintf('%s\n', err.message); % Display error message of the catch.
end

% If all imaginary parts are zero, remove imaginary part.
if sum(imag(data)) == 0, data = real(data); end




