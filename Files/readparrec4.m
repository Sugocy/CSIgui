function [imgmagn, imgphase, imgother, sliceinfo, examinfo] = readparrec4(filename, outputmode)
%%%% Description:                     Read Philips PAR-and REC-files v4.2 
%%% Creator: Ir. Q. van Houtum       Version: 4.0          Date: 2015-09
%%% --------------------------------------------------------------------
%%% Output variables:
%%%
%%% [magnitude images, phase images, other images, slice info, exam info] 
%%%         = readparrec(.*);
%%%
%%% Image matrices index =
%%%         width x height x slices x echoes x dynamics x cardiac phase 
%%%         and optional: scale(1:3).
%%%
%%% Input arguments:
%%%
%%% .* filename (String)         
%%%             / Including path, excluding .PAR/.REC, output set to float.
%%% .* outputmode (Integer)
%%%             / Scaled to (1)float, (2)display (3)REC-values and (4)all.
%%% .* filename + outputmode           
%%%             / See both input arguments
%%% .* no input                 
%%%             / Uidialog for file location and output set as float
%%%                
%%% From par-header (v4.2): 
%%%             Display scaling     -      Floating point scaling
%%%             DV = PV * RS + RI   -      FP = DV / (RS * SS) 
%%% ...with PV = Pixel Value from REC-file. Scale factors saved per slice.
%%%
%%%
%%% Sliceinfo is a struct with fields M, P and O representing slice
%%% information for each image-array. The cell-array, containing a struct
%%% with all the information, is indexed according to: 
%%%                 sliceinfo.M{ slice#, echo#, dyn#, card#);
%%% making it easy to find the data for each slice.

%% Initiation - Input arguments and file handling 
switch nargin 
    case 0
        [filename, filepath, findex] = uigetfile('*.PAR', 'Open PAR-file');    
        if findex == 0, return; end
        file1 = [filepath '\' filename(1:end-3) 'REC']; 
        file2 = [filepath '\' filename];
        outputmode = 1;
    case 1
        if  ischar(filename)
            outputmode = 1; file1 = [filename '.REC']; 
            file2 = [filename '.PAR'];
        elseif isfloat(filename)
            outputmode = filename;
            [filename, filepath, findex] = ...
                uigetfile('*.PAR', 'Open PAR-file'); 
            if findex == 0, return; end 
            file1 = [filepath '\' filename(1:end-3) 'REC'];
            file2 = [filepath '\' filename]; 
        end
        
    case 2
        file1 = [filename '.REC']; file2 = [filename '.PAR'];
end

%% Open PAR files
if exist('file1', 'var') && exist('file2', 'var')
    % --- Open file-ID handle --- %
    fidR = fopen(file1, 'r'); fidP = fopen(file2, 'r');
    
    % If fopen succeeded, continue:
    if fidR <0, warning('REC-file does not exist.'); 
        imgmagn = 0; imgphase = 0; imgother = 0; 
        sliceinfo = 0; examinfo = 0; return;
    else
        % Extract info from files using dir-api: only if file exists
        infoR = dir(file1); sizeR = infoR.bytes/1048576; % Megabytes (MB)
        infoP = dir(file2); sizeP = infoP.bytes/1048576; % Megabytes (MB)
    
        % Scan text from PAR-file and read data from REC file.
        if sizeP > 0,dataP = textscan(fidP,'%s','Delimiter','\n'); 
            dataP = char(dataP{1});
        end
        if sizeR > 0,dataR = fread(fidR, 'uint16'); end
        % Close file-handle ID
        fclose(fidR);   fclose(fidP);
    end
else, warning('Wrong input file-type, requires *par file'); return; 
end


%% Find starting index of general and image data in PAR
%if exist('dataP', 'var') == 1
for i = 1:size(dataP,1)
    templine = dataP(i,:);
    % if only spaces - eg. empty == start of data
    if sum(isspace(templine)) == size(dataP,2)
        if    ~exist('dataStartLine', 'var'), dataStartLine = i+1; 
        else, dataStopLine = i-1; 
        end
    % If dot as first char in line == general information
    elseif strcmp(templine(1),'.') == 1
        if    ~exist('infoStartLine', 'var'), infoStartLine = i;
        else, infoStopLine = i; 
        end
    end
end
% else % No data P created - abort
%     warning('No dataP variable found.'); 
%     imgmagn = 0; imgphase = 0; imgother = 0; sliceinfo = 0; examinfo = 0;
%     return;
% end


%% Get imagetype image dimensions from PAR
% Number of images and image types
Mcnt = 0; Pcnt = 0; Ocnt = 0;
nrslicem = 0;  nrdynm = 0; nrcardm = 0; nrechom = 0; 
nrslicep = 0;  nrdynp = 0; nrcardp = 0; nrechop = 0; 
nrsliceo = 0;  nrdyno = 0; nrcardo = 0; nrechoo = 0; 

% Loop each dataline which represents one image.
for i = dataStartLine:dataStopLine
    % Data line i:
    dataline = sscanf(dataP(i,:), '%f');
    
    % Image type: Magnitude, phase or other
    imgtype = dataline(5);                                   
    switch imgtype
        case 0, Mcnt = Mcnt + 1;                        %% Magnitude image 
            if (dataline(1) > nrslicem)                                     % Check slice index of M-image
                nrslicem = nrslicem+1; 
                matrszM  = [dataline(10) dataline(11)]; 
            end 
            if (dataline(2) > nrechom),nrechom = nrechom+1; end             % Check echoes index
            if (dataline(3) > nrdynm), nrdynm  = nrdynm +1; end             % Check dynamic index
            if (dataline(4) > nrcardm),nrcardm = nrcardm+1; end             % Check cardiac phase index
            
        case 3, Pcnt = Pcnt + 1;                       %% Phase images
            if (dataline(1) > nrslicep)                                     % Check slice index of P-image
                nrslicep =nrslicep+1; 
                matrszP = [dataline(10) dataline(11)]; 
            end
            if (dataline(2) > nrechop),nrechop = nrechop+1; end             % Check echoes index
            if (dataline(3) > nrdynp), nrdynp  = nrdynp +1; end             % Check dynamic index
            if (dataline(4) > nrcardp),nrcardp = nrcardp+1; end             % Check cardiac phase index 
            
        otherwise, Ocnt = Ocnt + 1;                    %% Other images
            if (dataline(1) > nrsliceo)                                     % Check slice index of P-image
                nrsliceo =nrsliceo+1; 
                matrszO = [dataline(10) dataline(11)]; 
            end
            if (dataline(2) > nrechoo),nrechoo = nrechoo+1; end             % Check echoes index
            if (dataline(3) > nrdyno), nrdyno  = nrdyno +1; end             % Check dynamic index
            if (dataline(4) > nrcardo),nrcardo = nrcardo+1; end             % Check cardiac phase index 
    end
end


%% Compare counted and stored data
% Each line will equal one image.
nrofimg = size(dataStartLine:dataStopLine,2); nrofcnt = (Ocnt+Pcnt+Mcnt);
if nrofimg ~= nrofcnt
    fprintf('Data-image counting error in PAR-file!\n'); 
    fprintf('%i %i %i\n',Ocnt, Pcnt, Mcnt);
    fprintf('%i', nrofimg);  return; 
end


%% Create storage variables for REC
if     outputmode ~= 4                                                     % Only 1 output type requested
    if Mcnt > 0, imgmagn  = NaN([matrszM nrslicem nrechom nrdynm nrcardm]); 
    else, imgmagn  = []; end   
    if Pcnt > 0, imgphase = NaN([matrszP nrslicep nrechop nrdynp nrcardp]); 
    else, imgphase = []; end
    if Ocnt > 0, imgother = NaN([matrszO nrsliceo nrechoo nrdyno nrcardo]); 
    else, imgother = []; end
elseif outputmode == 4                                                     % All 3 output types requested   
    if Mcnt > 0, imgmagn  = NaN([matrszM nrslicem nrechom nrdynm nrcardm 3]); 
    else, imgmagn   = []; end 
    if Pcnt > 0, imgphase = NaN([matrszP nrslicep nrechop nrdynp nrcardp 3]); 
    else, imgphase  = []; end
    if Ocnt > 0, imgother = NaN([matrszO nrsliceo nrechoo nrdyno nrcardo 3]); 
    else, imgother  = []; end
end

% To save slice specific information/properties structure
if Mcnt > 0, minfo = cell(nrslicem, nrechom, nrdynm, nrcardm); 
else, minfo = [];
end 
if Pcnt > 0, pinfo = cell(nrslicep, nrechop, nrdynp, nrcardp);
else, pinfo = [];
end 
if Ocnt > 0, oinfo = cell(nrsliceo, nrechoo, nrdyno, nrcardo);
else, oinfo = [];
end

%% Read image data from PAR and REC 
% i2     - Counts inside the REC file - HIGHLY IMPORTANT FOR THE ALGORITHM!
i2 = 0; 

% Start looping each line in PAR file and save the data from REC file
% accordingly: i1 = par-file line, i2 = rec-data count variable
for i = dataStartLine:dataStopLine
    % Par-info for image
    dataline = sscanf(dataP(i,:), '%f');
    % Slice info variable
    slinfo = struct;

                % ========================================== %
                % ----           Prepare REC data         ---%
                % ========================================== %
    
    % Scale factors for REC-data from PAR file.
    slinfo.RI = dataline(12);                                               % Rescale inter - px values on console
    slinfo.RS = dataline(13);                                               % Rescale slope - rescaling px values
    slinfo.SS = dataline(14);                                               % Scale slope   - px values as floating-pval
    
    % Image contrast window for display
    contrCenter = dataline(15); contrWidth  = dataline(16);
    slinfo.consolecontr   = ...
        [contrCenter-(contrWidth/2) contrCenter+(contrWidth/2)];
    
    % Load values from REC-file
    % Number of pixels/image extracted from a list
    rec_matrixsz = [dataline(10) dataline(11)]; 
    pxvalRec  = reshape(dataR( ((rec_matrixsz(1)*rec_matrixsz(2)*i2)+1) :...
        ((rec_matrixsz(1)*rec_matrixsz(2))*(i2+1)) ), rec_matrixsz(1:2));
    
    % Scaling of the REC-values
    % See formulas in PAR-header or readparrec help.
    if outputmode ~= 3                                                     % If REC-values are not requested output, scale to display values!
        if isfield(slinfo,{'RS', 'RI', 'SS'}) == 1                          % If scaling factors exist
            pxvalDisp  = (pxvalRec.*slinfo.RS) + slinfo.RI;                 % Calculate DV
            if (outputmode == 1) ||  (outputmode == 4)                      % If floating point output is requested, scale to floating values!
                pxvalFloat = pxvalDisp./...
                    (slinfo.RS*slinfo.SS);                                  % Calculate FP
            end
        end
    end
                % ========================================== %
                % --- Prepare image specific properties  --- %
                % ========================================== %
    
    % Save miscellaneous slice information 
    if     dataline(26) == 1, or = 'tra';
    elseif dataline(26) == 2, or = 'sag';
    elseif dataline(26) == 3, or = 'cor';
    else,                     or = 'Unknown';    
    end
    slinfo.acqori_tsc= or; % TRA/SAG/COR
    slinfo.angulation= [dataline(17) dataline(18) dataline(19)];            % ap,fh,rl
    slinfo.offcenter = [dataline(20) dataline(21) dataline(22)];            % ap,fh,rl   
    slinfo.resolution= [dataline(29) dataline(30) dataline(23)];            % pixel x,y,z
    slinfo.slthick   = dataline(23);
    slinfo.echonr    = dataline(2);  
    slinfo.scanperc  = dataline(9);
    slinfo.gapsize   = dataline(24); 
    slinfo.TE        = dataline(31); 
    slinfo.flipangle = dataline(36);
    slinfo.NSA       = dataline(35);
    slinfo.dynnr     = dataline(3);
    slinfo.cardnr    = dataline(4);
    slinfo.timeofscan= dataline(32);
    slinfo.timetrig  = dataline(33);

    
    % Ordering of structure fields
    slinfo  = orderfields(slinfo);
    
                % ========================================== %
                % --- Save all data to correct img array --- %
                % ========================================== %
    
    % --- Magnitude Images ----- %                                          % Magnitude images
    if  dataline(5) == 0,                                                   slinfo.image = 'M';                         

        % Store loaded magnitude REC-values:
        switch outputmode % X x Y x Slice x Echo x Dynamic x Card
            case 1, imgmagn(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalFloat;
            case 2, imgmagn(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalDisp;
            case 3, imgmagn(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalRec;
            case 4, imgmagn(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4),:) = cat(5, pxvalFloat, pxvalDisp, pxvalRec);
        end
        
        % Save sliceinfo to correct position in information cell-array
        minfo{dataline(1),dataline(2),dataline(3), dataline(4)} = slinfo;
        
        
    % --- Phase images -------- %                                           % Phase images
    elseif dataline(5) == 3,                                                slinfo.image = 'P';

        % Store loaded PHASE REC values:
        switch outputmode % X x Y x Slice x Echo x Dynamic x Card
            case 1, imgphase(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalFloat;
            case 2, imgphase(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalDisp;
            case 3, imgphase(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalRec;
            case 4, imgphase(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4),:) = cat(5, pxvalFloat, pxvalDisp, pxvalRec);
        end

        % Save sliceinfo to correct position in information cell-array
        pinfo{dataline(1),dataline(2),dataline(3), dataline(4)} = slinfo;
        
    % --- Other image ------- %                                             % Other images 
    elseif (dataline(5) ~= 3) && (dataline(5) ~= 0),                        slinfo.image = 'O';
        
        % Store loaded Other REC values:
        switch outputmode % X x Y x Slice x Echo x Dynamic x Card
            case 1, imgother(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalFloat;
            case 2, imgother(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalDisp;
            case 3, imgother(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4)) = pxvalRec;
            case 4, imgother(:,:,dataline(1),dataline(2),dataline(3),...
                    dataline(4),:) = cat(5, pxvalFloat, pxvalDisp, pxvalRec);
        end
        
        % Save sliceinfo to correct position in information cell-array
        oinfo{dataline(1),dataline(2),dataline(3), dataline(4)} = slinfo;
        
    end 
    
    % Increase image selection loop count, see pxvalREC of Rec-file.
    i2=i2+1;
end

%% Save all image specific data to output argument sliceinfo;
sliceinfo   = struct;
sliceinfo.M = minfo; sliceinfo.P = pinfo; sliceinfo.O = oinfo;



%% General information
% Store general exam information into structure examinfo
if nargout > 3, examinfo = struct;
for q = infoStartLine:infoStopLine
    % Get line and find start of data index after semicolon
    templine  = dataP(q,:); ind = strfind(templine,':')+4; 
    % Search for data/information fields
    % - Strings
    if strfind(templine, ' Patient name '),         tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.patientname = tmp; end
    if strfind(templine, ' Examination name '),     tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.examname = tmp; end
    if strfind(templine, ' Protocol name '),        tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.protname = tmp; end
    if strfind(templine, ' Examination date/time '),tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.examdatetime = tmp; end
    if strfind(templine, ' Patient position '),     tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.patientpos = tmp; end
    if strfind(templine, ' Preparation direction '),tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.prepdir = tmp; end
    if strfind(templine, ' Technique '),            tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.scandur = tmp; end
    if strfind(templine, ' Scan mode '),            tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.scanmode = tmp; end
    % - Integer and floating point + counted values.
    if strfind(templine, ' Acquisition nr '),       tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.acqnr = str2double(tmp); end
    if strfind(templine, ' of cardiac phases '),    tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.maxnrcard  = str2double(tmp);  examinfo.nrcardcnt = (nrcardm); end
    if strfind(templine, ' of slices/locations '),  tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.maxnrslice = str2double(tmp);  examinfo.nrslicecnt = (nrslicem); end
    if strfind(templine, ' of dynamics '),          tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.maxnrdyn   = str2double(tmp);  examinfo.nrdyncnt = (nrdynm); end
    if strfind(templine, ' of echoes '),            tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.maxnrechoes = str2double(tmp); end
    if strfind(templine, ' of mixes '),             tmp = templine(ind:end); tmp = strrep(tmp, '  ', ''); examinfo.maxnrmix = str2double(tmp); end
    if strfind(templine, ' Water Fat shift '),      tmp = templine(ind:end); tmp = sscanf(tmp, '%f'); examinfo.wfshift = tmp; end
    if strfind(templine, ' Scan resolution '),      tmp = templine(ind:end); tmp = sscanf(tmp, '%i'); examinfo.scanres = [tmp(1) tmp(2)]; end
    if strfind(templine, ' Repetition time [ms] '), tmp = templine(ind:end); tmp = sscanf(tmp, '%f'); examinfo.TR  = tmp; end
    if strfind(templine, ' FOV (ap,fh,rl) '),       tmp = templine(ind:end); tmp = sscanf(tmp, '%f'); examinfo.FOV = tmp; end
end
    examinfo = orderfields(examinfo);
end

% No output arguments: combine magnitude, phase and other image arrays
if nargout == 0
    imgmagn = {imgmagn; imgphase ; imgother};
    fprintf('No output arguments specified, all image types are combined.');
end


end
