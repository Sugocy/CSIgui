function csi_editSPAR(fpn, varargin)
%%%% Description:                                             Update SPAR.
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2017-07
%%% --------------------------------------------------------------------
%%% csi_editSPAR([filepath '\' filename])
%%% Input:
%%%         fpn     :       file name and path of SPAR-file, extension
%%%                         optional. If only filepath given, userinput by
%%%                         GUI will be enabled.
%%%
%%%         varargin:       Data which has to be written to parameters.
%%%                         Find parameters names by loading the spar-file
%%%                         using csi_readSPAR. Each field represents
%%%                         parameter which can be updated in SPAR.
%%%                         {'name' - 'value'}; If varargin{1} is an empty
%%%                         field, e.g. an odd sized varargin,
%%%                         userinput by GUI will be enabled.
%%%
%%% Contact: qhoutum2@umcutrecht.nl
%%%

user_mode = 1;

%% Handle input 1
[fp, fn, ~] = fileparts(fpn); 
if strcmp(fp(end) ,'\')
    fpn = [fp fn '.SPAR'];
else
    fpn = [fp '\' fn '.SPAR'];
end

% 1. Read SPAR file.
spar = csi_readSPAR(fpn); 
if ~isstruct(spar) 
    warning('SPAR-file does not exitst. Abort edit.'); return; 
end 


if nargin > 1
%% PROCESS INPUT PARAMETERS
   % Get varargin parameters
   userinp = varargin{1};
   
   % Set extra GUI base userinput on or off.
   if mod(size(userinp,2),2) == 1 % on
      user_mode = 1; userinp = userinp(2:end);    
   else                           % off
      user_mode = 0;
   end
   
   % Replace value with user input
   for qi = 1:2:size(userinp,2)
       attr = userinp{qi}; val = userinp{qi+1};
       % Set new value
       spar.(attr) = val;
   end

end

if nargin == 1 || user_mode == 1
%% GET USER INPUT

    % Ask qres inputs/time!
    fnames = fieldnames(spar); qres = 12;
    for qi = 1:ceil((size(fnames,1)/qres))

        if qi == ceil((size(fnames,1)/qres))    
            attr =fnames( ((qi-1)*qres)+1:end);
            val = getFields(spar, attr);
            val = cellfun(@num2str, val,'UniformOutput', 0);
            % DLG for user input
            ans_user{qi} = inputdlg(attr,'Input SPAR-parameters', 1,val);    
        else
            attr = fnames( ((qi-1)*qres)+1:(qi*qres));
            val = getFields(spar,attr);
            val = cellfun(@num2str, val,'UniformOutput', 0);
            % DLG for user input
            ans_user{qi} = inputdlg(attr,'Input SPAR-parameters', 1, val);
        end
        
        for ai = 1:size(attr,1)
            if isnan(str2double(ans_user{qi}{ai}))
                spar.(attr{ai}) = ans_user{qi}{ai};
            else
                spar.(attr{ai}) = str2double(ans_user{qi}{ai});
            end
        end
    end


end



%% UPDATE SPAR-INFO

% Read SPAR as lines
fid = fopen(fpn,'r'); 
if fid == -1, warning('SPAR-file not found!'); return; end
% Loop each line in SPAR-file
n = 1;  while ~feof(fid), spar_lines{n,1} = fgetl(fid); n=n+1; end
% Close file ID
fclose(fid);


fnames = fieldnames(spar); ind = NaN(1,size(fnames,1));
for fi = 1:size(fnames,1)
    attr = fnames{fi}; val  = spar.(attr);
    
    % Take care of one kown exception.
    if strcmp(attr, 'Spec_imageInPlaneTransf')
        attr = 'Spec.image in plane transf';
    end
    
    % Find attribute in SPAR
    ind_tmp = find(cellfun(@isempty, strfind(spar_lines, attr) )==0);
    if isempty(ind_tmp)
    % Name of parameter in SPAR different from attribute label
        fprintf('Parameter %s not found in SPAR, nor updated.\n', attr);
        ind_tmp = NaN;
    end
    
    % Found more lines where attr is found
    if size(ind_tmp,1) > 1
        for kk = 1:size(ind_tmp,1)
            if sum(ind_tmp(kk) == ind) == 0
                ind_tmp = ind_tmp(kk); break;
            end
        end
    end
    % Store index.
    ind(fi) = ind_tmp;
    
    % Update lines!
    if ~isnan(ind(fi))
        % Get original line
        tmp_line = spar_lines{ind(fi)};
        % Cut original line up to colon+1 and add new value
        
        cut_1 = tmp_line(1:strfind(tmp_line,':'));
        cut_2 = tmp_line(strfind(tmp_line,':')+1:end);
        if strcmp(cut_2(1),'[') && ~strcmp(val(1),'[') 
            tmp_line = [ cut_1 '[' val ']'];
        elseif ~isempty(val) && strcmp(val(1),'[') 
            tmp_line = [ cut_1     num2str(val) ];
        else
            tmp_line = [ cut_1 ' ' num2str(val) ];
        end
        

        % Write in spar_lines.
        spar_lines{ind(fi)} = tmp_line;
    end
end

%% REWRITE SPAR-FILE

% Read SPAR as lines
fid = fopen(fpn,'w'); 
if fid == -1, warning('SPAR-file not found!\n'); return; end
for kk = 1:size(spar_lines,1), fprintf(fid, '%s\r\n', spar_lines{kk}); end
fclose(fid); fprintf('SPAR-file rewritten.\n');













