function [data, varargout] = csi_readSDAT(fpn, dim)
%%%% Description:                                Read SDAT data.
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%% csi_readSDAT( filename, FID or SPECTRA, data-type)
%%% Input:
%%%         fpn     :       Path and filename to sdat-file
%%%         dim     :       Dimension indexing, if known. (Speeds up).
%%%                         Requires SPAR file of not given!
%%%
%%% Output:
%%%         data    :       SDAT data
%%%         vararg  :       SPAR info header as struct.
%%%
%%% Contact: qhoutum2@umcutrecht.nl
%%% Requires the VAX-read and write script-set for Matlab. 
    

% Process input file name
[fp,fn] = fileparts(fpn); fpn = [fp '\' fn]; % Excl. ext!

% Handle input
if nargin == 1
      % READ SPAR AS STRUCT.
      spar_info = csi_readSPAR([fpn '.SPAR']); 
      if ~isstruct(spar_info), data = NaN; varargout{1} = NaN; return; end
      
      % Find dimension paragraphs
      dimi = 1; go = 1;
      while go ~= 0
          if isfield(spar_info, sprintf('dim%i_pnts',dimi))
              % Get size of dimension: dimN_pnts.
              dim(dimi) = spar_info.(sprintf('dim%i_pnts',dimi));
              % Get label of dimension: dimN_direction.
              dim_label{dimi} = regexprep(...
                    spar_info.(sprintf('dim%i_direction',dimi)),'\W','');
              % Increase dim-iter
              dimi = dimi+1;
          else, go = 0; % quite while-loop.
          end
      end
      dim = num2cell(dim);
      
      % Set output labels in spar-info.
      forder = fieldnames(spar_info); spar_info.dim_labels = dim_label;
      spar_info = orderfields(spar_info, {'dim_labels', forder{:}});
      if nargout == 2, varargout{1} = spar_info; end
else
      dim = {dim}; if nargout == 2, varargout{1} = dim; end
end

% Get data length
fid = fopen([fpn '.SDAT'],'r','ieee-le');
% Abort if file not opened.
if fid == -1
    data = NaN; warning('SDAT-file does not exist.'); 
    return; 
end 
data_size = length(fread(fid));	fclose(fid);

% Read data into memory
fid  = fopen([fpn '.SDAT'],'r','ieee-le');
data = freadVAXG(fid , data_size, 'float32'); fclose(fid);

% Reshape and convert to complex data.
data = reshape(data,2,[]); data = complex(data(1,:),data(2,:));

% THIS IS BAD QUICK FIX CODE BAD BAD BAD
if isfield(spar_info,'spec_row_upper_val')
    disp badbadgoingbadcheckCSIreadsdat!
    dim{2} = spar_info.spec_row_upper_val;
end

% Reshape according array-dimensions from either SPAR or second input.
try data = squeeze(reshape(data,dim{:}));
catch err
    err
warning('Reshaping of SDAT-data unsuccesfull! Raw-indexing returned.');
end



