function status = csi_writeSDAT( fpn , data, data_type, varargin)
%%%% Description:                                Write CSI data as SPAR.
%%% Creator: Ir. Q. van Houtum       Version: 1.5          Date: 2017-08
%%% --------------------------------------------------------------------
%%% mrs_writeSDAT( filename, FID or SPECTRA, data-type)
%%% Input:
%%%         filename:       name with or without path to save SDAT file to.
%%%         data    :       FID or SPECTRA data. Default reads as FID. If 
%%%                         input are SPECTRA set type to 'spec'. 
%%%         type    :       Either 'fid' (default) or 'spec'.
%%%
%%%         varargin:       Cell array with label for the remaining
%%%                         index-dimensions excluding the time/freq-index.
%%%                         Used to store the dimensions in SPAR file
%%%                         correctly.
%%%
%%% CSI/MRS data is saved as FID in SDAT-files. SPAR file is created
%%% accordingly! Data size is set to parameters dimN_pnts, with N is 1 to
%%% number of index dimensions in data.
%%%
%%% Output:
%%%         status  :       #elements written to SDAT file. As the data is
%%%                         split into a real and imaginary part, count
%%%                         should be equal to numel(data)*2.
%%%
%%% Contact: qhoutum2@umcutrecht.nl


% If only file name and data are given SDAT related
if nargin == 2, data_type = 'fid'; end
% If extra data-labels are given SPAR related
if nargin == 4, data_dim_labels = varargin{1}; end


%% Write SDAT

% 1. Handle according input. Convert frequency data if necessary 
%    Send ND array of SPECTRA to return FIDS.
if strcmp(data_type, 'spec'), fid_of_spectra = csi_ifft(data);
else                          fid_of_spectra = data;
end

% 2. Create output file name and path
[fp,fn] = fileparts(fpn); fpn = [fp '\' fn]; % Excl. ext
    
% 3. Store size of data
data_size = size(fid_of_spectra);
% 4. Make list of data: loosing all dimensions
data = reshape( fid_of_spectra, 1, []);
% 5. Split imaginary and complex numbers numbers as seperate columns.
data =[ real(data); imag(data) ];
    
% 6. Write SDAT
file_id = fopen([fpn '.SDAT'],'w','ieee-le');         % file ID of SDAT file.
status  = fwriteVAXG(file_id, data, 'float32');       % Write as VAXG
fclose(file_id);                                      % Done.

% 7. Display info to user.
if numel(data) ~= status, warning('Saving SDAT-file failed!');
else, fprintf('SDAT-file written: %s\n',[fpn '.SDAT']);
end

% </end write SDAT>
% ======================================================================= %


%% WRITE SPAR
% This part writes an SPAR file with the SDAT file parameter lines 
% dim1_pnts allowing to reload and rearrange the data ranging 1:N.
if exist('data_dim_labels','var')
    csi_writeSPAR(fpn, data_size, data_dim_labels);
else
    csi_writeSPAR(fpn, data_size);
end
 
