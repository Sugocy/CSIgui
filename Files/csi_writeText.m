function csi_writeText(spec, fpn)
%%%% Description:                            Write spectrum to text file
%%% Creator: Ir. Q. van Houtum       Version: 1.0          Date: 2018-07
%%% --------------------------------------------------------------------
%%% csi_writeText( filename, FID or SPECTRA)
%%% Input:
%%%         spec    :       Either spectrum of fid data to be written.
%%%                         Expect the spectra at index 1.
%%%         fpn     :       Path and filename to text location
%%% 
%%% Contact: qhoutum2@umcutrecht.nl

if nargin == 1
    % Get filepath for destination
    [fn, fp, fi] = ...
        uiputfile({'*.txt', 'Text file'},'Save spectral data...');
    if fi == 0, return; end
    [~, fn,ext] = fileparts(fn);
else
    [fp, fn, ext] = fileparts(fpn);
end
if ~strcmp(fp(end),'\'), fp = [fp '\']; end

% Convert to a two column array
resh2exp = reshape(spec,[], 1); 
% Create complex and real part columns
dataR = real(resh2exp); dataI = imag(resh2exp);

format = '%32.32f %32.32f\n';
% Create 2 column array with real and imaginary part.
dataA = cat(2,dataR, dataI);


% Write array to file
fid = fopen([fp fn ext],'wt');
for ri = 1:size(dataA,1)
    fprintf(fid,format, dataA(ri,:)); 
end
fclose(fid);

% Also write away a size to reconstruct the data.
fid = fopen([fp fn '_ArraySize' ext],'wt');
fprintf(fid,'%32.32f ', size(spec)); 
fclose(fid);
