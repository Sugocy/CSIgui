function [val, foi] = getFieldValues(str, foi, case_sens, exact)
% Returns the values of [str]ucture.field name of interest or any field 
% that has the foi-string in it. Works with nested substructures.
%
% Input:
% str =         structure
% foi =         fieldname or string-in-fieldname of interest
% case_sens =   case sensitive on [0] or off [1], default is off.
% exact =       exact fieldname [1] or contains [0, default].
%
% Output:
% val =     cell-array with all values in the struct with "foi" in the 
%           fieldnames.
% foi =     all fieldnames it found containing foi including any 
%           substructures.
%
% Depencancy: findFieldnames(str, fname, case_sens, exact)
%
% quincyvanhoutum@gmail.com

if nargin < 3, case_sens = 1; end
if nargin < 4, exact = 0; end

% Get all fieldnames with string "foi" in it
foi = findFieldnames(str, foi, case_sens, exact);
foi_cut = cellfun(@strsplit, foi, repmat({'.'},size(foi)), 'Uniform', 0);

val = cell(1,numel(foi));
for kk = 1:numel(foi)
    nsubs = numel(foi_cut{kk});
    
    sub = str;
    for si = 1:nsubs-1, sub = sub.(foi_cut{kk}{si}); end
    val{kk} = sub.(foi_cut{kk}{nsubs});
end
