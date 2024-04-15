function flds = findFieldnames(struct_inp, query, case_sens)
% Find all fieldnames in a structure with multiple substructs and returning
% every fieldname containing the query.
%
% Input:
% struct_inp =  structure
% query =       string with (fieldname) string of interest.
% case_sens=    case sensitive search on [0] or off [1], default is off.
%
% Output:
% flds =        fields with "query" in it including any nested
%               substructures.
%
% quincyvanhoutum@gmail.com

if nargin < 3, case_sens = 1; end

% Fieldnames of structure
fn = fieldnames(struct_inp); 

% If any fieldnames contain query, add to output.
flds = fn(contains(fn, query,'IgnoreCase',case_sens)); 

% Loop each field to find nested structures and find query in those
% fieldnames.
for kk = 1:size(fn,1)    
    if isstruct( struct_inp.(fn{kk}) )
        % Field in struct_inp is a structure - read its fields
        tmp = findFieldnames(struct_inp.(fn{kk}), query, case_sens);
        
        if ~isempty(tmp)
            for ff = 1:size(tmp,1)
                flds = [flds; strjoin([fn{kk} tmp(ff)],'.')];
            end
        end
    end
end
