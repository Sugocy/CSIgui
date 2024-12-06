function val = getFields(str, foi)
% Given struct str and cell array with all fields of interest e.g. foi, it 
% will return ALL values for the given fieldnames if available in the
% struct.
%
% qhoutum2@umcutrecht.nl
val = {};
for fi = 1:size(foi,1)
    if isfield(str, foi{fi})
        val{fi,1} = str.(foi{fi});
    end
end

