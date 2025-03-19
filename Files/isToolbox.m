function tf = isToolbox(query)
% isToolbox returns 1 if true; query is an installed Toolbox, 0 if false;
% the query toolbox is not installed with Matlab.
%
% Case sensitive.

% Get all installed Toolboxes
v = ver; [installedToolboxes{1:length(v)}] = v(:).Name;
% Check if query is in installed Toolboxes.
tf = ismember(query, installedToolboxes);