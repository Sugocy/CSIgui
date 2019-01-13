function qry_ind = csi_findDimLabel(labels, query)
%%%% Description:               Returns index of query-strings in labels 
%%% Creator: Ir. Q. van Houtum       Version: 1.2          Date: 2017-07
%%% --------------------------------------------------------------------
%%%  
%%% Find index of specific dimension defined by query in data.raw from
%%% csi_loadData. Labels is a cell-string with a label per dimensions as
%%% returned from csi_loadData, data.labels
%%%
%%% Contact: qhoutum2@umcutrecht.nl

ksp_str = query; qry_ind = NaN(1,size(query,2));
for kk =1:size(query,2),
    if ~isempty(find(strcmp(labels, ksp_str{kk})))
        qry_ind(kk) = find(strcmp(labels, ksp_str{kk}));
    end 
end