function stats = csi_statistics_of_volume(data, filter, filter_val)
%%% Calculate some general statistics from a volume of data. Will return a
%%% structure with the following fields:
%%%
%%% mean, std, mode, freq, median, min, min_index, max, max_index.
%%%
%%% If filter and filter_value are given, will filter statistics accordingly.
%%% Example: given the SNR of each point in the data-array, will use only
%%% data with SNR >= filter_val.
%%%
%%% Will ingore NaN values.
%%%
%%% Dr. Quincy van Houtum, ELH-Institute Essen, v2 - 2024/06.
%%% quincyvanhoutum@gmail.com

% Process input
if nargin == 1, filter = ones(size(data)); filter_val = 1; end

% Output struct.
stats = struct;

% Set filter-out values to NaN.
data(filter < filter_val) = NaN;

% Mean +/- std
stats.mean = mean(data(:), 'omitnan'); 
stats.std = std(data(:), 'omitnan');

% Mode: rounded by 100th decimal
tmp = round(data.*100)./100; 
[stats.mode, stats.freq, C] = mode(tmp(~isnan(tmp))); 
if size(C{1},1) > 1
    stats.mode = C{1}'; 
    stats.freq = repmat(stats.freq,1,size(stats.mode,2)); 
end
clear('tmp');

% Median
stats.median = median(data(:), 'omitnan');

% Min/Max and location in array
tmp_ind = cell(1,numel(size(data)));

[stats.min, indMn] = min(data(:));
[tmp_ind{:}] = ind2sub(size(data),indMn);
stats.min_index = (tmp_ind);

[stats.max, indMx] = max(data(:));
[tmp_ind{:}] = ind2sub(size(data),indMx);
stats.max_index = (tmp_ind);

% Volume size
stats.size = size(data);

% Number of elements
stats.num_elements = numel(data);

% Number nans
stats.num_nan = sum(isnan(data(:)));

% Number of decimals
stats.num_dec_zeros = numzeros(stats.median);

