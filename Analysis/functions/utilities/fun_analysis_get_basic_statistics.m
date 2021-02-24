function stat_str = fun_analysis_get_basic_statistics(data, exactQ)
% fun_analysis_get_basic_statistics get the basic statistical properties
% about the input data
% Input: 
%   data: numerical vector. If it is not a vector, convert the data into
%   a vector automatically
% Output:
%   stat_str: structure with fields including num_data, min, max, etc. 
%    
% Implemented by Xiang Ji on 02/25/2019
%
if nargin < 2
    exactQ = true;
end

stat_str = fun_initialized_structure_array_with_fieldname_list({'num_data', 'min', 'max', 'mean', 'median', 'sum', 'std', 'range',...
    'hist_count', 'hist_edge', 'hist_bin_val', 'hist_probability', 'hist_pdf', 'hist_bin_size', 'hist_num_edge', 'hist_cumulate_count', ...
    'hist_cdf', 'prctile_th', 'prctile_val'});
if isempty(data)
    return;
end
if ~isvector(data)
    data = data(:);
end
% Remove NaN and Inf
data = data(isfinite(data));

prctile_th = [0, 0.1, 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.9, 100];

stat_str.num_data = double(numel(data));
stat_str.sum = double(sum(data));
stat_str.mean = double(stat_str.sum)/ double(stat_str.num_data);
stat_str.std = sqrt((double(sum(data.^2))/stat_str.num_data - stat_str.mean^2));

[stat_str.hist_count, stat_str.hist_edge] = histcounts(data);

hist_bin_size = stat_str.hist_edge(2:end) - stat_str.hist_edge(1:end-1);
stat_str.hist_bin_size = hist_bin_size;
stat_str.hist_num_edge = numel(stat_str.hist_edge);
stat_str.hist_cumulate_count = cumsum(stat_str.hist_count);
stat_str.hist_bin_val = movmean(double(stat_str.hist_edge), 2, 'Endpoints', 'discard');
stat_str.hist_probability = double(stat_str.hist_count) ./ stat_str.num_data;
stat_str.hist_pdf = stat_str.hist_probability ./ hist_bin_size;
stat_str.hist_cdf = stat_str.hist_cumulate_count ./stat_str.num_data;
if stat_str.num_data > 1e6 && stat_str.hist_num_edge > 10 && ~exactQ
    error('Legacy option');
    % If data size is large, use histogram to estiamte the max, min and median
    % to accelerate the computation.
%     stat_str.min = (stat_str.hist_edge(1) + stat_str.hist_edge(2))/2;
%     stat_str.max = (stat_str.hist_edge(end) + stat_str.hist_edge(end - 1))/2;
%     dist_to_median = stat_str.hist_cdf - 0.5;
%     [~, median_idx] = min(abs(dist_to_median));
%     stat_str.median = stat_str.hist_edge(median_idx);
else
    prctile_val = prctile(data, prctile_th);
    stat_str.min = prctile_val(1);
    stat_str.max = prctile_val(end);
    stat_str.median = double(median(data));
    stat_str.prctile_th = prctile_th;
    stat_str.prctile_val = prctile_val;
end
stat_str.range = stat_str.max - stat_str.min;
end