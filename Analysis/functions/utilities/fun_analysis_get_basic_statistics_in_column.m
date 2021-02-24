function column_stat_str = fun_analysis_get_basic_statistics_in_column(data_mat, stat_dim)
% fun_analysis_get_basic_statistics_in_column compute the statistics of
% each column in the input data matrix
% Input: 
%   data: numerical matrix
%   stat_dim: 1 or 2, index of dimension along which the statistics is
%   computed
%
% Output: 
%   in_bin_stat: strcutre array, each structure is defined in
%   fun_analysis_get_basic_statistics
%
% Implemented by Xiang Ji on 04/16/2020
if nargin < 2
    stat_dim = 1;
end

switch stat_dim
    case 1
        % Do nothing
    case 2
        data_mat = data_mat.';
    otherwise
        error('stat_dim should be either 1 or 2');
end

num_column = size(data_mat, 2);
for iter_bin = num_column : -1 : 1
    tmp_data = data_mat(:, iter_bin);
    in_bin_stat(iter_bin) = fun_analysis_get_basic_statistics(tmp_data);
end
scalar_field_name = {'num_data', 'min', 'max', 'mean', 'median', 'sum', 'std', ...
    'range', 'hist_num_edge'};
array_field_name = setdiff(fieldnames(in_bin_stat), scalar_field_name);
% array_field_name = {'hist_count', 'hist_edge', 'hist_bin_val', ...
%     'hist_probability', 'hist_pdf', 'hist_bin_size', ...
%     'hist_cumulate_count', 'hist_cdf', 'prctile_th', 'prctile_val'};
column_stat_str = struct;
for iter_field = 1 : numel(scalar_field_name)
    tmp_field_name = scalar_field_name{iter_field};
    tmp_vec = nan(1, num_column);
    for iter_col = 1 : num_column
        tmp_scalar = in_bin_stat(iter_col).(tmp_field_name);
        if ~isempty(tmp_scalar)
           tmp_vec(iter_col) = tmp_scalar;
        end
    end
    column_stat_str.(tmp_field_name) = tmp_vec;
end

for iter_field = 1 : numel(array_field_name)
    tmp_field_name = array_field_name{iter_field};
    column_stat_str.(tmp_field_name) = {in_bin_stat.(tmp_field_name)};
end
column_stat_str.prctile_th = column_stat_str.prctile_th{1};
column_stat_str.prctile_val = cat(1, column_stat_str.prctile_val{:});
end