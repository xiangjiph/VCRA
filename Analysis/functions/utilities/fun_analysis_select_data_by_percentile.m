function [data, selected_idx]= fun_analysis_select_data_by_percentile(data, ptl_low, ptl_high, remove_0_Q)
% fun_analysis_select_data_by_percentile remove nan in the data and select
% part of data according to the given percentile bound. 
%  
%
if nargin < 4
    remove_0_Q = false;
end
assert(ptl_low < ptl_high && ptl_high <= 100);
selected_Q = ~isnan(data);
if remove_0_Q
    selected_Q = selected_Q & (data ~= 0);
end
data = data(selected_Q);
selected_idx = find(selected_Q);
[data, sort_idx]= sort(data, 'ascend');
selected_idx = selected_idx(sort_idx);

num_data = numel(data);
idx_low = max(1, round(num_data * ptl_low / 100));
idx_high = min(num_data, round(num_data * ptl_high / 100));

data = data(idx_low : idx_high);
selected_idx = selected_idx(idx_low : idx_high);
end