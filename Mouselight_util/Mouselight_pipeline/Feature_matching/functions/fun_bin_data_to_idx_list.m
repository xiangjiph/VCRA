function [bin_cell_array, varargout] = fun_bin_data_to_idx_list(data, est_num_bin)%codegen
% fun_bin_data_to_idx_list bin the data according to their values and
% output the corresponding index list
% Input: 
%   data: numerical vector, assume to be integer
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: unique data value 

num_data = numel(data);
if ~issorted(data)
    [data, idx_list]= sort(data, 'ascend');
else
    idx_list = 1 : num_data;
end

if nargin < 2
    est_num_bin = min(2, round(data(end) - data(1)));
end
bin_size = 0;
bin_data = data(1);
bin_idx = zeros(1, round(num_data/2));
bin_value_list = zeros(est_num_bin,1);
bin_value_list(1) = data(1);
bin_cell_array = cell(est_num_bin,1);
num_bin = 0;
for idx = 1 : num_data
    tmp_data = data(idx);
    if tmp_data == bin_data
        bin_size = bin_size + 1;
        bin_idx(bin_size) = idx_list(idx);
    else
        num_bin = num_bin + 1;
        bin_cell_array{num_bin} = bin_idx(1 : bin_size);
        bin_data = tmp_data;
        bin_value_list(num_bin + 1) = bin_data;
        bin_idx(1) = idx_list(idx);
        bin_size = 1;
    end
end
num_bin = num_bin + 1;
bin_cell_array{num_bin} = bin_idx(1 : bin_size);
bin_cell_array(num_bin + 1 : end) = [];
bin_value_list = bin_value_list(1 : num_bin);
if nargout > 1
    varargout{1} = bin_value_list;
end
end

        