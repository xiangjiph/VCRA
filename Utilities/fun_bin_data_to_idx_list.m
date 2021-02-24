function [bin_cell_array, varargout] = fun_bin_data_to_idx_list(data) 
% fun_bin_data_to_idx_list bin the data according to their values and
% output the corresponding index list
% Input: 
%   data: numerical vector
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: unique data value 

if isempty(data)
    bin_cell_array = {};
    varargout{1} = [];
    return;
end
num_data = numel(data);
if ~issorted(data)
    [data, idx_list] = sort(data(:), 'ascend');
    % To be consistent with the previous behavior
    idx_list = idx_list.';
else
    idx_list = 1 : num_data;
end
est_num_bin = max(2, round(data(end) - data(1)) + 1);

if num_data < 65535
    idx_list = uint16(idx_list);
elseif num_data < 4294967295
    idx_list = uint32(idx_list);
end
%%
bin_value_list = zeros(est_num_bin,1);
bin_cell_array = cell(est_num_bin,1);
left_idx = 1;
bin_data = data(1);
num_bin = 0;
for right_idx = 1 : num_data
    tmp_data = data(right_idx);
    if tmp_data ~= bin_data
        num_bin = num_bin + 1;
        % copy 
        bin_cell_array{num_bin} = idx_list(left_idx : (right_idx - 1));
        bin_value_list(num_bin) = bin_data;
        % Update
        bin_data = tmp_data;
        % Re-initialize
        left_idx = right_idx;
    end
end
num_bin = num_bin + 1;
bin_cell_array{num_bin} = idx_list(left_idx : right_idx);
bin_cell_array(num_bin + 1 : end) = [];
bin_value_list(num_bin) = bin_data;
bin_value_list = bin_value_list(1 : num_bin);
if nargout > 1
    varargout{1} = bin_value_list;
end
%% Debug
% assert(sum(cellfun(@numel, bin_cell_array)) == num_data);
% bin_idx = cat(1, bin_cell_array{:});
% assert(numel(idx_list) == numel(bin_idx));
% assert(isempty(setdiff(bin_idx, idx_list)));
% assert(numel(unique(bin_idx)) == numel(bin_idx));
%% To do list:
% Test label2idx + histcount. 
end

        