function [bin_cell_array, varargout] = fun_bin_data_to_idx_list_by_edges(data, edge_list, include_endpointQ) 
% fun_bin_data_to_idx_list_by_edges bin the data according to their values
% and the edge and output the corresponding index list 
% Input:
%   data: numerical vector
%   edge_list: numerical vector, edge for binning
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: average of the bin value
if nargin < 3
    include_endpointQ = true;
end

if isempty(data)
    bin_cell_array = {};
    varargout{1} = [];
    return;
end
num_data = numel(data);
if ~issorted(data)
    [data, idx_list] = sort(data(:), 'ascend');
else
    idx_list = 1 : num_data;
end
num_bin = numel(edge_list) - 1;
edge_left = edge_list(1:end-1);
edge_right = edge_list(2:end);
edge_mean = (edge_left + edge_right) ./ 2;
if num_data < 65535
    idx_list = uint16(idx_list);
elseif num_data < 4294967295
    idx_list = uint32(idx_list);
end
%% Linear implementation
% Scan the sorted data once and copy the subarray once each bin. 
bin_cell_array = cell(num_bin,1);
current_bin_idx = 1;
tmp_left_idx = 1;
tmp_right_idx = 1;
while (tmp_right_idx <= num_data) && (current_bin_idx <= num_bin)
    % Include the data in the left bin, but not right bin
    % Search 
    while (tmp_right_idx <= num_data) && data(tmp_right_idx) < edge_right(current_bin_idx)
        if data(tmp_right_idx) >= edge_left(current_bin_idx)            
            tmp_right_idx = tmp_right_idx + 1;
        else
            tmp_left_idx = tmp_left_idx + 1;
            tmp_right_idx = tmp_left_idx;
        end
    end
    
    if tmp_right_idx <= num_data && current_bin_idx == num_bin 
        if include_endpointQ
            while tmp_right_idx <= num_data && data(tmp_right_idx) == edge_right(current_bin_idx)
                tmp_right_idx = tmp_right_idx + 1;
            end
        end
    end
    
    if (tmp_right_idx - tmp_left_idx) > 0
        bin_cell_array{current_bin_idx} = idx_list(tmp_left_idx : (tmp_right_idx - 1));
        tmp_left_idx = tmp_right_idx;
%         tmp_right_idx = tmp_right_idx + 1;
    end
    current_bin_idx = current_bin_idx + 1;
end

if nargout == 2
    varargout{1} = edge_mean;
elseif nargout == 3
    varargout{1} = edge_mean;
    data_cell = cell(size(bin_cell_array));
    for iter_cell = 1 : numel(bin_cell_array)
        data_cell{iter_cell} = data(bin_cell_array{iter_cell});
    end    
    varargout{2} = data_cell;
end
end     