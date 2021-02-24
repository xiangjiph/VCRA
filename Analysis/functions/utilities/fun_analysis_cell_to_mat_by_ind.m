function data_array = fun_analysis_cell_to_mat_by_ind(data_cell, cell_ind, array_size)

num_cell = numel(data_cell);
% Need to determine the content of the cell for initialization
is_valid_Q = ~cellfun(@isempty, data_cell);
valid_list_ind = find(is_valid_Q);
test_data = data_cell{valid_list_ind(1)};
issaclar_Q = isscalar(test_data);
if issaclar_Q
    data_array = nan(array_size);
else
    % Vector data. Put the data along the first indice
    data_array = nan(numel(test_data), array_size);
end

for iter_cell = 1 : num_cell
    if is_valid_Q(iter_cell)
        if isscalar_Q
            data_array(cell_ind(iter_cell)) = data_cell{iter_cell};
        else
            data_array(:, cell_ind(iter_cell)) = data_cell{iter_cell};
        end
    end
end
end