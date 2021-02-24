function data_array = fun_analysis_ai_str_cell_to_mat(ai_str_cell, field_name, array_ind, array_size)

is_not_empty_ind = find(~cellfun(@isempty, ai_str_cell));
if isempty(is_not_empty_ind)
    data_array = [];
    return;
end
array_ind = array_ind(is_not_empty_ind);

test_data = ai_str_cell{is_not_empty_ind(1)}.(field_name);

is_scalar_Q = isscalar(test_data);
is_vector_Q = isvector(test_data);

if is_scalar_Q
    data_array = nan(array_size);
elseif is_vector_Q
    data_array = nan([numel(test_data), array_size]);
else 
    data_array = cell(array_size);
end

for iter_str = 1 : numel(is_not_empty_ind)
    tmp_data = ai_str_cell{is_not_empty_ind(iter_str)}.(field_name);
    
    if is_scalar_Q
        data_array(array_ind(iter_str)) = tmp_data;
    elseif is_vector_Q
        data_array(:, array_ind(iter_str)) = tmp_data;
    else
        data_array{array_ind(iter_str)} = tmp_data;
    end
    
end

end