function feature_array = fun_merge_table_fields(input_table, field_name_list)
% fun_merge_structure_fields merges the field of the structure specified by
% the give field name array. Assume the number of elements in the merge
% direction are the same. 
% Input: 
%   input_str: structure with fields
%   field_name_list: cell array of strings, list of field name for merging
% Output: 
%   feature_array: N-by-M numerical array, where the (merge_dim) direction
%   is the data from different fields. 
num_field = numel(field_name_list);
field_idx = zeros(num_field, 1);
for iter_field = 1 : num_field
    tmp_field_name = field_name_list{iter_field};
    tmp_idx = find(strcmp(tmp_field_name, input_table.Properties.VariableNames));
    if isscalar(tmp_idx)
        field_idx(iter_field) = tmp_idx;
    elseif isempty(tmp_idx)
        error('Column name %s does not exist', tmp_field_name);
    else
        error('Column name %s is not unique', tmp_field_name);
    end
end
feature_array = table2array(input_table(:, field_idx));
end