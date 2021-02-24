function feature_array = fun_merge_structure_fields(input_str, field_name_list, merge_dim)
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
feature_array = cell(num_field, 1);
for iter_field = 1 : num_field
    tmp_field_name = field_name_list{iter_field};
    if isfield(input_str, tmp_field_name)
        feature_array{iter_field} = input_str.(tmp_field_name);
    else
        error('Field name %s does not exist', tmp_field_name);
    end
end
feature_array = cat(merge_dim, feature_array{:});
end