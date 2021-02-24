function str = fun_structure_field_slicing_by_index(str, selected_idx, slicing_dim, field_name)
% fun_structure_filed_indexing select the data in the filed values
% accoridng to the selectedQ. For scalar
% Input: 
%   str: structure with fields of size N-by-M
%   index_list: logical array of size N-by-1
% 
if nargin < 3
    slicing_dim = 1;
elseif nargin < 4
    field_name = fieldnames(str);
end

for iter_fn = 1 : numel(field_name)
    fn = field_name{iter_fn};
    field_value = str.(fn);
    if isscalar(field_value)
        continue
    else
        switch slicing_dim
            case 1
                str.(fn) = field_value(selected_idx, :);
            case 2
                str.(fn) = field_value(:, selected_idx);
            case 3
                str.(fn) = field_value(:, :, selected_idx);
            otherwise
                fprintf('Field %s will not be masked\n', fn);
        end
    end
end

end



