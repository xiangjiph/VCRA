function str = fun_structure_field_indexing(str, selectedQ, field_name)
% fun_structure_filed_indexing select the data in the filed values
% accoridng to the selectedQ. For scalar
% Input: 
%   str: structure with fields of size N-by-M
%   index_list: logical array of size N-by-1
% 
if ~islogical(selectedQ) || ~isvector(selectedQ)
    error('The indexing array must be a logical vector');
end
list_size = numel(selectedQ);
if nargin < 3
    field_name = fieldnames(str);
end
for iter_fn = 1 : numel(field_name)
    fn = field_name{iter_fn};
    field_value = str.(fn);
    if isscalar(field_value)
        continue
    else
        field_value_size = size(field_value);
        same_size_dim = find(field_value_size == list_size);
        if numel(same_size_dim) > 1
            same_size_dim = same_size_dim(1); 
        elseif isempty(same_size_dim)
            continue;
        end
        num_selected_elem = nnz(selectedQ);
        target_field_size = field_value_size;
        target_field_size(same_size_dim) = num_selected_elem;
        switch same_size_dim
            case 1
                str.(fn) = field_value(selectedQ, :);
            case 2
                str.(fn) = field_value(:, selectedQ);
            case 3
                str.(fn) = field_value(:, :, selectedQ);
            otherwise
                fprintf('Field %s will not be masked\n', fn);
        end
        str.(fn) = reshape(str.(fn), target_field_size);
    end
end

end



