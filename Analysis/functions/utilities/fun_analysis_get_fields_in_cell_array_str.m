function [layer_stat_map, varargout]= fun_analysis_get_fields_in_cell_array_str(input_cell_array, extract_filed_name, field_is_scalar_Q)
%
%
%
%
if nargin < 3
    field_is_scalar_Q = true;
end
missing_field_list = {};
layer_stat_map = [];
if ~isempty(input_cell_array)
    array_size = size(input_cell_array);
    num_cell = prod(array_size);
    valid_idx_1 = find(~cellfun(@isempty, input_cell_array), 1);
    % Assert each nonempty cell contains a structure of the same fields
    test_str = input_cell_array{valid_idx_1};
    common_field_name = fieldnames(test_str);
    num_fields = numel(common_field_name);
    layer_stat_map = fun_initialized_structure_array_with_fieldname_list(common_field_name);
    for iter_field = 1 : num_fields
        if field_is_scalar_Q
            layer_stat_map.(common_field_name{iter_field}) = nan(array_size);
        else
            layer_stat_map.(common_field_name{iter_field}) = cell(array_size);
        end
    end
    for iter_cell = 1 : num_cell
        tmp_str = input_cell_array{iter_cell};
        
        if ~isempty(tmp_str)
            tmp_field_name = fieldnames(tmp_str);
            if isempty(tmp_field_name)
                continue;
            else
                for iter_field = 1 : num_fields
                    tmp_fn = common_field_name{iter_field};
                    if ismember(tmp_fn, tmp_field_name)
                        % This is the rate-limiting step
                        tmp_field_value = tmp_str.(tmp_fn).(extract_filed_name);
                        if isscalar(tmp_field_value) && field_is_scalar_Q
                            % Only extract scalar field value.
                            layer_stat_map.(tmp_fn)(iter_cell) = tmp_field_value;
                        elseif ~isempty(tmp_field_value)
                            layer_stat_map.(tmp_fn){iter_cell} = tmp_field_value;
                        end
                    else
                        fprintf('The field %s does not exist. The cell indices is %d\n', tmp_fn, iter_cell);
                        missing_field_list{end+1} = tmp_fn;
                        % Field does not exist.
                    end
                end
            end
        end
    end
end
if nargout > 1
    varargout{1} = missing_field_list;
end
end