function data_str = fun_analysis_sc_structure_fields(data_str, field_name_cell, scale_factor)

if isempty(field_name_cell) || scale_factor == 1
    return;
end
num_corr_field = numel(field_name_cell);
for iter_field = 1 : num_corr_field
    tmp_field_name = field_name_cell{iter_field};
    tmp_field_value = data_str.(tmp_field_name);
    if isempty(tmp_field_value)
        continue;
    elseif isnumeric(tmp_field_value)
        data_str.(tmp_field_name) = tmp_field_value .* scale_factor;
    elseif iscell(tmp_field_value)
        tmp_num_cell = numel(tmp_field_value);
        if tmp_num_cell >= 1
            if isnumeric(tmp_field_value{1})
                for iter_cell = 1 : tmp_num_cell
                    tmp_field_value{iter_cell} = tmp_field_value{iter_cell} .* scale_factor;
                end
                data_str.(tmp_field_name) = tmp_field_value;
            else
                error('Not supported datatype');
            end
        end
    else
        error('Not supported datatype');
    end
end
end