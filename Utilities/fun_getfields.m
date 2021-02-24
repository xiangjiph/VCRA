function new_str = fun_getfields(data_str, field_name_cell)

new_str = struct;
num_field_name = numel(field_name_cell);

for iter_cell = 1 : num_field_name
   tmp_field_name = field_name_cell{iter_cell};
   new_str.(strrep(tmp_field_name, '.', '_')) = fun_getfield(data_str, tmp_field_name);
end
end