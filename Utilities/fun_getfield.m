function field_data = fun_getfield(field_data, field_name)

assert(isstruct(field_data) || istable(field_data) || isa(field_data, 'LinearModel'), 'The first input should be a structure or table.');
field_name_in_level = strsplit(field_name, '.');
num_field_level = numel(field_name_in_level);
for iter_level = 1 : num_field_level
    field_data = field_data.(field_name_in_level{iter_level});
end
end