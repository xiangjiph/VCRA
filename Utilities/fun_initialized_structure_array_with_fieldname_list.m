function str = fun_initialized_structure_array_with_fieldname_list(field_name_list, ini_val)
% fun_initialized_structure_array_with_fieldname_list generate a structure
% with given field name. The field value is initialized to be empty
% Input:
%   field_name_list: cell vector, each cell is a string array
% Output:
%   str: structure with fields
if nargin < 2
    ini_val = '[]';
end

num_field_name = numel(field_name_list);
arg_str = '';
for iter_fn = 1 : num_field_name
    if iter_fn ~= 1
        arg_str = arg_str + sprintf(", '%s', %s", field_name_list{iter_fn}, ini_val);
    else
        arg_str = arg_str + sprintf("'%s', %s", field_name_list{iter_fn}, ini_val);
    end
end
cmd_str = sprintf('struct(%s)', arg_str);
str = eval(cmd_str);
end