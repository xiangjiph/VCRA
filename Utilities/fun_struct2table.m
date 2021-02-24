function output_table = fun_struct2table(output_table)

assert(isstruct(output_table), 'The input should be a structure');
field_size = structfun(@(x) size(x, 1), output_table);
unique_field_size = unique(field_size);
assert(isscalar(unique_field_size), 'The number of column in fields of the structure are not the same');
if unique_field_size == 1
    output_table = struct2table(output_table, 'AsArray', true);
else
    output_table = struct2table(output_table);
end
end