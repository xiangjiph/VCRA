function str1 = fun_copy_fields_from_str2_to_str1(str1, str2)

fieldname2 = fieldnames(str2);
for iter_field = 1 : numel(fieldname2)
    str1.(fieldname2{iter_field}) = str2.(fieldname2{iter_field});
end
end
