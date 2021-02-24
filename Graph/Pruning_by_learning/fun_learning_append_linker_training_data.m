function str1 = fun_learning_append_linker_training_data(str1, str2)
% fun_learning_append_training_data append the list of str2 to the
% corresponding field of str1
% Input: 
%   str1, str2: structures defined in Annotation_pipeline_test_20190121
% Output:
%   str1: strucutre which merge the data, link_cc_ind and label in str1 and
%   str2( append str2 to str1)
if isfield(str1, 'data')
    str1.data = [str1.data; str2.data];
else
    error('Filed data does not exist');
end
if isfield(str1, 'raw_data')
    str1.raw_data = [str1.raw_data; str2.raw_data];
else
    error('Filed link_cc_ind does not exist');
end
if isfield(str1, 'label')
    fn_list = fieldnames(str1.label);
    for iter_fn = 1 : numel(fn_list)
        tmp_fn = fn_list{iter_fn};
        str1.label.(tmp_fn) = cat(1, str1.label.(tmp_fn), str2.label.(tmp_fn));
    end
else
    error('Filed link_cc_ind does not exist');
end
end