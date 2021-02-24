function selected_name = fun_learning_print_feature_name_list(data, idx)

data.feature_name = data.features.Properties.VariableNames;

if nargin < 2
    idx = 1 : numel(data.feature_name);
end

for iter_feature = 1 : numel(data.feature_name)
    fprintf('%d\t%s\n', iter_feature, data.feature_name{iter_feature});
end
selected_name = data.feature_name(idx);
end