function merged_data = fun_learning_merge_annotation_data(data_cell, object_class, data_class)
% fun_learning_merge_annotation_data merges the annotation data

if ~isa(data_cell, 'cell')
    data_cell = {data_cell};
end

num_block = numel(data_cell);
merged_data = struct;
merged_data.features = [];
merged_data.label = [];

for iter_block = 1 : num_block
    if isfield(data_cell{iter_block}, object_class)
        tmp_data = data_cell{iter_block}.(object_class).(data_class);
    else
        continue;
    end
    tmp_num_round = numel(tmp_data);
    for iter_round = 1 : tmp_num_round
        if isempty(merged_data.features)
            merged_data.features = tmp_data{iter_round}.features;
            merged_data.label = tmp_data{iter_round}.label;
        else
            var_name_1 = merged_data.features.Properties.VariableNames;
            var_name_2 = tmp_data{iter_round}.features.Properties.VariableNames;
            var_name_union = union(var_name_1, var_name_2);
            if numel(var_name_union) == numel(var_name_1) && numel(var_name_union) == numel(var_name_2)
                merged_data.features = cat(1, merged_data.features, tmp_data{iter_round}.features(:, var_name_1));
            else
               var_name_intersect = intersect(var_name_1, var_name_2);
               var_name_missing_in_1 = setdiff(var_name_2, var_name_intersect);
               var_name_missing_in_2 = setdiff(var_name_1, var_name_intersect);
               size_1 = size(merged_data.features);
               size_2 = size(tmp_data{iter_round}.features);
               % Use NAN to fill in the missing variables in both tables
               new_table = cat(1, merged_data.features(:, var_name_intersect), tmp_data{iter_round}.features(:, var_name_intersect));
               for iter_missing_1 = 1 : numel(var_name_missing_in_1)
                   new_table.(var_name_missing_in_1{iter_missing_1}) = cat(1, nan(size_1(1), 1), ...
                       table2array(tmp_data{iter_round}.features(:, var_name_missing_in_1{iter_missing_1})));
               end
               for iter_missing_2 = 1 : numel(var_name_missing_in_2)
                   new_table.(var_name_missing_in_2{iter_missing_2}) = cat(1, ...
                       table2array(merged_data.features(:, var_name_missing_in_2{iter_missing_2})), ...
                       nan(size_2(1), 1));
               end
               merged_data.features = new_table;
            end
            merged_data.label = cat(1, merged_data.label, tmp_data{iter_round}.label);
        end  
    end
end
end