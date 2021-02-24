classifier_folder_path = fullfile('Demo', 'Data');
if ~isfolder(classifier_folder_path)
    mkdir(classifier_folder_path);
end
%% Load training data
num_training_block = 1;
training_data = cell(num_training_block, 1);
for iter_block = 1 : num_training_block
    training_data{iter_block} = load('Demo_annotation_data.mat');
end
used_training_data_idx = 1:1;
num_block_used = numel(used_training_data_idx);
%% Merge data - normal link_ep1
link_ep1 = fun_learning_merge_annotation_data(training_data, 'link_ep1', 'normal');
linker = fun_learning_merge_annotation_data(training_data, 'linker', 'normal');
dim_short_link = fun_learning_merge_annotation_data(training_data, 'dim_short_link', 'normal');
%% Link with 1 endpoint
selected_feature_name = fun_learning_print_feature_name_list(link_ep1, [2, 4, 5, 6, 7, 10, 12, 21, 23, 25, 27, 29, 31, 33]);
classifier_link_ep1 = fun_learning_get_classifier(link_ep1, selected_feature_name, 'Bag');
% Save classifier
classifier_link_ep1.filepath = fullfile(classifier_folder_path, 'classifier_link_ep1_to_remove.mat');
save(classifier_link_ep1.filepath, '-struct', 'classifier_link_ep1');
%% Linker
linker_selected_feature_name = fun_learning_print_feature_name_list(linker, ...
    [3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 35, 36]);
classifier_linker = fun_learning_get_classifier(linker, linker_selected_feature_name, 'Bag');
classifier_linker.filepath = fullfile(classifier_folder_path, 'classifier_linker_to_remove.mat');
save(classifier_linker.filepath, '-struct', 'classifier_linker');
%% Dim short link
dim_short_link_feature_name = fun_learning_print_feature_name_list(dim_short_link, ...
    [2, 3, 4, 6, 7, 10, 12, 13, 14, 16, 17, 18, 22, 33, 36, 40, 41, 44, 45, 46, 49, 50, 51, 52, 56]);
classifier_dim_short_link = fun_learning_get_classifier(dim_short_link, dim_short_link_feature_name, 'Bag');
classifier_dim_short_link.filepath = fullfile(classifier_folder_path, 'classifier_dim_short_link_to_remove.mat');
save(classifier_dim_short_link.filepath, '-struct', 'classifier_dim_short_link');