set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack = 'ML_2018_08_15';
% stack = 'ML_2019_01_24';
stack = 'ML20200201';
grid_c_version = '240_cube_combined_5_o_1';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c.grid_ori;
script_root = DataManager.SCRIPT_PATH;
% classifier_folder_path = fullfile(script_root, 'Graph', 'Pruning_by_learning', sprintf('Annotation_%s',stack));
classifier_folder_path = fullfile(DataManager.fp_metadata_folder(dataset_name, stack), 'classifier');
if ~isfolder(classifier_folder_path)
    mkdir(classifier_folder_path);
end
%% Load training data
% combined_grid_idx_list = [219 220 306 338 873]; % For ML_2018_08_15
% combined_grid_idx_list = [591, 598, 801]; % For ML20190124
combined_grid_idx_list = [351, 520, 476, 925, 763];
num_training_block = numel(combined_grid_idx_list);
training_data = cell(num_training_block, 1);
for iter_block = 1 : num_training_block
    fn_training_data = DataManager.fp_training_data(dataset_name, stack, ...
        sprintf('%s_%d', grid_c_version, combined_grid_idx_list(iter_block)));
    training_data{iter_block} = load(fn_training_data);
end
used_training_data_idx = 1:5;
num_block_used = numel(used_training_data_idx);
%% Merge data - normal link_ep1
link_ep1 = fun_learning_merge_annotation_data(training_data, 'link_ep1', 'normal');
linker = fun_learning_merge_annotation_data(training_data, 'linker', 'normal');
dim_short_link = fun_learning_merge_annotation_data(training_data, 'dim_short_link', 'normal');
%% Link with 1 endpoint
% The annotation is are all valid - save the classififer
% selected_feature_name = fun_learning_print_feature_name_list(link_ep1, [2, 4, 7, 10, 12, 21, 23, 25, 27, 29, 31, 33]);

% selected_feature_name = fun_learning_print_feature_name_list(link_ep1, [2, 4, 5, 6, 7, 10, 12, 21, 23, 25, 27, 29, 31, 33]);
selected_feature_name = intersect(link_ep1.features.Properties.VariableNames, {'length', 'dt_max', 'dt_min', ...
    'dt_mean', 'dt_std', 'int_median', 'int_diff_ep2ep', 'dt_max_2_length', ...
    'dt_std_n', 'int_std_n', 'int_SNR', 'nearest_ep_dist', 'nearest_link_num_voxel',...
    'inner_product_epv_epv'});
classifier_link_ep1 = fun_learning_get_classifier(link_ep1, selected_feature_name, 'Bag');
plot_opt = struct;
plot_opt.FileType = 'png';
plot_opt.Output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Classifier');
plot_opt.FileName = ['Classifier_link_ep1_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image_in_several_formats(classifier_link_ep1.validation_set.stat.confusion_matrix_fig, plot_opt.FilePath);
% Save classifier
classifier_link_ep1.filepath = fullfile(classifier_folder_path, 'classifier_link_ep1_to_remove.mat');
save(classifier_link_ep1.filepath, '-struct', 'classifier_link_ep1');
%% Linker
% linker_selected_feature_name = fun_learning_print_feature_name_list(linker, [3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18, 20, 21, 23, 24, 25, 28, 36]);
% linker_selected_feature_name = fun_learning_print_feature_name_list(linker, ...
%     [3, 4, 6, 7, 8, 9, 10, 12, 13, 16, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 35, 36]);
linker_selected_feature_name = intersect(linker.features.Properties.VariableNames, ...
    {'int_o_mask_mean', 'int_o_mask_std', 'recon_int_std', 'recon_int_mean', ...
    'recon_int_median', 'recon_bg_int_mean', 'recon_bg_std', 'recon_SNR', ...
    'recon_radius', 'bg_std', 'link_ratio_o_mask', 'linker_SNR', 'linker_o_m_SNR', ...
    'length', 'ep2ep_dist', 'ep2ep_dist_2_length', 'dir_vec_1_to_2', ...
    'closest_skl_dist_2_closest_ep_dist_mean', 'closest_skl_dist_2_closest_ep_dist_std', ...
    'num_nearest_link_cc', 'num_nearest_dt_change_sign', 'num_seg_in_mask', ...
    'link_ori_num_voxel', 'cos_ep2bv_link', 'link_has_2ep_Q'});

classifier_linker = fun_learning_get_classifier(linker, linker_selected_feature_name, 'Bag');
plot_opt.FileName = ['Classifier_linker_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image_in_several_formats(classifier_linker.validation_set.stat.confusion_matrix_fig, plot_opt.FilePath);
classifier_linker.filepath = fullfile(classifier_folder_path, 'classifier_linker_to_remove.mat');
save(classifier_linker.filepath, '-struct', 'classifier_linker');
%% Dim short link
% dim_short_link_feature_name = fun_learning_print_feature_name_list(dim_short_link, [2, 10, 13, 14, 40, 49, 50, 51, 52, 56]);
% dim_short_link_feature_name = fun_learning_print_feature_name_list(dim_short_link, ...
%     [2, 3, 4, 6, 7, 10, 12, 13, 14, 16, 17, 18, 22, 33, 36, 40, 41, 44, 45, 46, 49, 50, 51, 52, 56]);
dim_short_link_feature_name = intersect(dim_short_link.features.Properties.VariableNames, ...
    {'length', 'ep2ep_dist', 'dt_ep1', 'dt_max', 'dt_min', 'dt_std', 'int_mean', ...
    'int_median', 'int_min', 'int_middle_point', 'int_min_idx_n', 'int_min_from_mid', ...
    'int_std', 'dt_e2e_2_length', 'dt_std_n', 'dt_ep_sum_2_ep_dist', ...
    'int_std_n', 'int_mid_2_max', 'int_min_2_max', 'int_mid_2_max_at_end', ...
    'int_diff_epavg_mid_n', 'mid_point_SNR', 'mask_SNR', 'skl_SNR', ...
    'ep2ep_vec_wrt_z'});
classifier_dim_short_link = fun_learning_get_classifier(dim_short_link, dim_short_link_feature_name, 'Bag');
plot_opt.FileName = ['Classifier_dim_short_link_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image_in_several_formats(classifier_dim_short_link.validation_set.stat.confusion_matrix_fig, plot_opt.FilePath);
classifier_dim_short_link.filepath = fullfile(classifier_folder_path, 'classifier_dim_short_link_to_remove.mat');
save(classifier_dim_short_link.filepath, '-struct', 'classifier_dim_short_link');
%% Dim short link - debug
% The original features were messed up in [219 220 306 338 873] for
% ML_2018_08_15. 
% The feature is computed in Recompute_dim_short_link_features.m
% dim_short_link_featur_name = fun_learning_print_feature_name_list(tmp.normal, [2, 10, 13, 14, 40, 49, 50, 51, 52, 56]);
% classifier_dim_short_link = fun_learning_get_classifier(tmp.normal, [2, 4, 5, 10, 13, 14, 49, 50, 40, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 47], 'Bag');
% classifier_dim_short_link = fun_learning_get_classifier(tmp.normal, dim_short_link_featur_name, 'Bag');
% plot_opt = struct;
% plot_opt.FileType = 'png';
% plot_opt.Output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Classifier');
% plot_opt.FileName = ['Classifier_dim_short_link_bag_confusion_matrix_random_sample.', plot_opt.FileType];
% plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
% fun_print_image(classifier_dim_short_link.validation_set.stat.confusion_matrix_fig, plot_opt);
% classifier_dim_short_link.filepath = fullfile(classifier_folder_path, 'classifier_dim_short_link_to_remove.mat');
% save(classifier_dim_short_link.filepath, '-struct', 'classifier_dim_short_link');
%% Combine training set of two brain
ML2018_dim_short_link_classifier = load('./Graph/Pruning_by_learning/Annotation_ML_2018_08_15/classifier_dim_short_link_to_remove.mat');

dim_short_link_c = dim_short_link;
dim_short_link_c.features = [dim_short_link_c.features; ML2018_dim_short_link_classifier.data.features];
dim_short_link_c.label = [dim_short_link_c.label; ML2018_dim_short_link_classifier.data.label];
% dim_short_link_feature_name = fun_learning_print_feature_name_list(dim_short_link, ...
%     [2, 3, 4, 6, 7, 10, 12, 13, 14, 16, 17, 18, 22, 33, 36, 40, 41, 44, 45, 46, 49, 50, 51, 52, 56]);
% dim_short_link_feature_name = fun_learning_print_feature_name_list(dim_short_link_c, [2, 10, 13, 14, 40, 49, 50, 51, 52, 56]);
dim_short_link_feature_name = intersect(dim_short_link.features.Properties.VariableNames, ...
    {'length', 'ep2ep_dist', 'dt_median', 'dt_max', 'dt_min', 'dt_std', 'int_mean', ...
    'int_median', 'int_min', 'int_middle_point', 'int_min_idx_n', 'int_min_from_mid', ...
    'int_std', 'dt_e2e_2_length', 'dt_std_n', 'dt_ep_sum_2_ep_dist', ...
    'int_std_n', 'int_mid_2_max', 'int_min_2_max', 'int_mid_2_max_at_end', ...
    'int_diff_epavg_mid_n', 'mid_point_SNR', 'mask_SNR', 'skl_SNR', ...
    'ep2ep_vec_wrt_z'});
classifier_dim_short_link = fun_learning_get_classifier(dim_short_link_c, dim_short_link_feature_name, 'Bag');
plot_opt.FileName = ['Classifier_dim_short_link_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image_in_several_formats(classifier_dim_short_link.validation_set.stat.confusion_matrix_fig, plot_opt.FilePath);
classifier_dim_short_link.filepath = fullfile(classifier_folder_path, 'classifier_dim_short_link_to_remove.mat');
save(classifier_dim_short_link.filepath, '-struct', 'classifier_dim_short_link');
%% Classifier for link with 1 endpoint (artefact)
link_ep1_artefact = fun_learning_merge_annotation_data(training_data, 'link_ep1', 'artefact');
% Add the true negative from the normal dataset
link_ep1_artefact.features = cat(1, link_ep1_artefact.features, link_ep1.features(~link_ep1.label, :));
link_ep1_artefact.features.ep_dir_vec_z = abs(link_ep1_artefact.features.ep_direction_vec(:, 3));
link_ep1_artefact.label = cat(1, link_ep1_artefact.label, link_ep1.label(~link_ep1.label));
% selected_link_ep1_af_name = fun_learning_print_feature_name_list(link_ep1_artefact, [2, 4, 7, 10, 12, 21, 23, 25, 27, 29, 31, 36]);
% selected_link_ep1_af_name = fun_learning_print_feature_name_list(link_ep1_artefact, [2, 4, 5, 6, 7, 10, 12, 21, 23, 25, 27, 29, 31, 33, 36]);
selected_link_ep1_af_name = intersect(link_ep1.features.Properties.VariableNames, {'length', 'dt_max', 'dt_min', ...
    'dt_mean', 'dt_std', 'int_median', 'int_diff_ep2ep', 'dt_max_2_length', ...
    'dt_std_n', 'int_std_n', 'int_SNR', 'nearest_ep_dist', 'nearest_link_num_voxel',...
    'inner_product_epv_epv'});
classifier_link_ep1_af = fun_learning_get_classifier(link_ep1_artefact, selected_link_ep1_af_name, 'Bag');
plot_opt = struct;
plot_opt.FileType = 'png';
plot_opt.Output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Classifier');
plot_opt.FileName = ['Classifier_link_ep1_artefact_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image_in_several_formats(classifier_link_ep1_af.validation_set.stat.confusion_matrix_fig, plot_opt.FilePath);
% Save classifier
classifier_link_ep1_af.filepath = fullfile(classifier_folder_path, 'classifier_link_ep1_artefact_to_remove.mat');
save(classifier_link_ep1_af.filepath, '-struct', 'classifier_link_ep1_af');