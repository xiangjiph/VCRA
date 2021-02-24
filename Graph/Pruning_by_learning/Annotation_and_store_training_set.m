% DataManager = FileManager;
% dataset_name = 'WholeBrain';
% stack = 'ML_2018_08_15';
% grid_c_version = '240_cube_combined_5_o_1';
% grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
% grid_info = grid_c.grid_ori;
% combined_grid_xyz_label = 151;
fn_training_data = DataManager.fp_training_data(dataset_name, stack, ...
    sprintf('%s_%d', grid_c_version, combined_grid_xyz_label));
classifier_root_path = fullfile(DataManager.SCRIPT_PATH, 'Graph', 'Pruning_by_learning', datestr(datetime('today')));
if ~isfolder(classifier_root_path)
    mkdir(classifier_root_path);
end
TD = struct;
TD.info = data_info;
template_training_data = struct;
template_training_data.features = [];
template_training_data.label = [];
template_training_data.raw_data = [];
% #######################
%% Link with one endpoint
% #######################
% Save annotation result
% annotate_link_ep1_label_1 = App_check_fake_gap_linker_result;
ep1_r1_str = fun_learning_get_training_data_from_annotation_label(annotate_link_ep1_label_1, tmp_ep1_features, int_link_0.ep1.link_cc_ind, false);

TD.ep1_normal_removed = [];
TD.ep1_artefact_removed = [];
TD.ep1_all_removed = [];

TD.ep1_normal_removed{end+1} = ep1_r1_str.normal;
TD.ep1_artefact_removed{end+1} = ep1_r1_str.artefact;
TD.ep1_all_removed{end+1} = ep1_r1_str.all;

% TD.ep1_normal_removed{end+1} = tmp.ep1_normal_removed{1};
% TD.ep1_artefact_removed{end+1} = tmp.ep1_broken_removed{1};
% TD.ep1_all_removed{end+1} = tmp.ep1_all_removed{1};
save(fn_training_data, '-struct', 'TD');
    %%  Test training
% Goal: low false positive rate ( no error delection ) is prefered. I.e. 
disp('Test performance of classifier and features');
% fun_learning_print_feature_name_list(TD.ep1_normal_removed{1})
test_ep1_feature = [2, 10, 21, 23, 27, 29, 31, 32, 33];
test_ep1_classifier = fun_learning_get_classifier(TD.ep1_normal_removed{1}, test_ep1_feature, 'AdaBoostM1');
disp(test_ep1_classifier.validation_set.stat);
    %%  Save/Load classifier
test_ep1_classifier.file_path = fullfile(classifier_root_path, 'classifier_link_ep1.mat');
save(test_ep1_classifier.file_path, '-struct', 'test_ep1_classifier');
%% Dim short link connecting capillaries
% To do
% 1. Specific function for candidate selection and feature computation
% 2. Need to incorporate the local mask information. Use distance transform
% to find the mask that associate to the link and compute its signal to
% noise level. 
% annotate_link_dim_short_label_1 = App_check_fake_gap_linker_result;
dim_short_str = fun_learning_get_training_data_from_annotation_label(annotate_link_dim_short_label_1, dim_short_link_features, vessel_graph_2.link.cc_ind(dim_short_link_label), false);

TD.link_dim_short_normal_removed = [];
TD.link_dim_short_normal_removed{end+1} = dim_short_str.normal;

TD.link_dim_short_all_removed = [];
TD.link_dim_short_all_removed{end+1} = dim_short_str.all;

TD.link_dim_short_artefact_removed = [];
TD.link_dim_short_artefact_removed{end+1} = dim_short_str.artefact;

save(fn_training_data, '-struct', 'TD');
    %%  Test training 
% fun_learning_print_feature_name_list(TD.link_dim_short_removed{1})
test_dim_link_feature = [13, 14, 37, 2, 18, 47, 49, 51, 52];
% test_dim_link_feature = [2, 4, 5, 10, 12, 13, 14, 15, 16, 18, 23, 24, 25, 26, 34, 37, 42, 45, 47];
disp('Test performance of classifier and features');
test_dim_link_classifier = fun_learning_get_classifier(TD.link_dim_short_removed{1}, test_dim_link_feature, 'Bag');
disp(test_dim_link_classifier.validation_set.stat);
% TD.link_dim_short_removed{1}.features( test_dim_link_classifier.validation_set.stat.false_positive_idx, :)
% debug_false_positive_link_label = dim_short_link_label(test_dim_link_classifier.validation_set.stat.false_positive_idx);
% tmp_debug = fun_graph_get_link_features(vessel_graph_2.link.cc_ind(debug_false_positive_link_label), vessel_image, vessel_mask_dt);
% App_check_fake_gap_linker(vessel_image, vessel_mask, vessel_skl_rc, vessel_graph_2.link.cc_ind(debug_false_positive_link_label))
% Conclusiong: rough estimation 
% TP 80% TN 99% FP 20% FN < 1% - because we have much more negative
% examples than positive examples. 
% histogram(table2array(dim_short_link_removed.features(~dim_short_link_removed.label, 49)))
    %%  Save/Load classifier
test_dim_link_classifier.file_path = fullfile(classifier_root_path, 'classifier_dim_short_link.mat');
save(test_dim_link_classifier.file_path, '-struct', 'test_dim_link_classifier');
%% Low signal to noise
% annotation_link_low_SNR_label_1 = App_check_fake_gap_linker_result;
low_snr_str = fun_learning_get_training_data_from_annotation_label(annotation_link_low_SNR_label_1, ...
    low_SNR_link_not_in_large_vessels_features, vessel_graph_2.link.cc_ind(low_SNR_link_not_in_large_vessels_label), false);

TD.link_low_SNR_thin_normal_removed = [];
TD.link_low_SNR_thin_all_removed = [];
TD.link_low_SNR_thin_artefact_removed = [];

TD.link_low_SNR_thin_all_removed{end+1} = low_snr_str.all;
TD.link_low_SNR_thin_normal_removed{end+1} = low_snr_str.normal;
TD.link_low_SNR_thin_artefact_removed{end+1} = low_snr_str.artefact;
save(fn_training_data, '-struct', 'TD');
    %% Test classifier
fun_learning_print_feature_name_list(TD.link_all_low_SNR_thin_removed{1})
test_link_low_SNR_thin_classifier = fun_learning_get_classifier(TD.link_all_low_SNR_thin_removed{1},...
    [2, 10, 13, 47, 49, 50], 'Bag');
disp(test_link_low_SNR_thin_classifier.validation_set.stat);
    %% Save classifier
test_link_low_SNR_thin_classifier.file_path = fullfile(classifier_root_path, 'classifier_low_SNR_thin_link.mat');
save(test_link_low_SNR_thin_classifier.file_path, '-struct', 'test_link_low_SNR_thin_classifier');
%% Link in short Loop
tmp_str = fun_learning_get_training_data_from_annotation_label(annotation_link_in_short_loops_label_1, link_in_short_loop_features, []);
TD.link_in_short_loop_normal_removed = [];
TD.link_in_short_loop_artefact_removed = [];
TD.link_in_short_loop_all_removed = [];
TD.link_in_short_loop_normal_removed{end+1} = tmp_str.normal;
TD.link_in_short_loop_artefact_removed{end+1} = tmp_str.artefact;
TD.link_in_short_loop_all_removed{end+1} = tmp_str.all;
save(fn_training_data, '-struct', 'TD');
    %%  Test training 
fun_learning_print_feature_name_list(TD.link_in_short_loop_all_removed{1})
test_link_in_short_loop_feature_idx = [2, 60, 58, 52, 48, 47, 45, 32];
% test_dim_link_feature = [2, 4, 5, 10, 12, 13, 14, 15, 16, 18, 23, 24, 25, 26, 34, 37, 42, 45, 47];
disp('Test performance of classifier and features');
test_short_loop_link_link_classifier = fun_learning_get_classifier(TD.link_in_short_loop_all_removed{1}, test_link_in_short_loop_feature_idx, 'Bag');
disp(test_short_loop_link_link_classifier.validation_set.stat);
% Performance: 
% For classifying all the annotated links, the performance of the current
% classifier is: 
% TP = 1, TN = 0.92, FP = 0, FN = 0.07 - good enough. 
% The false positive rate is the most important charactieristic here. 
    %% Save classifier
test_short_loop_link_link_classifier.file_path = fullfile(classifier_root_path, 'classifier_link_in_short_loop.mat');
save(test_short_loop_link_link_classifier.file_path, '-struct', 'test_short_loop_link_link_classifier');
%% Linker 
tmp_str = fun_learning_get_training_data_from_annotation_label(annotation_linker_label_1, linker_features, linker_str);
TD.linker_normal_removed = [];
TD.linker_artefact_removed = [];
TD.linker_all_removed = [];
TD.linker_normal_removed{end+1} = tmp_str.normal;
TD.linker_artefact_removed{end+1} = tmp_str.artefact;
TD.linker_all_removed{end+1} = tmp_str.all;
save(fn_training_data, '-struct', 'TD');
%% 
% Currently the positive sample is not enough. 
% fun_learning_print_feature_name_list(TD.liner_normal_removed{1})
test_linker_feature_idx = [7, 12, 22, 18, 26, 35];
disp('Test performance of classifier and features');
test_linker_classifier = fun_learning_get_classifier(TD.liner_normal_removed{1}, test_linker_feature_idx, 'Bag');
disp(test_linker_classifier.validation_set.stat);
    %% Save classifier
test_linker_classifier.file_path = fullfile(classifier_root_path, 'classifier_linker.mat');
save(test_linker_classifier.file_path, '-struct', 'test_linker_classifier');
    %% Visualization linker features
% figure;
% subplot(2,5,1)
% histogram(linker_features.recon_SNR(valid_linker_Q));
% hold on 
% histogram(linker_features.recon_SNR(annotation_linker_label_1.to_removeQ));
% subplot(2,5,2)
% histogram(linker_features.length(valid_linker_Q));
% hold on 
% histogram(linker_features.length(annotation_linker_label_1.to_removeQ));
% subplot(2,5,3)
% histogram(linker_features.closest_skl_dist_2_closest_ep_dist_mean(valid_linker_Q));
% hold on 
% histogram(linker_features.closest_skl_dist_2_closest_ep_dist_mean(annotation_linker_label_1.to_removeQ));
% subplot(2,5,4)
% histogram(linker_features.cos_ep2bv_link(valid_linker_Q));
% hold on 
% histogram(linker_features.cos_ep2bv_link(annotation_linker_label_1.to_removeQ));
% subplot(2,5,5)
% histogram(linker_features.link_ratio_o_mask(valid_linker_Q));
% hold on 
% histogram(linker_features.link_ratio_o_mask(annotation_linker_label_1.to_removeQ));
%% Link with one endpoint - round 2
% #######################
% Save annotation result
ep1_r2_str = fun_learning_get_training_data_from_annotation_label(annotation_link_ep1_label_2, link_ep1_feautres_2, int_link_2.ep1.link_cc_ind);
TD.ep1_normal_removed{end+1} = ep1_r2_str.normal;
TD.ep1_artefact_removed{end+1} = ep1_r2_str.artefact;
TD.ep1_all_removed{end+1} = ep1_r2_str.all;
save(fn_training_data, '-struct', 'TD');
%% Linker - round 2
tmp_str = fun_learning_get_training_data_from_annotation_label(annotation_linker_label_2, linker_features_2, linker_str_2);
TD.linker_normal_removed{end+1} = tmp_str.normal;
TD.linker_artefact_removed{end+1} = tmp_str.artefact;
TD.linker_all_removed{end+1} = tmp_str.all;
save(fn_training_data, '-struct', 'TD');
%% 