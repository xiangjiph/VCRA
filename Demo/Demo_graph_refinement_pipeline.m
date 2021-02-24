%% Graph refinement pipeline: 
% Required input: 
%   vessel_image: 3D numerical array
%   vessel_mask: 3D logical array (vessel voxel = true)
% Running precedure: 
% 1. Generate training data
%   Set opt.annotation_on_Q = true. GUI for annotation will pop up
%   automatically along the way of refinement. Click the "Save result"
%   button after finishing annotation. 
% 2. Train classifiers
%   Run Learning_train_graph_refinement_classifier.m 
%   We only used parts of the computed features for training classifiers at
%   the moment. 
% 3. Automatic graph refinement
%   Load trained classifiers from folder, set opt.annotation_on_Q = false. 
%   The function can called iteratively for multiple run of refinement. 
%
% Implemented by Xiang Ji (Kleinfeld Laboratory, Department of Physics, UC San Diego)
% Tested in MATLAB R2019b on Win 10 and CentOS 7.0
% If you have any question, please email xiangjiph@gmail.com
%%
addpath(genpath('.'));
clc;clear;close all;
test_data = load('Graph_refinemnt_test_data.mat');
%%
vessel_image = test_data.vessel_image;
vessel_mask = test_data.vessel_mask;
% volumeViewer(vessel_mask);
%% Convert mask to skeleton to graph 
vessel_skl = bwskel(vessel_mask);
vessel_graph = fun_skeleton_to_graph(vessel_skl);
vessel_mask_dt = bwdist(~vessel_mask);
%% Graph refinement setting - for annotation
opt = struct;
opt.output_graph_name = sprintf('%s_annotated', 'test_data_graph');
opt.output_skel_name = sprintf('%s_annotated', 'test_data');
% Classifier
opt.classifier = [];
% Parameter for automatic hair, self-loop and bilink loop removal
opt.max_self_loop_length = 30;
opt.internal_offset = 16;
opt.pruning_max_length = 2;
opt.max_bilink_loop_length = 15;
% Parameter for selecting candidate link for classifier
opt.select_dim_short_link = struct;
opt.select_dim_short_link.min_int_max = 22000;
opt.select_dim_short_link.dt_ep_sum_2_ep_dist_min = 0.15;
opt.select_dim_short_link.max_length = 25;

opt.select_short_loop.max_length = 150;
% Parameters for controling the conditions 
opt.annotation_on_Q = true; 
% Training data
TD = struct;
TD.info.filepath = fullfile('Demo', 'Data', 'Demo_annotation_data.mat');
opt.record.TD = TD;
%% Graph refinement 
[vessel_graph, opt.record] = fun_graph_refine_by_classifier(vessel_graph, vessel_image, vessel_mask_dt, opt);
%% Train Random forest classifiers for automatic graph refinement 
% Run Learning_train_graph_refinement_classifier.m 

% As in the test data, only one linker was annotated, there was no way to
% train the classifier. Nevertheless, the the script has been tested. 
% The following part uses the classifiers trained on our ex vivo images for
% illustration
%%
classifier_folder = fullfile('Demo', 'Data', 'Pretrained_classifier');
GR_classifier = struct;
GR_classifier.link_ep1 = load(fullfile(classifier_folder, 'classifier_link_ep1_to_remove.mat'), 'used_feature_name', 'classifier');
GR_classifier.dim_short_link = load(fullfile(classifier_folder, 'classifier_dim_short_link_to_remove.mat'), 'used_feature_name', 'classifier');
GR_classifier.linker = load(fullfile(classifier_folder, 'classifier_linker_to_remove.mat'), 'used_feature_name', 'classifier');
%% Graph refinement setting - for automatic refinement
opt = struct;
opt.output_graph_name = sprintf('%s_annotated', 'test_data_graph');
opt.output_skel_name = sprintf('%s_annotated', 'test_data');
% Classifier
opt.classifier = GR_classifier;
% Parameter for automatic hair, self-loop and bilink loop removal
opt.max_self_loop_length = 30;
opt.internal_offset = 16;
opt.pruning_max_length = 2;
opt.max_bilink_loop_length = 15;
% Parameter for selecting candidate link for classifier
opt.select_dim_short_link = struct;
opt.select_dim_short_link.min_int_max = 22000;
opt.select_dim_short_link.dt_ep_sum_2_ep_dist_min = 0.15;
opt.select_dim_short_link.max_length = 25;

opt.select_short_loop.max_length = 150;
% Parameters for controling the conditions 
opt.annotation_on_Q = false;
% Training data
opt.record.TD = [];
%% Automatic graph refinement
% Two iterations
for iter = 1 : 2
    [vessel_graph, opt.record] = fun_graph_refine_by_classifier(vessel_graph, vessel_image, vessel_mask_dt, opt);
end
fprintf('Finish automatic graph refinement\n');
fprintf('Number of deleted self-loops: %d\n', opt.record.num_self_loop);
fprintf('Number of deleted short bi-link loops: %d\n', opt.record.num_bilink_loop);
fprintf('Number of deleted false positive connections due to low z-resolution: %d\n', opt.record.num_dim_link);
fprintf('Number of deleted links with one internal endpoint: %d\n', opt.record.num_link_ep1 + opt.record.num_short_link_ep1);
fprintf('Number of linkers added: %d\n', numel(opt.record.added_linker_ind_w_ep));
fprintf('Number of remaining internal endpoints: %d\n', opt.record.int_link.num.ep); 