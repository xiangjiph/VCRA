clc;clear;close all;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
grid_version = '240_cube';
DataManager = FileManager;
%%
grid_sub = [17 13 33];
vessel_image = DataManager.load_block_data(dataset_name, stack, grid_version, grid_sub(1), ...
    grid_sub(2), grid_sub(3));
vessel_mask = DataManager.load_block_mask(dataset_name, stack, grid_version, ...
    grid_sub(1), grid_sub(2), grid_sub(3));
vessel_mask = fun_reconstruct_block_mask(vessel_mask);
volumeViewer(vessel_mask);
%% Convert mask to skeleton to graph 
vessel_skl = bwskel(vessel_mask);
vessel_graph = fun_skeleton_to_graph(vessel_skl);
vessel_mask_dt = bwdist(~vessel_mask);
%% Graph refinement setting
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
% Other parameters
opt.skel_recentering_Q = false;
% Training data
TD = struct;
TD.info.filepath = sprintf('Test_training_data.mat');
opt.record.TD = TD;
%% Graph refinement 
[vessel_graph, opt.record] = fun_graph_refine_by_classifier(vessel_graph, vessel_image, vessel_mask_dt, opt);
%%
test_data = struct;
test_data.vessel_image = vessel_image;
test_data.vessel_mask = vessel_mask;
test_data.vessel_graph = vessel_graph;
test_data.graph_refinement_opt = opt;
save('Graph_refinemnt_test_data.mat', '-struct', 'test_data');