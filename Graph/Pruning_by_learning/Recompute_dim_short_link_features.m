set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
grid_c_version = '240_cube_combined_5_o_1';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c.grid_ori;
gpuDevice(2);
%%
combined_grid_idx_list = [219 220 306 338 873];
num_training_block = numel(combined_grid_idx_list);
training_data = cell(num_training_block, 1);
for iter_block = 1 : num_training_block
    fn_training_data = DataManager.fp_training_data(dataset_name, stack, ...
        sprintf('%s_%d', grid_c_version, combined_grid_idx_list(iter_block)));
    training_data{iter_block} = load(fn_training_data);
end
used_training_data_idx = 1:5;
num_block_used = numel(used_training_data_idx);
%%
test_str = training_data{5};
data_info = test_str.info;
combined_grid_mmxx_grid = test_str.info.combined_grid_mmxx_grid;
%% Load data
disp('Load mask');
tic
vessel_mask = DataManager.load_blocks_files('mask', dataset_name, stack, grid_c.grid_ori.version, ...
    combined_grid_mmxx_grid(1):combined_grid_mmxx_grid(4), combined_grid_mmxx_grid(2):combined_grid_mmxx_grid(5), ...
    combined_grid_mmxx_grid(3):combined_grid_mmxx_grid(6), 'logical');
toc
disp('Load image');
tic
vessel_image = DataManager.load_blocks_files('image', dataset_name, stack, grid_c.grid_ori.version, ...
    combined_grid_mmxx_grid(1):combined_grid_mmxx_grid(4), combined_grid_mmxx_grid(2):combined_grid_mmxx_grid(5), ...
    combined_grid_mmxx_grid(3):combined_grid_mmxx_grid(6), 'uint16');
toc
%% Parameters
opt = struct;
% opt.classifier = GR_classifier;
opt.max_self_loop_length = 30;
opt.internal_offset = 16;
opt.pruning_max_length = 2;
opt.max_bilink_loop_length = 15;

opt.select_dim_short_link = struct;
opt.select_dim_short_link.min_int_max = 22000;
opt.select_dim_short_link.dt_ep_sum_2_ep_dist_min = 0.15;
opt.select_dim_short_link.max_length = 25;

opt.select_low_SNR_link.skl_SNR_max = 10;
opt.select_low_SNR_link.int_min_max = 40000;
opt.select_low_SNR_link.dt_ep1_max = 3;
opt.select_low_SNR_link.dt_ep2_max = 3;

opt.select_short_loop.max_length = 150;

opt.annotation_on_Q = true;
internal_offset = opt.internal_offset;
max_rc_int = 40000;
image_size = size(vessel_image);
%% Generate graph
tic
disp('Skeletonization'); % ~60 sec for 1 tile
vessel_skl = bwskel(vessel_mask);
toc
tic
disp('Skeleton to graph'); % 6 sec
vessel_graph = fun_skeleton_to_graph(vessel_skl);
toc
% Use the size of the node to determine the existance of holes in the mask 
% Find the position where the node size is greater than 10
tic
disp('Check the existance of holes in the mask');
large_node_label = find(vessel_graph.node.num_voxel_per_cc >= 5);
if ~isempty(large_node_label)
    large_node_cc_ind = vessel_graph.node.cc_ind(large_node_label);
    num_large_node_cc = numel(large_node_cc_ind);
    tmp_add_ind = cell(num_large_node_cc, 1);
    tmp_node_bbox_expand = 10;
    for iter_node = 1 : num_large_node_cc
        tmp_ind = large_node_cc_ind{iter_node};
        tmp_sub = fun_ind2sub(vessel_graph.num.mask_size, tmp_ind);
        tmp_bbox_min = min(tmp_sub);
        tmp_bbox_max = max(tmp_sub);
        tmp_bbox_expand = [max([1,1,1], tmp_bbox_min - tmp_node_bbox_expand), ...
            min(vessel_graph.num.mask_size, tmp_bbox_max + tmp_node_bbox_expand)];
        tmp_bbox_size = tmp_bbox_expand(4:6) - tmp_bbox_expand(1:3) + 1;
        tmp_bbox = vessel_mask(tmp_bbox_expand(1) : tmp_bbox_expand(4), ...
            tmp_bbox_expand(2) : tmp_bbox_expand(5), tmp_bbox_expand(3) : tmp_bbox_expand(6));
        tmp_bbox_fill_hole = imfill(tmp_bbox, 'holes');
        tmp_diff_local_ind = find(tmp_bbox_fill_hole & ~tmp_bbox);
        if ~isempty(tmp_diff_local_ind)
            tmp_diff_local_sub = fun_ind2sub(tmp_bbox_size, tmp_diff_local_ind);
            tmp_diff_global_sub = tmp_diff_local_sub + tmp_bbox_expand(1:3) - 1;
            tmp_diff_global_ind = sub2ind(vessel_graph.num.mask_size, tmp_diff_global_sub(:, 1), ...
                tmp_diff_global_sub(:,2), tmp_diff_global_sub(:,3));
            tmp_add_ind{iter_node} = tmp_diff_global_ind;
        end
    end
    tmp_add_ind = cat(1, tmp_add_ind{:});
    if ~isempty(tmp_add_ind)
        disp('Update the mask and the skeleton');
        vessel_mask(tmp_add_ind) = true;
        vessel_skl = bwskel(vessel_mask);
    end
end
toc
tic
disp('Distance transform'); % 85 sec for 1 tile
vessel_mask_dt = bwdist(~vessel_mask);
toc
disp('Recentering');
tic
rc_metric = min(max_rc_int, vessel_image);
rc_metric = rc_metric + cast(ceil(vessel_mask_dt), class(rc_metric));
% It's funny that the recentering remove so many short hairs....
vessel_skl_rc = fun_skeleton_recentering_within_mask(vessel_skl, rc_metric, vessel_mask);
toc
vessel_graph = fun_skeleton_to_graph(vessel_skl_rc);
%% Remove the internal endpoints and recompute the features 
vessel_graph_rf = fun_graph_delete_hairs_and_short_loops(vessel_graph, vessel_image, opt);
rm_link_ep1_cc = cat(1, test_str.link_ep1.normal{1}.raw_data(test_str.link_ep1.normal{1}.label), ...
    test_str.link_ep1.not_sure{1}.raw_data, test_str.link_ep1.artefact{1}.raw_data);

rm_link_ep1_label_str = fun_graph_get_label_of_cc(vessel_graph_rf.link.map_ind_2_label, rm_link_ep1_cc);
vessel_graph_rf = fun_graph_pruning_by_link_label(vessel_graph_rf, rm_link_ep1_label_str.unique_cc_label_list);

rm_dim_short_link_label_str = fun_graph_get_label_of_cc(vessel_graph_rf.link.map_ind_2_label, test_str.dim_short_link.all{1}.raw_data);
dim_short_link_label = rm_dim_short_link_label_str.unique_cc_label_list;

tmp_vessel_graph = fun_analysis_get_connectivity_graph(vessel_graph_rf);
tmp_dim_short_link_loop = fun_analysis_get_loops_in_graph_by_link_label(tmp_vessel_graph, dim_short_link_label, 'euclidean');


dim_short_link_features = fun_graph_get_link_features(vessel_graph_rf.link.cc_ind(dim_short_link_label), vessel_image, vessel_mask_dt);

assert(numel(dim_short_link_label) == numel(tmp_dim_short_link_loop.link_label), 'Number of input links are not the same as the output links. Debug');
dim_short_link_features.shortest_loop_length = tmp_dim_short_link_loop.loop_length;
dim_short_link_features.shortest_loop_geodesic_length = tmp_dim_short_link_loop.loop_geodesic_length;
dim_short_link_features.length_ratio_in_shortest_loop = dim_short_link_features.length ./ dim_short_link_features.shortest_loop_length;
dim_short_link_features.ep2ep_vec_wrt_z = dim_short_link_features.ep1_to_ep2_direction_vec(:,3);
%% Find the connected links and compute features
dim_short_link_features.neighbor_links_int_mean = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_int_median = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_int_min = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_int_std = nan(size(dim_short_link_label));

dim_short_link_features.neighbor_links_length_mean = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_length_max = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_length_min = nan(size(dim_short_link_label));
dim_short_link_features.neighbor_links_length_std = nan(size(dim_short_link_label));
tic
for iter_link = 1 : numel(dim_short_link_label)
    test_link_label = dim_short_link_label(iter_link);

    test_connected_node = vessel_graph_rf.link.connected_node_label(test_link_label, :);
    assert(all(test_connected_node), 'Link connects to at least one endpoint');
    neighbor_links_of_node_1 = vessel_graph_rf.node.connected_link_label{test_connected_node(1)};
    neighbor_links_of_node_1 = neighbor_links_of_node_1(neighbor_links_of_node_1 ~= test_link_label);

    neighbor_links_of_node_2 = vessel_graph_rf.node.connected_link_label{test_connected_node(2)};
    neighbor_links_of_node_2 = neighbor_links_of_node_2(neighbor_links_of_node_2 ~= test_link_label);

    % Determine the end of the link segments that connected to the node
    node_1_neighbor_ind = fun_graph_get_node_neighbor_voxel_ind(vessel_graph_rf.node.cc_ind{test_connected_node(1)}, image_size, 26);
    node_2_neighbor_ind = fun_graph_get_node_neighbor_voxel_ind(vessel_graph_rf.node.cc_ind{test_connected_node(2)}, image_size, 26);
    % Use the neighbor link voxel to find: 
    % 1. Basic statistics of the neighbor links 
    neighbor_links_cc_ind_node_1 = vessel_graph_rf.link.cc_ind(neighbor_links_of_node_1);
    neighbor_links_cc_ind_node_2 = vessel_graph_rf.link.cc_ind(neighbor_links_of_node_2);

    neighbor_links_length = cellfun(@numel, cat(1, neighbor_links_cc_ind_node_1, neighbor_links_cc_ind_node_1));
    neighbor_links_ind = cat(1, neighbor_links_cc_ind_node_1{:}, neighbor_links_cc_ind_node_2{:});
    neighbor_links_int = single(vessel_image(neighbor_links_ind));

    dim_short_link_features.neighbor_links_length_max(iter_link) = max(neighbor_links_length);
    dim_short_link_features.neighbor_links_length_mean(iter_link) = mean(neighbor_links_length);
    dim_short_link_features.neighbor_links_length_std(iter_link) = std(neighbor_links_length);
    dim_short_link_features.neighbor_links_length_min(iter_link) = min(neighbor_links_length);

    dim_short_link_features.neighbor_links_int_mean(iter_link) = mean(neighbor_links_int);
    dim_short_link_features.neighbor_links_int_median(iter_link) = median(neighbor_links_int);
    dim_short_link_features.neighbor_links_int_std(iter_link) = std(neighbor_links_int);
    dim_short_link_features.neighbor_links_int_min(iter_link) = min(neighbor_links_int);
end
toc
%% Test the classification with the current training set
tmp = fun_learning_get_training_data_from_annotation_label(test_str.dim_short_link.all{1}.label, dim_short_link_features, ...
    test_str.dim_short_link.all{1}.raw_data);
fun_learning_print_feature_name_list(tmp.normal)
classifier_dim_short_link = fun_learning_get_classifier(tmp.normal, [2, 4, 5, 10, 13, 14, 49, 50, 40, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 47], 'Bag');
% classifier_dim_short_link = fun_learning_get_classifier(tmp.normal, [2, 10, 13, 14, 49, 50, 40], 'Bag');
classifier_dim_short_link.validation_set.stat


%%
% 2. Fit the neighbor pixels with a line, check the correlation
% 3. Distance between two lines, also normalized by their sum of radius 
% 4. Look what other features were used in the paper. 