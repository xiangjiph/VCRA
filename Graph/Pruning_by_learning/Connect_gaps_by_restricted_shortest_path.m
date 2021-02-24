clc;clear;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'mouselight_1';
grid_c_version = '240_cube_combined_5';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
gpuDevice(2);
%% Load data
% fun_vis_combined_grid(dataset_name, stack, grid_c_version);
combined_grid_xyz_label = 173;
combined_grid_idx_1 = grid_c.bbox_grid_sub_list(combined_grid_xyz_label, 1);
combined_grid_idx_2 = grid_c.bbox_grid_sub_list(combined_grid_xyz_label, 2);
combined_grid_layer = grid_c.bbox_grid_sub_list(combined_grid_xyz_label, 3);
combined_grid_xy_ind = grid_c.bbox_xy_linear_idx_mat{combined_grid_layer}(combined_grid_idx_1, combined_grid_idx_2);
combined_grid_mmxx_grid = grid_c.bbox_xyz_mmxx_grid{combined_grid_layer}(combined_grid_xy_ind,:);
% Reduce the block size for local test
% combined_grid_mmxx_grid(1:3) = combined_grid_mmxx_grid(1:3) + 1;
% combined_grid_mmxx_grid(4:6) = combined_grid_mmxx_grid(4:6) - 1;
disp('Load segmentation');
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

fn_training_data_link_ep1 = DataManager.fp_training_data(dataset_name, stack, ...
    sprintf('link_ep1_%s_%d', grid_c_version, combined_grid_xyz_label));
fn_training_data_linker = DataManager.fp_training_data(dataset_name, stack, ...
    sprintf('linker_%s_%d', grid_c_version, combined_grid_xyz_label));
fn_training_data_fake_link = DataManager.fp_training_data(dataset_name, stack, ...
    sprintf('fake_link_%s_%d', grid_c_version, combined_grid_xyz_label));

fn_classifier_link_ep1 = DataManager.fp_classifier(dataset_name, stack, ...
    sprintf('link_ep1_%s_%d', grid_c_version, combined_grid_xyz_label));
fn_classifier_linker = DataManager.fp_classifier(dataset_name, stack, ...
    sprintf('linker_%s_%d', grid_c_version, combined_grid_xyz_label));
fn_classifier_fake_link = DataManager.fp_classifier(dataset_name, stack, ...
    sprintf('fake_link_%s_%d', grid_c_version, combined_grid_xyz_label));


%% Convert skeleton to graph 
tic
disp('Skeletonization'); % ~60 sec for 1 tile
vessel_skl = bwskel(vessel_mask);
toc
tic
disp('Distance transform'); % 85 sec for 1 tile
vessel_mask_dt = bwdist(~vessel_mask);
toc
disp('Skeleton recentering'); % ~22 sec 
tic
vessel_skl_rc = fun_skeleton_recentering(vessel_skl, vessel_image);

toc
tic
disp('Skeleton to graph'); % 6 sec
vessel_graph = fun_skeleton_to_graph(vessel_skl_rc);
toc
%% Pruning
% 1. Remove all the single voxel link with 1 endpint
internal_offset = 16;
% Remove 581 links here
vessel_graph_refine =  fun_graph_pruning_internal_short_hairs(vessel_graph, 2, internal_offset); 
int_link_0 = fun_graph_get_free_link(vessel_graph_refine, internal_offset);
% Remove links by pre-trained classifier
fake_link_classifier_str = load('./Graph/Pruning_by_learning/classifier_link_w_ep1_vol_length_gt_2_20190116.mat');
fake_link_classifier_str.training_features{end} = 'inner_product_epv2_dispv';
fake_link_classifier = fake_link_classifier_str.classifier;
% Get features
int_link_0.ep1.features = fun_graph_get_link_w_1_ep_features(int_link_0.ep1.link_cc_ind, vessel_image, vessel_mask_dt);
% Merge features used for classification
int_link_0.ep1.classification_feature = fun_merge_table_fields(int_link_0.ep1.features, fake_link_classifier_str.training_features);
% Classification
int_link_0.ep1.to_remove_Q = fake_link_classifier.predict(int_link_0.ep1.classification_feature);
% Check the removed link with internal endpoint 
% int_link_0.ep1.removed_link_idx = find(int_link_0.ep1.to_remove_Q);
% int_link_0.ep1.removed_link_length = int_link_0.ep1.length(int_link_0.ep1.removed_link_idx );
% [tmp_l, tmp_idx] = sort(int_link_0.ep1.removed_link_length , 'descend');
% App_check_fake_gap_linker(vessel_image, vessel_mask, vessel_skl_rc, int_link_0.ep1.link_cc_ind(int_link_0.ep1.removed_link_idx(tmp_idx)))
%% Remove fake links labeled by the classifier
% Remove 494 links here
vessel_graph_refine_2 = fun_graph_pruning_by_link_label(vessel_graph_refine, int_link_0.ep1.link_label(int_link_0.ep1.to_remove_Q));
% Update the skeleton 
vessel_skel_refine = false(vessel_graph_refine_2.num.mask_size);
vessel_skel_refine([vessel_graph_refine_2.link.pos_ind; vessel_graph_refine_2.node.pos_ind; ...
    vessel_graph_refine_2.isopoint.pos_ind]) = true;
vessel_mask_rc = vessel_mask | vessel_skel_refine;
int_link = fun_graph_get_free_link(vessel_graph_refine_2, internal_offset);
%% Connect gaps - Link endpoint to nearest skeleton voxels
%%  Find the target voxel by threshold relaxation
% Generate all the possible connection and then do a selection later. 
% 1. Gather voxel list for connecting 
%   a. Delete the duplicated pair
%   b. Find pairs for all the ep1, ep2, threshold relaxation of ep2 -
%   Classifier for removing fake ep2 before? Too few examples? 
num_pair = int_link.num.ep;
th_rlx_pair = struct;
th_rlx_pair.pair_ind = zeros(num_pair, 2);
th_rlx_pair.sub_1 = int_link.ep_sub;
th_rlx_pair.ind_1 = int_link.ep_ind;
th_rlx_pair.pair_ind(:,1) = th_rlx_pair.ind_1;
th_rlx_pair.link_label_1 = full(vessel_graph_refine_2.link.map_ind_2_label(th_rlx_pair.ind_1));
assert(all(th_rlx_pair.link_label_1~=0), 'Exist endpoint that is not a link voxel');
th_rlx_pair.link_cc_ind_1 = vessel_graph_refine_2.link.cc_ind(th_rlx_pair.link_label_1);

th_rlx_pair.ind_2 = zeros(num_pair, 1);
th_rlx_pair.sub_2 = zeros(3, num_pair);
th_rlx_pair.connecting_th = zeros(num_pair, 1);
th_rlx_pair.connecting_Z_score = zeros(num_pair, 1);
th_rlx_pair.pair_dist = zeros(num_pair, 1);

% Takes about 5 seconds
tic
for iter_link_idx = 1 : num_pair
    ep_1_sub = th_rlx_pair.sub_1 (iter_link_idx, :);
    target_voxel_str = fun_graph_get_nearest_skl_ind(vessel_image, vessel_mask_rc, vessel_skel_refine, ep_1_sub, int_link.ep_sub);
    % Record 
    if target_voxel_str.foundQ
        th_rlx_pair.pair_ind(iter_link_idx, 2) = target_voxel_str.global_ind;
        th_rlx_pair.sub_2(:, iter_link_idx) = target_voxel_str.global_sub';
        th_rlx_pair.connecting_th(iter_link_idx) = target_voxel_str.connected_th;
        th_rlx_pair.connecting_Z_score(iter_link_idx) = round((target_voxel_str.connected_th - target_voxel_str.bg_mean) / target_voxel_str.bg_std);
        th_rlx_pair.pair_dist(iter_link_idx) = target_voxel_str.dist_ep2target;
    end
end
toc
tmp_found_Q_list = th_rlx_pair.pair_ind(:, 2) > 0;
th_rlx_pair.sub_2 = th_rlx_pair.sub_2';
th_rlx_pair = fun_structure_field_indexing(th_rlx_pair, tmp_found_Q_list );
%%  Find the nearest endpoint pairs
% Since link with 2 endpoints are often short, finding the nearest endpoint
% should be done in three steps.
nearest_ep_pair = struct;
% ep1 to ep1
dist_ep1_to_ep1 = squareform(pdist(int_link.ep1.ep_sub));
search_range_ep1_to_ep1 = 60;% um
[mutual_nn_ep1_idx, mutual_nn_ep2_idx, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_ep1_to_ep1, true);
tmp_Q = mutual_nn_dist <= search_range_ep1_to_ep1;
mutual_nn_ep1_idx = mutual_nn_ep1_idx(tmp_Q);
mutual_nn_ep2_idx = mutual_nn_ep2_idx(tmp_Q);
mutual_nn_dist = mutual_nn_dist(tmp_Q);
nearest_ep_pair.ep1_ep1.ind_1 = int_link.ep1.ep_ind(mutual_nn_ep1_idx);
nearest_ep_pair.ep1_ep1.ind_2 = int_link.ep1.ep_ind(mutual_nn_ep2_idx);
nearest_ep_pair.ep1_ep1.pair_ind = cat(2, int_link.ep1.ep_ind(mutual_nn_ep1_idx), ...
    int_link.ep1.ep_ind(mutual_nn_ep2_idx));
nearest_ep_pair.ep1_ep1.pair_link_label = cat(2, nearest_ep_pair.ep1_ep1.ind_1, ...
    nearest_ep_pair.ep1_ep1.ind_2 );
nearest_ep_pair.ep1_ep1.pair_ep_dist = mutual_nn_dist;
% Delete duplicate:
tmp_pair_ind = sort(nearest_ep_pair.ep1_ep1.pair_ind, 2, 'ascend');
[tmp_pair_unique, tmp_unique_idx] = unique(tmp_pair_ind, 'rows', 'stable');
nearest_ep_pair.ep1_ep1.ind_1 = nearest_ep_pair.ep1_ep1.ind_1(tmp_unique_idx);
nearest_ep_pair.ep1_ep1.ind_2 = nearest_ep_pair.ep1_ep1.ind_2(tmp_unique_idx);
nearest_ep_pair.ep1_ep1.pair_ind = nearest_ep_pair.ep1_ep1.pair_ind(tmp_unique_idx, :);
nearest_ep_pair.ep1_ep1.pair_link_label = nearest_ep_pair.ep1_ep1.pair_link_label(tmp_unique_idx, :);
nearest_ep_pair.ep1_ep1.pair_ep_dist = nearest_ep_pair.ep1_ep1.pair_ep_dist(tmp_unique_idx);
% ep1 to ep2
dist_ep1_to_ep2 = pdist2(int_link.ep1.ep_sub, int_link.ep2.ep_sub);
search_range_ep1_to_ep2 = 25;% um
[mutual_nn_ep1_idx, mutual_nn_ep2_idx, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_ep1_to_ep2);
tmp_Q = mutual_nn_dist <= search_range_ep1_to_ep2;
mutual_nn_ep1_idx = mutual_nn_ep1_idx(tmp_Q);
mutual_nn_ep2_idx = mutual_nn_ep2_idx(tmp_Q);
mutual_nn_dist = mutual_nn_dist(tmp_Q);

nearest_ep_pair.ep1_ep1.ind_1 = int_link.ep1.ep_ind(mutual_nn_ep1_idx);
nearest_ep_pair.ep1_ep1.ind_2 = int_link.ep2.ep_ind(mutual_nn_ep2_idx);
nearest_ep_pair.ep1_ep2.pair_ind = cat(2, nearest_ep_pair.ep1_ep1.ind_1, nearest_ep_pair.ep1_ep1.ind_2);
nearest_ep_pair.ep1_ep2.pair_link_label = cat(2, int_link.ep1.link_label(mutual_nn_ep1_idx), int_link.ep2.ep_link_label(mutual_nn_ep2_idx));
nearest_ep_pair.ep1_ep2.pair_ep_dist = mutual_nn_dist;

% ep2 to ep2
dist_ep2_to_ep2 = squareform(pdist(int_link.ep2.ep_sub));
search_range_ep2_to_ep2 = 20;% um
tmp_max = max(dist_ep2_to_ep2(:));
for iter_ep2 = 1 : 2 : (int_link.ep2.num_cc) * 2
    dist_ep2_to_ep2(iter_ep2, iter_ep2+1) = tmp_max + 1;
    dist_ep2_to_ep2(iter_ep2+1, iter_ep2) = tmp_max + 1;
end
[mutual_nn_ep2_1_idx, mutual_nn_ep2_2_idx, mutual_nn_dist_ep2ep2] = fun_find_col_row_co_minimum(dist_ep2_to_ep2, true);
tmp_Q = mutual_nn_dist_ep2ep2 <= search_range_ep2_to_ep2;
mutual_nn_ep2_1_idx = mutual_nn_ep2_1_idx(tmp_Q);
mutual_nn_ep2_2_idx = mutual_nn_ep2_2_idx(tmp_Q);
mutual_nn_dist_ep2ep2 = mutual_nn_dist_ep2ep2(tmp_Q);

[~, tmp_unique_idx ] = unique(sort(cat(2, mutual_nn_ep2_1_idx, mutual_nn_ep2_2_idx), 2), 'row', 'stable');
nearest_ep_pair.ep2_ep2.pair_ind  = cat(2, int_link.ep2.ep_ind(mutual_nn_ep2_1_idx(tmp_unique_idx)), ...
    int_link.ep2.ep_ind(mutual_nn_ep2_2_idx(tmp_unique_idx)));
nearest_ep_pair.ep2_ep2.pair_link_label = cat(2, int_link.ep2.ep_link_label(mutual_nn_ep2_1_idx(tmp_unique_idx)),...
    int_link.ep2.ep_link_label(mutual_nn_ep2_2_idx(tmp_unique_idx)));
nearest_ep_pair.ep2_ep2.pair_ep_dist = mutual_nn_dist_ep2ep2(tmp_unique_idx);
%%  Combine the pair found by restricted shortest path with the endpoint pairs
pair_ind = cat(1, th_rlx_pair.pair_ind, nearest_ep_pair.ep1_ep1.pair_ind, ...
    nearest_ep_pair.ep1_ep2.pair_ind, nearest_ep_pair.ep2_ep2.pair_ind);
pair_unique = sort(pair_ind, 2, 'ascend');
pair_unique = unique(pair_unique, 'rows', 'stable');
sub_1_list = fun_ind2sub(image_size, pair_unique(:,1))';
sub_2_list = fun_ind2sub(image_size, pair_unique(:,2))';
% Generate linkers for all these pairs
clear linker_str
num_pair = size(pair_unique,1);
linker_found = false(num_pair, 1);
tic
for iter_pair = num_pair:-1:1
    linker_str(iter_pair) = fun_graph_connect_gap_p2p(vessel_image, vessel_mask_rc, ...
        sub_1_list(:, iter_pair), sub_2_list(:, iter_pair), 26);
%     linker_cell{iter_pair} = tmp_linker;
    linker_found(iter_pair) = linker_str(iter_pair).foundQ;
end
toc
linker_found_idx = find(linker_found);
%%  Compute features for the linker
linker_to_test = linker_str(linker_found);
pair_to_test = pair_unique(linker_found, :);
num_linkers = numel(linker_to_test);
training_data = struct;
% Intensity
training_data.data.int_mean = [linker_to_test.int_mean]';
training_data.data.int_std = [linker_to_test.int_std]';
training_data.data.int_om_mean = [linker_to_test.int_o_mask_mean]';
training_data.data.int_om_std = [linker_to_test.int_o_mask_std]';
training_data.data.int_med = [linker_to_test.int_med]';
% Signal to noise
training_data.data.connected_th_SNR = (([linker_to_test.connected_th] - [linker_to_test.bg_mean])./ [linker_to_test.bg_std])';
training_data.data.linker_SNR = (([linker_to_test.int_mean] - [linker_to_test.bg_mean])./ [linker_to_test.bg_std])';
training_data.data.linker_o_m_SNR = (([linker_to_test.int_o_mask_mean] - [linker_to_test.bg_mean])./ [linker_to_test.bg_std])';
% Out of the mask ratio
training_data.data.linker_o_m_ratio = [linker_to_test.link_ratio_o_mask]';
% Length
training_data.data.num_voxel = [linker_to_test.num_voxel]';
training_data.data.length = cellfun(@(x) fun_graph_ind_to_length(x), {linker_to_test.link_sub_w_ep})';
% Position and end to end distance
training_data.data.sub_1 = sub_1_list(:, valid_linker_idx)';
training_data.data.sub_2 = sub_2_list(:, valid_linker_idx)';
training_data.data.p2p_dist = sqrt(sum((training_data.data.sub_1 - training_data.data.sub_2).^2, 2));
% Ratio between the end-to-end distance and the linker length
training_data.data.p2p_dist_2_length = training_data.data.p2p_dist ./ training_data.data.length;
training_data.data.dir_vec_1_to_2 = training_data.data.sub_2 - training_data.data.sub_1;
training_data.data.dir_vec_1_to_2 = training_data.data.dir_vec_1_to_2 ./ vecnorm(training_data.data.dir_vec_1_to_2, 2, 2);
% Point 1 is an endpointQ:
training_data.data.p1_is_ep = full(vessel_graph_refine_2.endpoint.map_ind_2_label(pair_to_test(:,1)));
training_data.data.p2_is_ep = full(vessel_graph_refine_2.endpoint.map_ind_2_label(pair_to_test(:,2)));
training_data.data.p1_is_node = full(vessel_graph_refine_2.node.map_ind_2_label(pair_to_test(:,1)));
training_data.data.p2_is_node = full(vessel_graph_refine_2.node.map_ind_2_label(pair_to_test(:,2)));
training_data.data.p1_is_link = full(vessel_graph_refine_2.link.map_ind_2_label(pair_to_test(:,1)));
training_data.data.p2_is_link = full(vessel_graph_refine_2.link.map_ind_2_label(pair_to_test(:,2)));
training_data.data.p1_is_isopoint_Q = (~training_data.data.p1_is_ep) & (~training_data.data.p1_is_node) & (~training_data.data.p1_is_link);
training_data.data.p2_is_isopoint_Q = (~training_data.data.p2_is_ep) & (~training_data.data.p2_is_node) & (~training_data.data.p2_is_link);
% Endpoint to endpoint linker
training_data.data.ep2ep_Q = training_data.data.p1_is_ep & training_data.data.p2_is_ep;
training_data.data.ep2lk_Q = (training_data.data.p1_is_ep & training_data.data.p2_is_link) |...
    (training_data.data.p2_is_ep & training_data.data.p1_is_link);
% 
training_data.link_cc_ind_w_ep = {linker_to_test.link_ind_w_ep};
training_data.link_cc_sub_w_ep = {linker_to_test.link_sub_w_ep};
% Orientation of the linker
training_data.data.dir_vec_linker = cellfun(@fun_graph_compute_link_orientation, training_data.link_cc_sub_w_ep, 'UniformOutput', false);
training_data.data.dir_vec_linker = cat(2, training_data.data.dir_vec_linker{:})';
% Orientation of the linker outside the mask 
training_data.data.dir_vec_linker_o_m = cellfun(@fun_graph_compute_link_orientation, {linker_to_test.link_sub_o_m}, 'UniformOutput', false);
training_data.data.dir_vec_linker_o_m = cat(2, training_data.data.dir_vec_linker_o_m{:})';
% inner products
training_data.data.ip_vec_lkr_vec_lkr_om = abs(sum(training_data.data.dir_vec_linker_o_m .* training_data.data.dir_vec_linker, 2)); 
training_data.data.ip_vec_dir12_vec_lkr = abs(sum(training_data.data.dir_vec_1_to_2  .* training_data.data.dir_vec_linker, 2));
training_data.data.ip_vec_dir12_vec_lkr_om = abs(sum(training_data.data.dir_vec_1_to_2  .* training_data.data.dir_vec_linker_o_m, 2));
% Use the value in training_data.data.ip_vec_dir12_vec_lkr to replace
% the NAN in training_data.data.ip_vec_dir12_vec_lkr_om
tmpQ = isnan(training_data.data.ip_vec_dir12_vec_lkr_om);
training_data.data.ip_vec_dir12_vec_lkr_om(tmpQ) = training_data.data.ip_vec_dir12_vec_lkr(tmpQ);

training_data.data = struct2table(training_data.data);
training_data.feature_name = training_data.data.Properties.VariableNames;
training_data.info.dataset_name = dataset_name;
training_data.info.stack = stack;
training_data.info.grid_c_version = grid_c_version;
training_data.info.combined_grid_xyz_label = combined_grid_xyz_label;
training_data.info.note = [];
training_data.info.annotator = 'Xiang Ji';
%% Annotation
new_gap_linker_cc = {linker_str(linker_found).link_ind_w_ep};
App_check_fake_gap_linker(vessel_image, vessel_mask, vessel_skel_refine, new_gap_linker_cc)
annotation_ep1_to_nn = App_check_fake_gap_linker_result;
annotation_ep1_to_nn.to_remove_idx = find(annotation_ep1_to_nn.to_removeQ);
annotation_ep1_to_nn.keep_Q = (~ annotation_ep1_to_nn.not_sureQ )& (~ annotation_ep1_to_nn.to_removeQ);
%%  Train classifier 
classifier_str = struct;
classifier_str.training_data = training_data;

classifier_str.used_feature_idx = [1, 3, 6, 7, 8, 9, 10, 14, 15, 29, 30, 31];
classifier_str.used_feature_name = classifier_str.training_data.feature_name(classifier_str.used_feature_idx);
classifier_str.label.removeQ = annotation_ep1_to_nn.to_removeQ;
classifier_str.label.not_sureQ = annotation_ep1_to_nn.not_sureQ;
classifier_str.label.remove_or_not_sureQ = classifier_str.label.removeQ | classifier_str.label.not_sureQ;
classifier_str.num_example = size(classifier_str.training_data.data, 1);

tmp_training_set_Q = rand(classifier_str.num_example, 1) < 0.8;
classifier_str.training_set.data = table2array(classifier_str.training_data.data(tmp_training_set_Q, classifier_str.used_feature_idx));
classifier_str.training_set.label = classifier_str.label.remove_or_not_sureQ(tmp_training_set_Q);
classifier_str.validation_set.data = table2array(classifier_str.training_data.data(~tmp_training_set_Q, classifier_str.used_feature_idx));
classifier_str.validation_set.label = classifier_str.label.remove_or_not_sureQ(~tmp_training_set_Q);
classifier_str.classifier = fitcensemble(classifier_str.training_set.data, classifier_str.training_set.label, 'Method', 'AdaBoostM1');

test_label = classifier_str.validation_set.label;
predicted_label = classifier_str.classifier.predict(classifier_str.validation_set.data);
classification_stat.num_true_positive = nnz(test_label == true & predicted_label == true);
classification_stat.num_positive = nnz(predicted_label == true);
classification_stat.num_true_negative = nnz(test_label == false & predicted_label == false);
classification_stat.num_negative = nnz(predicted_label == false);
classification_stat.num_false_positive = nnz(test_label == false & predicted_label == true);
classification_stat.num_false_negative = nnz(test_label == true & predicted_label == false);

classification_stat.true_positive = classification_stat.num_true_positive / classification_stat.num_positive;
classification_stat.true_negative = classification_stat.num_true_negative / classification_stat.num_negative;
classification_stat.false_positive = classification_stat.num_false_positive / classification_stat.num_positive;
classification_stat.false_negative = classification_stat.num_false_negative / classification_stat.num_negative;

classifier_str.validation_set.stat = classification_stat;
% Save the classifier
% classifier_fp = fullfile(DataManager.SCRIPT_PATH, 'Graph', 'Pruning_by_learning', 'classifier_linker_to_remove_20190121.mat');
% save(classifier_fp, '-struct', 'classifier_str');
% Note: 
% 1. I am a bit worry about this part since the number of positive training
% set is limited and the accuracy of the classifier can vary quite a lot.
% 2. I need to implement functions for handeling annotation data. These
% data should be able to be integrated to improve the classifier. 
%% Add linker to the vessel graph 
linker_valid_idx = linker_found_idx(annotation_ep1_to_nn.keep_Q);
sub_1_list_valid = sub_1_list(:, linker_valid_idx)';
sub_2_list_valid = sub_2_list(:, linker_valid_idx)';
[sub_1_list_valid, sub_2_list_valid, kept_pair_Q ]= fun_graph_point_pairs_is_significantly_unique(sub_1_list_valid, sub_2_list_valid, 5);
linker_valid_idx = linker_valid_idx(kept_pair_Q);
linker_valid = linker_str(linker_valid_idx);
% Add the linker by modifying the skeleton: 
% Should implement a faster algorithm to operate on the graph directly 

vessel_skel_gap_fixed = vessel_skel_refine;
for iter_linker = 1 : numel(linker_valid)
    vessel_skel_gap_fixed(linker_valid(iter_linker).link_ind_w_ep) = true;
end
profile on 
vessel_graph_gap_fixed = fun_skeleton_to_graph(vessel_skel_gap_fixed);
profile off
profile viewer
vessel_graph_gap_fixed = fun_graph_pruning_internal_short_hairs(vessel_graph_gap_fixed, 2);
int_link_3 = fun_graph_get_free_link(vessel_graph_gap_fixed, internal_offset);
% Remove all the internal links with 2 endpoints - if a correct linker
% could be found, they would have been at least reduced to a link with 1
% endpoint. 
vessel_graph_gap_fixed = fun_graph_pruning_by_link_label(vessel_graph_gap_fixed, int_link_3.ep2.link_label);
vessel_skel_gap_fixed(cat(1, int_link_3.ep2.link_cc_ind{:})) = false;

annotation_ep1_round_2 = App_check_fake_gap_linker_result;
% Need to implement functions for automatically get the training dataset
% and also merge two training dataset 

fun_vis_link_surrounding_by_cc_ind(vessel_image, vessel_mask, vessel_skel_refine, int_link_3.ep1.link_cc_ind{2})


tmp_store.graph = vessel_graph_refine_2;
tmp_store.linker = linker_valid;
save('./data_for_adding_linker_to_graph.mat', '-struct', 'tmp_store');

%% Deal with the remaining internal endpoints
% 2. Linker selection - train classifier
% 3. Add the correct linker to the graph OR add to the skeleton and convert
% the corrected skeleton to graph 
% 4. Examine the low signal to noise level link, remove them - another
% classifier. 
% 5. Record the internal links with endpoint 
% 6. Use the remaining linker to modifiy the mask
% 7. Update the mask - save. Add skeleton information to the mask and
% record the corrected face direction of the mask cube
% 8. Network analysis
%%
%% Searching for possible endpoint to endpoint pair first
% The problem of this apporach is that the resulting classifier need to
% deal with multiple different cases and for ep2 to ep2 connection, the
% number of training examples is small. 
%% Collect endpoint pairs within a threshold values
% Distance between ep1 to ep1:
dist_ep1_to_ep1 = squareform(pdist(int_link.ep1.ep_sub));
search_range_ep1_to_ep1 = 50;% um
search_ep1_to_ep1_Q = (dist_ep1_to_ep1 <= search_range_ep1_to_ep1) & (dist_ep1_to_ep1 ~= 0); % Remove self-distance(0);
search_idx_ep1_to_ep1 = find(any(search_ep1_to_ep1_Q, 2));
search_pair_ep1_to_ep1 = cell(numel(search_idx_ep1_to_ep1),1);
search_pair_ep1_to_ep1_dist = cell(numel(search_idx_ep1_to_ep1),1);
for iter_ep = 1 : numel(search_idx_ep1_to_ep1)
    search_pair_ep1_to_ep1{iter_ep} = find(search_ep1_to_ep1_Q(search_idx_ep1_to_ep1(iter_ep), :))';
    search_pair_ep1_to_ep1_dist{iter_ep} = dist_ep1_to_ep1(search_idx_ep1_to_ep1(iter_ep),search_pair_ep1_to_ep1{iter_ep})';
end
search_pair_ep1_to_ep1_num = cellfun(@numel, search_pair_ep1_to_ep1);

ep1_ep1_idx_pair = cat(2, repelem(search_idx_ep1_to_ep1, search_pair_ep1_to_ep1_num), cat(1, search_pair_ep1_to_ep1{:}));
tmp_dist = cat(1, search_pair_ep1_to_ep1_dist{:});

ep1_ep1_ind_pair = int_link.ep1.ep_ind(ep1_ep1_idx_pair);
[nearest_ep_pair.ep1_ep1.pair_ind , tmp_unique_idx ] = unique(sort(ep1_ep1_ind_pair, 2), 'row', 'stable');
nearest_ep_pair.ep1_ep1.pair_link_label = int_link.ep1.link_label(ep1_ep1_idx_pair(tmp_unique_idx,:));
nearest_ep_pair.ep1_ep1.pair_ep_dist = tmp_dist(tmp_unique_idx);
%% ep1 to ep2 mutual nearest neighbor
dist_ep1_to_ep2 = pdist2(int_link.ep1.ep_sub, int_link.ep2.ep_sub);
search_range_ep1_to_ep2 = 25;% um
[mutual_nn_ep1_idx, mutual_nn_ep2_idx, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_ep1_to_ep2);
tmp_Q = mutual_nn_dist <= search_range_ep1_to_ep2;
mutual_nn_ep1_idx = mutual_nn_ep1_idx(tmp_Q);
mutual_nn_ep2_idx = mutual_nn_ep2_idx(tmp_Q);
mutual_nn_dist = mutual_nn_dist(tmp_Q);

nearest_ep_pair.ep1_ep2.pair_ind = cat(2, int_link.ep1.ep_ind(mutual_nn_ep1_idx), int_link.ep2.ep_ind(mutual_nn_ep2_idx));
nearest_ep_pair.ep1_ep2.pair_link_label = cat(2, int_link.ep1.link_label(mutual_nn_ep1_idx), int_link.ep2.ep_link_label(mutual_nn_ep2_idx));
nearest_ep_pair.ep1_ep2.pair_ep_dist = mutual_nn_dist;
%% ep2 to ep2 mutual nearest neighbor - exclude endpoint in the same
% connected component
dist_ep2_to_ep2 = squareform(pdist(int_link.ep2.ep_sub));
search_range_ep2_to_ep2 = 20;% um
tmp_max = max(dist_ep2_to_ep2(:));
for iter_ep2 = 1 : 2 : (int_link.ep2.num_cc) * 2
    dist_ep2_to_ep2(iter_ep2, iter_ep2+1) = tmp_max + 1;
    dist_ep2_to_ep2(iter_ep2+1, iter_ep2) = tmp_max + 1;
end
[mutual_nn_ep2_1_idx, mutual_nn_ep2_2_idx, mutual_nn_dist_ep2ep2] = fun_find_col_row_co_minimum(dist_ep2_to_ep2, true);
tmp_Q = mutual_nn_dist_ep2ep2 <= search_range_ep2_to_ep2;
mutual_nn_ep2_1_idx = mutual_nn_ep2_1_idx(tmp_Q);
mutual_nn_ep2_2_idx = mutual_nn_ep2_2_idx(tmp_Q);
mutual_nn_dist_ep2ep2 = mutual_nn_dist_ep2ep2(tmp_Q);

[~, tmp_unique_idx ] = unique(sort(cat(2, mutual_nn_ep2_1_idx, mutual_nn_ep2_2_idx), 2), 'row', 'stable');
nearest_ep_pair.ep2_ep2.pair_ind  = cat(2, int_link.ep2.ep_ind(mutual_nn_ep2_1_idx(tmp_unique_idx)), ...
    int_link.ep2.ep_ind(mutual_nn_ep2_2_idx(tmp_unique_idx)));
nearest_ep_pair.ep2_ep2.pair_link_label = cat(2, int_link.ep2.ep_link_label(mutual_nn_ep2_1_idx(tmp_unique_idx)),...
    int_link.ep2.ep_link_label(mutual_nn_ep2_2_idx(tmp_unique_idx)));
nearest_ep_pair.ep2_ep2.pair_ep_dist = mutual_nn_dist_ep2ep2(tmp_unique_idx);
%% Generating linkers for post-selection by threshold relaxation guided shortest path
% vessel_mask_rc = vessel_mask | vessel_skl_rc;
nearest_ep_pair.ep1_ep1.num_pair = size(nearest_ep_pair.ep1_ep1.pair_ind,1);
nearest_ep_pair.ep1_ep1.linker_str = cell(nearest_ep_pair.ep1_ep1.num_pair, 1);
nearest_ep_pair.ep1_ep1.linker_foundQ = false(nearest_ep_pair.ep1_ep1.num_pair, 1);
image_size = size(vessel_image);
tic
for iter_pair = 1 : nearest_ep_pair.ep1_ep1.num_pair % Takes 228sec for 500 links
    fprintf('Finished %d/%d\n', iter_pair, nearest_ep_pair.ep1_ep1.num_pair);
    test_ind_1 = nearest_ep_pair.ep1_ep1.pair_ind(iter_pair,1);
    test_ind_2 = nearest_ep_pair.ep1_ep1.pair_ind(iter_pair,2);
    test_sub_1 = fun_ind2sub(image_size, test_ind_1);
    test_sub_2 = fun_ind2sub(image_size, test_ind_2);
    connecting_link_str = fun_graph_connect_gap_p2p(vessel_image, vessel_mask_rc, test_sub_1, test_sub_2, 26);
    if connecting_link_str.foundQ
        nearest_ep_pair.ep1_ep1.linker_foundQ(iter_pair) = true;
        nearest_ep_pair.ep1_ep1.linker_str{iter_pair} = connecting_link_str;
    end
end
toc

nearest_ep_pair.ep1_ep2.num_pair = size(nearest_ep_pair.ep1_ep2.pair_ind,1);
nearest_ep_pair.ep1_ep2.linker_str = cell(nearest_ep_pair.ep1_ep2.num_pair, 1);
nearest_ep_pair.ep1_ep2.linker_foundQ = false(nearest_ep_pair.ep1_ep2.num_pair, 1);
tic
for iter_pair = 1 : nearest_ep_pair.ep1_ep2.num_pair % Takes 0.38sec for 14 links
    fprintf('Finished %d/%d\n', iter_pair, nearest_ep_pair.ep1_ep2.num_pair);
    test_ind_1 = nearest_ep_pair.ep1_ep2.pair_ind(iter_pair,1);
    test_ind_2 = nearest_ep_pair.ep1_ep2.pair_ind(iter_pair,2);
    test_sub_1 = fun_ind2sub(image_size, test_ind_1);
    test_sub_2 = fun_ind2sub(image_size, test_ind_2);
    connecting_link_str = fun_graph_connect_gap_p2p(vessel_image, vessel_mask_rc, test_sub_1, test_sub_2, 26);
    if connecting_link_str.foundQ
        nearest_ep_pair.ep1_ep2.linker_foundQ(iter_pair) = true;
        nearest_ep_pair.ep1_ep2.linker_str{iter_pair} = connecting_link_str;
    end
end
toc

nearest_ep_pair.ep2_ep2.num_pair = size(nearest_ep_pair.ep2_ep2.pair_ind,1);
nearest_ep_pair.ep2_ep2.linker_str = cell(nearest_ep_pair.ep2_ep2.num_pair, 1);
nearest_ep_pair.ep2_ep2.linker_foundQ = false(nearest_ep_pair.ep2_ep2.num_pair, 1);
tic
for iter_pair = 1 : nearest_ep_pair.ep2_ep2.num_pair % Takes 0.09sec for 1 links
    fprintf('Finished %d/%d\n', iter_pair, nearest_ep_pair.ep2_ep2.num_pair);
    test_ind_1 = nearest_ep_pair.ep2_ep2.pair_ind(iter_pair,1);
    test_ind_2 = nearest_ep_pair.ep2_ep2.pair_ind(iter_pair,2);
    test_sub_1 = fun_ind2sub(image_size, test_ind_1);
    test_sub_2 = fun_ind2sub(image_size, test_ind_2);
    connecting_link_str = fun_graph_connect_gap_p2p(vessel_image, vessel_mask_rc, test_sub_1, test_sub_2, 26);
    if connecting_link_str.foundQ
        nearest_ep_pair.ep2_ep2.linker_foundQ(iter_pair) = true;
        nearest_ep_pair.ep2_ep2.linker_str{iter_pair} = connecting_link_str;
    end
end
toc
% vis_local_bbox_mmll = fun_vis_link_surrounding_by_cc_ind(vessel_image, vessel_mask, vessel_skl_rc, matched_pair.ep1_ep1.linker_str{end-1}.link_ind_w_ep,  50);
%% Annotate ep1 to ep1 
nearest_ep_pair.ep1_ep1.valid_pair_idx = find(~cellfun(@isempty, nearest_ep_pair.ep1_ep1.linker_str));
nearest_ep_pair.ep1_ep1.num_valid_pair = numel(nearest_ep_pair.ep1_ep1.valid_pair_idx);
nearest_ep_pair.ep1_ep1.linker_cc = cell(nearest_ep_pair.ep1_ep1.num_valid_pair, 1);
nearest_ep_pair.ep1_ep1.valid_linker_num_voxel = zeros(nearest_ep_pair.ep1_ep1.num_valid_pair, 1);
for iter_pair = 1 : nearest_ep_pair.ep1_ep1.num_valid_pair 
    tmp_str = nearest_ep_pair.ep1_ep1.linker_str{nearest_ep_pair.ep1_ep1.valid_pair_idx(iter_pair)};
    nearest_ep_pair.ep1_ep1.linker_cc{iter_pair} = tmp_str.link_ind_w_ep;
    nearest_ep_pair.ep1_ep1.valid_linker_num_voxel(iter_pair) = numel(tmp_str.link_ind_w_ep);
end
App_check_fake_gap_linker(vessel_image, vessel_mask, vessel_skl_refined, nearest_ep_pair.ep1_ep1.linker_cc)