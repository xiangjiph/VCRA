clc;clear;close all;
DataManager = FileManager;
dataset_name = 'DKLab';
stack = 'Blinder_Tsai';
folder_fp = DataManager.fp_graph_folder(dataset_name, stack);
itk_folder_name = 'itk';
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
    'combined_stat');
data_name = sprintf('%s_%s', dataset_name, stack);
%% Load graphs
vg_cell = {};

updated_fp = fullfile(DataManager.fp_analysis_data_folder('WholeBrain', 'ML_2018_08_15'), ...
    'Test_graph_w_avl_1.mat');
vg_cell{end+1} = DataManager.load_data(updated_fp);

bt_graph_list = {'au', 'co', 'db', 'av'};
for i = 1 : numel(bt_graph_list)
    tmp_fp = fullfile(folder_fp, sprintf('%s_%s_%s_graph.mat', dataset_name, stack, ...
        bt_graph_list{i}));
    tmp_g = DataManager.load_data(tmp_fp);
    vg_cell{end+1} = tmp_g;
end
%% Parameters
flow_parameters = struct;
% Pressure in Pascal. From Blinder 2013 Nat Neuroscience
flow_parameters.vein_endpoint_pressure = 0 * 133.322; 
flow_parameters.pressure_diff_artery_to_vein = 50 * 133.322; 

flow_parameters.artery_endpoint_pressure = flow_parameters.vein_endpoint_pressure + ...
    flow_parameters.pressure_diff_artery_to_vein; 

flow_parameters.capillary_label = 1;
flow_parameters.artery_label = 2;
flow_parameters.vein_label = 3;

flow_parameters.voxel_size_m = 1e-6;

% 
opt_graph_feature = struct;
computeQ = struct;
computeQ.basic_link_feature = true;
computeQ.basic_node_feature = true;
computeQ.basic_reconstruction = false;
computeQ.node_path_to_nearest_neighbor = true;
computeQ.link_shortest_path = true;
computeQ.dist_tissue_2_vessel = false;
computeQ.dimension = false;
computeQ.max_z_proj = false;
computeQ.capillary_branching_order = true;
computeQ.capillary_artery_order = NaN; % Placeholder for furture development 
computeQ.capillary_vein_order = NaN; % Placeholder for furture development

computeQ.noncapillary_DT = false; % Pick the vessels are are not capillaries, analyze the distance transform properties. 
computeQ.sgl_seg_rm_ptb = false;
computeQ.dist_to_brain_surface = false; % Compute the distance between the vessels and the surface of the brain. 

computeQ.ep2ep_anisotropy = false;
computeQ.volume_weighted_ep2ep_anisotropy = false;
computeQ.capillary_ep2ep_anisotropy = false;
computeQ.capillary_volume_weighted_ep2ep_anisotropy = false;    

opt_graph_feature.computeQ = computeQ;
opt_graph_feature.vis_dim_fit = false;
opt_graph_feature.dim_fit_cutoff_length = 25;
opt_graph_feature.recon_max_error_rate = 0.1;
opt_graph_feature.merge_neighbor_nodes_Q = false;
opt_graph_feature.max_merge_length = 5;
opt_graph_feature.always_merge_max_length = 1;
opt_graph_feature.DT_scale_factor = 0.5;
opt_graph_feature.nonmax_win_size_um = 40;

opt_graph_feature.capillary_max_radius = 3.5; % For saperating capillary vs large vessels for local statistics. 

opt_graph_feature.brain_mask_fp = [];
opt_graph_feature.brain_mask_downsample_rate = 16;

min_dist_to_boundary = 50;
cap_max_r_um = 3.5;
%% Compute static structural feaures
num_graphs = numel(vg_cell);
vg_s_cell = cell(size(vg_cell));
for i = 1 : num_graphs
    t_tic = tic;
    tmp_vg = vg_cell{i};
    [tmp_vg_s, recon_mask] = fun_analysis_compute_graph_features_by_labels(tmp_vg,...
        [], [], opt_graph_feature);

    tmp_vg_s.link.features.type = nan(tmp_vg_s.link.num_cc, 1);
    for i_l = 1 : tmp_vg_s.link.num_cc
        tmp_ind = tmp_vg_s.link.cc_ind{i_l};
        tmp_label = unique(full(tmp_vg_s.label(tmp_ind)));
        if ~isscalar(tmp_label)
            fprintf('The link has mixed voxel label\n');
            % Check if more than 1 label is greater than 1
            tmp_label = tmp_label(tmp_label > 1);
            assert(isscalar(tmp_label), 'The link contians voxel labeled as artery and vein');
        end
        tmp_vg_s.link.features.type(i_l) = tmp_label;
    end
    % Artery and vein branch order:
    tmp_cap_Q = tmp_vg_s.link.features.dt_median <= cap_max_r_um;    
    tmp_vg_s.link.features.artery_branch_order = fun_analysis_get_capillary_order_to_labeled_vessels_masked(tmp_vg_s, ...
        tmp_vg_s.link.features.type, 1, 2, tmp_cap_Q);
    tmp_vg_s.link.features.vein_branch_order = fun_analysis_get_capillary_order_to_labeled_vessels_masked(tmp_vg_s, ...
        tmp_vg_s.link.features.type, 1, 3, tmp_cap_Q);

    int_bbox_mmxx = [ones(1,3) * (min_dist_to_boundary + 1), ...
        tmp_vg_s.num.mask_size - min_dist_to_boundary + 1];
    is_internal_link_Q = fun_voxel_sub_in_bbox_mmxx_Q(tmp_vg_s.link.features.mid_sub,...
        int_bbox_mmxx);
    tmp_vg_s.link.features.is_internal_Q = is_internal_link_Q;
    vg_s_cell{i} = tmp_vg_s;
    fprintf("Finish processing graph %d. Elapsed time is %.2f\n", i, toc(t_tic));
end
%% Artery to vein: number of branch orders
max_branch_order = 10;
[bov_of_boa1_c, boa_of_bov1_c, loop_num_edge_c, boa_c, bov_c] = deal(cell(num_graphs, 1));

for i = 1 : num_graphs
    tmp_vg_s = vg_s_cell{i};
    is_internal_cap_Q = (tmp_vg_s.link.features.type == 1 &...
        tmp_vg_s.link.features.is_internal_Q & ...
        tmp_vg_s.link.features.dt_median <= cap_max_r_um);
    int_cap_boa = tmp_vg_s.link.features.artery_branch_order(is_internal_cap_Q);
    int_cap_bov = tmp_vg_s.link.features.vein_branch_order(is_internal_cap_Q);

    int_cap_bov_of_boa1 = int_cap_bov(int_cap_boa == 1);
    int_cap_bov_of_boa1 = int_cap_bov_of_boa1(isfinite(int_cap_bov_of_boa1 ));
    int_cap_boa_of_bov1 = int_cap_boa(int_cap_bov == 1);
    int_cap_boa_of_bov1 = int_cap_boa_of_bov1(isfinite(int_cap_boa_of_bov1 ));
    
    bov_of_boa1_c{i} = int_cap_bov_of_boa1;
    boa_of_bov1_c{i} = int_cap_boa_of_bov1;
    loop_num_edge_c{i} = tmp_vg_s.link.features.shortest_loop_geodesic_length(is_internal_cap_Q);
    boa_c{i} = int_cap_boa;
    bov_c{i} = int_cap_bov;
end

cb_boa1 = cat(1, bov_of_boa1_c{:});
cb_boa1_sz = cellfun(@numel, bov_of_boa1_c);
cb_boa1_l = repelem((1:num_graphs)', cb_boa1_sz , 1);
cb_boa1_stat = fun_analysis_get_basic_statistics(cb_boa1);

cb_bov1 = cat(1, boa_of_bov1_c{:});
cb_bov1_sz = cellfun(@numel, boa_of_bov1_c);
cb_bov1_l = repelem((1:num_graphs)', cb_bov1_sz, 1);
cb_bov1_stat = fun_analysis_get_basic_statistics(cb_bov1);

fig_hdl = figure;
fig_hdl.Position(4) = fig_hdl.Position(4) * 2;
ax1 = subplot(2,1,1);
boxplot(cb_boa1, cb_boa1_l);
ax1.YLabel.String = 'BOV of BOA1';
ax1.YLim(1) = 0;
grid(ax1, 'on');
ax1.XLabel.String = 'Graph index';
yyaxis(ax1, 'right');
plot(ax1, cb_boa1_sz);
ax1.YLabel.String = 'Number of edges';
ax1.Title.String = sprintf('sz: %d, mean %.1f, med %d', cb_boa1_stat.num_data, ...
    cb_boa1_stat.mean, cb_boa1_stat.median);

ax2 = subplot(2, 1, 2);
boxplot(cb_bov1, cb_bov1_l);
ax2.YLabel.String = 'BOA of BOV1';
grid(ax2, 'on');
ax2.YLim(1) = 0;
ax2.XLabel.String = 'Graph index';
yyaxis(ax2, 'right');
plot(ax2, cb_bov1_sz);
ax2.YLabel.String = 'Number of edges';
ax2.Title.String = sprintf('sz: %d, mean %.1f, med %d', cb_bov1_stat.num_data, ...
    cb_bov1_stat.mean, cb_bov1_stat.median);

fig_fp = fullfile(vis_folder, sprintf('%s_%s_artery_vein_branch_order_dist.png', data_name, 'combined'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Number of edges in shortest loop vs branch order
max_branch_order = 8;
cb_loop_num_edge = cat(1, loop_num_edge_c{:});
cb_boa = cat(1, boa_c{:});
cb_bov = cat(1, bov_c{:});

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax1 = subplot(1,2,1);

tmp_selected_Q = (cb_boa >= 1 & cb_boa <= max_branch_order);
tmp_y = cb_loop_num_edge(tmp_selected_Q);
tmp_x = cb_boa(tmp_selected_Q);
[tmp_x_idx, tmp_x_v] = fun_bin_data_to_idx_list(tmp_x);
tmp_y_cell = fun_bin_data_to_cells_by_ind(tmp_y, tmp_x_idx);
tmp_y_stat = fun_analysis_get_basic_statistics_in_bins(tmp_y, tmp_x_idx);
tmp_x_num = cellfun(@numel, tmp_x_idx);

bp1 = boxplot(ax1, tmp_y, tmp_x);
ax1.XLabel.String = 'BOA';
ax1.YLabel.String = 'Shortest loop length';
ax1.Title.String = sprintf('# %d', nnz(tmp_selected_Q));
ax1.YLim = [0, 25];
hold(ax1, 'on');
plot(ax1, tmp_x_v, [tmp_y_stat.mean], 'k.');
grid(ax1);
yyaxis(ax1, 'right');
plt1 = plot(ax1, tmp_x_v, tmp_x_num);
ax1.YLabel.String = 'Number of edges';
ax1.YLim(1) = 0;

ax2 = subplot(1,2,2);
tmp_selected_Q = (cb_bov >= 1 & cb_bov <= max_branch_order);
tmp_y = cb_loop_num_edge(tmp_selected_Q);
tmp_x = cb_bov(tmp_selected_Q);
[tmp_x_idx, tmp_x_v] = fun_bin_data_to_idx_list(tmp_x);
tmp_y_stat = fun_analysis_get_basic_statistics_in_bins(tmp_y, tmp_x_idx);
tmp_x_num = cellfun(@numel, tmp_x_idx);

boxplot(ax2, tmp_y, tmp_x);
ax2.XLabel.String = 'BOV';
ax2.YLabel.String = 'Shortest loop length';
ax2.Title.String = sprintf('# %d', nnz(tmp_selected_Q));
ax2.YLim = [0, 25];
hold(ax2, 'on');
plot(ax2, tmp_x_v, [tmp_y_stat.mean], 'k.');
grid(ax2);
yyaxis(ax2, 'right');
plt2 = plot(ax2, tmp_x_v, tmp_x_num);
ax2.YLabel.String = 'Number of edges';
ax2.YLim(1) = 0;

fig_fp = fullfile(vis_folder, sprintf('%s_%s_shortest_loop_length_vs_BOs_max_%d.png',...
    data_name, 'combined', max_branch_order));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
