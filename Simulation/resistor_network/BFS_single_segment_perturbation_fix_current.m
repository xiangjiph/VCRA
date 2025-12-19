clc;clear;close all;
DataManager = FileManager();
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
grid_version = '240_cube';
vessel_graph_0 = load('WholeBrain_ML_2018_08_15_240_cube_combined_5_o_1_578_annotated_20210305.mat');
grid_c_version = '240_cube_combined_5_o_1';
skl_version = '240_cube_rec';

vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'BFS_single_segment_perturbation');
data_name = sprintf('%s_%s_c5o1_%d_BFS_SSP', dataset_name, stack, ...
    vessel_graph_0.info.combined_grid_xyz_label);
itk_folder_name = sprintf('%s_annotated', grid_c_version);
min_dist_to_boundary = 50;
%% Simulation parameters
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
%% Preprocessing the graph
% Get the largest graph connected component 
% This needs to be done first to keep the edge label consistent
[vessel_graph, valid_node_label] = fun_simulation_blood_flow_preprocessing_graph(vessel_graph_0, ...
    flow_parameters);
%% Add static structural features
vessel_graph_analyzed.link.features = fun_analysis_get_link_features_by_label(...
    vessel_graph, [], {'geometry'});
% Artery and vein branch order: 
vessel_graph.link.features.artery_branch_order = fun_analysis_get_capillary_order_to_labeled_vessels(vessel_graph, ...
    vessel_graph.link.features.type, 1, 2);
% Reverse branch order - start from capilalries
vessel_graph.link.features.artery_branch_order_r = fun_analysis_get_capillary_order_to_labeled_vessels(vessel_graph, ...
    vessel_graph.link.features.type, 2, 1);

vessel_graph.link.features.vein_branch_order = fun_analysis_get_capillary_order_to_labeled_vessels(vessel_graph, ...
    vessel_graph.link.features.type, 1, 3);
vessel_graph.link.features.vein_branch_order_r = fun_analysis_get_capillary_order_to_labeled_vessels(vessel_graph, ...
    vessel_graph.link.features.type, 3, 1);

int_bbox_mmxx = [ones(1,3) * (min_dist_to_boundary + 1), ...
    vessel_graph.num.mask_size - min_dist_to_boundary + 1];
is_internal_link_Q = fun_voxel_sub_in_bbox_mmxx_Q(vessel_graph_analyzed.link.features.mid_sub,...
    int_bbox_mmxx);
%%
tic
resistor_network = ResistorNetwork(vessel_graph, valid_node_label, flow_parameters);
rnssp = RNSSP(resistor_network);
rnssp.compute_pseudoinv_UWL();
toc
%% Self-modulation factor, faa vs vessel segment length
CuAu_list = rnssp.get_fss_factors(rnssp.InvL, rnssp.RN.link_conductance);
% Stat
vsl_len_um = vessel_graph.link.features.length(resistor_network.link_label) * 1e6;
length_bin = 0 : 10 : 450;
ratio_bin = 1e-8 : 0.05 : 1.05;
% Overlay the stat curve
[link_bin_idx, link_bin_val] = fun_bin_data_to_idx_list_by_edges(vsl_len_um, length_bin);
is_valid_bin_Q = cellfun(@numel, link_bin_idx) > 10;
link_bin_val = link_bin_val(is_valid_bin_Q);
fss_in_bin_stat = fun_analysis_get_basic_statistics_in_bins(CuAu_list, link_bin_idx);
fss_in_bin_med = cat(1, fss_in_bin_stat(is_valid_bin_Q).prctile_val);
fss_in_bin_med = fss_in_bin_med(:, 7:9);

fig_hdl = figure;
fig_hdl.Position(3:4) = [600, 450];
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, vsl_len_um, CuAu_list, ...
    'XBinEdges', length_bin, 'YBinEdges', ratio_bin, 'DisplayStyle', 'tile')
hold(ax_hdl, 'on');
errorbar(ax_hdl, link_bin_val, fss_in_bin_med(:, 2), ...
    fss_in_bin_med(:, 2) - fss_in_bin_med(:, 1), ...
    fss_in_bin_med(:, 3) - fss_in_bin_med(:, 2), 'LineWidth', 1.5);
ax_hdl.ColorScale = 'log';
ax_hdl.XLabel.String = 'Segment Length (\mum)';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$f^{\alpha}_{\alpha}$';
cbar_hdl = colorbar(ax_hdl);
cbar_hdl.Label.String = 'Number of segments';
cbar_hdl.Limits(1) = 1;
ax_hdl.FontSize = 12;
% Endpoint veins: 
non_ep_vein_label = setdiff(find(vessel_graph.link.features.type == 3), ...
    resistor_network.vein_ep_link_label);
non_ep_vein_idx = resistor_network.map_int_link_label_to_idx(non_ep_vein_label);
non_ep_vein_idx_valid_Q = ~isnan(non_ep_vein_idx);
non_ep_vein_label = non_ep_vein_label(non_ep_vein_idx_valid_Q);
non_ep_vein_idx = non_ep_vein_idx(non_ep_vein_idx_valid_Q);

vein_len = vessel_graph.link.features.length(non_ep_vein_label);
vein_med_radius = vessel_graph.link.features.radius(non_ep_vein_label);
vein_r = CuAu_list(non_ep_vein_idx);
hold(ax_hdl, 'on');
sca_vein_hdl = scatter(ax_hdl, vein_len * 1e6, vein_r, 'b.');
% Show arteries specifically
non_ep_artery_label = setdiff(find(vessel_graph.link.features.type == 2), ...
    resistor_network.artery_ep_link_label);
non_ep_artery_idx = resistor_network.map_int_link_label_to_idx(non_ep_artery_label);
non_ep_artery_idx_valid_Q = ~isnan(non_ep_artery_idx);
non_ep_artery_label = non_ep_artery_label(non_ep_artery_idx_valid_Q);
non_ep_artery_idx = non_ep_artery_idx(non_ep_artery_idx_valid_Q);

art_len = vessel_graph.link.features.length(non_ep_artery_label);
art_med_radius = vessel_graph.link.features.radius(non_ep_artery_label);
art_r = CuAu_list(non_ep_artery_idx);
hold(ax_hdl, 'on');
sca_art_hdl = scatter(ax_hdl, art_len * 1e6, art_r, 'r.');

legend(ax_hdl, [sca_vein_hdl, sca_art_hdl], {'Vein', 'Artery'});

fig_fp = fullfile(vis_folder, sprintf('%s_fix_flux_f_aa.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%  Histogram of self-modulation factor for arteries, veins and capillaries
% ratio_bin = 1e-8 : 0.05 : 1.05;
hist_plot_setting = {'Normalization', 'pdf', 'FaceAlpha', 0.5};
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, art_r, ratio_bin, hist_plot_setting{:}, 'FaceColor', '#D95319');
hold(ax_hdl, 'on');
histogram(ax_hdl, vein_r, ratio_bin, hist_plot_setting{:}, 'FaceColor', '#0072BD');
is_capillary_Q = true(size(CuAu_list));
is_capillary_Q([non_ep_artery_idx; non_ep_vein_idx]) = false;
hold(ax_hdl, 'on');
histogram(ax_hdl, CuAu_list(is_capillary_Q), ratio_bin, hist_plot_setting{:}, 'FaceColor', '#77AC30');
legend(ax_hdl, {'Artery', 'Vein', 'Capillary'}, 'Location', 'best');
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.XLabel.String = '$f^{\alpha}_{\alpha}$';
ax_hdl.YLabel.String = 'PDF';

fig_fp = fullfile(vis_folder, sprintf('%s_fix_flux_f_aa_hist_vs_vsl_type.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
% Arteries are less "rigid" than veins. 
%% Effect on neighboring vessels: 
% Construct link linear response matrix f^mn_ab
Fmn_ab = rnssp.construct_Fmn_ab(rnssp.InvL);

% edge pair-wise distance 
edge_mid_dist = squareform(single(pdist(vessel_graph_analyzed.link.features.mid_sub))); 

% edge geodesic distance 
edge_g_dist = zeros(size(edge_mid_dist), 'uint8');
max_g_dist = 8;
% should it be directional? but we don't know the direction - without
% simulation 
t_tic = tic;
for iter_el = 1 : numel(rnssp.RN.link_label)
    tmp_edge = rnssp.RN.link_label(iter_el);
    tmp_nb_el = fun_graph_get_neighbor_link_by_order(vessel_graph, tmp_edge, max_g_dist);
    tmp_num_nb_e = cellfun(@numel, tmp_nb_el);
    tmp_nb_e_o = repelem((1 : max_g_dist)', tmp_num_nb_e, 1);
    tmp_nb_el = cat(1, tmp_nb_el{:});
    % for consistancy, each column is for each edge
    edge_g_dist(tmp_nb_el, tmp_edge) = tmp_nb_e_o;
end
fprintf('Finish computing edge geodesic distances. Elapsed time is %.2f seconds\n', toc(t_tic));
%% Fba / (1 - faa) vs euclidean distance between edges middle point 
faa = diag(Fmn_ab);
fba_2_1mfaa = Fmn_ab ./ (1 - faa.');
fba_2_1mfaa(1 : size(fba_2_1mfaa, 1) + 1 : end) = 0;

tmp_gt1_ind = find(fba_2_1mfaa > 1.01);
tmp_gt1_sub = fun_ind2sub(size(fba_2_1mfaa), tmp_gt1_ind);
tmp_gt1_edge_a = unique(tmp_gt1_sub(:, 2), 'stable');
tmp_gt1_edge_a_faa = faa(tmp_gt1_edge_a);
for i = 1 : numel(tmp_gt1_edge_a)
    tmp_el = tmp_gt1_edge_a(i);
    if abs(faa(tmp_el) - 1) < 1e-6 
        fba_2_1mfaa(tmp_el, :) = 0;
        fba_2_1mfaa(:, tmp_el) = 0;
    end
end
fba_2_1mfaa = abs(single(fba_2_1mfaa));

max_dist = 1000;
% dist_bin = 0 : 25 : 500;
dist_bin = 10 .^ linspace(0, log10(max_dist), 20);
fnmab_bin = 10 .^ linspace(-6, 0, 20);
%% fba / (1 - faa) vs euclidean distance for different vessel types
loglog_fit_min_x = 50;
[fig_hdl] = rnssp.plot_fba_2_1mfaa_vs_euclidean_dist_by_vsl_type(edge_mid_dist, fba_2_1mfaa, ...
                dist_bin, fnmab_bin, vessel_graph.link.features.type, ...
                is_internal_link_Q, loglog_fit_min_x);

fig_fp = fullfile(vis_folder, sprintf('%s_fix_flux_fba21mfaa_vs_edist_hist2.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% fba / (1 - faa) vs geodesic distance for different vessel types
plt_x_bin = 0.5 : (max_g_dist + 0.5);
[fig_hdl] = rnssp.plot_fba_2_1mfaa_vs_geodesic_dist_by_vsl_type(single(edge_g_dist), fba_2_1mfaa, ...
                plt_x_bin, fnmab_bin, vessel_graph.link.features.type, ...
                is_internal_link_Q);

fig_fp = fullfile(vis_folder, sprintf('%s_fix_flux_fba21mfaa_vs_gdist_hist2.png', data_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
