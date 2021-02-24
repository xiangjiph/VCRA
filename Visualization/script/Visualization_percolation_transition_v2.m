clc;clear;close all;
DataManager = FileManager;
vessel_graph = load('test_graph.mat');
vessel_graph_wo_ep = fun_graph_pruning_internal_short_hairs(vessel_graph, inf, 0);
[up_l_cc_label, up_n_cc_label] = fun_analysis_get_link_graph_cc_label(vessel_graph_wo_ep);
up_cc_node_ind = fun_bin_data_to_idx_list(up_n_cc_label);
up_cc_num_node = cellfun(@numel, up_cc_node_ind);

graph_ud_str = fun_analysis_get_connectivity_graph(vessel_graph_wo_ep);
graph_ud = graph_ud_str.graph_w;

graph_degree = degree(graph_ud);
phi_c = (mean(graph_degree)) / (mean(graph_degree.^2) - mean(graph_degree));
%%
percolation_p =  0.05 : 0.1 : 0.95;
percolation_str = fun_analysis_percolation_transition_by_bond_removal(graph_ud, percolation_p,  ...
    50, false);
%% Reconstruction of the unperturbed vascular network
vasc_mask = fun_graph_reconstruction(vessel_graph_wo_ep);
vasc_mask_center = crop_center_box(vasc_mask, ceil(vessel_graph.num.mask_size / 2), 120);
% fig_hdl = figure;
% v_hdl_up = volshow(vasc_mask_center);
% % volumeViewer(vasc_mask);
%% Remove half of the links and visualize the connected component in the reconstruction
vessel_skel_ind = cat(1, vessel_graph_wo_ep.link.pos_ind, vessel_graph_wo_ep.node.pos_ind);
max_est_num_cc = 256 * 10;
% label_color_map = jet(max_est_num_cc);
% Randomize the color map
% label_color_map = label_color_map(randperm(max_est_num_cc), :);
label_color_map = [1, 0, 0];

rm_link_ratio_list = 0.0 : 0.1 : 0.9;
for iter_idx = 1 : numel(rm_link_ratio_list)
    rm_link_ratio = rm_link_ratio_list(iter_idx);
    % rm_link_ratio = 0.50;
    num_rm_link = round(vessel_graph_wo_ep.link.num_cc * rm_link_ratio);
    rm_link_label = randsample(vessel_graph_wo_ep.link.num_cc, num_rm_link, false);
    rm_link_ind = cat(1, vessel_graph_wo_ep.link.cc_ind{rm_link_label});
    
    vessek_skel_ind_p = setdiff(vessel_skel_ind, rm_link_ind);
    
    vessel_graph_p = fun_skeleton_to_graph(vessek_skel_ind_p, vessel_graph_wo_ep.num.mask_size);
    
    %% get node graph cc label
    [p_l_cc_label, p_n_cc_label] = fun_analysis_get_link_graph_cc_label(vessel_graph_p);
    % Randomize the cc label for visualization
    p_l_cc_size = fun_bin_data_to_idx_list(p_l_cc_label);
    p_l_cc_size = cellfun(@numel, p_l_cc_size);
    [~, p_l_cc_size_ind] = sort(p_l_cc_size, 'descend');
    random_cc_label_map = nan(size(p_l_cc_size_ind));
    random_cc_label_map(p_l_cc_size_ind) = 1 : numel(p_l_cc_size);
%     p_num_cc = max(p_l_cc_label);
%     random_cc_label_map = randperm(p_num_cc);
    p_l_cc_label = random_cc_label_map(p_l_cc_label);
    p_n_cc_label = random_cc_label_map(p_n_cc_label);
    
    recon_ind = cat(1, vessel_graph_p.link.pos_ind, vessel_graph_p.node.pos_ind);
    recon_r = full(vessel_graph.radius(recon_ind));
    recon_label = cat(1, repelem(p_l_cc_label, vessel_graph_p.link.num_voxel_per_cc, 1), ...
        repelem(p_n_cc_label, vessel_graph_p.node.num_voxel_per_cc, 1));
    
    cc_node_ind = fun_bin_data_to_idx_list(p_n_cc_label);
    cc_num_node = cellfun(@numel, cc_node_ind);
    recon_mask = fun_skeleton_reconstruction_label_aprox(...
        recon_ind, recon_r, recon_label, vessel_graph.num.mask_size, 0.01);
    
    recon_mask_center = crop_center_box(recon_mask, ceil(vessel_graph.num.mask_size/2), 120);
    % Find the largest connected component in the cube
    [recon_mask_label_bin, recon_maks_bin_value]= fun_bin_data_to_idx_list(recon_mask_center(recon_mask_center ~= 0));
    num_pixel_in_cc = cellfun(@numel, recon_mask_label_bin);
    [num_pixel_in_cc, tmp_sort_idx] = sort(num_pixel_in_cc, 'descend');
    recon_mask_center = recon_mask_center == recon_maks_bin_value(tmp_sort_idx(1));
    
    [ax_hdl, fig_hdl] = fun_vis_labeled_array(recon_mask_center, label_color_map);
    ax_hdl.XAxis.Visible = 'off';
    ax_hdl.YAxis.Visible = 'off';
    ax_hdl.ZAxis.Visible = 'off';
    ax_hdl.Title.String = [];
    ax_hdl.Color = 'none';
    ax_hdl.Box = 'on';
    vis_fp = fullfile(DataManager.fp_visualization_folder(...
        vessel_graph.info.dataset_name, vessel_graph.info.stack), 'paper', ...
        sprintf('%s_%s_percolation_visualization_pc_%03d_largest_cc.png', ...
            vessel_graph.info.dataset_name, vessel_graph.info.stack, ...
            round((1 - rm_link_ratio) * 100)));
    fun_print_image_in_several_formats(fig_hdl, vis_fp);
    delete(fig_hdl);       
end