function [input_graph, varargout] = fun_analysis_compute_graph_features_by_labels(input_graph, link_label_list, node_label_list, opt)
% A wrap up function for computing link and node features

%% Parameters and initialization
persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
if nargin < 2
    link_label_list = [];
    node_label_list = [];
end
if nargin < 4
    opt = DataManager.load_task_parameters('default', 'internal_subgrid_analysis.xml');
else
    computeQ = opt.computeQ;
    vis_dimension_fit = opt.vis_dim_fit;
end

if isfield(opt, 'overwrite_computed_featureQ') && opt.overwrite_computed_featureQ
    % Potential bug may exist, since there might be some dependence in the
    % order of feature calculated. 
    % Enable partially update the graph and the subgrid features
    computeQ = update_computeQ(input_graph, computeQ);
end

if isfield(opt, 'recon_max_error_rate')
    recon_max_error_rate = opt.recon_max_error_rate;
else
    recon_max_error_rate = 0.1;
end
if isempty(link_label_list)
    link_label_list = (1 : input_graph.link.num_cc).';
end
if isempty(node_label_list)
    node_label_list = (1 : input_graph.node.num_cc).';
end
input_graph.node.computed_features_label = node_label_list;
input_graph.link.computed_features_label = link_label_list;

computeQ.graph = computeQ.link_shortest_path || computeQ.node_path_to_nearest_neighbor;
computeQ.recon_mask = computeQ.dimension || computeQ.dist_tissue_2_vessel || ...
    computeQ.basic_reconstruction;
image_size = input_graph.num.mask_size;
num_node = numel(node_label_list);
num_link = numel(link_label_list);
all_vessel_radius = [];

input_graph.stat = struct;
%% Reconstruction
if computeQ.recon_mask && ~computeQ.dist_tissue_2_vessel
    vessel_recon_mask = fun_graph_to_reconstruction_mask(input_graph, false, recon_max_error_rate);
else
    vessel_recon_mask_label = fun_graph_to_reconstruction_mask(input_graph, true, recon_max_error_rate, 'link');
    vessel_recon_mask = (vessel_recon_mask_label ~= 0);    
end
if nargout == 2 && exist('vessel_recon_mask', 'var')
    varargout{1} = vessel_recon_mask;
end
if computeQ.basic_reconstruction
    input_graph.num.Mask_volume = nnz(vessel_recon_mask);  
    input_graph.num.Block_volume = prod(image_size);
    input_graph.num.Mask_volume_ratio = input_graph.num.Mask_volume/ input_graph.num.Block_volume;
end
%% Debug load image
% data_info = input_graph.info;
% block_im = DataManager.load_blocks_files('image', data_info.dataset_name, data_info.stack, '240_cube', ...
%     data_info.combined_grid_mmxx_grid(1) : data_info.combined_grid_mmxx_grid(4), ...
%     data_info.combined_grid_mmxx_grid(2) : data_info.combined_grid_mmxx_grid(5), ...
%     data_info.combined_grid_mmxx_grid(3) : data_info.combined_grid_mmxx_grid(6));
% 
% block_seg_mask = DataManager.load_blocks_files('mask', data_info.dataset_name, data_info.stack, '240_cube', ...
%     data_info.combined_grid_mmxx_grid(1) : data_info.combined_grid_mmxx_grid(4), ...
%     data_info.combined_grid_mmxx_grid(2) : data_info.combined_grid_mmxx_grid(5), ...
%     data_info.combined_grid_mmxx_grid(3) : data_info.combined_grid_mmxx_grid(6));
% 
% vis_mask = uint8(block_seg_mask);
% vis_mask(vessel_recon_mask) = 2;
% vis_mask(input_graph.link.pos_ind) = 3;
% vis_mask(input_graph.node.pos_ind) = 4;
% DataManager.visualize_itksnap(block_im, vis_mask)
%% Fractal dimension
if computeQ.dimension
    if isfield(opt, 'dim_fit_cutoff_length')
        dimension_fitting_cutoff_length = opt.dim_fit_cutoff_length; %um
    else
        dimension_fitting_cutoff_length = 25;
    end
    dimension_properties = fun_analysis_get_fractal_dimension(vessel_recon_mask, dimension_fitting_cutoff_length, vis_dimension_fit);
    input_graph.stat.dimension = dimension_properties; 
end
clearvars vessel_recon_mask
%% Link properties
if computeQ.basic_link_feature
    input_graph.link.features = fun_analysis_get_link_features_by_label(input_graph, link_label_list, {'geometry', 'dt'});
    % Deal with 0 radius link
    tmp_0_r_link_label = find(input_graph.link.features.dt_median == 0);
    if ~isempty(tmp_0_r_link_label)
        for iter_link_idx = tmp_0_r_link_label
            tmp_link_label = link_label_list(iter_link_idx);
            tmp_node_label = input_graph.link.connected_node_label(tmp_link_label, :);
            if tmp_node_label(1) ~= 0
                tmp_r1 = full(mean(input_graph.radius(input_graph.node.cc_ind{tmp_node_label(1)})));
            else
                tmp_r1 = 0;
            end
            if tmp_node_label(2) ~= 0
                tmp_r2 = full(mean(input_graph.radius(input_graph.node.cc_ind{tmp_node_label(2)})));
            else
                tmp_r2 = 0;
            end
            if tmp_r1 && tmp_r2
                tmp_r = (tmp_r1 + tmp_r2)/2;
            elseif tmp_r1 || tmp_r2
                tmp_r = max(tmp_r1, tmp_r2);
            else
                warning('The radius of both the node are 0. Link label %d\n', iter_link_idx);
            end
            input_graph.link.features.dt_median(iter_link_idx) = tmp_r;    
        end
    end
    % These numbers only make sense if the features of all the links are
    % computed. 
    input_graph.num.total_link_length = sum(input_graph.link.features.length);
    input_graph.num.total_surface_area = sum(input_graph.link.features.surface_area);
end
% Correlation between orientation and radius
% tmp_valid_Q = all(isfinite(input_graph.link.features.ep1_to_ep2_direction_vec), 2) & ...
%     isfinite(input_graph.link.features.dt_median) & input_graph.link.features.dt_median <= 3.5;
% [tmp1, tmp2] = corr(input_graph.link.features.dt_median(tmp_valid_Q), ...
%     abs(input_graph.link.features.ep1_to_ep2_direction_vec(tmp_valid_Q, 3)))
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, input_graph.link.features.dt_median(tmp_valid_Q), ...
%     input_graph.link.features.ep1_to_ep2_direction_vec(tmp_valid_Q, 3), 'DisplayStyle', 'tile')
% cbar = colorbar(ax_hdl);
%% Capillary branching order from major vessels
if isfield(computeQ, 'capillary_branching_order') && computeQ.capillary_branching_order
    if isempty(all_vessel_radius)
        all_vessel_radius = fun_analysis_get_cc_median_radius(input_graph.link.cc_ind, input_graph.radius);
    end
    assert(numel(all_vessel_radius) == input_graph.link.num_cc, 'The number of element in the radius vector is different from the number of links in the graph');
    is_large_vessel_Q = all_vessel_radius > opt.capillary_max_radius;
    capillary_order = fun_analysis_get_capillary_order_to_labeled_vessels(input_graph, is_large_vessel_Q, false, true);
    input_graph.link.features.is_large_vessel_Q = is_large_vessel_Q(link_label_list);
    input_graph.link.features.capillary_branching_order = capillary_order(link_label_list);
end
%% Node properties
if computeQ.basic_node_feature
    input_graph.node.features = fun_analysis_get_node_cc_features_by_label(input_graph, node_label_list);
    branching_geometric_features = fun_analysis_get_node_branching_geometric_features_by_label(input_graph, node_label_list, link_label_list);
    input_graph.node.features = fun_copy_fields_from_str2_to_str1(input_graph.node.features, branching_geometric_features);
end
%% Geodesic and physical path length to the neighbor nodes
if computeQ.graph
    input_graph.graph = fun_analysis_get_connectivity_graph(input_graph);
end
if computeQ.node_path_to_nearest_neighbor
    path_to_neighbor_str = fun_analysis_get_path_to_neighbors_by_node_label(input_graph.graph.graph_w, node_label_list, input_graph.node.features.neighbor_node_label);
    input_graph.node.features = fun_copy_fields_from_str2_to_str1(input_graph.node.features, path_to_neighbor_str);
end
%% Compute the loop properties of links
if computeQ.link_shortest_path
    shortest_loop_str = fun_analysis_get_loops_in_graph_by_link_label(input_graph.graph, link_label_list, 'euclidean');
    
    map_link_label_2_list_idx = zeros(input_graph.link.num_cc, 1);
    map_link_label_2_list_idx(link_label_list) = 1 : num_link;    
    
    input_graph.link.features.shortest_loop_length = nan(num_link, 1);
    input_graph.link.features.shortest_loop_length(map_link_label_2_list_idx(shortest_loop_str.link_label)) = shortest_loop_str.loop_length;
    
    input_graph.link.features.shortest_loop_geodesic_length = nan(num_link, 1);
    input_graph.link.features.shortest_loop_geodesic_length(map_link_label_2_list_idx(shortest_loop_str.link_label)) = shortest_loop_str.loop_geodesic_length;
    
    input_graph.link.features.shortest_loop_node_label = cell(num_link, 1);
    input_graph.link.features.shortest_loop_node_label(map_link_label_2_list_idx(shortest_loop_str.link_label)) = shortest_loop_str.loop_node_label;
    
    input_graph.link.features.shortest_loop_link_label = cell(num_link, 1);
    input_graph.link.features.shortest_loop_link_label(map_link_label_2_list_idx(shortest_loop_str.link_label)) = shortest_loop_str.loop_link_label;
end
%% Distance transform properties
if isfield(computeQ, 'sgl_seg_rm_ptb') && computeQ.sgl_seg_rm_ptb
    dt_pur_str = fun_analysis_reconstruction_space_link_perturbation(...
        vessel_recon_mask_label, input_graph, input_graph.graph.graph_uw, ...
        link_label_list, opt.DT_scale_factor, opt.nonmax_win_size_um);
    
    input_graph.link.features = cat(2, input_graph.link.features, dt_pur_str);
    input_graph.link.features.nearest_tissue_radius = sqrt((input_graph.link.features.nearest_tissue_volume + ...
        input_graph.link.features.volume) ./ (pi * input_graph.link.features.length)) - ...
        input_graph.link.features.dt_median;
    
elseif computeQ.dist_tissue_2_vessel
    if isempty(link_label_list) || isempty(node_label_list)
        [input_graph.link.features.nearest_tissue_volume, ...
            input_graph.link.features.nearest_tissue_dt_mean, ...
            input_graph.link.features.nearest_tissue_dt_median,...
            input_graph.link.features.nearest_tissue_dt_max] = deal([]);
    else
        recon_prop = fun_analysis_reconstruction_space_properties(vessel_recon_mask_label, opt.DT_scale_factor, false);
        input_graph.stat.dist_tissue_2_vessel = recon_prop.global_stat;

        [~, idx_in_all, idx_in_list] = intersect(recon_prop.closest_link_label, link_label_list, 'stable');
        % Record global statistics
        input_graph.link.features.nearest_tissue_volume = nan(num_link, 1);
        input_graph.link.features.nearest_tissue_volume(idx_in_list) = recon_prop.link_nearest_tissue_volume(idx_in_all);
        input_graph.link.features.nearest_tissue_dt_mean = nan(num_link, 1);
        input_graph.link.features.nearest_tissue_dt_mean(idx_in_list) = recon_prop.link_nearest_tissue_dt_mean(idx_in_all);
        input_graph.link.features.nearest_tissue_dt_median = nan(num_link, 1);
        input_graph.link.features.nearest_tissue_dt_median(idx_in_list) = recon_prop.link_nearest_tissue_dt_median(idx_in_all);
        input_graph.link.features.nearest_tissue_dt_max = nan(num_link, 1);
        input_graph.link.features.nearest_tissue_dt_max(idx_in_list) = recon_prop.link_nearest_tissue_dt_max(idx_in_all);
        
        input_graph.link.features.nearest_tissue_radius = sqrt((input_graph.link.features.nearest_tissue_volume + ...
            input_graph.link.features.volume) ./ (pi * input_graph.link.features.length)) - ...
            input_graph.link.features.dt_median;
    end
end
%% Distance transform properties of large vessels
if isfield(computeQ, 'noncapillary_DT') && computeQ.noncapillary_DT
    large_vessel_recon_prop = fun_analysis_reconstruction_noncapillary_space_properties(input_graph, opt.capillary_max_radius, opt.DT_scale_factor, recon_max_error_rate);
    
    large_vessel_recon_prop = large_vessel_recon_prop(link_label_list, :);
    input_graph.link.features = cat(2, input_graph.link.features, large_vessel_recon_prop);
    
    input_graph.link.features.noncapillary_nearest_tissue_radius = ...
        sqrt((input_graph.link.features.noncapillary_nearest_tissue_dt_volume + ...
            input_graph.link.features.volume) ./ (pi * input_graph.link.features.length)) - ...
            input_graph.link.features.dt_median;
end
%% Distance between vessel and the surface of the brain
persistent brain_mask_dt brain_mask_ds_rate brain_mask_fp;
if computeQ.dist_to_brain_surface
    if isempty(brain_mask_dt) || ~strcmp(brain_mask_fp, opt.brain_mask_fp)
        if isfield(opt, 'brain_mask')
            brain_mask_dt = opt.brain_mask;    
        else
            error('Missing brain_mask field in opt structure');
            % But it seems that the file cannot be read by several
            % processes simultaneously.           
        end
        brain_mask_fp = opt.brain_mask_fp;        
        brain_mask_dt = bwdist(~brain_mask_dt);
        brain_mask_ds_rate = round(input_graph.info.dataset_size ./ size(brain_mask_dt));
        brain_mask_ds_rate = unique(brain_mask_ds_rate);
        assert(isscalar(brain_mask_ds_rate), 'Nonisotropic downsampling rate');
        brain_mask_dt = brain_mask_dt .* brain_mask_ds_rate;
    end
    % Crop the brain mask according to the input_graph bounding box
    region_bbox_mmxx_ds = round(input_graph.info.bbox_mmxx ./ brain_mask_ds_rate);
    region_bbox_mmxx_ds(1:3) = max(region_bbox_mmxx_ds(1:3), 1);
    region_bbox_mmxx_ds(4:6) = min(region_bbox_mmxx_ds(4:6), size(brain_mask_dt));
    region_bbox_mmll_ds = region_bbox_mmxx_ds;
    region_bbox_mmll_ds(4:6) = region_bbox_mmll_ds(4:6) - region_bbox_mmll_ds(1:3) + 1;
    
    region_mask_dt = crop_bbox3(brain_mask_dt, region_bbox_mmll_ds);
    % Compute the distance between links and the surface of the brain. 
    link_cc = input_graph.link.cc_ind(link_label_list);
    link_cc = fun_analysis_get_scaled_cc_ind(link_cc, image_size, 1 / brain_mask_ds_rate, region_bbox_mmll_ds(4:6));
    input_graph.link.features.dist_to_brain_surface_mean = fun_analysis_get_cc_distance_to_brain_surface(link_cc, region_mask_dt, @mean);
    
    node_cc = input_graph.node.cc_ind(node_label_list);
    node_cc = fun_analysis_get_scaled_cc_ind(node_cc, image_size, 1 / brain_mask_ds_rate, region_bbox_mmll_ds(4:6));    
    input_graph.node.features.dist_to_brain_surface_mean = fun_analysis_get_cc_distance_to_brain_surface(node_cc, region_mask_dt, @mean);
end
%% Save the information
input_graph.info.feature_computedQ = computeQ;
end

function computeQ = update_computeQ(input_graph, computeQ)
if isfield(input_graph, 'info') && isfield(input_graph.info, 'feature_computedQ')
   computed_Q = input_graph.info.feature_computedQ;
   computed_feature_name = fieldnames(computed_Q);
   to_compute_feature_name = fieldnames(computeQ);
   field_to_turn_off = intersect(computed_feature_name, to_compute_feature_name);
   for iter_field = 1 : numel(field_to_turn_off)
       computeQ.(field_to_turn_off{iter_field}) = false;
   end    
end
end