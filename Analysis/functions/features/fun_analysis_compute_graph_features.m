function [input_graph, varargout]= fun_analysis_compute_graph_features(input_graph, opt)
% A wrap up function for computing link and node features
%% Parameters
if nargin < 2
    computeQ = struct;
    computeQ.basic_link_feature = true;
    computeQ.basic_node_feature = true;
    computeQ.basic_reconstruction = true;
    computeQ.node_path_to_nearest_neighbor = false;
    computeQ.link_shortest_path = false;
    computeQ.dist_tissue_2_vessel = true;
    computeQ.dimension = false;
    vis_dimension_fit = false;
    
else
    computeQ = opt.computeQ;
    vis_dimension_fit = opt.vis_dim_fit;
    
end
computeQ.graph = computeQ.link_shortest_path || computeQ.node_path_to_nearest_neighbor;
computeQ.recon_mask = computeQ.dimension || computeQ.dist_tissue_2_vessel || ...
    computeQ.basic_reconstruction;
image_size = input_graph.num.mask_size;
recon_cc = cat(1, input_graph.link.cc_ind, input_graph.node.cc_ind);
recon_r = fun_graph_get_cc_radius(recon_cc, input_graph.radius, 'all');
recon_r = max(2, recon_r);
recon_ind = cat(1, input_graph.link.pos_ind, input_graph.node.pos_ind);
recon_label = cat(1, input_graph.link.label, -input_graph.node.label);
input_graph.stat = struct;

if computeQ.recon_mask
    vessel_recon_mask = fun_skeleton_reconstruction(recon_ind, recon_r, image_size);
    if nargout == 2
        varargout{1} = vessel_recon_mask;
    end
end

if computeQ.basic_reconstruction
    input_graph.num.Mask_volume = nnz(vessel_recon_mask);  
    input_graph.num.Block_volume = prod(image_size);
    input_graph.num.Mask_volume_ratio = input_graph.num.Mask_volume/ input_graph.num.Block_volume;
end%% Node properties
if computeQ.basic_node_feature
    input_graph.node.features = fun_analysis_get_node_cc_features(input_graph);
end
%% Link properties
if computeQ.basic_link_feature
    input_graph.link.features = fun_analysis_get_link_features(input_graph, {'geometry', 'dt'});
    % Deal with 0 radius link
    tmp_0_r_link_label = find(input_graph.link.features.dt_median == 0);
    if ~isempty(tmp_0_r_link_label)
        for iter_link_label = tmp_0_r_link_label
            tmp_node_label = input_graph.link.connected_node_label(iter_link_label, :);
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
                warning('The radius of both the node are 0. Link label %d\n', iter_link_label);
            end
            input_graph.link.features.dt_median(iter_link_label) = tmp_r;    
        end
    end
    input_graph.num.total_link_length = sum(input_graph.link.features.length);
    input_graph.num.total_surface_area = sum(input_graph.link.features.surface_area);
end
%% Geodesic and physical path length to the neighbor vessels
if computeQ.graph
    input_graph.graph = fun_analysis_get_connectivity_graph(input_graph);
end
if computeQ.node_path_to_nearest_neighbor
    path_to_neighbor_str = fun_analysis_get_path_to_neighbors(input_graph.graph.graph_w, input_graph.node.features.neighbor_node_label);
    % Copy the field to input graph
    input_graph.node.features.path_to_neighbor_length = path_to_neighbor_str.path_to_neighbor_length;
    input_graph.node.features.path_to_neighbor_geodesic_length = path_to_neighbor_str.path_to_neighbor_geodesic_length;
    input_graph.node.features.path_to_nearest_neighbor = path_to_neighbor_str.path_to_nearest_neighbor;
    input_graph.node.features.path_to_nearest_neighbor_length = path_to_neighbor_str.path_to_nearest_neighbor_length;
    input_graph.node.features.path_to_nearest_neighbor_geodesic_length = path_to_neighbor_str.path_to_nearest_neighbor_geodesic_length;
end
%% Compute the loop properties of links
if computeQ.link_shortest_path
    shortest_loop_str = fun_analysis_get_loops_in_graph(input_graph.graph, 'euclidean');
    input_graph.link.features.shortest_loop_length = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.shortest_loop_length(shortest_loop_str.link_label) = shortest_loop_str.loop_length;
    
    input_graph.link.features.shortest_loop_geodesic_length = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.shortest_loop_geodesic_length(shortest_loop_str.link_label) = shortest_loop_str.loop_geodesic_length;
    
    input_graph.link.features.shortest_loop_node_label = cell(input_graph.link.num_cc, 1);
    input_graph.link.features.shortest_loop_node_label(shortest_loop_str.link_label) = shortest_loop_str.loop_node_label;
    
    input_graph.link.features.shortest_loop_link_label = cell(input_graph.link.num_cc, 1);
    input_graph.link.features.shortest_loop_link_label(shortest_loop_str.link_label) = shortest_loop_str.loop_link_label;
end
%% Distance transform properties
if computeQ.dist_tissue_2_vessel
    vessel_recon_mask_label = fun_skeleton_reconstruction_label(recon_ind, recon_r, recon_label, image_size);
    recon_prop = fun_analysis_reonstruction_space_properties(vessel_recon_mask, vessel_recon_mask_label, false);
    input_graph.stat.dist_tissue_2_vessel = recon_prop.global_stat;
    % Record global statistics
    input_graph.link.features.nearest_tissue_volume = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.nearest_tissue_volume(recon_prop.closest_link_label) = recon_prop.link_nearest_tissue_volume;
    input_graph.link.features.nearest_tissue_dt_mean = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.nearest_tissue_dt_mean(recon_prop.closest_link_label) = recon_prop.link_nearest_tissue_dt_mean;
    input_graph.link.features.nearest_tissue_dt_median = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.nearest_tissue_dt_median(recon_prop.closest_link_label) = recon_prop.link_nearest_tissue_dt_median;
    input_graph.link.features.nearest_tissue_dt_max = nan(input_graph.link.num_cc, 1);
    input_graph.link.features.nearest_tissue_dt_max(recon_prop.closest_link_label) = recon_prop.link_nearest_tissue_dt_max;

    input_graph.node.features.nearest_tissue_volume = nan(input_graph.node.num_cc, 1);
    input_graph.node.features.nearest_tissue_volume(recon_prop.closest_node_label) = recon_prop.node_nearest_tissue_volume;
    input_graph.node.features.nearest_tissue_dt_mean = nan(input_graph.node.num_cc, 1);
    input_graph.node.features.nearest_tissue_dt_mean(recon_prop.closest_node_label) = recon_prop.node_nearest_tissue_dt_mean;
    input_graph.node.features.nearest_tissue_dt_median = nan(input_graph.node.num_cc, 1);
    input_graph.node.features.nearest_tissue_dt_median(recon_prop.closest_node_label) = recon_prop.node_nearest_tissue_dt_median;
    input_graph.node.features.nearest_tissue_dt_max = nan(input_graph.node.num_cc, 1);
    input_graph.node.features.nearest_tissue_dt_max(recon_prop.closest_node_label) = recon_prop.node_nearest_tissue_dt_max;
end
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
end