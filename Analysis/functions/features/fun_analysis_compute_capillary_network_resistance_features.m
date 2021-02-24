function network_resist_str = fun_analysis_compute_capillary_network_resistance_features(vessel_graph)

%% Parameters
internal_offset = 240;
fitting_start_dist = 100;
min_num_datt_point_for_fitting = 100;
num_bin_for_fitting = 100;
%% Initialize the output structure
persistent str_template; 
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'fitting',...
        'fitting_measures', 'internal_offset', 'fitting_start_dist', ...
        'network_bbox_size', 'valid_node_sub', 'r_pair'});
end
network_resist_str = str_template;
%% Get capillary network
vessel_med_r = cellfun(@(x) median(full(vessel_graph.radius(x))), vessel_graph.link.cc_ind);
% Remove all the links of radius greater than a thresold to get the
% capillary network
while any(vessel_med_r == 0)
    cap_cc_Q = (vessel_med_r > 0 & vessel_med_r <= capillary_max_r);
    cap_ind = cat(1, vessel_graph.link.cc_ind{cap_cc_Q}, vessel_graph.node.cc_ind{:});
    vessel_graph = fun_skeleton_to_graph(cap_ind, test_im_size);
    vessel_graph = fun_graph_add_radius(vessel_graph, vessel_graph.radius, min_r);
    vessel_med_r = cellfun(@(x) median(full(vessel_graph.radius(x))), vessel_graph.link.cc_ind);
end
%% Construct the capillary graph and compute the network resistance
vessel_graph.link.features = fun_simulation_blood_flow_compute_network_features(vessel_graph, 1e-6);
% Get the largest connected component in the graph
[vessel_graph.node.graph_cc_label, vessel_graph.link.graph_cc_label ] = ...
    fun_graph_get_connected_component_label_for_nodes(vessel_graph);
[node_cc_label_list_idx, node_cc_label_unique] = fun_bin_data_to_idx_list(vessel_graph.node.graph_cc_label);
[~, largest_cc_ind] = max(cellfun(@numel, node_cc_label_list_idx));
% Construct the conductance matrix and compute the node to node network
% resistance.
valid_node_label = find(vessel_graph.node.graph_cc_label == node_cc_label_unique(largest_cc_ind));
[conductance_cap, ~] = fun_simulation_blood_flow_get_conductance_matrix(vessel_graph, valid_node_label);
r_pair_cap = fun_analysis_compute_node_to_node_network_resistance(conductance_cap, true);
% Compute the center of mass of the node
valid_node_sub = cellfun(@(x) mean(fun_ind2sub(vessel_graph.num.mask_size, x), 1), vessel_graph.node.cc_ind(valid_node_label), 'UniformOutput', false);
valid_node_sub = cat(1, valid_node_sub{:});
% Select the node that are far from the boundary

is_far_from_boundary_Q = fun_voxel_sub_in_bbox_mmxx_Q(valid_node_sub, [internal_offset, internal_offset, internal_offset, ...
    vessel_graph.num.mask_size - internal_offset + 1]);

if ~any(is_far_from_boundary_Q) 
    return; 
end

valid_node_sub = valid_node_sub(is_far_from_boundary_Q, :);
valid_node_pdist = squareform(pdist(valid_node_sub));
r_pair_cap = r_pair_cap(is_far_from_boundary_Q, is_far_from_boundary_Q);
% Bin node wise network resistance by node to node distance
bin_dist_data  = valid_node_pdist;
bin_resistance_data = r_pair_cap;
is_not_diagonal = ~eye(size(bin_dist_data), 'logical');
bin_resistance_data = bin_resistance_data(is_not_diagonal);
bin_dist_data = bin_dist_data(is_not_diagonal);
dist_edge = linspace(min(bin_dist_data (:)), max(bin_dist_data (:)), num_bin_for_fitting);
[dist_ind_cell] = fun_bin_data_to_idx_list_by_edges(bin_dist_data , dist_edge);
num_point_in_cell = cellfun(@numel, dist_ind_cell);
valid_bin_Q = num_point_in_cell > min_num_datt_point_for_fitting;
dist_ind_cell = dist_ind_cell(valid_bin_Q);
resistance_binned_cell = fun_bin_data_to_cells_by_ind(bin_resistance_data, dist_ind_cell);
avg_binned_resistance = cellfun(@mean, resistance_binned_cell);
avg_dist = movmean(dist_edge, 2, 'Endpoints', 'discard');
avg_dist = avg_dist(valid_bin_Q);
% Linear fit
used_for_fitting_Q = avg_dist > fitting_start_dist;
avg_resist_fit = avg_binned_resistance(used_for_fitting_Q);
avg_dist_fit = avg_dist(used_for_fitting_Q)';
[fit_obj, gof_obj] = fit(avg_dist_fit, avg_resist_fit, 'poly1');
%% Save result
network_resist_str.fitting = fit_obj;
network_resist_str.fitting_measures = gof_obj;
network_resist_str.internal_offset = internal_offset;
network_resist_str.fitting_start_dist = fitting_start_dist;
network_resist_str.valid_node_sub = valid_node_sub;
network_resist_str.r_pair = r_pair_cap;
network_resist_str.network_bbox_size = vessel_graph.num.mask_size;
end