function [vessel_graph, valid_node_label] = fun_simulation_blood_flow_preprocessing_graph(vessel_graph, ...
    flow_parameters)


vessel_graph = fun_graph_sort_link_cc_ind_from_node_1_to_2(vessel_graph);
vessel_graph.endpoint.type = full(vessel_graph.label(vessel_graph.endpoint.pos_ind));
% Purning without removing open-ended large vessels
large_vessel_endpoint_ind = vessel_graph.endpoint.pos_ind(...
    full(vessel_graph.label(vessel_graph.endpoint.pos_ind)) == flow_parameters.artery_label | ...
    full(vessel_graph.label(vessel_graph.endpoint.pos_ind)) == flow_parameters.vein_label);
vessel_graph = fun_graph_pruning_internal_short_hairs_masked(vessel_graph, ...
    inf, 0, large_vessel_endpoint_ind);
% Some features have been computed 
vessel_graph.node.features = fun_analysis_get_node_cc_features(vessel_graph);

vessel_graph.link.features = fun_simulation_blood_flow_compute_network_features(vessel_graph,...
    flow_parameters.voxel_size_m);

% Find the connected components in the graph 
[vessel_graph.node.graph_cc_label, vessel_graph.link.graph_cc_label] = ...
    fun_graph_get_connected_component_label_for_nodes(vessel_graph);
% Find the connected components with both artery and veins
valid_cc_label = fun_simulation_blood_flow_get_valid_graph_cc_label(...
    vessel_graph.link.graph_cc_label, vessel_graph.link.features.type);
%%
assert(isscalar(valid_cc_label)); % Can be relaxed later
vessel_graph.link.features.no_endpoint_Q = all(vessel_graph.link.connected_node_label > 0, 2);
vessel_graph.node.features.is_valid_Q = false(vessel_graph.node.num_cc, 1);
for iter_node = 1 : vessel_graph.node.num_cc
    if vessel_graph.node.graph_cc_label(iter_node) == valid_cc_label
        tmp_link_label = vessel_graph.node.connected_link_label{iter_node};
        % Check if the link is in the connected component and does not have an
        % endpoint
        tmp_link_node = vessel_graph.link.connected_node_label(tmp_link_label, :);
        tmp_link_no_ep_Q = all(tmp_link_node > 0, 2);
        tmp_num_valid_link = nnz(tmp_link_no_ep_Q);
        if tmp_num_valid_link > 1
            vessel_graph.node.features.is_valid_Q(iter_node) = true;
        end
    end
end
valid_node_label = find(vessel_graph.node.features.is_valid_Q);
end