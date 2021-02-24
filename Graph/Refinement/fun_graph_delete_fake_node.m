function vessel_graph = fun_graph_delete_fake_node(vessel_graph)
% fun_graph_delete_fake_node convert node of degree 1 to endpoint and node
% of degree 2 to link voxels, delete node of degree 0
% 
% Implemented by Xiang Ji on 02/23/2019
node_label_to_delete = find(vessel_graph.node.num_link == 0);
node_label_to_link = find(vessel_graph.node.num_link == 2);
node_label_to_ep = find(vessel_graph.node.num_link == 1);
if ~isempty(node_label_to_link)
    for iter_node = 1 : numel(node_label_to_link)
        [vessel_graph, ~]= fun_graph_pruning_convert_node_to_link(vessel_graph, node_label_to_link(iter_node));
    end
end
if ~isempty(node_label_to_ep)
    for iter_node = 1 : numel(node_label_to_ep)
        [vessel_graph, ~]= fun_graph_pruning_convert_node_to_endpoint(vessel_graph, node_label_to_ep(iter_node));
    end
end

if ~isempty(node_label_to_delete)
    for iter_node = 1 : numel(node_label_to_delete)
        vessel_graph.node.cc_ind{node_label_to_delete(iter_node)} = [];
        vessel_graph.node.connected_link_label{node_label_to_delete(iter_node)} = [];
    end
end
vessel_graph = fun_graph_relabel(vessel_graph);
end