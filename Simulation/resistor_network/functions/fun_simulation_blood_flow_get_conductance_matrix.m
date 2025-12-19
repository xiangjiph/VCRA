function [conductance_mat, varargout]= fun_simulation_blood_flow_get_conductance_matrix(vessel_graph, node_label)

if nargin < 2
    node_label = 1 : vessel_graph.node.num_cc;
end
% Construct the conductance matrix without source term
num_valid_node = numel(node_label);
map_node_label_to_list_idx = nan(vessel_graph.node.num_cc, 1);
map_node_label_to_list_idx(node_label) = 1 : num_valid_node;

connected_link = cat(1, vessel_graph.node.connected_link_label{node_label});
connected_link = unique(connected_link, 'stable');

link_node_pair = vessel_graph.link.connected_node_label(connected_link, :);
% The link should not have endpoint, neither is it a self-loop
is_valid_pair_Q = ~any(link_node_pair == 0, 2) & (link_node_pair(:, 1) ~= link_node_pair(:, 2)); 

link_node_pair = link_node_pair(is_valid_pair_Q, :);
valid_link_label = connected_link(is_valid_pair_Q);
% Get the corresponding conductance
link_conductance = vessel_graph.link.features.conductance(connected_link);
link_conductance = link_conductance(is_valid_pair_Q);
% Convert the node pair to indices in the sparse matrix
node_pair_idx = map_node_label_to_list_idx(link_node_pair);
% Remove the links that are effectively link with endpoint in the graph. If
% a node is connected to three links, among which, two are the links with
% endpoint, then the node become an endpoint in the graph and thus has
% been removed previously. 
is_valid_pair_Q = all(~isnan(node_pair_idx), 2);
valid_link_label = valid_link_label(is_valid_pair_Q);
link_conductance = link_conductance(is_valid_pair_Q);
node_pair_idx = node_pair_idx(is_valid_pair_Q, :);
% Add non-diagonal term - conductance of multiple links between two nodes
% are added automatically in the sparse function. 
conductance_mat = sparse(node_pair_idx(:, 1), node_pair_idx(:, 2), ...
    link_conductance, num_valid_node, num_valid_node);

conductance_mat = conductance_mat + conductance_mat.';
% Add the diagonal term
diagonal_term = full(sum(conductance_mat, 2));
diagonal_term = sparse(1:num_valid_node, 1:num_valid_node, diagonal_term, ...
    num_valid_node, num_valid_node);
conductance_mat = diagonal_term - conductance_mat;

if nargout > 1
    varargout{1} = valid_link_label;
end
end