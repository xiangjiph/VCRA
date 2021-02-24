function [node_cc_label, varargout]= fun_analysis_get_node_graph_cc_label(vessel_graph)
% fun_analysis_get_link_graph_cc_label gives the node in the same
% graph connected component the same label
% Input: 
%   vessel_graph: output of fun_skeleton_to_graph
% Output: 
%   node_graph_cc_label: N-by-1 numerical vector, where N is the number of
%   node in the graph
%   num_cc: numerical scalar, number of graph connected components 
num_node = vessel_graph.node.num_cc;
node_cc_label = zeros(num_node, 1);
num_labeled_node = 0;
tmp_cc_label = 0;
while num_labeled_node < num_node
   tmp_current_node_label = find(node_cc_label == 0, 1);
   tmp_cc_label = tmp_cc_label + 1;
   while ~isempty(tmp_current_node_label)
       node_cc_label(tmp_current_node_label) = tmp_cc_label;
       num_labeled_node = num_labeled_node + numel(tmp_current_node_label);
       tmp_connected_link = cat(1, vessel_graph.node.connected_link_label{tmp_current_node_label});
       tmp_connected_link = unique(tmp_connected_link(tmp_connected_link > 0));
       tmp_connected_node = cat(1, vessel_graph.link.connected_node_label(tmp_connected_link, :));
       tmp_connected_node = unique(tmp_connected_node(tmp_connected_node > 0));
       tmp_not_labeled_Q = (node_cc_label(tmp_connected_node) == 0);
       tmp_current_node_label = tmp_connected_node(tmp_not_labeled_Q);
   end
end
if nargout > 1
    varargout{1} = tmp_cc_label;
end
assert(all(node_cc_label > 0), 'Exist unlabeled nodes');
end