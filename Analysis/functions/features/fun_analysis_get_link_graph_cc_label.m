function [link_graph_cc_label, varargout]= fun_analysis_get_link_graph_cc_label(vessel_graph)
% fun_analysis_get_link_graph_cc_label gives the link in the same
% graph connected component the same label
% Input: 
%   vessel_graph: output of fun_skeleton_to_graph
% Output: 
%   link_graph_cc_label: N-by-1 numerical vector, where N is the number of
%   link in the graph
% 
num_link = vessel_graph.link.num_cc;
num_labeled_link = 0;
tmp_cc_label = 0;
link_graph_cc_label = zeros(num_link, 1);
node_graph_cc_label = zeros(vessel_graph.node.num_cc, 1);

while num_labeled_link < num_link
    tmp_current_link_label = find(link_graph_cc_label == 0, 1);
    
    tmp_cc_label = tmp_cc_label + 1;
    while ~isempty(tmp_current_link_label)
        link_graph_cc_label(tmp_current_link_label) = tmp_cc_label;
        tmp_connected_node = cat(2, vessel_graph.link.connected_node_label(tmp_current_link_label, :));
        tmp_connected_node = unique(tmp_connected_node(tmp_connected_node > 0));
        tmp_unlabeled_node_Q = (node_graph_cc_label(tmp_connected_node) == 0);
        node_graph_cc_label(tmp_connected_node(tmp_unlabeled_node_Q)) = tmp_cc_label;        
        
        tmp_connected_link = cat(1, vessel_graph.node.connected_link_label{tmp_connected_node});
        tmp_not_labeled_Q = (link_graph_cc_label(tmp_connected_link) == 0);
        tmp_current_link_label = tmp_connected_link(tmp_not_labeled_Q);        
    end
    num_labeled_link = nnz(link_graph_cc_label);
end
assert(all(link_graph_cc_label > 0), 'Exist unlabeled links');
if nargout > 1
    is_isolated_node_Q = node_graph_cc_label <= 0;
    if any(is_isolated_node_Q)
        fprintf('Exist unlabeled nodes. Label as NaN.\n');
        node_max_label = max(node_graph_cc_label);
        node_graph_cc_label(is_isolated_node_Q) = (node_max_label + 1) : ...
            1 : (node_max_label + nnz(is_isolated_node_Q));
    end
    varargout{1} = node_graph_cc_label;
end

end