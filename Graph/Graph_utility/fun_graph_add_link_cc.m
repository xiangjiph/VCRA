function new_graph = fun_graph_add_link_cc(old_graph, link_cc_ind)



%% Approach One: recompute the entire graph 
if isempty(link_cc_ind)
    new_graph = old_graph;
    return;
end
voxel_list = unique(cat(1, old_graph.node.cc_ind{:}, old_graph.link.cc_ind{:}, ...
    old_graph.endpoint.pos_ind, old_graph.isopoint.pos_ind, ...
    link_cc_ind{:}));
new_graph = fun_skeleton_to_graph(voxel_list, old_graph.num.mask_size);
%% Approach Two: direct operation on the graph 
% Case: 
% 1. Endpoint to endpoint: 
%   a. Merge two link with single endpoints - delete 1 link, modify the
%   link cc list, zero one of the link label, reduce the number of links,
%   re-generate the map ...
%   b. Delete two endpoints
% 2. Endpoint to link
%   a. Merge original link with tne linker. 
%   b. Create a node at the other end of the linker - a new 3 degree link -
%   split the original link, add a new link label, record the new
%   connection
%   c. Delete one endpoint
% 3. Endpoint to node
%   a. Merge original link with the linker
%   b. The degree of node + 1, new connection, 
% To summarize, don't do it now...
% num_add_link = numel(link_cc_ind);
% ep_1_ind = zeros(num_add_link, 1);
% ep_2_ind = zeros(num_add_link, 1);
% 
% for iter_ep = 1 : num_add_link
%     tmp_cc = link_cc_ind{iter_ep};
%     tmp_e1_ind = tmp_cc(1);
%     tmp_e2_ind = tmp_cc(end);
%     % Add single linker to the graph 
%     
%     
%     ep_1_ind(iter_ep) = tmp_cc(1);
%     ep_2_ind(iter_ep) = tmp_cc(end);
% end
% e1_is_endpointQ = full(old_graph.endpoint.map_ind_2_label(ep_1_ind)) > 0;
% e1_is_linkQ = full(old_graph.link.map_ind_2_label(ep_1_ind)) > 0;
% e1_is_nodeQ = full(old_graph.node.map_ind_2_label(ep_1_ind)) > 0;
% e1_is_isopointQ = (~e1_is_endpointQ) & (~e1_is_linkQ) & (~e1_is_nodeQ);
% 
% e2_is_endpointQ = full(old_graph.endpoint.map_ind_2_label(ep_2_ind)) > 0;
% e2_is_linkQ = full(old_graph.link.map_ind_2_label(ep_2_ind)) > 0;
% e2_is_nodeQ = full(old_graph.node.map_ind_2_label(ep_2_ind)) > 0;
% e2_is_isopointQ = (~e2_is_endpointQ) & (~e2_is_linkQ) & (~e2_is_nodeQ);
end