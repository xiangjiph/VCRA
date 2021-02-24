function vessel_graph = fun_graph_delete_new_links_with_endpoint(vessel_graph, old_map_ep_ind_2_label, max_rm_num_voxel, verboseQ)
% fun_graph_delete_new_links_with_endpoint finds the links with endpoints
% due to graph annotation and remove these link. 

if nargin < 3
    max_rm_num_voxel = 5;
    verboseQ = false;
elseif nargin < 4
    verboseQ = false;
end
exist_in_old_Q = full(old_map_ep_ind_2_label(vessel_graph.endpoint.pos_ind));
new_endpoint_ind = vessel_graph.endpoint.pos_ind(~exist_in_old_Q);
if isempty(new_endpoint_ind)
    return;
end
new_endpoint_link_label = full(vessel_graph.link.map_ind_2_label(new_endpoint_ind));
assert(all(new_endpoint_link_label > 0), 'Exist endpoint not in the link');
% Find the unique indices, since some of the link can have two endpoints 
new_endpoint_link_label = unique(new_endpoint_link_label); 
new_endpoint_link_num_voxel = vessel_graph.link.num_voxel_per_cc(new_endpoint_link_label);
link_label_to_remove = new_endpoint_link_label(new_endpoint_link_num_voxel <= max_rm_num_voxel);
[vessel_graph, new_endpoint_ind] = fun_graph_delete_internal_links(vessel_graph, link_label_to_remove);

if verboseQ
    fprintf('Number of links to be deleted: %d\n', numel(link_label_to_remove));
    fprintf('Length of the links with endpoints are: \n');
    disp(new_endpoint_link_num_voxel');
    fprintf('Number of newly created endpoints: %d\n', numel(new_endpoint_ind));
end
end