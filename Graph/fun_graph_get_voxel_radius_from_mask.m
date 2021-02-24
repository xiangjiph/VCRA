function graph_str = fun_graph_get_voxel_radius_from_mask(graph_str, mask, distance_type)
% graph_get_voxel_radius_from_mask extract the radius for each skeleton
% voxel from the distance transform of the mask(the same mask for comptuing
% the skeleton and then the graph). 
% Input:
%   graph_str: struct output by skeleton3D_to_graph
%   mask: mask for computng the distance and previously the skeleton
%   distance_type: type of distance used for distance transform. 
% Output: 
%   graph_str: a list of radius for each voxel of the links and
%   nodes(link.link_radius, node.node_radius) and their median value
%   (link_med_radius, node_med_radius) are added to the original graph_str
if nargin < 3
    distance_type = 'euclidean';
end


mask_distance_transform = bwdist(~mask,distance_type);
% Link
for link_idx = 1 : length(graph_str.link)
    link_voxel_list = graph_str.link(link_idx).link_voxels;
    radius_list = zeros(graph_str.link(link_idx).link_length,1);
    for voxel_idx = 1 : graph_str.link(link_idx).link_length
        radius_list(voxel_idx) = mask_distance_transform(link_voxel_list(voxel_idx));
    end
    graph_str.link(link_idx).link_radius = radius_list;
    graph_str.link(link_idx).link_med_radius = median(radius_list);
end
% Node
for node_idx = 1 : length(graph_str.node)
    node_voxel_list = graph_str.node(node_idx).idx;
    radius_list = zeros(length(node_voxel_list),1);
    for voxel_idx = 1 : length(node_voxel_list)
        radius_list(voxel_idx) = mask_distance_transform(node_voxel_list(voxel_idx));
    end
    graph_str.node(node_idx).node_radius = radius_list;
    graph_str.node(node_idx).node_med_radius = median(radius_list);
end
end
