function skl_mask = fun_graph_delete_node_voxel(skl_mask)

neighbor_count = convn(single(skl_mask > 0), ones(3,3,3, 'single'), 'same');
is_node_Q = neighbor_count > 3;
skl_mask(is_node_Q) = 0;
end