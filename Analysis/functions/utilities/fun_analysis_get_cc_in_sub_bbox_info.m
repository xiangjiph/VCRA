function inbbox_str = fun_analysis_get_cc_in_sub_bbox_info(grid_c_info, vessel_graph, grid_c_label)


grid_array_sub = grid_c_info.bbox_grid_sub_list(grid_c_label, :);
grid_array_ind = sub2ind(grid_c_info.grid_size, grid_array_sub(1), grid_array_sub(2), grid_array_sub(3));
tmp_internal_subgrid_bbox_mmxx_dataset = grid_c_info.internal_subgrid_bbox_mmxx{grid_array_ind};
grid_c_mmxx_pixel = grid_c_info.bbox_xyz_mmxx_pixel_list(grid_c_label, :);

tmp_internal_subgrid_bbox_mmxx_grid_c = tmp_internal_subgrid_bbox_mmxx_dataset - ...
    [grid_c_mmxx_pixel(1:3) - 1, grid_c_mmxx_pixel(1:3) - 1];
assert(~iscolumn(tmp_internal_subgrid_bbox_mmxx_grid_c), 'The internal subgrid bbox mmxx list is a column vector');
num_internal_subgrid = size(tmp_internal_subgrid_bbox_mmxx_grid_c, 1);

link_label_in_bbox_cell = cell(num_internal_subgrid, 1);
node_label_in_bbox_cell = cell(num_internal_subgrid, 1);
link_voxel_in_bbox_ratio_cell = cell(num_internal_subgrid, 1);
for iter_int_grid = 1 : num_internal_subgrid
    tmp_int_bbox = tmp_internal_subgrid_bbox_mmxx_grid_c(iter_int_grid, :);
    % All the links that pass through the bounding box
    tmp_link_sub = fun_ind2sub(vessel_graph.num.mask_size, vessel_graph.link.pos_ind);
    tmp_node_sub = fun_ind2sub(vessel_graph.num.mask_size, vessel_graph.node.pos_ind);
    tmp_link_voxel_in_bbox_Q = all(bsxfun(@ge, tmp_link_sub, tmp_int_bbox(1:3)), 2) & ...
        all(bsxfun(@le, tmp_link_sub, tmp_int_bbox(4:6)), 2);
    tmp_node_voxel_in_bbox_Q = all(bsxfun(@ge, tmp_node_sub, tmp_int_bbox(1:3)), 2) & ...
        all(bsxfun(@le, tmp_node_sub, tmp_int_bbox(4:6)), 2);
    tmp_voxel_link_label = full(vessel_graph.link.map_ind_2_label(...
        vessel_graph.link.pos_ind(tmp_link_voxel_in_bbox_Q)));
    tmp_voxel_node_label = full(vessel_graph.node.map_ind_2_label(...
        vessel_graph.node.pos_ind(tmp_node_voxel_in_bbox_Q)));
    
    [tmp_voxel_link_label_cell, tmp_link_in_bbox_label] = fun_bin_data_to_idx_list(tmp_voxel_link_label);
    
    tmp_node_in_bbox_label = unique(tmp_voxel_node_label);
    
    tmp_link_num_voxel = vessel_graph.link.num_voxel_per_cc(tmp_link_in_bbox_label);
    tmp_link_num_voxel_in_bbox = cellfun(@numel, tmp_voxel_link_label_cell);
    tmp_link_in_bbox_ratio = (tmp_link_num_voxel_in_bbox ./ tmp_link_num_voxel);
    
    link_label_in_bbox_cell{iter_int_grid} = tmp_link_in_bbox_label;
    node_label_in_bbox_cell{iter_int_grid} = tmp_node_in_bbox_label;
    link_voxel_in_bbox_ratio_cell{iter_int_grid} = tmp_link_in_bbox_ratio;
end
inbbox_str.link_label_in_bbox_cell = link_label_in_bbox_cell;
inbbox_str.node_label_in_bbox_cell = node_label_in_bbox_cell;
inbbox_str.link_voxel_in_bbox_ratio_cell = link_voxel_in_bbox_ratio_cell;

internal_link_lable = unique(cat(1, link_label_in_bbox_cell{:}));
num_internal_link = numel(internal_link_lable);
internal_node_label = unique(cat(1, node_label_in_bbox_cell{:}));
num_internal_node = numel(internal_node_label);

link_label_2_list_ind = sparse(internal_link_lable, ones(num_internal_link, 1), ...
    1 : num_internal_link);
node_label_2_list_ind = sparse(internal_node_label, ones(num_internal_node, 1), ...
    1 : num_internal_node);

inbbox_str.internal_link_label = internal_link_lable;
inbbox_str.num_internal_link = num_internal_link;
inbbox_str.internal_node_label = internal_node_label;
inbbox_str.num_internal_node = num_internal_node;
inbbox_str.link_label_2_list_ind = link_label_2_list_ind;
inbbox_str.node_label_2_list_ind = node_label_2_list_ind;
end