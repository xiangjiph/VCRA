vessel_graph = load('Demo\Data\Demo_vessel_graph.mat');
% Get all the voxel indices
vxl_ind = cat(1, vessel_graph.link.cc_ind{:}, vessel_graph.node.cc_ind{:}, ...
    vessel_graph.isoloop.cc_ind{:}, vessel_graph.isopoint.pos_ind);
num_voxel = numel(vxl_ind);
% Convert the voxel indices to subscripts in 3D array
% The voxel resolution of the vessel graph is 1 x 1 x 1 micrometer^3. The
% size of the brain is ~ 10 mm on each dimension. Therefore, uint16 is
% sufficient for the voxel coordiantes
vxl_sub = zeros([num_voxel, 3], 'uint16');
[vxl_sub(:, 1), vxl_sub(:, 2), vxl_sub(:, 3)] = ind2sub(vessel_graph.num.mask_size, vxl_ind);
% Get the neighboring voxel list indices (26-neighbors)
% Pad the array 
mask_size_pad = vessel_graph.num.mask_size + 2;
vxl_ind_pad = sub2ind(mask_size_pad, vxl_sub(:, 1) + 1, vxl_sub(:, 2) + 1,...
    vxl_sub(:, 3) + 1);
vxl_neighbor_ind = sparse(vxl_ind_pad, 1, 1 : num_voxel, ...
    prod(mask_size_pad), 1);

vxl_neighbor_ind_add = fun_skeleton_neighbor_add_coeff_3D(mask_size_pad, 26, true);
connected_vxl_list_idx = full(vxl_neighbor_ind(bsxfun(@plus, vxl_ind_pad', vxl_neighbor_ind_add)));
has_neighbor_Q = connected_vxl_list_idx > 0;
% Count the number of neighboring voxels in 26-neighbors for each voxel
num_neighbor_vxl = sum(has_neighbor_Q, 1);
% Get the connected voxel pairs: N-by-2 matrix, where N is the number of
% connected voxel pairs
voxel_connection_list = cat(2, repelem([1 : num_voxel]', num_neighbor_vxl', 1), ...
    connected_vxl_list_idx(has_neighbor_Q));
