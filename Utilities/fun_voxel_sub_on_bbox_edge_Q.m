function on_edge_Q = fun_voxel_sub_on_bbox_edge_Q(voxel_sub, bbox_mmxx)
% 
if isempty(voxel_sub)
    on_edge_Q = [];
    return;
end
assert((size(voxel_sub, 2) == 3) && (numel(bbox_mmxx) == 6));
on_edge_Q = any(bsxfun(@eq, voxel_sub, bbox_mmxx(1:3)), 2) |...
    any(bsxfun(@eq, voxel_sub, bbox_mmxx(4:6)), 2);
end