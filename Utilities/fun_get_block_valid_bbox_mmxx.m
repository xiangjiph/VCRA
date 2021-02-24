function bbox_mmxx = fun_get_block_valid_bbox_mmxx(grid_sub, grid_valid_array, block_size, overlap)


if ~islogical(grid_valid_array)
    grid_valid_array = logical(grid_valid_array);
end
assert(isscalar(overlap), 'The overlap should be a scalar');
grid_size = size(grid_valid_array);
valid_array_pad = padarray(grid_valid_array, [1,1,1], 0, 'both');

ind_add_pad = fun_skeleton_neighbor_add_coeff_3D(grid_size + 2, 6, true);
grid_ind_pad = sub2ind(grid_size + 2, grid_sub(1) + 1, grid_sub(2) + 1, ...
    grid_sub(3) + 1);
has_valid_neighbor_Q = valid_array_pad(bsxfun(@plus, grid_ind_pad, ind_add_pad))';
% The order of the adder is [top, left, back, front, right, bottom] and we
% want [back, left, top, front, right, bottom];
has_valid_neighbor_Q = has_valid_neighbor_Q([3,2,1,4,5,6]);

half_overlap = overlap / 2;
assert(mod(half_overlap, 1) == 0, 'Half of the overlap is not an integer');

bbox_mmxx = [1, 1, 1, block_size];
bbox_mmxx(1:3) = bbox_mmxx(1:3) + has_valid_neighbor_Q(1:3) .* (half_overlap);
bbox_mmxx(4:6) = bbox_mmxx(4:6) - has_valid_neighbor_Q(4:6) .* (half_overlap);
end