function in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(sub, bbox_mmxx)
% Test if the voxel is in the bounding box specified by the maximum and
% minimum value of the coordinates
% Input: 
%   sub: N-by-M numerical array, where N is the number of voxels
%   bbox_mmxx: 
%   bbox_mmxx: bounding box [min_1, min_2, min_3, max_1, max_2, max_3]
% Output: 
%   in_bbox_Q: true is the voxel is inside the bounding box
num_dim = numel(bbox_mmxx) / 2;
if ~all(bbox_mmxx(1:num_dim) <= bbox_mmxx(num_dim+1 : 2 * num_dim))
    warning('bbox_mmxx has negative box length');
    in_bbox_Q = false(size(sub, 1), 1);
end
assert(size(sub, 2) == num_dim, ...
    'The number of column in the coordinate array should be consistent with the bounding box parameter');
in_bbox_Q = all(bsxfun(@ge, sub, bbox_mmxx(1:num_dim)), 2) & ...
    all(bsxfun(@le, sub, bbox_mmxx((num_dim + 1):end)), 2);
end