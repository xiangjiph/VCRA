function voxel_sub = fun_stitching_merge_surface_voxels(voxel_sub, merge_block_size)
% fun_stitching_merge_surface_voxels merges the input voxel subscript list
% by computing the average position of the voxels within the block, whose
% size is specified by merge_block_size
% Input: 
%   voxel_sub: N-by-3, coordinate of the voxel position in 3D space
%   merge_block_size: 1-by-3 numerical vector, size of the block for merging. 
% Output: 
%   voxel_sub: N'-by-3 numerical array after merging
sub_min = min(voxel_sub(:, 1:3), [], 1);
sub_max = max(voxel_sub(:, 1:3), [], 1);
image_size = sub_max - sub_min + 1;
downsampled_image_size = ceil(image_size ./ merge_block_size);

num_edge_voxel = zeros(downsampled_image_size);

mean_value = zeros([size(voxel_sub,2), downsampled_image_size]);
% mean_edge_pos_2 = zeros(downsampled_image_size);
% mean_edge_pos_3 = zeros(downsampled_image_size);
bbox_sub = ceil((1 + bsxfun(@minus, voxel_sub(:, 1:3), sub_min))./ merge_block_size);
for iter1 = 1 : size(bbox_sub, 1)
    num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + 1;
    mean_value(:, bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
        mean_value(:, bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, :)';
%     mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 2);
%     mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 3);
end
voxel_sub = reshape(permute(mean_value, [2,3,4,1]), prod(downsampled_image_size), [])./ max(1, num_edge_voxel(:));
voxel_sub = voxel_sub(all(voxel_sub(:, 1:3) > 0, 2),:);
end