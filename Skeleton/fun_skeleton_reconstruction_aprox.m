function mask = fun_skeleton_reconstruction_aprox(skeleton_voxel, skeleton_radius, block_size, max_error_rate)
% fun_skeleton_reconstruction_aprox perform multi-scale reconstruction. The
% accuracy of the reconstruction is controled by the accuracy_goal. 
% Input: 
%   skeleton_voxel: numerical vector, skeleton voxel indices
%   skeleton_radius: numerical vector, radius of the skeleton voxel, or,
%   distance between the skeleton and the nearest vessel mask edge
%   block_size: size of the block, for converting the skeleton voxel
%   indices into subscripts. 
%   max_error_rate: numerical scalar less than 1.Default value is 0.1. The
%   relative error for vessels of difference size
% Output: 
%   mask: 3D logical array of size block_size
%
% Implemented by Xiang Ji on 04/05/2019
if nargin < 4
    max_error_rate = 0.1; %
end
mask = false(block_size);
if isempty(skeleton_voxel)
    return
end
% Classify skeleton voxel
max_radius = max(skeleton_radius);
max_downsample_ratio = max(0, floor(log2(max_radius * max_error_rate)));

downsample_ratio = 2 .^ (0 : max_downsample_ratio);
downsample_max_scale = downsample_ratio./ max_error_rate;
downsample_max_scale(end) = ceil(max_radius);
downsample_min_scale = [0, downsample_max_scale(1:end-1)];
num_scale = numel(downsample_ratio);
for iter_scale = 1 : num_scale
    tmp_scale_ratio = downsample_ratio(iter_scale);
    tmp_Q = (downsample_min_scale(iter_scale) < skeleton_radius & ...
        skeleton_radius <= downsample_max_scale(iter_scale));
    tmp_r = skeleton_radius(tmp_Q);
    tmp_ind = skeleton_voxel(tmp_Q);
    tmp_block_size = round(block_size ./ tmp_scale_ratio);
    if tmp_scale_ratio > 1
        % Recompute the indices of the skeleton pixel and radius
        tmp_r = round(tmp_r ./ tmp_scale_ratio);
        tmp_sub = fun_ind2sub(block_size, tmp_ind);
        tmp_sub_rz = round(tmp_sub ./ tmp_scale_ratio);
        % Check in the bounding box
        tmp_valid_Q = all(bsxfun(@ge, tmp_sub_rz, 1), 2) & ...
            all(bsxfun(@le, tmp_sub_rz, tmp_block_size), 2);
        tmp_sub_rz = tmp_sub_rz(tmp_valid_Q, :);
        tmp_r = tmp_r(tmp_valid_Q, :);
        tmp_ind = sub2ind(tmp_block_size, tmp_sub_rz(:, 1), tmp_sub_rz(:, 2), ...
            tmp_sub_rz(:, 3));
        [tmp_ind, tmp_unique_idx, ~] = unique(tmp_ind, 'stable');
        tmp_r = tmp_r(tmp_unique_idx);
    end
    tmp_mask = fun_skeleton_reconstruction(tmp_ind, tmp_r, tmp_block_size);
    if tmp_scale_ratio > 1
        tmp_mask = imresize3(uint8(tmp_mask), block_size, 'Method', 'nearest') > 0;
    end
    mask = mask | tmp_mask;
end
end

