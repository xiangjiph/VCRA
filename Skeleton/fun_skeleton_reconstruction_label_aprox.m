function mask_label = fun_skeleton_reconstruction_label_aprox(skeleton_voxel, skeleton_radius, skeleton_label, block_size, max_error_rate)
% fun_skeleton_reconstruction_label reconstruct the vessel labeled array. 
% Input: 
%   skeleton_voxel: numerical vector, linear indices location of the vessel
%   centerline
%   skeleton_radius: numerical vector, radius of the vessel center line
%   voxle 
%   skeleton_lable: lable of the centerline voxel, can be the label of the
%   link and node
%   block_size: 3-by-1 numerical vector, size of the mask for
%   reconstruction 
%   max_error_rate: numerical scalar less than 1.Default value is 0.1. The
%   relative error for vessels of difference size
% Output: 
%   mask_label: reconstruction of the vessel mask, with each voxel labeled
%   with the label of the skeleton voxel from which it is reconstructed. 
%
% Implemented by Xiang Ji on 02/20/2019   
if nargin < 5
    max_error_rate = 0.1;
end

label_max = max(abs(skeleton_label));
if label_max < intmax('int8')
    mask_label = zeros(block_size, 'int8');
elseif label_max < intmax('int16')
    mask_label = zeros(block_size, 'int16');
elseif label_max < intmax('int32')
    mask_label = zeros(block_size, 'int32');
else
    mask_label = zero(block_size, 'int64');
end
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
% Reconstruct from larger scale to smaller scale
for iter_scale = num_scale : -1 : 1
    tmp_scale_ratio = downsample_ratio(iter_scale);
    tmp_Q = (downsample_min_scale(iter_scale) < skeleton_radius & ...
        skeleton_radius <= downsample_max_scale(iter_scale));
    tmp_r = skeleton_radius(tmp_Q);
    tmp_ind = skeleton_voxel(tmp_Q);
    tmp_label = skeleton_label(tmp_Q);
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
        tmp_label = tmp_label(tmp_valid_Q, :);
        tmp_ind = sub2ind(tmp_block_size, tmp_sub_rz(:, 1), tmp_sub_rz(:, 2), ...
            tmp_sub_rz(:, 3));
        [tmp_ind, tmp_unique_idx, ~] = unique(tmp_ind, 'stable');
        tmp_r = tmp_r(tmp_unique_idx);
        tmp_label = tmp_label(tmp_unique_idx);
    end
    if ~isempty(tmp_ind)
        tmp_mask = fun_skeleton_reconstruction_label(tmp_ind, tmp_r, tmp_label, tmp_block_size);
        if tmp_scale_ratio > 1
            tmp_mask = imresize3(tmp_mask, block_size, 'Method', 'nearest');
        end
        % Overwrite mask label
        tmp_ind = find(tmp_mask ~=0);
        mask_label(tmp_ind) = tmp_mask(tmp_ind);    
    end
end
end

