function local_mask = fun_vis_voxel_reconstruction_local(voxel_sub, voxel_radius, label)
% This function find the local bounding box for the input voxels,
% reconstruct the vessel mask for these voxels and display in the local
% coordinate. 
% Input: 
%   voxel_sub: N-by-3 numerical array, coordiante of the voxels
%   voxel_radius: N-by-1 numerical vector, radius of the voxels
% Output: 
%   local_mask: local reconstruction of the input voxels
%
% Implemented by Xiang Ji on 07/15/2019

% Determine the bounding box for viusalization
if nargin < 3
    label = [];
end
max_r = ceil(max(voxel_radius, [], 'omitnan'));
bbox_min = min(voxel_sub - max_r, [], 1);
bbox_max = max(voxel_sub + max_r, [], 1);
bbox_size = bbox_max - bbox_min + 1;
voxel_sub_local = voxel_sub - bbox_min + 1;
voxel_ind_local = sub2ind(bbox_size, voxel_sub_local(:, 1), voxel_sub_local(:, 2), ...
    voxel_sub_local(:, 3));
if isempty(label) || numel(label) ~= numel(voxel_radius)
    local_mask = fun_skeleton_reconstruction(voxel_ind_local, voxel_radius, bbox_size);
else
    local_mask = fun_skeleton_reconstruction_label(voxel_ind_local, voxel_radius, label, bbox_size);
    % Shift the mask label for visualization
%     unique_label = unique(label, 'sorted');
end
figure;
vol_hdl = volshow(local_mask);
vol_hdl.Colormap = jet(256);
end