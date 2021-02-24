function new_link_str = fun_graph_split_node_get_new_link(cc_ind_1, cc_ind_2, node_ind, vessel_mask_dt, vessel_skl, num_voxel_fit)
% fun_graph_delete_node_get_new_link find the shortest path that starts at
% one point, go through a point and end at the other point. 
% Input: 
%   midpoint_str: output of fun_graph_delete_node_get_link_pair_mid_point
%   vessel_image: 3D numerical array
%   vessel_mask: 3D logical array, the mask of the vessels
%   vessel_skl: 3D logical array, the skeleton of the vessels
% Output: 
%   new_link_str: strcture with fields: 
%       link_sub_w_ep: subscript of link voxel list, including two
%       endpoints. 
%       link_ind_w_ep: index of the link voxel list, including two
%       endpoints. 
%       num_voxel: number of voxels in the link voxel list
%       connectivity: connectivity used for find the shortest path, default
%       value is 26;
%
% Implemented by Xiang Ji on 02/16/2019
%
% Modified by Xiang Ji on 02/17/2019
% 1. Merge with fun_graph_delete_node_get_link_pair_mid_point
if nargin < 6
    num_voxel_fit = 8;
end

image_size = size(vessel_mask_dt);

cc_sub_1 = fun_ind2sub(image_size, cc_ind_1);
cc_sub_2 = fun_ind2sub(image_size, cc_ind_2);
node_sub = fun_ind2sub(image_size, node_ind);
if ~isvector(node_sub)
    node_sub = mean(node_sub, 1);
end
node_dt = max(vessel_mask_dt(node_ind));

num_voxel_removed = round(max(node_dt/sqrt(2)));
% Remove the voxels near the node
num_cc_1 = size(cc_sub_1, 1);
kept_cc_sub_1 = cc_sub_1(min(num_voxel_removed + 1, num_cc_1):num_cc_1, :);
fit_sub_1 = kept_cc_sub_1(1 : min(num_voxel_fit, size(kept_cc_sub_1,1)), :);

num_cc_2 = size(cc_sub_2, 1);
kept_cc_sub_2 = cc_sub_2(min(num_voxel_removed + 1, num_cc_2):num_cc_2, :);
fit_sub_2 = kept_cc_sub_2(1 : min(num_voxel_fit, size(kept_cc_sub_2,1)), :);
%% Find the closest point on the linear fit of the nearby paired skeleton to the node
% 3D linear regression 
fit_sub = [fit_sub_1; fit_sub_2];
fit_sub_mean = mean(fit_sub, 1);
fit_sub_c = fit_sub - fit_sub_mean;
cov_mat = (fit_sub_c' * fit_sub_c)./(size(fit_sub_c, 1) - 1);
% Diagonal the covariant matrix and get the largest eigenvalue and its
% eigenvector
[eig_vec, eig_val] = eig(cov_mat, 'vector');
[eig_val, tmp_idx] = sort(eig_val, 'descend');
eig_vec = eig_vec(:, tmp_idx);

% Closest point on the fitted line to the node position 
fit_vec = eig_vec(:,1);
sub_closest_to_node = fit_sub_mean' + ((node_sub - fit_sub_mean) * fit_vec) * fit_vec;

% Save information 
midpoint_str = struct;
midpoint_str.pos = sub_closest_to_node';
midpoint_str.sub = round(sub_closest_to_node');
midpoint_str.dist_midpoint2node = sqrt(sum((node_sub - midpoint_str.pos).^2));
midpoint_str.skl_sub_mean = fit_sub_mean;
midpoint_str.cov = cov_mat;
midpoint_str.cov_eig_vec = eig_vec;
midpoint_str.cov_eig_val = eig_val;
midpoint_str.fit_cov = fit_vec;
% Estimate the standard deviation in the radial direction 
midpoint_str.std_r = sqrt(eig_val(2) + eig_val(3));
%% Get the new link 
ep1_sub = kept_cc_sub_1(1, :);
ep2_sub = kept_cc_sub_2(1, :);
mid_sub = round(sub_closest_to_node');
dist_mid2ep1 = sqrt(sum((mid_sub - ep1_sub).^2));
dist_mid2ep2 = sqrt(sum((mid_sub - ep2_sub).^2));
max_dist_mid2ep = max(dist_mid2ep1, dist_mid2ep2);
bbox_r = max(2, round(max_dist_mid2ep * 1.5));

[local_mask, local_bbox_mmxx] = crop_center_box(vessel_mask_dt, mid_sub, bbox_r);
local_mask = local_mask > 0;
local_bbox_size = local_bbox_mmxx(4:6) - local_bbox_mmxx(1:3) + 1;
ep1_local_sub = ep1_sub - local_bbox_mmxx(1:3) + 1;
ep2_local_sub = ep2_sub - local_bbox_mmxx(1:3) + 1;
mid_local_sub = mid_sub - local_bbox_mmxx(1:3) + 1;

ep1_local_ind = sub2ind(local_bbox_size, ep1_local_sub(1), ep1_local_sub(2), ep1_local_sub(3));
ep2_local_ind = sub2ind(local_bbox_size, ep2_local_sub(1), ep2_local_sub(2), ep2_local_sub(3));
mid_local_ind = sub2ind(local_bbox_size, mid_local_sub(1), mid_local_sub(2), mid_local_sub(3));

line_mask_mid2ep1 = fun_graph_get_p2p_cylinder_mask(local_bbox_size, mid_local_sub, ep1_local_sub, 2);
line_mask_mid2ep2 = fun_graph_get_p2p_cylinder_mask(local_bbox_size, mid_local_sub, ep2_local_sub, 2);
search_mask = imclose(line_mask_mid2ep1 | line_mask_mid2ep2, strel('cube', 2)) & local_mask;
% Remove all the nearby skeleton voxels form the search mask, add the three
% points back 
local_skl = crop_center_box(vessel_skl, mid_sub, bbox_r);
search_mask = search_mask & (~local_skl);
search_mask(ep1_local_ind) = true;
search_mask(ep2_local_ind) = true;
search_mask(mid_local_ind) = true;
% Check connectivity
search_mask_label = bwlabeln(search_mask);
assert(search_mask_label(ep1_local_ind) == search_mask_label(ep2_local_ind) && ...
    search_mask_label(ep2_local_ind) == search_mask_label(mid_local_ind) && ...
    search_mask_label(mid_local_ind) ~= 0,...
    'Either one of the points to be link is out of the search mask or they are not in the same connected component');
% Find the shortest path
search_metric = fun_graph_get_shortest_path_distance_metric(ones(local_bbox_size), search_mask);

connectivity = 26;
shortest_path_str = fun_graph_shortest_path_through_3pts(search_metric, ep1_local_sub, mid_local_sub, ep2_local_sub, connectivity);
%% Save information
% Convert the path to the global cooredinate
linker_sub_global = shortest_path_str.link_sub_w_ep + local_bbox_mmxx(1:3) - 1;
linker_ind_global = sub2ind(image_size, linker_sub_global(:,1), linker_sub_global(:,2), linker_sub_global(:,3));
kept_cc_ind_1 = sub2ind(image_size, kept_cc_sub_1(:,1), kept_cc_sub_1(:,2), kept_cc_sub_1(:,3));
kept_cc_ind_2 = sub2ind(image_size, kept_cc_sub_2(:,1), kept_cc_sub_2(:,2), kept_cc_sub_2(:,3));
new_link_str.link_cc_ind = cat(1, flip(kept_cc_ind_1), linker_ind_global(2:end-1), ...
    kept_cc_ind_2);
new_link_str.node_cc_ind = node_ind;
new_link_str.connectivity = connectivity;
new_link_str.midpoint_info = midpoint_str;
end