function gap_linking_str = fun_graph_get_linker_ep2ep_th_rlx(vessel_image, vessel_mask, vessel_skl, p_1_pos, p_2_pos, connectivity)
% fun_graph_connect_gap_ep2ep finds the link that connect the two given
% endpoitns
% Input: 
%   vessel_image: 3D numerical array, image of the vessel 
%   vessel_mask: 3D logical array, mask of the vessel 
%   ep_1_pos: numerical scalar or 3-by-1 numerical array, subscript of
%   point 1 in vessel_image. If scalar, assumed to be the voxel index in
%   the image, converted to subscript automatically
%   ep_2_pos: same as above, for the secdong point
%   connectivity: 8, 18 or 26. Number of neighbors for searching the
%   shortest path
% Output: 
%   connect_link_str: structure with fields: 
%       link_ind_w_ep: numerical vector, indices of link voxels in the metric array,
%       including the two endpoints
%       total_dist: numerical scalar, sum of the metric value of the found link 
%       voxel
%       link_sub_wo_ep: N-by-3 numerical array, subscript of the link
%       voxels, without the subscript of the two input endpoints
%       connectivity: same as the input
%       (The following fields are computed iff vessel_mask is given)
%       mask_mean: mean of image value of the voxel in the mask
%       mask_std: standard deviation of image value of the voxel outside the mask
%       bg_mean: mean of image value of the voxel in the mask
%       bg_std: standard deviation of image value of the voxel outside the mask
%       
% Author: Xiang Ji
% Date: 01, 09, 2019
% Modified on 03/06/2019
% 1. Fix a bug when more than 1 link cc happen to be in the mask... This
% can happen when a link is incorrectly deleted in the graph refinement
% process....
% Modified on 03/14/2019
% 1. Better estimation of the signal to noise leel by reconstructing the
% linker mask 
% To do list: 
% 1. Use input parser. Add the option for not computing the features
if nargin < 6
    connectivity = 26;
end
% min_cc_size = 0;
image_size = size(vessel_image);
if isscalar(p_1_pos)
    p_1_pos = fun_ind2sub(image_size, p_1_pos);
end
if isscalar(p_2_pos)
    p_2_pos = fun_ind2sub(image_size, p_2_pos);
end
if iscolumn(p_1_pos)
    p_1_pos = p_1_pos';
end
if iscolumn(p_2_pos)
    p_2_pos = p_2_pos';
end
% straight_link_mask_r = 3;
r_ep_dt = 5;

% The middle position of two endpoints
eps_mid_sub = round((p_1_pos + p_2_pos)./2);
% Distance between the middle point and the two end point - fix this value
% for now
eps_mid_to_ep_r = [sqrt(sum((p_1_pos - eps_mid_sub).^2)), sqrt(sum((p_2_pos - eps_mid_sub).^2))];
eps_mid_to_ep_r_max = max(5, max(eps_mid_to_ep_r));
local_center_bbox_r = round(2 * eps_mid_to_ep_r_max );
ep_in_diff_cc_Q = false;

persistent str_template
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'foundQ',...
        'ep_1_sub', 'ep_1_ind', 'ep_2_sub', 'ep_2_ind', 'search_bbox_r', 'search_bbox_mmxx',...
        'threshold_list', 'bg_mean', 'bg_std', 'connected_th', ...
        'search_full_image_Q', 'link_sub_w_ep', 'link_ind_w_ep', 'num_voxel', ...
        'int', 'int_mean', 'int_std', 'int_med', 'connectivity', 'out_mask_Q', ...
        'int_o_mask', 'int_o_mask_mean', 'int_o_mask_std', 'link_ratio_o_mask', ...
        'link_sub_o_m', 'link_ind_o_m', 'recon_radius', 'recon_int_mean', 'recon_int_std',...
        'recon_int_median', 'recon_bg_int_mean', 'recon_bg_std', 'recon_SNR', 'recon_mask_ind', ...
        'recon_mask_om_ratio', 'closest_skl_ind', 'closest_skl_dt_2_r', 'closest_skl_dist_2_closest_ep_dist'});
end
gap_linking_str = str_template;
gap_linking_str.foundQ = false;
gap_linking_str.ep_1_sub = p_1_pos;
gap_linking_str.ep_2_sub = p_2_pos;
gap_linking_str.ep_1_ind = sub2ind(image_size, p_1_pos(1), p_1_pos(2), p_1_pos(3));
gap_linking_str.ep_2_ind = sub2ind(image_size, p_2_pos(1), p_2_pos(2), p_2_pos(3));

while ~ep_in_diff_cc_Q 
    gap_linking_str.search_bbox_r = local_center_bbox_r;    
    [local_int, local_bbox_mmxx] = crop_center_box(vessel_image, eps_mid_sub, local_center_bbox_r);
    gap_linking_str.search_bbox_mmxx = local_bbox_mmxx;
   
    local_mask = crop_center_box(vessel_mask, eps_mid_sub, local_center_bbox_r);
    local_int = single(local_int);
    local_mask_size = size(local_mask);
    % Position of the endpoint in the local coordinate
    ep_local_sub_1 = p_1_pos - local_bbox_mmxx(1:3) + 1;
    ep_local_sub_2 = p_2_pos - local_bbox_mmxx(1:3) + 1;
    if all(ep_local_sub_1 > 0) && all(ep_local_sub_1 <= local_mask_size)
        ep_local_ind_1 = sub2ind(local_mask_size, ep_local_sub_1(1), ep_local_sub_1(2), ...
            ep_local_sub_1(3));
    else
%         disp('Two endpoints are always in the same connected component. Terminate searching');
        return;
    end
    if all(ep_local_sub_2 > 0) && all(ep_local_sub_2 <= local_mask_size)
        ep_local_ind_2 = sub2ind(local_mask_size, ep_local_sub_2(1), ep_local_sub_2(2), ...
            ep_local_sub_2(3));
    else
%         disp('Two endpoints are always in the same connected component. Terminate searching');
        return;
    end
    local_cc = bwconncomp(local_mask, 26);
    local_mask_label = labelmatrix(local_cc);
    ep_cc_label_1 = local_mask_label(ep_local_ind_1);
    ep_cc_label_2 = local_mask_label(ep_local_ind_2);
    assert(ep_cc_label_1 >0 && ep_cc_label_2>0, 'The skeleton voxel is outside the vessel mask');
    if ep_cc_label_1 == ep_cc_label_2
        if local_center_bbox_r >= eps_mid_to_ep_r_max
%             disp('Two endpoints belong to the same connected component. Reduce the mask size');
            local_center_bbox_r = local_center_bbox_r  - max(1, min(round(local_center_bbox_r * 0.1), 5));
            if local_center_bbox_r <= 1
                return;
            end
        else
%             disp('Two endpoints are always in the same connected component. Terminate searching');
            return;
        end
    else
        ep_in_diff_cc_Q = true;
    end
end
%% Get search mask by distance transform
local_skl = crop_center_box(vessel_skl, eps_mid_sub, local_center_bbox_r);
local_skl_cc = bwconncomp(local_skl);
local_skl_label = labelmatrix(local_skl_cc);

ep1_cc_mask = (local_mask_label == ep_cc_label_1);
ep1_cc_link_mask = false(local_mask_size);
tmp_label_1 = local_skl_label(ep_local_ind_1);
ep1_cc_link_mask(local_skl_cc.PixelIdxList{tmp_label_1}) = true;
ep1_cc_link_mask = ep1_cc_link_mask & ep1_cc_mask;

ep2_cc_mask = (local_mask_label == ep_cc_label_2);
ep2_cc_link_mask = false(local_mask_size);
tmp_label_2 = local_skl_label(ep_local_ind_2);
ep2_cc_link_mask(local_skl_cc.PixelIdxList{tmp_label_2}) = true;
ep2_cc_link_mask = ep2_cc_link_mask & ep2_cc_mask;

dt_search_mask_ep1 = fun_graph_get_endpoint_link_mask(ep_local_ind_1, ep1_cc_link_mask, r_ep_dt);
dt_search_mask_ep2 = fun_graph_get_endpoint_link_mask(ep_local_ind_2, ep2_cc_link_mask, r_ep_dt);

dt_search_mask = dt_search_mask_ep1 & dt_search_mask_ep2;
% DataManager = FileManager;
% vis_mask = zeros(local_mask_size, 'uint8');
% vis_mask(local_mask) = 1;
% vis_mask(ep1_cc_link_mask) = 2;
% vis_mask(ep2_cc_link_mask) = 3;
% DataManager.visualize_itksnap(local_int, vis_mask);
%% Find connecting mask by threshold relaxation 
local_bg_mean = mean(local_int(~local_mask));
local_bg_std = std(local_int(~local_mask));
local_vessel_int_mean = mean(local_int(local_mask));
% Generate threshold list
num_step = round((local_vessel_int_mean - local_bg_mean)/local_bg_std);
th_list = local_bg_mean + (num_step : -1 : -3) * local_bg_std;
gap_linking_str.threshold_list = th_list;
gap_linking_str.bg_mean = local_bg_mean;
gap_linking_str.bg_std=  local_bg_std;
connected_mask = [];
for iter_th = 1 : length(th_list)
    test_th = th_list(iter_th);
    test_mask = (local_int > test_th & dt_search_mask);
    if ~test_mask(ep_local_ind_1) || ~test_mask(ep_local_ind_2)
        % Endpoint not in the mask
        continue;
    end
    test_mask_labeled = bwlabeln(test_mask, 26);
    if test_mask_labeled(ep_local_ind_1) == test_mask_labeled(ep_local_ind_2)
        % Connected. Use this mask
        connected_mask = (test_mask_labeled == test_mask_labeled(ep_local_ind_1));
        gap_linking_str.foundQ = true;
        break;        
    end
end
if isempty(connected_mask)
%     disp('No connected mask found. Terminate.');
    return;
end

gap_linking_str.connected_th = test_th;
if test_th == th_list(end)
    gap_linking_str.search_full_image_Q = true;
else
    gap_linking_str.search_full_image_Q = false;
end
%% Find matched boundary voxel position 
% Voxels in the origianl masks that connected to the connected mask: 
linker_mask_dilate = imdilate(connected_mask & ~ ep1_cc_mask & ~ep2_cc_mask, strel('cube', 3));
% Boundary voxels of the mask that the endpoint is in:
% ep1_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask_ep1 & ep1_cc_mask;
% ep2_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask_ep2 & ep2_cc_mask;
ep1_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask & ep1_cc_mask;
ep2_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask & ep2_cc_mask;
% Find the nearest boundary voxel pair: 
ep1_cc_boundary_ind = find(ep1_cc_mask_boundary_valid);
ep1_cc_boundary_sub = fun_ind2sub(local_mask_size, ep1_cc_boundary_ind);

ep2_cc_boundary_ind = find(ep2_cc_mask_boundary_valid);
ep2_cc_boundary_sub  = fun_ind2sub(local_mask_size, ep2_cc_boundary_ind);
% distance between distance points: 
boundary_voxel_dist = pdist2(ep1_cc_boundary_sub, ep2_cc_boundary_sub);
% Find the minimum distance: 
[~, tmp_min_ind] = min(boundary_voxel_dist(:));
[ep1_cc_bp_idx, ep2_cc_bp_idx] = ind2sub(size(boundary_voxel_dist), tmp_min_ind);

ep1_bv_sub = ep1_cc_boundary_sub(ep1_cc_bp_idx, :);
ep2_bv_sub = ep2_cc_boundary_sub(ep2_cc_bp_idx, :);
%% Find the shortest path 
% Exclude skeleton voxels except for the two endpoints
skl_wo_ep_mask = ep1_cc_link_mask & ep2_cc_link_mask;
skl_wo_ep_mask(ep_local_ind_1) = false;
skl_wo_ep_mask(ep_local_ind_2) = false;
connected_mask = imclose(connected_mask, strel('sphere', 2));
% % Add the connecting straight line between two endpoint
% straight_line_mask = fun_graph_get_p2p_cylinder_mask(size(local_mask), ep_local_sub_1, ep_local_sub_2, straight_link_mask_r);
search_mask = dt_search_mask & (ep1_cc_mask | ep2_cc_mask | connected_mask) & ~skl_wo_ep_mask;
search_metric = fun_graph_get_shortest_path_distance_metric(local_int, search_mask);
shortest_path_str = fun_graph_shortest_path_through_4pts(search_metric, ep_local_sub_1, ep1_bv_sub, ep2_bv_sub, ep_local_sub_2, connectivity);

% fun_vis_link_surrounding_by_cc_ind(local_int, local_mask_label == 2, local_skl, shortest_path_str.link_ind_w_ep, 60)
%% Save information
% Convert to the global image coordinate (with endpoints)
gap_linking_str.link_sub_w_ep = shortest_path_str.link_sub_w_ep + local_bbox_mmxx(1:3) - 1;
gap_linking_str.link_ind_w_ep = sub2ind(image_size, gap_linking_str.link_sub_w_ep(:,1), ...
    gap_linking_str.link_sub_w_ep(:,2), gap_linking_str.link_sub_w_ep(:,3));

gap_linking_str.num_voxel = numel(gap_linking_str.link_ind_w_ep);
% Average intensity of the path found 
gap_linking_str.int = local_int(shortest_path_str.link_ind_w_ep);
gap_linking_str.int_mean = mean(gap_linking_str.int);
gap_linking_str.int_std = std(gap_linking_str.int);
gap_linking_str.int_med = median(gap_linking_str.int);
gap_linking_str.connectivity = connectivity;

% Get the linker outside the mask
out_of_mask_link_ind_Q = ~local_mask(shortest_path_str.link_ind_w_ep);
out_of_mask_link_ind = shortest_path_str.link_ind_w_ep(out_of_mask_link_ind_Q);
gap_linking_str.out_mask_Q = out_of_mask_link_ind_Q;
gap_linking_str.link_sub_o_m = gap_linking_str.link_sub_w_ep(out_of_mask_link_ind_Q,:);
gap_linking_str.link_ind_o_m = gap_linking_str.link_ind_w_ep(out_of_mask_link_ind_Q);

gap_linking_str.int_o_mask = local_int(out_of_mask_link_ind);
gap_linking_str.int_o_mask_mean = mean(gap_linking_str.int_o_mask);
gap_linking_str.int_o_mask_std = std(gap_linking_str.int_o_mask);
gap_linking_str.link_ratio_o_mask = numel(out_of_mask_link_ind) / gap_linking_str.num_voxel;
%% Reconstruct the mask and estimate the signal to noise level more accurately. 
% Estimate the radius 
sample_radius = 5;
sample_tissue_dist_to_linker = 2;
%   Get the radius of the voxel near the two endpoints
skl_ind = find(local_skl);
skl_sub = fun_ind2sub(local_mask_size, skl_ind);
% Skeleton voxel closed to first endpoint & is in the same connected component
ep_2_skl_dist = pdist2(p_1_pos, skl_sub);
is_ep_neighbor_skl_Q = ep_2_skl_dist <= sample_radius;
ep_neighbor_skl_ind = skl_ind(is_ep_neighbor_skl_Q);
ep_neighbor_skl_ind = ep_neighbor_skl_ind(ep1_cc_link_mask(ep_neighbor_skl_ind));
ep_neighbor_radius_med_1 = double(max(1, round(median(local_skl(ep_neighbor_skl_ind)))));
% Skeleton voxel closed to the second endpoint & is in the same connected
% component
ep_2_skl_dist = pdist2(p_2_pos, skl_sub);
is_ep_neighbor_skl_Q = ep_2_skl_dist <= sample_radius;
ep_neighbor_skl_ind = skl_ind(is_ep_neighbor_skl_Q);
ep_neighbor_skl_ind = ep_neighbor_skl_ind(ep2_cc_link_mask(ep_neighbor_skl_ind));
ep_neighbor_radius_med_2 = double(max(1, round(median(local_skl(ep_neighbor_skl_ind)))));

ep_neighbor_radius_med = min(ep_neighbor_radius_med_1, ep_neighbor_radius_med_2);

% Reconstruct the linker mask 
linker_recon_mask = fun_skeleton_reconstruction_single_strel(shortest_path_str.link_ind_w_ep, ...
    local_mask_size, strel('sphere', ep_neighbor_radius_med).Neighborhood);
linker_recon_mask_num_voxel_0 = nnz(linker_recon_mask);
linker_recon_mask = linker_recon_mask & ~local_mask;
linker_recon_mask_num_voxel_1 = nnz(linker_recon_mask);
linker_recon_mask_nonoverlapping_ratio = linker_recon_mask_num_voxel_1 / linker_recon_mask_num_voxel_0;

linker_recon_mask_dilated = imdilate(linker_recon_mask, strel('sphere', sample_tissue_dist_to_linker)) & ~local_mask;
linker_neighbor_background = ~linker_recon_mask & linker_recon_mask_dilated;
% Compute background statistics
neighbor_bg_int = local_int(linker_neighbor_background);
neighbor_bg_int_mean = mean(neighbor_bg_int);
neighbor_bg_int_std = std(neighbor_bg_int);

linker_recon_mask_int = local_int(linker_recon_mask);
linker_recon_int_mean = mean(linker_recon_mask_int);
linker_recon_int_median = median(linker_recon_mask_int);
linker_recon_int_std = std(linker_recon_mask_int);
linker_recon_SNR = (linker_recon_int_mean - neighbor_bg_int_mean)/neighbor_bg_int_std;

linker_recon_ind = find(linker_recon_mask);
linker_sub = fun_ind2sub(local_mask_size, linker_recon_ind);
linker_global_sub = linker_sub + local_bbox_mmxx(1:3) - 1;
linker_global_ind = sub2ind(image_size, linker_global_sub(:,1), linker_global_sub(:,2), ...
    linker_global_sub(:,3));

gap_linking_str.recon_radius = ep_neighbor_radius_med;
gap_linking_str.recon_int_mean = linker_recon_int_mean;
gap_linking_str.recon_int_std = linker_recon_int_std;
gap_linking_str.recon_int_median = linker_recon_int_median;

gap_linking_str.recon_bg_int_mean = neighbor_bg_int_mean;
gap_linking_str.recon_bg_std = neighbor_bg_int_std;
gap_linking_str.recon_SNR = linker_recon_SNR;

gap_linking_str.recon_mask_ind = {linker_global_ind};
gap_linking_str.recon_mask_om_ratio = linker_recon_mask_nonoverlapping_ratio;
%% Distance between the linker skeleton voxle and the nearest existing voxel vs their radius
dist_linker_2_existing_skl = pdist2(shortest_path_str.link_sub_w_ep, skl_sub);
[min_dist, min_skl_idx] = min(dist_linker_2_existing_skl, [], 2);
min_dist_skl_ind = skl_ind(min_skl_idx);
min_dist_skl_sub = fun_ind2sub(local_mask_size, min_dist_skl_ind);
min_dist_skl_sub_global = min_dist_skl_sub + local_bbox_mmxx(1:3) - 1;
min_dist_skl_ind_global = sub2ind(image_size, min_dist_skl_sub_global(:,1), ...
    min_dist_skl_sub_global(:,2), min_dist_skl_sub_global(:,3));

min_dist_skl_dt = local_skl(min_dist_skl_ind);
min_dist_2_skl_r = min_dist ./ min_dist_skl_dt;

dist_linker_2_ep1 = pdist2(shortest_path_str.link_sub_w_ep, shortest_path_str.link_sub_w_ep(1, :));
dist_linker_2_ep2 = pdist2(shortest_path_str.link_sub_w_ep, shortest_path_str.link_sub_w_ep(end, :));
dist_linker_2_nearest_ep = min(dist_linker_2_ep1, dist_linker_2_ep2);

min_dist_2_skl_vs_min_dist_2_ep = min_dist ./ dist_linker_2_nearest_ep;

gap_linking_str.closest_skl_ind = {min_dist_skl_ind_global};
gap_linking_str.closest_skl_dt_2_r = min_dist_2_skl_r;
gap_linking_str.closest_skl_dist_2_closest_ep_dist = min_dist_2_skl_vs_min_dist_2_ep;
end

function sub = fun_ind2sub(block_size, ind)
if numel(block_size) == 2
    sub = zeros(numel(ind), 2);
    [sub(:,1), sub(:,2)] = ind2sub(block_size, ind);
elseif numel(block_size) == 3
    sub = zeros(numel(ind), 3);
    [sub(:,1), sub(:,2), sub(:,3)] = ind2sub(block_size, ind);
else
    error('Unsupported block size');
end
end