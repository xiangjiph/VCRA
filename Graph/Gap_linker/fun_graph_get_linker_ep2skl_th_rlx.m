function gap_linking_str = fun_graph_get_linker_ep2skl_th_rlx(vessel_image, vessel_mask, vessel_skl, ep_sub, ep_all_sub)
% fun_graph_get_linker_ep2skl use threshold relaxation and shortest path on
% designed metric to find a 3D line segment that connect the inut endpoint
% to a skeleton voxel in the nearby connected component. 
% Input: 
%   vessel_image: 3D numerical array
%   vessel_mask: 3D logical array, mask of the vessel
%   vessel_skl: 3D logical array, the skeleton of the vessel
%   ep_1_sub: 1-by-3 vector, subscript of the endpoint
%   ep_all_sub: N-by-3 numerical array, subscripts of the endpoints in the
%   vessel skeleton. Optional input, can be [].
% Output: 
%   gap_linking_str: structure with fields defined below.
%
% Written by Xiang Ji, UC San Diego
% Date: 02/09/2019
% Modified by Xiang Ji, 03/14/2019
% 1. Search for the path between connected components first, then find the
% voxel where the linker enters the target mask, and search for the closest
% skeleton of the entering voxel. 
% 2. Break skeleton into 
% Modified by Xiang Ji, 03/26/2019
% 1. Fix a bug that target boundary voxel and connected target voxel are
% not in the same connected components and therefore no target skeleton voxel can
% be found. 
if nargin < 5
    ep_all_sub = [];
end
image_size = double(size(vessel_image));
if isscalar(ep_sub)
    warning('Input endpoint subscript is a scalar. Assumed to be index and convert it to subscript automatically');
    ep_sub = fun_ind2sub(image_size, ep_sub);
end

local_bbox_r = 25;
local_bbox_r_max = 40;
max_use_nearby_endpoint_dist = 5;
r_link_cc_dt = 5;
connectivity = 26;
sample_radius = 5;
sample_tissue_dist_to_linker = 2;
imclose_th_mask = false;
% break_skl_cc = true;
persistent str_template;
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
gap_linking_str.ep_1_sub = ep_sub;
gap_linking_str.ep_1_ind = sub2ind(image_size, ep_sub(1), ep_sub(2), ep_sub(3));
%% Search for the local connected components other than the one contains the
% endpoint. 
cc_contain_skl_Q = false;
bbox_expanded_Q = true;
while ~cc_contain_skl_Q
    % Crop local image, mask and skeleton
    [local_int, local_bbox_mmxx] = crop_center_box(vessel_image, ep_sub, local_bbox_r);
    local_int = single(local_int);
    local_mask = crop_center_box(vessel_mask, ep_sub, local_bbox_r);
    local_skl = crop_center_box(vessel_skl, ep_sub, local_bbox_r);
    % Update local mask
    local_mask = local_mask | local_skl;
    local_mask_size = size(local_mask);
    % Position of the endpoint in the local coordinate
    ep_local_sub = ep_sub - local_bbox_mmxx(1:3) + 1;
    ep_local_ind = sub2ind(local_mask_size, ep_local_sub(1), ep_local_sub(2), ...
        ep_local_sub(3));
%     if break_skl_cc
%         local_skl = fun_graph_delete_node_voxel(local_skl);
%     end
%     local_skl_cc = bwconncomp(local_skl);
%     local_skl_labeled = labelmatrix(local_skl_cc);
    
    % Connected components in the local mask
    local_cc = bwconncomp(local_mask, connectivity);
    local_mask_label = labelmatrix(local_cc);
%     ep_skl_label = local_skl_labeled(ep_local_ind);
    ep_cc_label = local_mask_label(ep_local_ind);
    
    ep_cc_mask = (local_mask_label == ep_cc_label);
    
%     ep_cc_link_mask = false(local_mask_size);
%     ep_cc_link_mask(local_skl_cc.PixelIdxList{ep_skl_label}) = true;
    ep_cc_link_mask = ep_cc_mask & local_skl;
    dt_search_mask_ep = fun_graph_get_endpoint_link_mask(ep_local_ind, ep_cc_link_mask, r_link_cc_dt);
    
    assert(ep_cc_label>0, 'The skeleton voxel is outside the vessel mask');
    local_cc_selected = false(local_cc.NumObjects, 1);
    for iter_cc = 1 : local_cc.NumObjects
        if iter_cc ~= ep_cc_label
            tmp_ind = local_cc.PixelIdxList{iter_cc};
            tmp_ind_in_dt_mask = tmp_ind(dt_search_mask_ep(tmp_ind));
            if ~isempty(tmp_ind_in_dt_mask)
                if any(local_skl(tmp_ind_in_dt_mask))
                    % The cc in the dt_mask contains skeleton voxels
                    cc_contain_skl_Q = true;
                    local_cc_selected(iter_cc) = true;
                    local_cc.PixelIdxList{iter_cc} = tmp_ind_in_dt_mask;
                end
            end
        end
    end
    
    if ~cc_contain_skl_Q
        if (local_bbox_r <= local_bbox_r_max)
            local_bbox_r = local_bbox_r + 5;
%             disp('Increase the size of the local mask.');
        else
%             fprintf('No nearby connected components.\n');
            return;
        end
    elseif bbox_expanded_Q
        % Convert to false, expand the bounding box to include the nearest
        % skeleton. The estimation is based on the distance transform of
        % the skeleton voxel
        cc_contain_skl_Q = false;
        local_bbox_r = local_bbox_r + double(ceil(max(local_skl(:))));
        bbox_expanded_Q = false;
    else
        local_cc_selected_label = find(local_cc_selected);
        valid_cc_mask_label = zeros(local_mask_size, 'like', local_mask_label);
        for iter_cc = local_cc_selected_label'
            valid_cc_mask_label(local_cc.PixelIdxList{iter_cc}) = iter_cc;
        end
    end
end
% fun_vis_link_surrounding_by_cc_ind(local_int, logical(ep_cc_mask), local_skl, linker_bv2bv, 60)

gap_linking_str.search_bbox_r = local_bbox_r;
gap_linking_str.search_bbox_mmxx = local_bbox_mmxx;
% Compute local background average intensty and its standard deviation.
local_bg_mean = mean(local_int(~local_mask));
local_bg_std = std(local_int(~local_mask));
%% Threshold relaxation
% Estimate the compute threshold list
cc_ep_int_mean = mean(local_int(local_mask_label == ep_cc_label));
num_step = round((cc_ep_int_mean - local_bg_mean)/local_bg_std);
th_list = local_bg_mean + (num_step : -1 : 0) * local_bg_std;
gap_linking_str.threshold_list = th_list;
gap_linking_str.bg_mean = local_bg_mean;
gap_linking_str.bg_std=  local_bg_std;
% Iteratively reduce the threshold until the the input endpoint is in the
% same connected component with a skeleton voxel that is originally in the
% other connected component
connected_cc_label = [];
for iter_th = 1 : length(th_list)
    test_th = th_list(iter_th);
    test_mask = (local_int > test_th) & dt_search_mask_ep;
    if imclose_th_mask
        test_mask = imclose(test_mask, strel('cube', 2));
    end
    if ~test_mask(ep_local_ind)
        % Endpoint not in the mask
        continue;
    end
    % Tiny image close? - Will it cause problems? 
%     test_mask = imclose(test_mask, strel('cube', 3));    
    test_mask_labeled = bwlabeln(test_mask, connectivity);
    test_ep_label = test_mask_labeled(ep_local_ind);
    % Test if any of the other connected component in the original mask is
    % connected
    for iter_cc = local_cc_selected_label'
        if any(test_mask_labeled(local_cc.PixelIdxList{iter_cc}) == test_ep_label)
            % Connected. 
            connected_cc_label = iter_cc;
            connected_mask = (test_mask_labeled == test_ep_label);
            % Notice that threshold relaxation can merge more than 1
            % connected component in the original mask. Since we break the
            % loop here, then getting the target_cc_mask, we cannot just
            % use the cc label for masking. 
            break;
        end
    end
    if ~isempty(connected_cc_label)
        skl_in_target_cc = connected_mask & local_skl & ~ep_cc_link_mask;
        if any(skl_in_target_cc(:))
            gap_linking_str.foundQ = true;
            target_cc_label = unique(valid_cc_mask_label(skl_in_target_cc));
            break;
        end
    end
end
% Note that this threshold can be further refined by binary search between
% the found threshold and the last unconnected threhsold
if isempty(connected_cc_label) || ~gap_linking_str.foundQ
%    disp('No label found'); 
   return;    
end
gap_linking_str.connected_th = test_th;
if test_th <= local_bg_mean
    gap_linking_str.search_full_image_Q = true;
else
    gap_linking_str.search_full_image_Q = false;
end
%% Find voxel triples for connection: 
% 1. Get all the link voxels in the mask that the endpoint is in 
% 2. Get the position of all the boundary voxel of the mask that the
% endpoint is in 
% 3. Use distance transform to select the mask boundary voxels that are
% closer to the endpoint than to the other link voxels in the mask
% 4. Among these boundary voxels, find the one that is closest to the
% boundary voxel of the other connected component. 
% 5. Notice that for very large vessels, the distance between the
% centerline and the boundary vessels can be pretty large. 
% tmp = fun_graph_get_closest_geodesic_voxel_pair(linker_mask_dilate & dt_search_mask_ep, ...
%     ep_cc_boundary_ind, target_cc_boundary_ind);dis
% each connected component in target cc mask should contains the skeleton
% voxels. 
target_cc_mask = false(local_mask_size);
for iter_label = 1 : numel(target_cc_label)
    target_cc_mask(valid_cc_mask_label == target_cc_label(iter_label)) = true;
end
target_cc_mask = (target_cc_mask & connected_mask);
% Voxels in the origianl masks that connected to the connected mask: 
linker_mask_dilate = imdilate(connected_mask & ~ ep_cc_mask & ~valid_cc_mask_label, strel('cube', 3));
% Boundary voxels of the mask that the endpoint is in:
ep_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask_ep & ep_cc_mask;
target_cc_mask_boundary_valid = linker_mask_dilate & dt_search_mask_ep & target_cc_mask;
% Find the nearest boundary voxel pair: 
ep_cc_boundary_ind = find(ep_cc_mask_boundary_valid);
ep_cc_boundary_sub = fun_ind2sub(local_mask_size, ep_cc_boundary_ind);

target_cc_boundary_ind = find(target_cc_mask_boundary_valid);
target_cc_boundary_sub  = fun_ind2sub(local_mask_size, target_cc_boundary_ind);
% distance between distance points: 
boundary_voxel_dist = pdist2(ep_cc_boundary_sub, target_cc_boundary_sub);
% Find the minimum distance: 
[~, tmp_min_ind] = min(boundary_voxel_dist(:));
[ep_cc_bp_idx, target_cc_bp_idx] = ind2sub(size(boundary_voxel_dist), tmp_min_ind);

% ep_bv_sub = ep_cc_boundary_sub(ep_cc_bp_idx, :);
ep_bv_ind = ep_cc_boundary_ind(ep_cc_bp_idx);
% tg_bv_sub = target_cc_boundary_sub(target_cc_bp_idx, :);
tg_bv_ind = target_cc_boundary_ind(target_cc_bp_idx);

skl_in_target_cc = target_cc_mask & local_skl & ~ep_cc_link_mask & connected_mask;
skl_in_target_cc_sub = fun_ind2sub(local_mask_size, find(skl_in_target_cc));
if isempty(skl_in_target_cc_sub)
    % This can happen if there are two potential cc to connect to,
    % threshold relaxation connect this two cc 
    gap_linking_str.foundQ = false;
    return;
end
%% Find the shortest path
% Avoid the shortest path visit any of the centerline skeleton except for
% the two endpoints of the linker
target_cc_link_mask = target_cc_mask & local_skl;
% 
search_mask = dt_search_mask_ep & (ep_cc_mask | target_cc_mask | connected_mask);
search_metric = fun_graph_get_shortest_path_distance_metric(local_int, search_mask);
search_graph = fun_graph_get_search_graph(search_metric, 26);
linker_ep2bv = shortestpath(search_graph, ep_local_ind, ep_bv_ind, 'Method', 'positive');
assert(~isempty(linker_ep2bv), 'Cannot find path from the endpoint to the mask boundary');
linker_bv2bv = shortestpath(search_graph, ep_bv_ind, tg_bv_ind, 'Method', 'positive');
assert(~isempty(linker_bv2bv), 'Cannot find path connecting two boundary voxels in different connected components');
% Find the first voxel in the target mask
linker_bv2bv_is_in_target_mask_Q = target_cc_mask(linker_bv2bv);
bv_idx_first_in_target_mask = find(linker_bv2bv_is_in_target_mask_Q, 1);
linker_bv2bv = linker_bv2bv(1:bv_idx_first_in_target_mask);
tg_bv_ind_new = linker_bv2bv(bv_idx_first_in_target_mask);
tg_bv_sub_new = fun_ind2sub(local_mask_size, tg_bv_ind_new);
%% Find the nearest skeleton voxel
tg_bv_skl_dist = pdist2(tg_bv_sub_new, skl_in_target_cc_sub);
[~, min_dist_skl_idx] = min(tg_bv_skl_dist);
min_dist_skl_local_sub = skl_in_target_cc_sub(min_dist_skl_idx, :);
target_voxel_sub = min_dist_skl_local_sub;
% Check if nearby endpoint is in the ep_all_sub list
% Local position of the endpoints
if ~isempty(ep_all_sub)
    ep_all_local_sub_Q = all(bsxfun(@ge, ep_all_sub, local_bbox_mmxx(1:3)), 2) &...
        all(bsxfun(@le, ep_all_sub, local_bbox_mmxx(4:6)),2);
    ep_all_local_sub = bsxfun(@minus, ep_all_sub, local_bbox_mmxx(1:3) - 1);
    ep_all_local_sub = ep_all_local_sub(ep_all_local_sub_Q,:); 
    if ~isempty(ep_all_local_sub)
        % Distance between the endpoints in the local bounding box and the skeleton
        % voxels that closest to the endpoint
        dist_skl_to_eps = pdist2(min_dist_skl_local_sub, ep_all_local_sub);
        [min_dist_skl_to_eps, tmp_idx] = min(dist_skl_to_eps);
        closest_ep_label = test_mask_labeled(ep_all_local_sub(tmp_idx,1), ep_all_local_sub(tmp_idx,2), ...
            ep_all_local_sub(tmp_idx,3));
        if (min_dist_skl_to_eps <= max_use_nearby_endpoint_dist) && (connected_cc_label == closest_ep_label)
            tmp_new_sub = ep_all_local_sub(tmp_idx, :);
            if search_mask(tmp_new_sub(1), tmp_new_sub(2), tmp_new_sub(3))
                % Then endpoint is very close to the found skeleton voxle and
                % they belong to the same connected components, and the new
                % node is in the search mask 
                target_voxel_sub = ep_all_local_sub(tmp_idx, :);            
            end
        end
    end    
end
target_voxel_ind = sub2ind(local_mask_size, target_voxel_sub(1), target_voxel_sub(2), target_voxel_sub(3));
linker_bv2ep = shortestpath(search_graph, tg_bv_ind_new, target_voxel_ind, 'Method', 'positive');
assert(~isempty(linker_bv2ep), 'Cannot find path connecting the boundary voxel and the target skeleton voxel');
shortest_path_str.link_ind_w_ep = cat(2, linker_ep2bv, linker_bv2bv(2:end-1), linker_bv2ep)';
shortest_path_str.link_sub_w_ep = fun_ind2sub(local_mask_size, shortest_path_str.link_ind_w_ep);
shortest_path_str.num_target = 1;
shortest_path_str.connectivity = connectivity;
shortest_path_str.total_dist = sum(search_metric(shortest_path_str.link_ind_w_ep));
%% Save information
gap_linking_str.ep_2_sub = target_voxel_sub + local_bbox_mmxx(1:3) - 1;
gap_linking_str.ep_2_ind = sub2ind(image_size, gap_linking_str.ep_2_sub(1), ...
    gap_linking_str.ep_2_sub(2), gap_linking_str.ep_2_sub(3));

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
%   Get the radius of the voxel near the two endpoints
skl_ind = find(local_skl);
skl_sub = fun_ind2sub(local_mask_size, skl_ind);
% Skeleton voxel closed to first endpoint & is in the same connected component
ep_2_skl_dist = pdist2(ep_local_sub, skl_sub);
is_ep_neighbor_skl_Q = ep_2_skl_dist <= sample_radius;
ep_neighbor_skl_ind = skl_ind(is_ep_neighbor_skl_Q);
ep_neighbor_skl_ind = ep_neighbor_skl_ind(ep_cc_link_mask(ep_neighbor_skl_ind));
ep_neighbor_radius_med_1 = double(max(1, round(median(local_skl(ep_neighbor_skl_ind)))));
% Skeleton voxel closed to the second endpoint & is in the same connected
% component
ep_2_skl_dist = pdist2(target_voxel_sub, skl_sub);
is_ep_neighbor_skl_Q = ep_2_skl_dist <= sample_radius;
ep_neighbor_skl_ind = skl_ind(is_ep_neighbor_skl_Q);
ep_neighbor_skl_ind = ep_neighbor_skl_ind(target_cc_link_mask(ep_neighbor_skl_ind));
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

% Save information
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
%% Distance between the linker skeleton voxel and the nearest existing voxel vs their radius
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
%% For visualization and debug
% fun_vis_link_surrounding_by_cc_ind(local_int, local_mask, local_skl, linker_bv2bv, 60)
end