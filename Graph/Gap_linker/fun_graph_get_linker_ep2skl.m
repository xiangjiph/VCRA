function gap_linking_str = fun_graph_get_linker_ep2skl(vessel_image, vessel_mask, vessel_skl, ep_sub, ep_all_sub)
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
if nargin < 5
    ep_all_sub = [];
end
image_size = double(size(vessel_image));
if isscalar(ep_sub)
    warning('Input endpoint subscript is a scalar. Assumed to be index and convert it to subscript automatically');
    ep_sub = fun_ind2sub(image_size, ep_sub);
end

local_bbox_r = 30;
local_bbox_r_max = 50;
max_use_nearby_endpoint_dist = 5;
r_link_cc_dt = 5;
connectivity = 26;

gap_linking_str = fun_initialized_structure_array_with_fieldname_list({'foundQ',...
    'ep_1_sub', 'ep_1_ind', 'ep_2_sub', 'ep_2_ind', 'search_bbox_r', 'search_bbox_mmxx',...
    'threshold_list', 'bg_mean', 'bg_std', 'connected_th', ...
    'search_full_image_Q', 'link_sub_w_ep', 'link_ind_w_ep', 'num_voxel', ...
    'int', 'int_mean', 'int_std', 'int_med', 'connectivity', 'out_mask_Q', ...
    'int_o_mask', 'int_o_mask_mean', 'int_o_mask_std', 'link_ratio_o_mask', ...
    'link_sub_o_m', 'link_ind_o_m' });
gap_linking_str.foundQ = false;
gap_linking_str.ep_1_sub = ep_sub;
gap_linking_str.ep_1_ind = sub2ind(image_size, ep_sub(1), ep_sub(2), ep_sub(3));
%% Search for the local connected components other than the one contains the
% endpoint. 
cc_contain_skl_Q = false;
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
    % Connected components in the local mask
    local_cc = bwconncomp(local_mask, connectivity);
    local_mask_label = labelmatrix(local_cc);
    ep_cc_label = local_mask_label(ep_local_ind);
    
    ep_cc_mask = (local_mask_label == ep_cc_label);
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
        if local_bbox_r <= local_bbox_r_max 
            local_bbox_r = local_bbox_r + 5;
            disp('Increase the size of the local mask.');
        else
            fprintf('No nearby connected components.\n');
            return;
        end
    else
        local_cc_selected_label = find(local_cc_selected);
        valid_cc_mask_label = zeros(local_mask_size, 'like', local_mask_label);
        for iter_cc = local_cc_selected_label'
            valid_cc_mask_label(local_cc.PixelIdxList{iter_cc}) = iter_cc;
        end
    end
end
gap_linking_str.foundQ = true;
gap_linking_str.search_bbox_r = local_bbox_r;
gap_linking_str.search_bbox_mmxx = local_bbox_mmxx;
gap_linking_str.bg_mean = mean(local_int(~local_mask));
gap_linking_str.bg_std =  std(local_int(~local_mask));
%% Find the nearest endpoint directly 
% boundary voxels of the valid cc: 
valid_cc_mask = valid_cc_mask_label > 0;
valid_cc_boundary = (valid_cc_mask) & (~imerode(valid_cc_mask, strel('cube', 3)));
valid_cc_ind = find(valid_cc_boundary);
valid_cc_sub = fun_ind2sub(local_mask_size, valid_cc_ind);

ep_cc_boundary = (ep_cc_mask) & (~imerode(ep_cc_mask, strel('cube', 3)));
ep_cc_boundary = ep_cc_boundary & dt_search_mask_ep;
ep_cc_ind = find(ep_cc_boundary);
ep_cc_sub = fun_ind2sub(local_mask_size, ep_cc_ind);

% distance between distance points: 
boundary_voxel_dist = pdist2(ep_cc_sub, valid_cc_sub);
% Find the minimum distance: 
[~, tmp_min_ind] = min(boundary_voxel_dist(:));
[ep_cc_bp_idx, target_cc_bp_idx] = ind2sub(size(boundary_voxel_dist), tmp_min_ind);

ep_bv_sub = ep_cc_sub(ep_cc_bp_idx, :);
tg_bv_sub = valid_cc_sub(target_cc_bp_idx, :);
tg_bv_ind = valid_cc_ind(target_cc_bp_idx);
tg_cc_label = valid_cc_mask_label(tg_bv_ind);
tg_cc_mask = (valid_cc_mask_label == tg_cc_label);
%% Find the nearest skeleton voxel
local_skl_masked = tg_cc_mask & local_skl;
local_skl_sub = fun_ind2sub(local_mask_size, find(local_skl_masked));
tg_bv_skl_dist = pdist2(tg_bv_sub, local_skl_sub);
[~, min_dist_skl_idx] = min(tg_bv_skl_dist);
min_dist_skl_local_sub = local_skl_sub(min_dist_skl_idx, :);
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
        closest_ep_label = valid_cc_mask_label(ep_all_local_sub(tmp_idx,1), ep_all_local_sub(tmp_idx,2), ...
            ep_all_local_sub(tmp_idx,3));
        if (min_dist_skl_to_eps <= max_use_nearby_endpoint_dist) && (tg_cc_label == closest_ep_label) && ...
                ep_all_local_sub(tmp_idx,1) ~= ep_local_sub(1) && ...
                ep_all_local_sub(tmp_idx,2) ~= ep_local_sub(2) && ...
                ep_all_local_sub(tmp_idx,3) ~= ep_local_sub(3)
            if min_dist_skl_to_eps > eps
                fprintf('Change the closest point to the endpoint\n');
            end
            % Then endpoint is very close to the found skeleton voxle and
            % they belong to the same connected components
            target_voxel_sub = ep_all_local_sub(tmp_idx, :);            
        end
    end    
end
%% Find the shortest path
% Avoid the shortest path visit any of the centerline skeleton except for
% the two endpoints of the linker
ep2ep_dist = sqrt(sum((ep_local_sub - target_voxel_sub).^2));
target_skl_mask = tg_cc_mask & local_skl;
target_skl_mask(target_voxel_sub(1), target_voxel_sub(2), target_voxel_sub(3)) = false;
connected_mask = fun_graph_get_p2p_cylinder_mask(local_mask_size, ep_local_sub, target_voxel_sub, ep2ep_dist);
search_mask = dt_search_mask_ep & ~target_skl_mask & (ep_cc_mask | tg_cc_mask | connected_mask);
search_metric = fun_graph_get_shortest_path_distance_metric(local_int, search_mask);
shortest_path_str = fun_graph_shortest_path_through_4pts(search_metric, ep_local_sub, ep_bv_sub, tg_bv_sub, target_voxel_sub, connectivity);
fun_vis_link_surrounding_by_cc_ind(local_int, local_mask, local_skl, shortest_path_str.link_ind_w_ep)
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


end