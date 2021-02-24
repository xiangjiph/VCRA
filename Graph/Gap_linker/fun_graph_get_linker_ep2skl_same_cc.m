function gap_linking_str = fun_graph_get_linker_ep2skl_same_cc(vessel_image, vessel_mask, vessel_skl, ep_sub, ep_all_sub)
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

local_bbox_r = 20;
local_bbox_r_max = 50;
r_link_cc_dt = 5;
connectivity = 26;
max_use_nearby_endpoint_dist = 5;
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
target_skl_voxel_found_Q = false;
while ~target_skl_voxel_found_Q
    % Crop local image, mask and skeleton
    [local_int, local_bbox_mmxx] = crop_center_box(vessel_image, ep_sub, local_bbox_r);
    local_mask_size = local_bbox_mmxx(4:6) - local_bbox_mmxx(1:3) + 1;
    ep_local_sub = ep_sub - local_bbox_mmxx(1:3) + 1;
    ep_local_ind = sub2ind(local_mask_size, ep_local_sub(1), ep_local_sub(2), ...
        ep_local_sub(3));
    
    local_int = single(local_int);

    local_skl = crop_center_box(vessel_skl, ep_sub, local_bbox_r);
    local_skl_cc = bwconncomp(local_skl);
    if local_skl_cc.NumObjects >= 2
        local_mask = crop_center_box(vessel_mask, ep_sub, local_bbox_r);
        local_mask = local_mask | local_skl;
        local_mask_label = bwlabeln(local_mask);
        
        local_skl_label = labelmatrix(local_skl_cc);
        ep_skl_label = local_skl_label(ep_local_ind);
        ep_skl_mask = (local_skl_label == ep_skl_label);
        
        ep_mask_label = local_mask_label(ep_local_ind);
        
        is_neighbor_cc_Q = true(local_skl_cc.NumObjects, 1);
        is_neighbor_cc_Q(ep_skl_label) = false;
        neighbor_skl_ind = cat(1, local_skl_cc.PixelIdxList{is_neighbor_cc_Q});
        neighbor_skl_mask_cc = local_mask_label(neighbor_skl_ind);
        is_in_the_same_cc_Q = (neighbor_skl_mask_cc == ep_mask_label);
        
        dt_search_mask_ep = fun_graph_get_endpoint_link_mask(ep_local_ind, ep_skl_mask, r_link_cc_dt);
        is_in_dt_search_mask_Q = dt_search_mask_ep(neighbor_skl_ind);
        
        valid_neighbor_skl_ind = neighbor_skl_ind(is_in_the_same_cc_Q & is_in_dt_search_mask_Q);
        if ~isempty(valid_neighbor_skl_ind)
            % Compute the pairwise distance between the endpoint and the
            % neighboring skeleton voxels
            valid_neighbor_skl_sub = fun_ind2sub(local_mask_size, valid_neighbor_skl_ind);
            dist_ep2skl = pdist2(ep_local_sub, valid_neighbor_skl_sub);
            [~, min_neighbor_idx] = min(dist_ep2skl);
            nearest_skl_ind = valid_neighbor_skl_ind(min_neighbor_idx);
            target_voxel_sub = valid_neighbor_skl_sub(min_neighbor_idx, :);
            target_skl_voxel_found_Q = true;
        end
    end
    
    if ~target_skl_voxel_found_Q 
        if (local_bbox_r <= local_bbox_r_max)
            local_bbox_r = local_bbox_r + 5; 
        else
            return;
        end
    end
end
%% Shift the target skl voxel to the nearest endpoint
if ~isempty(ep_all_sub)
    ep_all_local_sub_Q = all(bsxfun(@ge, ep_all_sub, local_bbox_mmxx(1:3)), 2) &...
        all(bsxfun(@le, ep_all_sub, local_bbox_mmxx(4:6)),2);
    if any(ep_all_local_sub_Q)
        ep_all_local_sub = bsxfun(@minus, ep_all_sub(ep_all_local_sub_Q, :), local_bbox_mmxx(1:3) - 1);
        % Distance between the endpoints in the local bounding box and the skeleton
        % voxels that closest to the endpoint
        dist_skl_to_eps = pdist2(target_voxel_sub, ep_all_local_sub);
        [min_dist_skl_to_eps, tmp_idx] = min(dist_skl_to_eps);
        closest_ep_label = local_mask_label(ep_all_local_sub(tmp_idx,1), ep_all_local_sub(tmp_idx,2), ...
            ep_all_local_sub(tmp_idx,3));
        if (min_dist_skl_to_eps <= max_use_nearby_endpoint_dist) && (ep_mask_label == closest_ep_label)
            % Then endpoint is very close to the found skeleton voxle and
            % they belong to the same connected components
            target_voxel_sub = ep_all_local_sub(tmp_idx, :);
        end
    end    
end
%%
gap_linking_str.search_bbox_r = local_bbox_r;
gap_linking_str.search_bbox_mmxx = local_bbox_mmxx;
%% Threshold relaxation
% Compute local background average intensty and its standard deviation.
local_bg_mean = mean(local_int(~local_mask));
local_bg_std = std(local_int(~local_mask));
% % Estimate the compute threshold list
% cc_ep_int_mean = mean(local_int(local_mask_label == ep_cc_label));
% num_step = round((cc_ep_int_mean - local_bg_mean)/local_bg_std);
% th_list = local_bg_mean + (num_step : -1 : 0) * local_bg_std;
% gap_linking_str.threshold_list = th_list;
gap_linking_str.bg_mean = local_bg_mean;
gap_linking_str.bg_std = local_bg_std;
% % Iteratively reduce the threshold until the the input endpoint is in the
% % same connected component with a skeleton voxel that is originally in the
% % other connected component
% connected_cc_label = [];
% for iter_th = 1 : length(th_list)
%     test_th = th_list(iter_th);
%     test_mask = (local_int > test_th) & dt_search_mask_ep;
%     if ~test_mask(ep_local_ind)
%         % Endpoint not in the mask
%         continue;
%     end
%     % Tiny image close? - Will it cause problems? 
% %     test_mask = imclose(test_mask, strel('cube', 3));    
%     test_mask_labeled = bwlabeln(test_mask, connectivity);
%     test_ep_label = test_mask_labeled(ep_local_ind);
%     % Test if any of the other connected component in the original mask is
%     % connected
%     for iter_cc = local_cc_selected_label'
%         if any(test_mask_labeled(local_cc.PixelIdxList{iter_cc}) == test_ep_label)
%             % Connected. 
%             connected_cc_label = iter_cc;
%             connected_mask = (test_mask_labeled == test_ep_label);
%             % Notice that threshold relaxation can merge more than 1
%             % connected component in the original mask. Since we break the
%             % loop here, then getting the target_cc_mask, we cannot just
%             % use the cc label for masking. 
%             break;
%         end
%     end
%     if ~isempty(connected_cc_label)
%         local_skl_masked = connected_mask & local_skl & ~ep_cc_link_mask;
%         if any(local_skl_masked(:))
%             gap_linking_str.foundQ = true;
%             break;
%         end
%     end
% end
% % Note that this threshold can be further refined by binary search between
% % the found threshold and the last unconnected threhsold
% if isempty(connected_cc_label) || ~gap_linking_str.foundQ
% %    disp('No label found'); 
%    return;    
% end
% gap_linking_str.connected_th = test_th;
% if test_th <= local_bg_mean
%     gap_linking_str.search_full_image_Q = true;
% else
%     gap_linking_str.search_full_image_Q = false;
% end
%% Find the shortest path
% Avoid the shortest path visit any of the centerline skeleton except for
% the two endpoints of the linker
target_skl_mask = local_skl;
target_skl_mask([nearest_skl_ind, ep_local_ind]) = false;
% 
search_mask = dt_search_mask_ep & ~target_skl_mask & (local_mask_label == ep_skl_label);

% dist_metric = ones(local_mask_size);
search_metric = fun_graph_get_shortest_path_distance_metric(local_int, search_mask);
shortest_path_str = fun_graph_shortest_path_between_points(search_metric, ep_local_sub, target_voxel_sub, connectivity);
% fun_vis_link_surrounding_by_cc_ind(local_int, local_mask, local_skl, shortest_path_str.link_ind_w_ep, 60)
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