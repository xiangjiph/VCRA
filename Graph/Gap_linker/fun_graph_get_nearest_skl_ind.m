function target_voxel_str = fun_graph_get_nearest_skl_ind(vessel_image, vessel_mask, vessel_skl, ep_sub, ep_all_sub, return_local_arrayQ)
% fun_graph_get_nearest_skl_ind use threshold relaxation to find the
% closest skeleton voxels that locates in a connected component other than
% the connected component that ep_1_sub belong to. 
% Input: 
%   vessel_image: 3D numerical array
%   vessel_mask: 3D logical array, mask of the vessel
%   vessel_skl: 3D logical array, the skeleton of the vessel
%   ep_1_sub: 1-by-3 vector, subscript of the endpoint
% Output: 
%   target_voxel_str: structure with fields



if nargin < 5
    ep_all_sub = [];
    return_local_arrayQ = false;
elseif nargin < 6
    return_local_arrayQ = false;
end
image_size = double(size(vessel_image));
if isscalar(ep_sub)
    warning('Input endpoint subscript is a scalar. Assumed to be index and convert it to subscript automatically');
    ep_sub = fun_ind2sub(image_size, ep_sub);
end

local_bbox_r = 20;
local_bbox_r_max = 40;
min_cc_size = 250; % If too small, might find connected components with no skeleton voxels inside. 
max_use_nearby_endpoint_dist = 5;

local_cc_to_connect_num = 0;
target_voxel_str = fun_initialized_structure_array_with_fieldname_list({'foundQ', ...
    'ep_global_sub', 'ep_global_ind', 'search_bbox_r', 'search_bbox_mmxx', 'ep_local_sub', ...
   'threshold_list', 'bg_mean', 'bg_std', 'connected_th', 'search_full_image_Q', ...
   'dist_ep2target', 'local_sub', 'target_global_sub', 'target_global_ind', 'local_int', ...
   'local_mask', 'connected_mask' });
target_voxel_str.foundQ = false;
target_voxel_str.ep_global_sub = ep_sub;
target_voxel_str.ep_global_ind = sub2ind(image_size, ep_sub(1), ep_sub(2), ep_sub(3));
%% Search for the local connected components other than the one contains the
% endpoint. 
while ~local_cc_to_connect_num
    
    target_voxel_str.search_bbox_r = local_bbox_r;
    % Crop local image, mask and skeleton
    [local_int, local_bbox_mmxx] = crop_center_box(vessel_image, ep_sub, local_bbox_r);
    target_voxel_str.search_bbox_mmxx = local_bbox_mmxx;
    local_int = single(local_int);
    local_mask = crop_center_box(vessel_mask, ep_sub, local_bbox_r);
    local_skl = crop_center_box(vessel_skl, ep_sub, local_bbox_r);
    % Update local mask
    local_mask = local_mask | local_skl;
    
    % Compute local background average intensty and its standard deviation.
    local_bg_mean = mean(local_int(~local_mask));
    local_bg_std = std(local_int(~local_mask));
    local_mask_size = size(local_mask);
    % Position of the endpoint in the local coordinate
    ep_local_sub = ep_sub - local_bbox_mmxx(1:3) + 1;
    target_voxel_str.ep_local_sub = ep_local_sub;
    ep_local_ind = sub2ind(local_mask_size, ep_local_sub(1), ep_local_sub(2), ...
        ep_local_sub(3));
    % Connected components in the local mask
    local_cc = bwconncomp(local_mask, 26);
    local_mask_label = labelmatrix(local_cc);
    ep_cc_label = local_mask_label(ep_local_ind);
    assert(ep_cc_label>0, 'The skeleton voxel is outside the vessel mask');
    local_cc_size = cellfun(@numel, local_cc.PixelIdxList);
    % Remoce small connecterd component
    local_cc_selected = local_cc_size > min_cc_size;
    % Get connected component other than the one that contaisn the endpoint
    local_cc_selected(ep_cc_label) = false;
    local_cc_to_connect = local_cc.PixelIdxList(local_cc_selected);
    local_cc_to_connect_num = numel(local_cc_to_connect);
    if local_cc_to_connect_num ~= 0
        % Check if the cc contains the skeleton 
        cc_contain_skl_Q = any(local_skl( cat(1, local_cc_to_connect{:})));
    else
        cc_contain_skl_Q = false;
    end
    if local_cc_to_connect_num == 0 || ~cc_contain_skl_Q
        if local_bbox_r <= local_bbox_r_max 
            local_bbox_r = local_bbox_r + 5;
            disp('Increase the size of the local mask.');
        else
            disp('No nearby connected components.');
            return;
        end
    end
end
local_search_mask = (local_mask_label > 0) & (local_mask_label ~= ep_cc_label);
%% Threshold relaxation
% Estimate the compute threshold list
cc_ep_int_mean = mean(local_int(local_mask_label == ep_cc_label));
num_step = round((cc_ep_int_mean - local_bg_mean)/local_bg_std);
th_list = local_bg_mean + (num_step : -1 : -3) * local_bg_std;
target_voxel_str.threshold_list = th_list;
target_voxel_str.bg_mean = local_bg_mean;
target_voxel_str.bg_std=  local_bg_std;
% Iteratively reduce the threshold until the the input endpoint is in the
% same connected component with a skeleton voxel that is originally in the
% other connected component
connected_cc_label = [];
for iter_th = 1 : length(th_list)
    test_th = th_list(iter_th);
    test_mask = local_int > test_th;
    if ~test_mask(ep_local_ind)
        % Endpoint not in the mask
        continue;
    end
    test_mask_labeled = bwlabeln(test_mask, 26);
    test_ep_label = test_mask_labeled(ep_local_ind);
    % Test if any of the other connected component in the original mask is
    % connected
    for iter_cc = 1 : local_cc_to_connect_num
        if any(test_mask_labeled(local_cc_to_connect{iter_cc}) == test_ep_label)
            % Connected. 
            connected_cc_label = iter_cc;
            connected_mask = (test_mask_labeled == test_ep_label);
            break;
        end
    end    
    if ~isempty(connected_cc_label)
        local_skl_masked = connected_mask & local_search_mask & local_skl;
        local_skl_sub = fun_ind2sub(local_mask_size, find(local_skl_masked));
        if ~isempty(local_skl_sub)
            target_voxel_str.foundQ = true;
            break;
        end
    end
end
% Note that this threshold can be further refined by binary search between
% the found threshold and the last unconnected threhsold
if isempty(connected_cc_label) || isempty(local_skl_sub)
   disp('No label found'); 
   return;    
end
target_voxel_str.connected_th = test_th;
if test_th == th_list(end)
    target_voxel_str.search_full_image_Q = true;
else
    target_voxel_str.search_full_image_Q = false;
end
%% Find the nearest skeleton voxel

ep_skel_dist = pdist2(ep_local_sub, local_skl_sub);
[target_voxel_str.dist_ep2target, min_dist_skl_idx] = min(ep_skel_dist);
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
        closest_ep_label = test_mask_labeled(ep_all_local_sub(tmp_idx,1), ep_all_local_sub(tmp_idx,2), ...
            ep_all_local_sub(tmp_idx,3));
        if (min_dist_skl_to_eps <= max_use_nearby_endpoint_dist) && (connected_cc_label == closest_ep_label)
            if min_dist_skl_to_eps > eps
                fprintf('Change the closest point to the endpoint\n');
            end
            % Then endpoint is very close to the found skeleton voxle and
            % they belong to the same connected components
            target_voxel_sub = ep_all_local_sub(tmp_idx, :);            
        end
    end    
end

target_voxel_str.local_sub = target_voxel_sub;
target_voxel_str.target_global_sub = target_voxel_sub + local_bbox_mmxx(1:3) - 1;
target_voxel_str.target_global_ind = sub2ind(image_size, target_voxel_str.target_global_sub(1), target_voxel_str.target_global_sub(2), ...
    target_voxel_str.target_global_sub(3));

if return_local_arrayQ
    target_voxel_str.local_int = local_int;
    target_voxel_str.local_mask = local_mask;
    target_voxel_str.connected_mask = connected_mask;
end

end