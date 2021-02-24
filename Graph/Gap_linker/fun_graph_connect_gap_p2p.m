function gap_linking_str = fun_graph_connect_gap_p2p(vessel_image, vessel_mask, p_1_pos, p_2_pos, connectivity)
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
if nargin < 4
    connectivity = 26;
    vessel_mask = [];
elseif nargin < 5
    vessel_mask = [];
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
straight_link_mask_r = 3;


% The middle position of two endpoints
eps_mid_sub = round((p_1_pos + p_2_pos)./2);
% Distance between the middle point and the two end point - fix this value
% for now
eps_mid_to_ep_r = [sqrt(sum((p_1_pos - eps_mid_sub).^2)), sqrt(sum((p_2_pos - eps_mid_sub).^2))];
eps_mid_to_ep_r_max = max(eps_mid_to_ep_r);
local_center_bbox_r = round(2 * eps_mid_to_ep_r_max );
ep_in_diff_cc_Q = false;

gap_linking_str = struct('foundQ', false ,'ep_1_sub', [], 'ep_2_sub', [], 'search_bbox_r',[], 'search_bbox_mmxx', [],  ...
    'threshold_list', [], 'bg_mean', [], 'bg_std', [], 'connected_th', [], 'search_full_image_Q', [], ...
    'link_sub_w_ep', [], 'link_ind_w_ep', [], 'num_voxel', [], 'int', [], 'int_mean', [], 'int_std', [],  'int_med',  [], ...
    'connectivity', [], 'out_mask_Q', [], 'int_o_mask', [], 'int_o_mask_mean', [], ...
    'int_o_mask_std', [], 'link_ratio_o_mask', [], 'link_sub_o_m', [], 'link_ind_o_m', []);
gap_linking_str.ep_1_sub = p_1_pos;
gap_linking_str.ep_2_sub = p_2_pos;
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
        disp('Two endpoints are always in the same connected component. Terminate searching');
        return;
    end
    if all(ep_local_sub_2 > 0) && all(ep_local_sub_2 <= local_mask_size)
        ep_local_ind_2 = sub2ind(local_mask_size, ep_local_sub_2(1), ep_local_sub_2(2), ...
            ep_local_sub_2(3));
    else
        disp('Two endpoints are always in the same connected component. Terminate searching');
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
        else
            disp('Two endpoints are always in the same connected component. Terminate searching');
            return;
        end
    else
        ep_in_diff_cc_Q = true;
    end
end
%% Construct the adjacent matrix
%%  Compute local background average intensty and its standard deviation.
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
    test_mask = local_int > test_th;
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
connected_mask = imclose(connected_mask, strel('sphere', 2));
if isempty(connected_mask)
    disp('No connected mask found. Terminate.');
    return;
end
gap_linking_str.connected_th = test_th;
if test_th == th_list(end)
    gap_linking_str.search_full_image_Q = true;
else
    gap_linking_str.search_full_image_Q = false;
end
% Add the connecting straight line between two endpoint
straight_line_mask = fun_graph_get_p2p_cylinder_mask(size(local_mask), ep_local_sub_1, ep_local_sub_2, straight_link_mask_r);
% volumeViewer(straight_line_mask | connected_mask)
% Find the shortest path in the connected mask
% Construct distance metric
search_metric = fun_graph_get_shortest_path_distance_metric(local_int, connected_mask | straight_line_mask);
%%
% Search the shortest path between points with the given distance metric
shortest_path_str = fun_graph_shortest_path_between_points(search_metric, ep_local_sub_1, ep_local_sub_2, connectivity);

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