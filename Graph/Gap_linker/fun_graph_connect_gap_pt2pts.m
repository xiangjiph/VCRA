function link_str = fun_graph_connect_gap_pt2pts(vessel_image, ep_1_sub, ep_2_sub, connectivity, vessel_mask)
% fun_graph_connect_gap_ep2ep finds the link that connect the two given
% endpoitns
% Input: 
%   vessel_image: 3D numerical array, image of the vessel 
%   ep_1_sub: 1-by-3 numerical array, subscript of point 1 in vessel_image
%   ep_2_sub: N-by-3 numerical array, subscript of point 2 in vessel_image
%   connectivity: 8, 18 or 26. Number of neighbors for searching the
%   shortest path
%   vessel_mask: 3D logical array, mask of the vessel 
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
image_size = size(vessel_image);
% The middle position of two endpoints
eps_mid_sub = round((ep_1_sub + ep_2_sub)./2);
% Distance between the middle point and the two end point
eps_mid_to_ep_r = [sqrt(sum((ep_1_sub - eps_mid_sub).^2)), sqrt(sum((ep_2_sub - eps_mid_sub).^2))];
local_center_bbox_r = round(2 * max(eps_mid_to_ep_r) );
[local_int, local_bbox_mmxx] = crop_center_box(vessel_image, eps_mid_sub, local_center_bbox_r);
local_int = single(local_int);
% Mask for the valid voxels to search - Can be modified later
valid_mask = local_int > 0;

% Define distance metric
local_int_max = max(local_int(:));
metric = local_int_max - local_int;
metric(~valid_mask) = inf;
% Subscript in the cropped image stack
ep_1_local_sub = ep_1_sub - local_bbox_mmxx(1:3) + 1;
ep_2_local_sub = ep_2_sub - local_bbox_mmxx(1:3) + 1;
% Search the shortest path between points with the given distance metric
shortest_path_str = fun_graph_shortest_path_between_points(metric, ep_1_local_sub, ep_2_local_sub, connectivity );

% Convert to the global image coordinate (with endpoints)
link_str = struct;
link_str.link_sub_w_ep = shortest_path_str.link_sub_w_ep + local_bbox_mmxx(1:3) - 1;
link_str.link_ind_w_ep = sub2ind(image_size, link_str.link_sub_w_ep(:,1), ...
    link_str.link_sub_w_ep(:,2), link_str.link_sub_w_ep(:,3));
% Average intensity of the path found 
link_str.avg_int = local_int_max - shortest_path_str.total_dist / numel(shortest_path_str.link_ind_w_ep);
link_str.link_int = vessel_image(link_str.link_ind_w_ep);
link_str.connectivity = connectivity;
if ~isempty(vessel_mask)
    local_mask = crop_center_box(vessel_mask, eps_mid_sub, local_center_bbox_r);
    link_str.mask_mean = mean(local_int(local_mask));
    link_str.mask_std = std(local_int(local_mask));
    link_str.bg_mean = mean(local_int(~local_mask));
    link_str.bg_std = std(local_int(~local_mask));
end