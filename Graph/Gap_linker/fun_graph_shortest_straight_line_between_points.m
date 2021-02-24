function p2p_straight_linker_str = fun_graph_shortest_straight_line_between_points(p_1_pos, p_2_pos, image_size)

%% Initialize
p2p_straight_linker_str = fun_initialized_structure_array_with_fieldname_list({'foundQ',...
    'ep_1_sub', 'ep_1_ind', 'ep_2_sub', 'ep_2_ind', 'search_bbox_mmxx',...
    'link_sub_w_ep', 'link_ind_w_ep', 'num_voxel'});
p2p_straight_linker_str.foundQ = false;
%% Preprocessing
if isscalar(p_1_pos)
    p_1_pos = fun_ind2sub(mask_size, p_1_pos);
end
if isscalar(p_2_pos)
    p_2_pos = fun_ind2sub(mask_size, p_2_pos);
end
if iscolumn(p_1_pos)
    p_1_pos = p_1_pos';
end
if iscolumn(p_2_pos)
    p_2_pos = p_2_pos';
end
p2p_straight_linker_str.ep_1_sub = p_1_pos;
p2p_straight_linker_str.ep_2_sub = p_2_pos;
p2p_straight_linker_str.ep_1_ind = sub2ind(image_size, p_1_pos(1), p_1_pos(2), p_1_pos(3));
p2p_straight_linker_str.ep_2_ind = sub2ind(image_size, p_2_pos(1), p_2_pos(2), p_2_pos(3));
%% Find the path
% Find the bounding box of the two points 
bbox_mm = min(p_1_pos, p_2_pos) - 1;
bbox_xx = max(p_1_pos, p_2_pos) + 1;

bbox_size = bbox_xx - bbox_mm + 1;

local_sub_1 = p_1_pos - bbox_mm + 1;
local_sub_2 = p_2_pos - bbox_mm + 1;

% Construct the straight line mask 
straight_line_mask = single(fun_graph_get_p2p_cylinder_mask(bbox_size, ...
    local_sub_1, local_sub_2, 1));
% Set the voxels not in the mask to inf for finding the shortest path 
straight_line_mask(~straight_line_mask) = inf;
% Find the shortest path 
shortest_path_str = fun_graph_shortest_path_between_points(straight_line_mask, ...
    local_sub_1, local_sub_2, 26);
if isempty(shortest_path_str.link_ind_w_ep)
    fprintf('Fail to find the link that connect the two points\n');
    return;
else
    p2p_straight_linker_str.foundQ = true;
end
p2p_straight_linker_str.search_bbox_mmxx = [bbox_mm, bbox_xx];
p2p_straight_linker_str.link_sub_w_ep = shortest_path_str.link_sub_w_ep + bbox_mm - 1;
p2p_straight_linker_str.link_ind_w_ep = sub2ind(image_size, p2p_straight_linker_str.link_sub_w_ep(:, 1), ...
    p2p_straight_linker_str.link_sub_w_ep(:, 2), p2p_straight_linker_str.link_sub_w_ep(:, 3));
p2p_straight_linker_str.num_voxel = numel(p2p_straight_linker_str.link_ind_w_ep);
end