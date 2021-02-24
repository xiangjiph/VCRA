function gap_linking_str = fun_graph_connect_gap_ep2skl(vessel_image, vessel_mask, vessel_skl, ep_1_sub, ep_all_sub)
% fun_graph_connect_gap_ep2skl search the shortest path betwwen the
% endpoint and the skeleton within the masked found by threshold
% relaxation.
% Input: 
%   vessel_image: 3D numercal array, image
%   vessel_mask: 3D logical array, mask of the vessel
%   vessel_skl: 3D logical array, mask of the skeleton
%   ep_1_sub: 1x3 numerical vector, subscript of the endpoint in vessel mask 
% Output: 
%   gat_linking_str: structure with fields:
%       foundQ: true if a shortest path is found
%       threshold_list: numerical vector, threshold candidate of threshold
%       relaxation
%       bg_mean: average intensity in the background of the local image
%       bg_std: standard deviation of the intensity in the background of
%       the local iamge
%       
% To improve: 
% 1. High false connection rate: ~30%. Multiple reasons for the false
% connection: 
%   a). Endpoint to endpoint connection: One of the endpoints has found a
%   linker, but the linker does not end at the other endpoint. - Maybe need
%   to search for the endpoint first, before searching for the link voxel. 
%   b). Endpoint to endpoint connection: too far away, more than 35 um
%   apart, get connected to the nearby endpoints - Trade off... Increasing
%   the searching range also increase the false positive rate. 
%   c). Isolated endpoints: I don't know where they should be connected or
%   not. 
%   d). Real fake connection from fake link hairs. 
if nargin < 5
    ep_all_sub = [];
end
gap_linking_str = struct;
target_voxel_str = fun_graph_get_nearest_skl_ind(vessel_image, vessel_mask, vessel_skl, ep_1_sub, ep_all_sub);
% Copy the fields 
tmp_fn = fieldnames(target_voxel_str);
for iter_fn = 1 : numel(tmp_fn)
    if ~any(strcmp(tmp_fn{iter_fn} , {'local_int', 'local_mask', 'local_sub', 'ep_local_sub', 'global_sub', 'global_ind'}))
        gap_linking_str.(tmp_fn{iter_fn}) = target_voxel_str.(tmp_fn{iter_fn});
    end
end

if target_voxel_str.foundQ
    image_size = size(vessel_image);
    local_int = target_voxel_str.local_int;
    local_mask = target_voxel_str.local_mask;
    % Construct distance metric
    search_metric = fun_graph_get_shortest_path_distance_metric(local_int, target_voxel_str.connected_mask);
    % Search for the shortest path 
    local_ep_link_str = fun_graph_shortest_path_between_points(search_metric, target_voxel_str.ep_local_sub, target_voxel_str.local_sub, 26);
    % Out of the mask voxels: 
    out_of_mask_link_ind_Q = ~local_mask(local_ep_link_str.link_ind_w_ep);
    
    if isa(local_ep_link_str.link_sub_w_ep, 'cell')
        [~, max_avg_idx] = max(local_int_max - (local_ep_link_str.total_dist./cellfun(@numel, local_ep_link_str.link_ind_w_ep)));
        max_avg_link_local_sub = local_ep_link_str.link_sub_w_ep{max_avg_idx};
        max_avg_link_local_ind = local_ep_link_str.link_ind_w_ep{max_avg_idx};
    else
        max_avg_link_local_sub = local_ep_link_str.link_sub_w_ep;
        max_avg_link_local_ind = local_ep_link_str.link_ind_w_ep;
    end

    max_avg_link_global_sub = round(max_avg_link_local_sub + target_voxel_str.search_bbox_mmxx(1:3) - 1);
    max_avg_link_global_ind = sub2ind(image_size, max_avg_link_global_sub(:,1), ...
        max_avg_link_global_sub(:,2), max_avg_link_global_sub(:,3));
    gap_linking_str.num_voxel = numel(max_avg_link_local_ind);

    gap_linking_str.int = local_int(max_avg_link_local_ind);
    gap_linking_str.int_std = std(gap_linking_str.int);
    gap_linking_str.int_mean = mean(gap_linking_str.int);
    gap_linking_str.int_med = median(gap_linking_str.int);
    gap_linking_str.link_ind_w_ep = max_avg_link_global_ind;
    gap_linking_str.link_sub_w_ep = max_avg_link_global_sub;
    
    
    gap_linking_str.link_ind_o_mask = max_avg_link_global_ind(out_of_mask_link_ind_Q);
    gap_linking_str.link_sub_o_mask = max_avg_link_global_sub(out_of_mask_link_ind_Q,:);
    % Get the linker outside the mask
    out_of_mask_link_ind = max_avg_link_local_ind(out_of_mask_link_ind_Q);
    gap_linking_str.out_mask_Q = out_of_mask_link_ind_Q;
    gap_linking_str.int_o_mask = local_int(max_avg_link_local_ind(out_of_mask_link_ind_Q));
    gap_linking_str.int_o_mask_mean = mean(gap_linking_str.int_o_mask);
    gap_linking_str.int_o_mask_std = std(gap_linking_str.int_o_mask);
    gap_linking_str.link_ratio_o_mask = numel(out_of_mask_link_ind) / gap_linking_str.num_voxel;
end
end