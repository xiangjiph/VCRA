function linker_features = fun_graph_get_linker_features(linker_str_array, vessel_graph)
% fun_graph_get_linker_features computes features for fake linker
% classifier from linker_str
% Input: 
%   linker_str_array: structure array computed by fun_graph_connect_gap_p2p
% Output: 
%   linker_features: linker features
% 
% Author: Xiang Ji
% Date: 01/21/2019
linker_features = struct;
% Original features
feature_name = {'int_mean', 'int_std', 'int_o_mask_mean', 'int_o_mask_std', ...
    'int_med', 'recon_int_std', 'recon_int_mean', 'recon_int_median', ...
    'recon_bg_int_mean', 'recon_bg_std', 'recon_mask_om_ratio', 'recon_SNR',...
    'recon_radius', 'num_voxel', 'bg_mean', 'bg_std', 'connected_th', 'link_ratio_o_mask'};
for iter_feature = 1 : numel(feature_name)
    linker_features.(feature_name{iter_feature}) = [linker_str_array.(feature_name{iter_feature})]';
end
% Normalize intensity by maximum possible value: 
% int_max = intmax('uint16');
% linker_features.int_mean = linker_features.int_mean./int_max;
% linker_features.int_med = linker_features.int_med ./ int_max;
%% Signal to noise
% skeleton
linker_features.connected_th_SNR = (linker_features.connected_th - linker_features.bg_mean) ./ (linker_features.bg_std);
linker_features.linker_SNR = ( linker_features.int_mean - linker_features.bg_mean )./ linker_features.bg_std;
linker_features.linker_o_m_SNR = (linker_features.int_o_mask_mean - linker_features.bg_mean)./ linker_features.bg_std;
%% Linker geometry
% Length
linker_features.length = cellfun(@(x) fun_graph_sub_to_length(x), {linker_str_array.link_sub_w_ep})';
% Position and end to end distance
sub_1 = cat(1, linker_str_array.ep_1_sub);
sub_2 = cat(1, linker_str_array.ep_2_sub);
linker_features.ep2ep_dist = sqrt(sum((sub_1 - sub_2).^2, 2));
% Ratio between the end-to-end distance and the linker length
linker_features.ep2ep_dist_2_length = linker_features.ep2ep_dist ./ linker_features.length;
linker_features.dir_vec_1_to_2 = sub_2 - sub_1;
linker_features.dir_vec_1_to_2 = linker_features.dir_vec_1_to_2 ./ vecnorm(linker_features.dir_vec_1_to_2, 2, 2);
%% Closest voxel statistics
linker_features.closest_skl_dist_2_closest_ep_dist_mean = cellfun(@(x) mean(x, 'omitnan'), {linker_str_array.closest_skl_dist_2_closest_ep_dist});
linker_features.closest_skl_dist_2_closest_ep_dist_std = cellfun(@(x) std(x, 'omitnan'), {linker_str_array.closest_skl_dist_2_closest_ep_dist});
%% Include links 
num_linker = numel(linker_str_array);

linker_features.closest_skl_dist_2_closest_ep_dist_mean = nan(num_linker, 1);
linker_features.closest_skl_dist_2_closest_ep_dist_std = nan(num_linker, 1);
linker_features.num_nearest_link_cc = nan(num_linker, 1);
linker_features.num_nearest_dt_change_sign = nan(num_linker, 1);
linker_features.num_seg_in_mask = nan(num_linker, 1);
linker_features.out_of_mask_ratio = nan(num_linker, 1);

linker_features.link_ori_num_voxel = nan(num_linker, 1);
linker_features.link_ep_orientation = nan(3, num_linker);
linker_features.linker_orientation = nan(3, num_linker);

pca_max_num_voxel = 10;
image_size = vessel_graph.num.mask_size;
linker_ep1_link_2epQ = false(num_linker, 1);
for iter_linker = 1 : num_linker
    tmp_linker_ind = linker_str_array(iter_linker).link_ind_w_ep;
    tmp_linker_omQ = linker_str_array(iter_linker).out_mask_Q;
    tmp_linker_skl_dist_2_ep_dist = linker_str_array(iter_linker).closest_skl_dist_2_closest_ep_dist;
    tmp_linker_skl_dist_2_radius = linker_str_array(iter_linker).closest_skl_dt_2_r;
    tmp_linker_skl_dist_2_radius_l_1 = tmp_linker_skl_dist_2_radius < 1;
    % Number of times when the linker is within the reconstruction mask of
    % its nearby skeleton voxels. 
    linker_features.num_seg_in_mask(iter_linker) = sum(abs(diff(tmp_linker_skl_dist_2_radius_l_1)));
    % Number of times when the distance to the closest voxels change sign. 
    linker_features.num_nearest_dt_change_sign(iter_linker) = nnz(diff(diff(tmp_linker_skl_dist_2_radius)>0));
    
    tmp_linker_ep1_link_label = full(vessel_graph.link.map_ind_2_label(linker_str_array(iter_linker).ep_1_ind));
    linker_ep1_link_2epQ(iter_linker) = ~any(vessel_graph.link.connected_node_label(tmp_linker_ep1_link_label, :));
    tmp_linker_ep1_link_ind = vessel_graph.link.cc_ind{tmp_linker_ep1_link_label};
    tmp_linker_ep1_link_num_voxel = numel(tmp_linker_ep1_link_ind);
    linker_features.link_ori_num_voxel(iter_linker) = tmp_linker_ep1_link_num_voxel;
    linker_features.closest_skl_dist_2_closest_ep_dist_mean(iter_linker) = mean(tmp_linker_skl_dist_2_ep_dist, 'omitnan');
    linker_features.closest_skl_dist_2_closest_ep_dist_std(iter_linker) = std(tmp_linker_skl_dist_2_ep_dist, 'omitnan');
    
    tmp_linker_nearest_skl_ind = linker_str_array(iter_linker).closest_skl_ind{:};
    linker_features.num_nearest_link_cc(iter_linker) = numel(unique(full(vessel_graph.link.map_ind_2_label(tmp_linker_nearest_skl_ind))));
    % Orientation of the link
    % Make sure the endpoint is at the end of the ind list
    if tmp_linker_ep1_link_ind(1) == tmp_linker_ind(1)
        tmp_linker_ep1_link_ind = flip(tmp_linker_ep1_link_ind);
    end    
    tmp_linker_ep1_link_sub = fun_ind2sub(image_size, tmp_linker_ep1_link_ind);
    if size(tmp_linker_ep1_link_sub,1) > 2
        
        tmp_1_sub = tmp_linker_ep1_link_sub(max(1, (end-min(pca_max_num_voxel, tmp_linker_ep1_link_num_voxel) + 1)):end,:);
        [ep_1_vec] = pca(tmp_linker_ep1_link_sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
        % make sure the vector is point along the direction of the vessel: 
        ep_voxel_vec = tmp_1_sub(end, :) - tmp_1_sub(1, :);
        tmp_inner_produce = ep_voxel_vec * ep_1_vec;
        if tmp_inner_produce < 0
            linker_features.link_ep_orientation(:, iter_linker) = - ep_1_vec;
        else
            linker_features.link_ep_orientation(:, iter_linker) = ep_1_vec;
        end
    end
    % Direction from the first endpoint to the point before entering the
    % target cc mask
    tmp_last_om_idx = find(tmp_linker_omQ, 1, 'last');
    tmp_last_om_sub = linker_str_array(iter_linker).link_sub_w_ep(tmp_last_om_idx, :);
    tmp_ep_sub = linker_str_array(iter_linker).ep_1_sub;
    tmp_ep2bv_vec = tmp_last_om_sub - tmp_ep_sub;
    tmp_ep2bv_vec_norm = sqrt(sum(tmp_ep2bv_vec.^2));
    tmp_ep2bv_vec = tmp_ep2bv_vec ./ tmp_ep2bv_vec_norm;
    linker_features.linker_orientation(:, iter_linker) = tmp_ep2bv_vec;
    linker_features.out_of_mask_ratio(iter_linker) = mean(tmp_linker_omQ);
end
linker_features.link_ep_orientation = linker_features.link_ep_orientation';
linker_features.linker_orientation = linker_features.linker_orientation';
linker_features.cos_ep2bv_link = sum(linker_features.link_ep_orientation .* linker_features.linker_orientation, 2);
linker_features.link_has_2ep_Q = linker_ep1_link_2epQ;
%% Train classifier
if num_linker > 1
    linker_features = struct2table(linker_features);
else
    linker_features = struct2table(linker_features, 'AsArray', true);
end
end
% function vec = fun_graph_compute_link_orientation(sub)
% % fun_graph_compute_link_orientation compute the orientation of the link
% % connected component by PCA. 
% % Input: 
% %   sub: N-by-3 numerical array, subscript of the link voxels 
% % Output:
% %   vec: 3-by-1 array, the normalized eigenvector correspond to the greatest
% %   eigenvalue.
% 
% if size(sub, 1) < 3
%     vec = [nan;nan;nan];
% else
%     [vec] = pca(sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
% end
% 
% end