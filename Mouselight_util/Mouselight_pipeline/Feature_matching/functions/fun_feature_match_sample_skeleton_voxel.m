function [des_1_sub, des_2_sub_shift, des_1_list_ind, des_2_list_ind] = fun_feature_match_sample_skeleton_voxel(des_1_sub, ...
    des_2_sub_shift, des_1_label, des_2_label, opt)
% fun_feature_match_sample_descriptor is a helper function for
% searchpair_vessel. It samples the descriptors from two point clouds to
% reduce the computation load of the coherent point drift
%

if isempty(des_1_sub) || isempty(des_2_sub_shift)
    return;
end

if isfield(opt, 'total_num_descriptor')
    total_num_descriptor = opt.total_num_descriptor;
else
    total_num_descriptor = 1e4;
end
if isfield(opt, 'rm_descriptor_by_pdist2_Q')
    rm_descriptor_by_pdist2_Q = opt.rm_descriptor_by_pdist2_Q;
    max_disp_pixel_yxz = opt.max_disp_pixel_yxz;
else
    rm_descriptor_by_pdist2_Q = false;
end
total_num_descriptor_neighbor_search = sqrt(100 * 1024 * 1024 * 1024 / 4);
%% Remove voxels in one set that are far away from any other voxels in the other set
if rm_descriptor_by_pdist2_Q
    if numel(des_1_label) > total_num_descriptor_neighbor_search || ...
            numel(des_2_label) > total_num_descriptor_neighbor_search
        warning('The neighbor searching require more than 100GB of memory. Terminated');
        des_1_list_ind = (1 : size(des_1_sub, 1))';
        des_2_list_ind = (1 : size(des_2_sub, 1))';
    else
        tmp_pdist = pdist2(single(des_1_sub(:,1)), single(des_2_sub_shift(:,1)));
        %     tmp_pdist3 = (tmp_pdist.^2) ./3;
        tmp_pdist_reasonable = tmp_pdist < max_disp_pixel_yxz(1);
        tmp_pdist = pdist2(single(des_1_sub(:,2)), single(des_2_sub_shift(:,2)));
        %     tmp_pdist3 = tmp_pdist3 + (tmp_pdist.^2) ./3;
        tmp_pdist_reasonable = tmp_pdist_reasonable & tmp_pdist < max_disp_pixel_yxz(2);
        tmp_pdist = pdist2(single(des_1_sub(:,3)), single(des_2_sub_shift(:,3)));
        %     tmp_pdist3 = sqrt(tmp_pdist3 + (tmp_pdist.^2));
        tmp_pdist_reasonable = tmp_pdist_reasonable & tmp_pdist < max_disp_pixel_yxz(3);
        des_1_list_ind = find(any(tmp_pdist_reasonable, 2));
        des_2_list_ind = find(any(tmp_pdist_reasonable, 1)');
        clearvars tmp_pdist tmp_pdist_reasonable
    end
else
    des_1_list_ind = (1 : size(des_1_sub, 1))';
    des_2_list_ind = (1 : size(des_2_sub, 1))';
end
if isempty(des_1_list_ind) || isempty(des_2_list_ind)
    des_1_sub = des_1_sub(sample_ind_1, :);
    des_2_sub_shift = des_2_sub_shift(sample_ind_2, :);
    return;
end
%% Uniformally sample the voxels in each connected component
num_voxel_1 = size(des_1_sub, 1);
num_voxel_2 = size(des_2_sub_shift, 1);

% sample_step = min(floor([num_voxel_1, num_voxel_2] ./ total_num_descriptor));
sample_step = 2;
if num_voxel_1 > total_num_descriptor && num_voxel_2 > total_num_descriptor
    des_1_list_ind = des_1_list_ind(1 : sample_step : end);
    des_2_list_ind = des_2_list_ind(1 : sample_step : end);
end
des_1_sub = des_1_sub(des_1_list_ind, :);
des_1_label = des_1_label(des_1_list_ind);
des_2_sub_shift = des_2_sub_shift(des_2_list_ind, :);
des_2_label = des_2_label(des_2_list_ind);
%% Keep the long connected components
if size(des_1_sub, 1) > total_num_descriptor
    tmp_idx_list_1 = fun_bin_data_to_idx_list(des_1_label);
    tmp_seg_length_1 = cellfun(@numel, tmp_idx_list_1);
    [tmp_seg_length_1, tmp_seg_idx] = sort(tmp_seg_length_1, 'descend');
    tmp_idx_list_1 = tmp_idx_list_1(tmp_seg_idx);
    tmp_cumsum_voxel = cumsum(tmp_seg_length_1);
    [~, cutoff_idx] = min(abs(tmp_cumsum_voxel - total_num_descriptor));
    des_1_selected_ind = cat(2, tmp_idx_list_1{1:cutoff_idx});
    des_1_sub = des_1_sub(des_1_selected_ind, :);
    des_1_list_ind = des_1_list_ind(des_1_selected_ind);
end

if size(des_2_sub_shift, 1) > total_num_descriptor
    tmp_idx_list_2 = fun_bin_data_to_idx_list(des_2_label);
    tmp_seg_length_2 = cellfun(@numel, tmp_idx_list_2);
    [tmp_seg_length_2, tmp_seg_idx] = sort(tmp_seg_length_2, 'descend');
    tmp_idx_list_2 = tmp_idx_list_2(tmp_seg_idx);
    tmp_cumsum_voxel = cumsum(tmp_seg_length_2);
    [~, cutoff_idx] = min(abs(tmp_cumsum_voxel - total_num_descriptor));
    des_2_selected_ind = cat(2, tmp_idx_list_2{1:cutoff_idx});
    des_2_sub_shift = des_2_sub_shift(des_2_selected_ind, :);
    des_2_list_ind = des_2_list_ind(des_2_selected_ind);
end
%% Other trial on sampling:
%     cc_idx_1 = fun_bin_data_to_idx_list(des_1_label);
%     num_voxel_1 = numel(des_1_label);
%     num_cc_voxel_1 = cellfun(@length, cc_idx_1);
%     cc_idx_2 = fun_bin_data_to_idx_list(des_2_label);
%     num_voxel_2 = numel(des_2_label);
%     num_cc_voxel_2 = cellfun(@length, cc_idx_2);
%     % Map from idx to label
%     idx_2_label_1 = repelem(1:numel(num_cc_voxel_1), num_cc_voxel_1);
%     idx_2_label_2 = repelem(1:numel(num_cc_voxel_2), num_cc_voxel_2);
%
%     tmp_cc_label = 10;
%     tmp_cc_idx = cc_idx_1{tmp_cc_label};
%     tmp_nn_idx = des_1_nn_idx(tmp_cc_idx);
%     tmp_nn_idx_unique = unique(tmp_nn_idx);
%     tmp_matched_ratio = numel(tmp_nn_idx_unique)/numel(tmp_nn_idx);
%     tmp_nn_label = unique(idx_2_label_2(tmp_nn_idx_unique));
%     tmp_nn_dist = des_1_nn_dist(tmp_cc_idx);
%     % visualize the nearest neighbor pair
%
%     tmp_vis_sub_1 = des_1_sub(tmp_cc_idx,:);
%     tmp_vis_sub_2 = des_2_sub_shift(cat(2, cc_idx_2{tmp_nn_label(1)}),:);
% %     figure;
%     clf
%     scatter3(tmp_vis_sub_1(:,1), tmp_vis_sub_1(:,2), tmp_vis_sub_1(:,3));
%     hold on
%     scatter3(tmp_vis_sub_2(:,1), tmp_vis_sub_2(:,2), tmp_vis_sub_2(:,3));
%     legend('Tile 1', 'Tile 2');
%
% Maybe better solution: chop skeletons into short pieces and sample uniformly in
% the overlapping space.
end