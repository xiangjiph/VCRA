function est_r = fun_radius_estimation_by_fitting(tile_image, tile_sub, est_r_um, ...
    est_edge_int_n, ori_vec, tile_vxl_size_um, local_max_search_Q)

if nargin < 7
    local_max_search_Q = false;
end
vis_Q = false;
est_r = nan;
if ~isfinite(est_r_um) || est_r_um <= 0
    return;
end
%% Find the local intensity maximum and update the tile_sub position 
if local_max_search_Q
%     num_dims = numel(tile_sub);    
    local_int_max_search_r = max([1,1,1], min([3, 3, 1], ceil(est_r_um ./ tile_vxl_size_um)));
    [local_im, local_bbox_mmxx] = crop_center_box(tile_image, tile_sub, local_int_max_search_r);
    [local_int_max, local_max_ind] = max(local_im(:)); %#ok<ASGLU>
    local_max_sub = fun_ind2sub(size(local_im), local_max_ind);
    tile_sub = local_max_sub + local_bbox_mmxx(1:3) - 1;
end
%% Get the intensity profile along the radial direction 
int_profile_str_xy = fun_radius_estimation_get_xy_int_profile(tile_image, tile_sub, ...
    est_r_um, ori_vec, tile_vxl_size_um, vis_Q);
%% Fitting to determine the position of the boundary at sub-voxel resolution 
zero_ind = (numel(int_profile_str_xy.r) + 1) / 2;
neg_int_smooth = int_profile_str_xy.int_val_n_smooth(1 : zero_ind);
neg_last_le_th_ind = find(neg_int_smooth <= est_edge_int_n, 1, 'last');
neg_first_ge_th_ind = find(neg_int_smooth >= est_edge_int_n, 1, 'first');
if neg_last_le_th_ind == neg_first_ge_th_ind
    neg_r_th = int_profile_str_xy.r(neg_last_le_th_ind);
else
    neg_r_th = int_profile_str_xy.r(neg_last_le_th_ind) + ...
        (est_edge_int_n - neg_int_smooth(neg_last_le_th_ind)) / (neg_int_smooth(neg_first_ge_th_ind) - neg_int_smooth(neg_last_le_th_ind)) * ...
        int_profile_str_xy.dr;
end

pos_int_smooth = int_profile_str_xy.int_val_n_smooth(zero_ind : end);
pos_last_gt_th_ind = find(pos_int_smooth >= est_edge_int_n, 1, 'last');
pos_first_le_th_ind = find(pos_int_smooth <= est_edge_int_n, 1, 'first');
if pos_last_gt_th_ind == pos_first_le_th_ind
    pos_r_th = int_profile_str_xy.r(pos_last_gt_th_ind + zero_ind - 1);
else
    pos_r_th = int_profile_str_xy.r(pos_last_gt_th_ind + zero_ind - 1) + ...
        (pos_int_smooth(pos_last_gt_th_ind) - est_edge_int_n) / (pos_int_smooth(pos_last_gt_th_ind) - pos_int_smooth(pos_first_le_th_ind)) * ...
        int_profile_str_xy.dr;
end
est_r = (pos_r_th - neg_r_th) / 2;
if isempty(est_r)
    est_r = nan;
end
%% visualization
if vis_Q
    fig_hdl = figure;
    ax_hdl = subplot(1,2,1);
    scatter(ax_hdl, int_profile_str_xy.r, int_profile_str_xy.int_val_n);
    hold(ax_hdl, 'on');
    scatter(ax_hdl, int_profile_str_xy.r, int_profile_str_xy.int_val_n_smooth);
    
    local_im = crop_center_box(tile_image, tile_sub, [20. 20, 0]);
    ax_hdl_2 = subplot(1,2,2);
    imagesc(ax_hdl_2, local_im);
    colormap(ax_hdl_2, 'gray');
    ax_hdl_2.DataAspectRatio = [1,1,1];
end
end