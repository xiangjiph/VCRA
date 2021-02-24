function dt_ori_str = fun_analysis_cube_dt_grad_pca(brain_mask_dt, bbox_mmxx)

bbox_mmxx(:, 4:6) = min(bbox_mmxx(:, 4:6), size(brain_mask_dt));
num_bbox = size(bbox_mmxx, 1);
bbox_dt_ori = nan(3, num_bbox);
bbox_dt_ori_s = nan(3, num_bbox);
bbox_dist_2_brain_surface = nan(num_bbox, 1);
for iter_bbox = 1 : num_bbox
    tmp_bbox_mmll = bbox_mmxx(iter_bbox, :);
    tmp_bbox_mmll(4:6) = tmp_bbox_mmll(4:6) - tmp_bbox_mmll(1:3) + 1;
    tmp_dt_array = crop_bbox3(brain_mask_dt, tmp_bbox_mmll);
    tmp_dt_array_in_brain_Q = (tmp_dt_array > 0);
    if any(tmp_dt_array_in_brain_Q)
        bbox_dist_2_brain_surface(iter_bbox) = mean(tmp_dt_array(tmp_dt_array_in_brain_Q));
        [tmp_dt_1, tmp_dt_2, tmp_dt_3] = fun_gradient3D(tmp_dt_array, [1,2,3]);
        tmp_grad_vec = cat(2, tmp_dt_1(:), tmp_dt_2(:), tmp_dt_3(:));
        tmp_grad_vec = tmp_grad_vec(any(tmp_grad_vec ~= 0, 2), :);
        [tmp_u, tmp_s, ~] = svd( (tmp_grad_vec.' * tmp_grad_vec), 'econ');
        bbox_dt_ori(:, iter_bbox) = tmp_u(:, 1);
        bbox_dt_ori_s(:, iter_bbox) = tmp_s([1,5,9]);
    end
end
dt_ori_str.ori_eig_vec = bbox_dt_ori.';
dt_ori_str.ori_eig_val = bbox_dt_ori_s.';
dt_ori_str.ori_eig_ratio = dt_ori_str.ori_eig_val ./ sum(dt_ori_str.ori_eig_val, 2);
dt_ori_str.avg_dist_2_surface = bbox_dist_2_brain_surface;
end