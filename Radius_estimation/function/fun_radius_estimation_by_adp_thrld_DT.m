function [est_r, varargout] = fun_radius_estimation_by_adp_thrld_DT(tile_image, seg_sub_in_tile_im, ...
    r_rough_um, seg_ori_vec, r_ori_to_edge_int_n, tile_voxel_size, num_iteration)
%% Initialization
persistent imclose_strel skel_sample_sub
vis_Q = false;
est_r = nan;
if nargout > 1
    varargout{1} = seg_sub_in_tile_im;
end
if iscolumnvector(seg_sub_in_tile_im)
    seg_sub_in_tile_im = seg_sub_in_tile_im';
end
if ~isfinite(r_rough_um) || r_rough_um <= 0
    return;
end
seg_sub_bbox_mm = min(seg_sub_in_tile_im, [], 1);
seg_sub_bbox_xx = max(seg_sub_in_tile_im, [], 1);
num_seg_voxels = size(seg_sub_in_tile_im, 1);
seg_mid_sub = seg_sub_in_tile_im(ceil(num_seg_voxels/2), :);
rendered_im_size = size(tile_image);
assert(all(seg_sub_bbox_xx <= rendered_im_size), 'Segment voxel coordinate out of range');
%% Parameters
search_expand_coeff = 1.5;
exp_times = 0;
max_exp_times = 3;
bg_ptl = 0.1;
min_est_bg_std = 200;
max_est_bg_est = 1000;
min_crop_window_r = 3;
skel_int_sample_r = 2;
imclose_cube_size = 3;

avg_voxel_size_xy = (tile_voxel_size(1) + tile_voxel_size(2))/2;

if isempty(imclose_strel)
    imclose_strel = strel('cube', imclose_cube_size);
end
%% Calculation that need to be done once for each segments
% Estimate the tilt angle
% ori_vec = fun_radius_estimation_get_segment_orientation_vector(seg_sub_in_tile_im);
ori_vec_z = abs(seg_ori_vec(3));
% Estimate the skeleton central intensity
if isempty(skel_sample_sub)
    skel_sample_strel = strel('disk', skel_int_sample_r); % disk of radius 2 pixels on xy direction
    skel_sample_strel = skel_sample_strel.Neighborhood;
    skel_sample_sub = fun_ind2sub(size(skel_sample_strel), find(skel_sample_strel));
    skel_sample_sub = skel_sample_sub - skel_sample_sub(ceil(size(skel_sample_sub, 1) / 2), :);
end

seg_neighbor_int_max = nan(num_seg_voxels, 1);
% For each skeleton voxels, sample its neighbors pixel intensity within
% the upsampling error range
for iter_vxl = 1 : num_seg_voxels
    tmp_vxl_sub = seg_sub_in_tile_im(iter_vxl, :);
    tmp_vxl_sub_xy = tmp_vxl_sub(1:2);
    assert(isrow(tmp_vxl_sub_xy), 'tmp_vxl_sub_xy is not a row vector');
    tmp_vxl_neighbor_xy = bsxfun(@plus, tmp_vxl_sub_xy, skel_sample_sub);
    tmp_valid_neighbor_Q = all(bsxfun(@ge, tmp_vxl_neighbor_xy, [1,1]), 2) & ...
        all(bsxfun(@le, tmp_vxl_neighbor_xy, rendered_im_size(1:2)), 2);
    tmp_vxl_neighbor_xy = tmp_vxl_neighbor_xy(tmp_valid_neighbor_Q, :);
    tmp_vxl_neighbor_ind = sub2ind(rendered_im_size, tmp_vxl_neighbor_xy(:, 1), ...
        tmp_vxl_neighbor_xy(:, 2), repelem(tmp_vxl_sub(3), nnz(tmp_valid_neighbor_Q), 1));
    tmp_vlx_neighbor_int = tile_image(tmp_vxl_neighbor_ind);
    seg_neighbor_int_max(iter_vxl) = max(tmp_vlx_neighbor_int);
end
skel_seg_max_int = mean(seg_neighbor_int_max);
%% Calculation that need iterative update
est_r = r_rough_um;
window_crop_r = r_rough_um;
update_window_Q = true;
est_r_record = nan(num_iteration, 1);

for iter_est = 1 : num_iteration
    if ~isfinite(est_r)
        if exp_times < max_exp_times
            est_r = r_rough_um;
            update_window_Q = true;
            exp_times = exp_times + 1;
            search_expand_coeff = search_expand_coeff * 1.5;
            window_crop_r = r_rough_um; % Good? 
        else
            break;
        end
    elseif est_r == 0
        est_r = nan;
        break;
    elseif window_crop_r < est_r
        window_crop_r = est_r;
        update_window_Q = true;
    end
    if update_window_Q
        update_window_Q = false;
        %% Determine the cropping boudning box
        im_crop_size_um = max(min_crop_window_r, ceil(window_crop_r .* search_expand_coeff));
        im_crop_r_pxl = ceil(im_crop_size_um ./ tile_voxel_size);
        im_crop_r_pxl(3) = 1;
        im_crop_bbox_mm = max([1,1,1], seg_sub_bbox_mm - im_crop_r_pxl);
        im_crop_bbox_xx = min(rendered_im_size, seg_sub_bbox_xx + im_crop_r_pxl);
        im_crop_bbox_ll = im_crop_bbox_xx - im_crop_bbox_mm + 1;
        im_crop_bbox_mmll = [im_crop_bbox_mm, im_crop_bbox_ll];
        im_crop_bbox_mmxx = [im_crop_bbox_mm, im_crop_bbox_xx];
        num_local_im_voxel = prod(im_crop_bbox_mmxx(4:6) - im_crop_bbox_mmxx(1:3) + 1);
        % Crop the image and cast to single precision data type
        local_render_im = single(crop_bbox3(tile_image, im_crop_bbox_mmll));
        %% Background estimation
        local_render_im_int_sort = sort(local_render_im(:), 'ascend');
        local_render_im_int_min = local_render_im_int_sort(round(num_local_im_voxel * bg_ptl));
        local_render_im_int_std = std(local_render_im_int_sort(1 : round(num_local_im_voxel * bg_ptl)));
        local_render_im_int_std = max(local_render_im_int_std, min_est_bg_std);
        if local_render_im_int_std > max_est_bg_est
%             fprintf('Image with high background noise level\n');
            local_render_im_int_std = max_est_bg_est;
        end
        local_render_bg_max = local_render_im_int_min + 2 * local_render_im_int_std;
    end
    %% Estimate the edge intensity for thresholding
    min_edge_ratio = r_ori_to_edge_int_n(ori_vec_z, double(est_r));
    if min_edge_ratio > 1
        warning('The normalized edge intensity w.r.t. center intensity is greater than 1');
    end
    local_th_int = (skel_seg_max_int - local_render_im_int_min) * min_edge_ratio + local_render_im_int_min;
    if local_th_int < local_render_bg_max
        %         fprintf('The local intensity threshold is lower than the estimated noise level\n');
        local_th_int = local_render_bg_max;
    end
    %% Threshold the image
    local_dt_max_search_r = max([1,1,1], min([4,4,1], ceil(est_r ./ tile_voxel_size)));
    local_render_sub = seg_mid_sub - im_crop_bbox_mmxx(1:3) + 1;
    dt_search_z_min = max(1, local_render_sub(3) - local_dt_max_search_r(3));
    dt_search_z_max = min(im_crop_bbox_ll(3), local_render_sub(3) + local_dt_max_search_r(3));
    
    dt_search_im = local_render_im(:, :, dt_search_z_min : dt_search_z_max);
    
    dt_search_mask = dt_search_im >= local_th_int;
    % Imclose to smooth the mask - morphological close is necessary for
    % improving accuracy
    dt_search_mask = imclose(dt_search_mask, imclose_strel);
    %% Estimate radius from 2D distance transform
    local_render_mask_dt = zeros(size(dt_search_mask), 'single');
    for iter_sec = 1 : size(dt_search_mask, 3)
        local_render_mask_dt(:, :, iter_sec) = bwdist(~dt_search_mask(:, :, iter_sec));
    end
    % Find the largest DT in the neighbor of the input voxel
    dt_search_sub = local_render_sub;
    dt_search_sub(3) = dt_search_sub(3) - dt_search_z_min + 1;
        
    [r_est_neighbor, r_est_neighbor_bbox_mmxx]= crop_center_box(local_render_mask_dt, dt_search_sub, local_dt_max_search_r);
    % Scale back to radius in um
    [est_r, est_r_max_ind] = max(r_est_neighbor, [], 'all', 'linear');
    est_r = double(est_r) .* avg_voxel_size_xy;    
    est_r_max_sub = fun_ind2sub(r_est_neighbor_bbox_mmxx(4:6) - r_est_neighbor_bbox_mmxx(1:3) + 1, ...
        est_r_max_ind);
    est_r_max_sub = est_r_max_sub + r_est_neighbor_bbox_mmxx(1:3) - 1;
    est_r_max_sub(3) = est_r_max_sub(3) + dt_search_z_min - 1;
    est_r_max_sub = est_r_max_sub + im_crop_bbox_mmll(1:3) - 1;
    assert(isscalar(est_r), 'Estimated radius is not a scalar');
    
    if isfinite(est_r)
        est_r_record(iter_est) = est_r;
    end
        
    rep_ind = find(est_r_record == est_r);
    
    if iter_est >= 2
        if any(abs(diff(est_r_record)) > 0.25) && ori_vec_z < sind(60) && est_r > 1
            fprintf('Debug\n');
            keyboard;
       end        
    end
    
    if numel(rep_ind) >= 2
        assert(issorted(rep_ind, 'ascend'), 'rep_ind is not in ascending order');
        last_rep_ind = (rep_ind(end-1) + 1) : 1 : rep_ind(end);
        est_r = mean(est_r_record(last_rep_ind));
        break;
    end
end
if ~isfinite(est_r) || est_r <= 0
    est_r = nan;
    return;
else
    if nargout > 1
        varargout{1} = est_r_max_sub;
    end
end

if vis_Q
    %%
    [r_est_neighbor, r_est_bbox_mmxx] = crop_center_box(local_render_mask_dt, dt_search_sub, local_dt_max_search_r);
    [r_max, r_max_ind] = max(r_est_neighbor(:));
    r_max_sub = fun_ind2sub(r_est_bbox_mmxx(4:6) - r_est_bbox_mmxx(1:3) + 1, ...
        r_max_ind);
    r_max_local_sub = r_est_bbox_mmxx(1:3) + r_max_sub - 1;
    
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
    ax_1 = subplot(1,3,1);
    imagesc(ax_1, local_render_im(:, :, local_render_sub(3)));
    hold(ax_1, 'on');
    scatter(ax_1, local_render_sub(2), local_render_sub(1), 50, 'LineWidth', 4, 'Marker', '+');
    colormap(ax_1, 'jet');
    ax_1.Title.String = sprintf('Full resolution image\nSegment elevation angle %d\\circ', round(asind(ori_vec_z)));
    ax_1.DataAspectRatio = [1,1,1];
    cbar_hdl_1 = colorbar(ax_1);
    cbar_hdl_1.Label.String = 'Intensity';
    cbar_hdl_1.FontSize = 14;
    ax_1.XAxis.Visible = 'off';
    ax_1.YAxis.Visible = 'off';
    ax_1.FontSize = 14;
    ax_2 = subplot(1,3,2);
    imagesc(ax_2, local_render_mask_dt(:, :, r_max_sub(3)) > 0);
    hold(ax_2, 'on');
    scatter(ax_2, r_max_local_sub(2), r_max_local_sub(1), 50, 'LineWidth', 4, 'Marker', '+');
    ax_2.DataAspectRatio = [1,1,1];
    ax_2.Title.String = sprintf('Thresholded image\nThreshold %d', round(local_th_int));
    ax_2.XAxis.Visible = 'off';
    ax_2.YAxis.Visible = 'off';
    ax_2.FontSize = 14;
    ax_2.Colormap(1, :) = [0, 0, 0];
    ax_2.Colormap(end, :) = [1, 0, 0];
    cbar_hdl_2 = colorbar(ax_2);
    cbar_hdl_2.Label.String = 'Pixel value';
    cbar_hdl_2.FontSize = 14;
    
    ax_3 = subplot(1,3,3);
    imagesc(ax_3, local_render_mask_dt(:, :, r_max_sub(3)) .* avg_voxel_size_xy);
    colormap(ax_3, 'jet');
    hold(ax_3, 'on');
    scatter(ax_3, r_max_local_sub(2), r_max_local_sub(1), 50, 'LineWidth', 4, 'Marker', '+');
    ax_3.DataAspectRatio = [1,1,1];
    ax_3.Title.String = sprintf('Distance transform\nEstimated radius %.2f \\mum', est_r);
    cbar_hdl_3 = colorbar(ax_3);
    cbar_hdl_3.Label.String = 'Distance to edge (\mum)';
    cbar_hdl_3.FontSize = 14;
    ax_3.XAxis.Visible = 'off';
    ax_3.YAxis.Visible = 'off';
    ax_3.FontSize = 14;
    %%
    dataset_name = 'WholeBrain';
    stack = 'ML_2018_08_15';
    DataManager = FileManager;
    vis_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Radius_estimation', ...
        sprintf('%s_%s_RE_r_%d_nm_ele_agl_%d_deg_iter_%d.png', dataset_name, stack, round(est_r * 1000), round(asind(ori_vec_z)), iter_est));
    fun_print_image_in_several_formats(fig_hdl, vis_fp);
end
end