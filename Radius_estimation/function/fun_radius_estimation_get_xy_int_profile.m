function int_profile_str = fun_radius_estimation_get_xy_int_profile(tile_image, vxl_sub, est_r_um, seg_ori_vec, voxel_size_um, vis_Q)


if nargin < 6
    vis_Q = false;
end
if iscolumnvector(vxl_sub)
    vxl_sub = vxl_sub';
end
% Convert the vessel segment orientation vector to the unit vector along
% the radial direction
if numel(seg_ori_vec) == 3
    est_radial_ori = seg_ori_vec(1:2);
end
if any(est_radial_ori ~= 0)
    % Get the vector along the radial direction
    est_radial_ori = est_radial_ori([2, 1]);
    est_radial_ori(1) = -est_radial_ori(1);
    est_radial_ori = est_radial_ori ./ sqrt(sum(est_radial_ori.^2));
else
    est_radial_ori = [0, 1];
end
if isrowvector(est_radial_ori)
    est_radial_ori = est_radial_ori';
end
%% Parameters
sample_max_ratio = 12.8;
search_range_coeff = 2;
min_search_range_pxl = 15;
% Searching range of narest intensity profile near the radial direction. 
% +/- 30 degree 
rot_theta_list = pi/36 * [-6 : 1 : 6];
%% Get intensity interpolation
num_dims = numel(vxl_sub);
est_r_pxl = est_r_um ./ voxel_size_um(1:2);
est_r_pxl_min = min(est_r_pxl);
xy_search_range = max(min_search_range_pxl, round(search_range_coeff * est_r_pxl));
xy_search_range_min = min(xy_search_range);
switch num_dims
    case 2
        [test_plane_im, test_plane_bbox_mmxx] = crop_center_box(tile_image, vxl_sub, xy_search_range);
    case 3
        [test_plane_im, test_plane_bbox_mmxx] = crop_center_box(tile_image, vxl_sub, [xy_search_range, 0]);
    otherwise
        error('Unexpected error');
end
test_plane_im = double(test_plane_im);
[test_plane_sub_1, test_plane_sub_2] = ndgrid(test_plane_bbox_mmxx(1):test_plane_bbox_mmxx(num_dims+1), ...
    test_plane_bbox_mmxx(2):test_plane_bbox_mmxx(num_dims+2));
test_plane_int_interpolation = griddedInterpolant(test_plane_sub_1, test_plane_sub_2, test_plane_im, 'linear', 'nearest');
%% Determine the radial direction on xy-plane

num_theta = numel(rot_theta_list);

min_int_profile_std = inf;
% min_std_theta_ind = 0;
for iter_theta = 1 : num_theta
    rot_theta = rot_theta_list(iter_theta);
    rot_mat = [cos(rot_theta), -sin(rot_theta); sin(rot_theta), cos(rot_theta)];
    test_radial_ori = (rot_mat * est_radial_ori)';
    
    % int_profile_xy_pos = test_voxel_local_max_sub(1:2) + est_radial_ori .* ...
    %     [-xy_search_range_min : 1 : xy_search_range_min]';
    int_profile_xy_pos = vxl_sub(1:2) + test_radial_ori .* ...
        [-xy_search_range_min : 1 : xy_search_range_min]';
    
    int_profile_xy_pos_r = int_profile_xy_pos - vxl_sub(1:2);
    int_profile_xy_pos_r = int_profile_xy_pos_r .* voxel_size_um(1:2);
    tmp_int_profile_xy_r_abs = sqrt(sum(int_profile_xy_pos_r.^2, 2));
    tmp_ind = 1 : round((numel(tmp_int_profile_xy_r_abs)-1) / 2);
    tmp_int_profile_xy_r = tmp_int_profile_xy_r_abs;
    tmp_int_profile_xy_r(tmp_ind) = - tmp_int_profile_xy_r(tmp_ind);    
    
    int_profile_xy_val = test_plane_int_interpolation(int_profile_xy_pos(:, 1), ...
        int_profile_xy_pos(:, 2));
    tmp_int_profile_xy_val_min = min(int_profile_xy_val);
    [tmp_int_profile_xy_val_max0, int_profile_xy_val_max_ind] = max(int_profile_xy_val);
    
    tmp_int_profile_xy_val_n = (int_profile_xy_val - tmp_int_profile_xy_val_min) ./ (tmp_int_profile_xy_val_max0 - tmp_int_profile_xy_val_min);
    % Average over pixels near the maximum.
    max_sample_ind_sub = max(1, int_profile_xy_val_max_ind - round(est_r_pxl_min/sample_max_ratio));
    max_sample_ind_sup = min(numel(int_profile_xy_val), int_profile_xy_val_max_ind + round(est_r_pxl_min/sample_max_ratio));
    
    tmp_int_profile_xy_val_max = mean(int_profile_xy_val(max_sample_ind_sub : max_sample_ind_sup));
    
    int_profile_xy_val = min(tmp_int_profile_xy_val_max, int_profile_xy_val);
    % Normalization
    if tmp_int_profile_xy_val_max == tmp_int_profile_xy_val_min
        tmp_int_profile_xy_val_n_smooth = ones(size(int_profile_xy_val));
    else
        tmp_int_profile_xy_val_n_smooth = (int_profile_xy_val - tmp_int_profile_xy_val_min) ./ (tmp_int_profile_xy_val_max - tmp_int_profile_xy_val_min);
    end
    
    % Approximate the profile as gaussian, compute the standard deviation
    weighted_avg_mean =  (tmp_int_profile_xy_val_n_smooth' * tmp_int_profile_xy_r) / sum(tmp_int_profile_xy_val_n_smooth);
    weighted_std =  (tmp_int_profile_xy_val_n_smooth' * (tmp_int_profile_xy_r - weighted_avg_mean).^2) / sum(tmp_int_profile_xy_val_n_smooth);
    
    assert(isfinite(weighted_std), 'Intensity profile std is not finite.');    
    if weighted_std < min_int_profile_std
        min_int_profile_theta = rot_theta;
        min_int_profile_mean = weighted_avg_mean;
        min_int_profile_std = weighted_std;
        int_profile_xy_val_n_smooth = tmp_int_profile_xy_val_n_smooth;
        int_profile_xy_val_n = tmp_int_profile_xy_val_n;
        int_profile_xy_r_abs = tmp_int_profile_xy_r_abs;
        int_profile_xy_r = tmp_int_profile_xy_r;
        int_profile_xy_val_max0 = tmp_int_profile_xy_val_max0;
        int_profile_xy_val_max = tmp_int_profile_xy_val_max;
        int_profile_xy_val_min = tmp_int_profile_xy_val_min;
%         min_std_theta_ind = iter_theta;
        %% Visualization
        if vis_Q
            tmp_xy_pos = int_profile_xy_pos;
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,1];
            ax_hdl = subplot(1,2,1);
            imagesc(ax_hdl, test_plane_im);
            hold(ax_hdl, 'on');
            line_sub = tmp_xy_pos - test_plane_bbox_mmxx(1:2);
            scatter(ax_hdl, line_sub(:, 2), line_sub(:, 1));
            ax_hdl_2 = subplot(1,2,2);
            scatter(ax_hdl_2, tmp_int_profile_xy_r, tmp_int_profile_xy_val_n_smooth);
            hold(ax_hdl_2, 'on');
            plot(ax_hdl_2, tmp_int_profile_xy_r, exp(- (tmp_int_profile_xy_r - weighted_avg_mean).^2 ./ (2 * weighted_std^2)));
            ax_hdl_2.Title.String = sprintf('Rotated angle %.2f rad', rot_theta);
        end
    end
end

%%
int_profile_str = struct;
int_profile_str.rotate_angle_theta = min_int_profile_theta;
int_profile_str.rotate_angle_deg = min_int_profile_theta * 180 / pi;
int_profile_str.gaussian_mean = min_int_profile_mean;
int_profile_str.gaussian_std = min_int_profile_std;
int_profile_str.int_val_n = int_profile_xy_val_n;
int_profile_str.int_val_n_smooth = int_profile_xy_val_n_smooth;
int_profile_str.r_abs = int_profile_xy_r_abs;
int_profile_str.r = int_profile_xy_r;
int_profile_str.dr = abs(int_profile_str.r_abs(1) - int_profile_str.r_abs(2));
int_profile_str.max_int = int_profile_xy_val_max0;
int_profile_str.max_int_smooth = int_profile_xy_val_max;
int_profile_str.min_int = int_profile_xy_val_min;

end