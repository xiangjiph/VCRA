function int_profile_str = fun_radius_estimation_get_z_int_profile(tile_image, local_max_sub, ...
    est_r, est_radial_ori, local_min_int, voxel_size)

if iscolumnvector(local_max_sub)
    local_max_sub = local_max_sub';
end

sample_max_ratio = inf;
z_search_range_coeff = 2.5;
z_extended_valid_l_um = 4;
z_extended_valid_l_pxl = round(z_extended_valid_l_um / voxel_size(3));
%% Get intensity interpolation
est_r_pxl = est_r ./ voxel_size(3);
min_search_range_pxl = z_extended_valid_l_pxl * 1.5; % extra valid range in the simulation 
z_search_half_range = max(min_search_range_pxl, round(z_search_range_coeff * est_r_pxl + z_extended_valid_l_pxl));
z_search_upper_lim = min(size(tile_image, 3), local_max_sub(3) + z_search_half_range);
z_search_lower_lim = max(1, local_max_sub(3) - z_search_half_range);
z_search_z_list = (z_search_lower_lim : 1 : z_search_upper_lim)';
z_search_r_list = (z_search_z_list - local_max_sub(3)) * voxel_size(3);

z_int_profile = double(squeeze(tile_image(local_max_sub(1), local_max_sub(2), ...
    z_search_z_list)));
% Interpolation 
z_int_profile_interpolation = griddedInterpolant(z_search_r_list, z_int_profile, 'linear', 'nearest');
z_search_r_list = (z_search_r_list(1) : (voxel_size(3)/2) : z_search_r_list(end))';
z_int_profile = z_int_profile_interpolation(z_search_r_list);
%% Smooth the profile before normalization
[int_profile_z_val_max0, int_profile_xy_val_max_ind] = max(z_int_profile);
if int_profile_z_val_max0 <= local_min_int
    warning('Maximum intensity along z-direction is smaller than the estimated local minimum intensity on xy-direction');
    local_min_int = min(z_int_profile);
end
assert(isfinite(local_min_int) && local_min_int >= 0);
z_int_profile_n = (z_int_profile - local_min_int) ./ (int_profile_z_val_max0 - local_min_int);

int_profile_z_val_max = mean(z_int_profile(max(1, int_profile_xy_val_max_ind - round(est_r_pxl/sample_max_ratio)) : ...
    min(numel(z_int_profile), int_profile_xy_val_max_ind + round(est_r_pxl/sample_max_ratio))));

z_int_profile = min(int_profile_z_val_max, z_int_profile);
% Normalization
z_int_profile_n_smooth = (z_int_profile - local_min_int) ./ (int_profile_z_val_max - local_min_int);

int_profile_str = struct;
int_profile_str.int_val_n = z_int_profile_n;
int_profile_str.int_val_n_smooth = z_int_profile_n_smooth;
int_profile_str.r = z_search_r_list;
int_profile_str.r_abs =  abs(z_search_r_list);
int_profile_str.dr = abs(z_search_r_list(1) - z_search_r_list(2));
int_profile_str.max_int = int_profile_z_val_max0;
int_profile_str.max_int_smooth = int_profile_z_val_max;
int_profile_str.ori_z = abs(est_radial_ori);
end