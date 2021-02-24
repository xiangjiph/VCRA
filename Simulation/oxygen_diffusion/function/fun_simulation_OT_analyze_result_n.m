function [oxy_result_str, varargout] = fun_simulation_OT_analyze_result_n(dt_array, pO2_array, ...
    internal_offset, local_max_window_size)
oxy_result_str = struct;

mask_size = size(pO2_array);
% Determine the interal region bounding box
int_bbox = [ones(1, 3) .* internal_offset, mask_size - internal_offset];
int_bbox(4:6) = int_bbox(4:6) - int_bbox(1:3) + 1;

% tmp_mask = crop_bbox3(vis_data_str.vessel_mask_rz, int_bbox);
vsl_mask = (dt_array == 0);
int_dt = crop_bbox3(dt_array, int_bbox);
int_pO2 = crop_bbox3(pO2_array, int_bbox);
int_vsl_mask = int_dt == 0;
int_dt = int_dt(~int_vsl_mask);
int_pO2 = int_pO2(~int_vsl_mask);

bin_x_edge = 0.5 : ceil(max(int_dt));
oxy_result_str.pO2_stat_in_bin = fun_analysis_get_y_stat_in_x_bin(int_dt, int_pO2, bin_x_edge);
%% Individual histogram
oxy_result_str.pO2_stat = fun_analysis_get_basic_statistics(int_pO2);
oxy_result_str.dt_stat = fun_analysis_get_basic_statistics(int_dt);
%% Finding the local DT maxima
is_dt_lm_Q = fun_array_local_maximum(dt_array, local_max_window_size);
is_dt_lm_Q = is_dt_lm_Q & ~vsl_mask;
dt_lm_ind = find(is_dt_lm_Q);
dt_lm_sub = fun_ind2sub(mask_size, dt_lm_ind);

is_not_near_boundary_Q = all(dt_lm_sub > internal_offset, 2) & ...
    all(dt_lm_sub < mask_size - internal_offset, 2);
dt_lm_sub = dt_lm_sub(is_not_near_boundary_Q, :);
dt_lm_ind = dt_lm_ind(is_not_near_boundary_Q);

dt_lm_v = dt_array(dt_lm_ind);
dt_lm_pO2_v = pO2_array(dt_lm_ind);
%% Find local oxygen concentration minimum
is_pO2_lm_Q = fun_array_local_maximum(- pO2_array, local_max_window_size);
is_pO2_lm_Q = is_pO2_lm_Q & ~vsl_mask;
pO2_lm_ind = find(is_pO2_lm_Q);
pO2_lm_sub = fun_ind2sub(mask_size, pO2_lm_ind);

is_not_near_boundary_Q = all(pO2_lm_sub > internal_offset, 2) & ...
    all(pO2_lm_sub < mask_size - internal_offset, 2);
pO2_lm_sub = pO2_lm_sub(is_not_near_boundary_Q, :);
pO2_lm_ind = pO2_lm_ind(is_not_near_boundary_Q);

pO2_lm_v = pO2_array(pO2_lm_ind);
pO2_lm_dt_v = dt_array(pO2_lm_ind);
%% Compute the distance between paired DT maximum and pO2 minimum
pdist_oxy_2_dt = pdist2(pO2_lm_sub, dt_lm_sub);
[list_idx_1, list_idx_2, min_dist] = fun_find_col_row_co_minimum(pdist_oxy_2_dt);
%% Compute the percentile of pO2
[oxy_result_str.hist2.count, oxy_result_str.hist2.dt_edge, oxy_result_str.hist2.pO2_edge] = ...
    histcounts2(int_dt, int_pO2);
%% Output result
oxy_result_str.internal_offset = internal_offset;
oxy_result_str.lm_wd_size = local_max_window_size;

oxy_result_str.dt_lm.dt_v = dt_lm_v;
oxy_result_str.dt_lm.pO2_v = dt_lm_pO2_v;
oxy_result_str.dt_lm.sub = dt_lm_sub;
oxy_result_str.dt_lm.ind = dt_lm_ind;

oxy_result_str.pO2_lm.dt_v = pO2_lm_dt_v;
oxy_result_str.pO2_lm.pO2_v = pO2_lm_v;
oxy_result_str.pO2_lm.sub = pO2_lm_sub;
oxy_result_str.pO2_lm.ind = pO2_lm_ind;

oxy_result_str.dist_lms.paired_pO2_list_idx = list_idx_1;
oxy_result_str.dist_lms.paired_dt_list_idx = list_idx_2;
oxy_result_str.dist_lms.paired_dist = min_dist;
if nargout > 1
    varargout{1} = int_dt;
    varargout{2} = int_pO2;
end
end