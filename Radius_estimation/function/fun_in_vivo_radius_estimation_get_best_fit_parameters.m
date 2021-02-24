function est_psf_str = fun_in_vivo_radius_estimation_get_best_fit_parameters(simu_prof_str, ...
    est_r, est_z, exp_profile_xy, exp_profile_z)

est_psf_str = struct;

est_z = abs(est_z);

[~, vessel_radius_ind] = min(abs(simu_prof_str.vsl_rds_list_um - est_r));
[~, vessel_ori_ind] = min(abs(simu_prof_str.vsl_ori_z_comp - est_z));
num_para_combinations = numel(simu_prof_str.edge_int_str);
para_array_size = size(simu_prof_str.edge_int_str);
%% Algorithm 2 (Vectorized)
% Expand the interpolation range for cross-correlation fitting
% Only use the intensity simulation within 1.5 r

% This value was 1.5 for ML20180815 and ML20190124. Reduce to 1 because of
% the prevalent line-shift artefact in ML20200201
interpolation_selection_r_coeff_xy = 1.5; 
interpolation_selection_r_coeff_z = 2;
z_extended_valid_range_um = 4;
psf_int_xy_selected_Q = (exp_profile_xy.r_abs <= est_r * interpolation_selection_r_coeff_xy);
psf_int_z_selected_Q = (exp_profile_z.r_abs <= (est_r * interpolation_selection_r_coeff_z + z_extended_valid_range_um/2));

psf_int_xy = exp_profile_xy.r(psf_int_xy_selected_Q);
psf_int_xy_initial_shift_ind = round((psf_int_xy(1) - exp_profile_xy.r(1)) / exp_profile_xy.dr);

psf_int_z = exp_profile_z.r(psf_int_z_selected_Q);
psf_int_z_initial_shift_ind = round((psf_int_z(1) - exp_profile_z.r(1)) / exp_profile_z.dr);

psf_int_dist_xy = nan(numel(psf_int_xy), num_para_combinations);
psf_int_dist_z = nan(numel(psf_int_z), num_para_combinations);
% Flip the indices for correlation here
psf_int_xy = flip(psf_int_xy);
psf_int_z = flip(psf_int_z);
for iter_psf = 1 : num_para_combinations
    tmp_psf_int_r_dist_interp = simu_prof_str.edge_int_str{iter_psf}.normalized_radial_int_dist_itp{vessel_ori_ind, vessel_radius_ind};
    if ~isempty(tmp_psf_int_r_dist_interp)
        psf_int_dist_xy(:, iter_psf) = tmp_psf_int_r_dist_interp(abs(psf_int_xy));
    end
    tmp_psf_int_r_dist_vrt_interp = simu_prof_str.edge_int_str{iter_psf}.normalized_vertical_radial_int_dist_itp{vessel_ori_ind, vessel_radius_ind};
    if ~isempty(tmp_psf_int_r_dist_vrt_interp)
        psf_int_dist_z(:, iter_psf) = tmp_psf_int_r_dist_vrt_interp(abs(psf_int_z));
    end
end

% Compute the cross-correlation between the simulated intensity profile and
% the measured intensity profile
corr_xy = convn(exp_profile_xy.int_val_n_smooth, psf_int_dist_xy, 'full');
corr_xy = corr_xy ./ sqrt(sum(psf_int_dist_xy.^2, 1) .* sum(exp_profile_xy.int_val_n_smooth .^ 2, 1));

corr_z = convn(exp_profile_z.int_val_n_smooth, psf_int_dist_z, 'full');
corr_z = corr_z ./ sqrt(sum(psf_int_dist_z .^ 2, 1) .* sum(exp_profile_z.int_val_n_smooth .^ 2, 1));
% Find the largest correlation score for each psf
[corr_xy_psf_max, corr_xy_psf_max_ind] = max(corr_xy, [], 1);
[corr_z_psf_max, corr_z_psf_max_ind] = max(corr_z, [], 1);
% Find the largest correlation for both xy and z direction 
if est_z <= sqrt(2)/2
    corr_z_weight = 1 - est_z;
    psf_corr_mat = (corr_xy_psf_max + corr_z_psf_max * corr_z_weight) ./ (1 + corr_z_weight);
else
    psf_corr_mat = corr_xy_psf_max;
end
[psf_corr_mat_max, psf_corr_max_ind] = max(psf_corr_mat);

psf_corr_mat_xy = reshape(corr_xy_psf_max, para_array_size);
psf_xcorr_max_xy = corr_xy_psf_max(psf_corr_max_ind);
psf_xcorr_max_profile_xy = psf_int_dist_xy(:, psf_corr_max_ind);
psf_xcorr_max_ind_xy = corr_xy_psf_max_ind(psf_corr_max_ind);
psf_xcorr_max_ind_xy = psf_xcorr_max_ind_xy - psf_int_xy_initial_shift_ind;

psf_corr_mat_z = reshape(corr_z_psf_max, para_array_size);
psf_xcorr_max_z = corr_z_psf_max(psf_corr_max_ind);
psf_xcorr_max_profile_z = psf_int_dist_z(:, psf_corr_max_ind);
psf_xcorr_max_ind_z = corr_z_psf_max_ind(psf_corr_max_ind);
psf_xcorr_max_ind_z = psf_xcorr_max_ind_z - psf_int_z_initial_shift_ind;

[best_fit_para_1_ind, best_fit_para_2_ind] = ind2sub(para_array_size, psf_corr_max_ind);
%% Save data
est_psf_str.xy_corr_mat = psf_corr_mat_xy;
est_psf_str.xy_corr_max = psf_xcorr_max_xy;
est_psf_str.xy_corr_max_ind = psf_xcorr_max_ind_xy;
est_psf_str.xy_max_corr_int_profile = psf_xcorr_max_profile_xy;
est_psf_str.xy_max_corr_int_r = psf_int_xy;
est_psf_str.xy_max_corr_shift_ind = (psf_xcorr_max_ind_xy - numel(psf_int_xy));
est_psf_str.xy_max_corr_shift_um = est_psf_str.xy_max_corr_shift_ind * abs(exp_profile_xy.r_abs(2) - exp_profile_xy.r_abs(1));

est_psf_str.z_corr_mat = psf_corr_mat_z;
est_psf_str.z_corr_max = psf_xcorr_max_z;
est_psf_str.z_corr_max_ind = psf_xcorr_max_ind_z;
est_psf_str.z_max_corr_int_profile = psf_xcorr_max_profile_z;
est_psf_str.z_max_corr_int_r = psf_int_z;
est_psf_str.z_max_corr_shift_ind = (psf_xcorr_max_ind_z - numel(psf_int_z));
est_psf_str.z_max_corr_shift_um = est_psf_str.z_max_corr_shift_ind * abs(exp_profile_z.r_abs(2) - exp_profile_z.r_abs(1));

est_psf_str.corr_mat = reshape(psf_corr_mat, para_array_size);
est_psf_str.corr_max = psf_corr_mat_max;

est_psf_str.best_fit_parameter = [simu_prof_str.rbc_int_n_list(best_fit_para_1_ind), ...
    simu_prof_str.vsl_bc_lyr_thkn_um_list(best_fit_para_2_ind)];
est_psf_str.best_fit_rbc_int_n_ind = best_fit_para_1_ind;
est_psf_str.best_fit_blt_ind = best_fit_para_2_ind;
%% Visualization
% fprintf('Best fit PSF FWHM (%.2f, %.2f, %.2f). Corr max = %f\n', est_psf_str.best_fit_FWHM, est_psf_str.corr_max);
% fig_hdl = figure;
% ax_hdl_1 = subplot(1,2,1);
% plot(ax_hdl_1, exp_profile_xy.r, exp_profile_xy.int_val_n_smooth);
% hold(ax_hdl_1, 'on');
% % plot(ax_hdl_1, psf_int_xy + est_psf_str.xy_max_corr_shift_um, psf_int_dist_xy(:, psf_corr_max_ind));
% plot(ax_hdl_1, psf_int_xy + est_psf_str.xy_max_corr_shift_um, psf_int_dist_xy(:, 628));
% % 
% ax_hdl_2 = subplot(1,2,2);
% plot(ax_hdl_2, exp_profile_z.r, exp_profile_z.int_val_n_smooth);
% hold(ax_hdl_2, 'on');
% plot(ax_hdl_2, psf_int_z + est_psf_str.z_max_corr_shift_um, psf_xcorr_max_profile_z);
end