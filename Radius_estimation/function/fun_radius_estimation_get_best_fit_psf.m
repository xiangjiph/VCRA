function est_psf_str = fun_radius_estimation_get_best_fit_psf(psf_est_int, ...
    est_r, est_z, exp_profile_xy, exp_profile_z)

est_psf_str = struct;

est_z = abs(est_z);

[~, vessel_radius_ind] = min(abs(psf_est_int.vessel_radius_list_um - est_r));
[~, vessel_ori_ind] = min(abs(psf_est_int.vessel_ori_z_comp - est_z));
num_psf = numel(psf_est_int.psf_FWHM);
psf_array_size = size(psf_est_int.psf_FWHM);
%% Algorithm 1
% [psf_corr_mat_xy, psf_corr_mat_z, psf_corr_mat] = deal(nan(size(psf_est_int.psf_FWHM)));
% psf_corr_max = 0;
% 
% psf_xcorr_max_profile_xy = [];
% psf_xcorr_max_xy = 0;
% psf_xcorr_max_ind_xy = 0;
% 
% psf_xcorr_max_profile_z = [];
% psf_xcorr_max_z = 0;
% psf_xcorr_max_ind_z = 0;
% 
% tic
% for iter_psf = 1 : num_psf
%     % xy-direction
%     tmp_psf_int_r_dist =  psf_est_int.edge_int_str{iter_psf}.normalized_radial_int_dist{vessel_ori_ind, vessel_radius_ind};
%     tmp_psf_int_r_dist_interp = griddedInterpolant(tmp_psf_int_r_dist(:, 1), tmp_psf_int_r_dist(:, 2), 'linear', 'nearest');
%     tmp_int_profile = tmp_psf_int_r_dist_interp(exp_profile_xy.r_abs);
%     
%     tmp_int_xcorr = xcorr(tmp_int_profile, exp_profile_xy.int_val_n_smooth, 'normalized');
%     [tmp_max_xcorr_xy, tmp_max_ind_xy] = max(tmp_int_xcorr);
%     psf_corr_mat_xy(iter_psf) = tmp_max_xcorr_xy;
%     
%     % z_direction
%     tmp_psf_int_r_dist_vrt = psf_est_int.edge_int_str{iter_psf}.normalized_vertical_radial_int_dist{vessel_ori_ind, vessel_radius_ind};
%     tmp_psf_int_r_dist_vrt_interp = griddedInterpolant(tmp_psf_int_r_dist_vrt(:, 1), ...
%         tmp_psf_int_r_dist_vrt(:, 2), 'linear', 'nearest');
%     tmp_int_vrt_profile = tmp_psf_int_r_dist_vrt_interp(exp_profile_z.r_abs);
%     
%     tmp_int_vrt_xcorr = xcorr(tmp_int_vrt_profile, exp_profile_z.int_val_n_smooth, 'normalized');
%     [tmp_max_xcorr_z, tmp_max_ind_z] = max(tmp_int_vrt_xcorr);
%     psf_corr_mat_z(iter_psf) = tmp_max_xcorr_z;
%     
%     if est_z <= sin(pi/4)
% %         psf_corr_mat(iter_psf) = (tmp_max_xcorr_xy + tmp_max_xcorr_z)/2;
%         psf_corr_mat(iter_psf) = (tmp_max_xcorr_xy + tmp_max_xcorr_z * est_z)/ (1 + est_z);
%     else
%         psf_corr_mat(iter_psf) = tmp_max_xcorr_xy;
%     end
%     
%     if psf_corr_mat(iter_psf) > psf_corr_max
%         psf_corr_max = psf_corr_mat(iter_psf);
%         
%         psf_xcorr_max_xy = tmp_max_xcorr_xy;
%         psf_xcorr_max_profile_xy = tmp_int_profile;
%         psf_xcorr_max_ind_xy = tmp_max_ind_xy;
%         
%         psf_xcorr_max_z = tmp_max_xcorr_z;
%         psf_xcorr_max_profile_z = tmp_int_vrt_profile;
%         psf_xcorr_max_ind_z = tmp_max_ind_z;
%     end    
% end
% % Select the best fit
% [psf_corr_mat_max, psf_corr_max_ind] = max(psf_corr_mat(:));
% [psf_corr_max_sub_1, psf_corr_max_sub_2] = ind2sub(size(psf_corr_mat), psf_corr_max_ind);
% toc
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


psf_int_dist_xy = nan(numel(psf_int_xy), num_psf);
psf_int_dist_z = nan(numel(psf_int_z), num_psf);
% Flip the indices for correlation here
psf_int_xy = flip(psf_int_xy);
psf_int_z = flip(psf_int_z);
for iter_psf = 1 : num_psf
    tmp_psf_int_r_dist_interp =  psf_est_int.edge_int_str{iter_psf}.normalized_radial_int_dist_itp{vessel_ori_ind, vessel_radius_ind};
%     psf_int_dist_xy(:, iter_psf) = tmp_psf_int_r_dist_interp(exp_profile_xy.r_abs);
    psf_int_dist_xy(:, iter_psf) = tmp_psf_int_r_dist_interp(abs(psf_int_xy));
    
    tmp_psf_int_r_dist_vrt_interp = psf_est_int.edge_int_str{iter_psf}.normalized_vertical_radial_int_dist_itp{vessel_ori_ind, vessel_radius_ind};
%     psf_int_dist_z(:, iter_psf) = tmp_psf_int_r_dist_vrt_interp(exp_profile_z.r_abs);    
    psf_int_dist_z(:, iter_psf) = tmp_psf_int_r_dist_vrt_interp(abs(psf_int_z));
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
if est_z <= 1
    corr_z_weight = 1 - est_z;
    psf_corr_mat = (corr_xy_psf_max + corr_z_psf_max * corr_z_weight) ./ (1 + corr_z_weight);
else
    psf_corr_mat = corr_xy_psf_max;
end
[psf_corr_mat_max, psf_corr_max_ind] = max(psf_corr_mat);

psf_corr_mat_xy = reshape(corr_xy_psf_max, psf_array_size);
psf_xcorr_max_xy = corr_xy_psf_max(psf_corr_max_ind);
psf_xcorr_max_profile_xy = psf_int_dist_xy(:, psf_corr_max_ind);
psf_xcorr_max_ind_xy = corr_xy_psf_max_ind(psf_corr_max_ind);
psf_xcorr_max_ind_xy = psf_xcorr_max_ind_xy - psf_int_xy_initial_shift_ind;

psf_corr_mat_z = reshape(corr_z_psf_max, psf_array_size);
psf_xcorr_max_z = corr_z_psf_max(psf_corr_max_ind);
psf_xcorr_max_profile_z = psf_int_dist_z(:, psf_corr_max_ind);
psf_xcorr_max_ind_z = corr_z_psf_max_ind(psf_corr_max_ind);
psf_xcorr_max_ind_z = psf_xcorr_max_ind_z - psf_int_z_initial_shift_ind;

[psf_corr_max_sub_1, psf_corr_max_sub_2] = ind2sub(psf_array_size, psf_corr_max_ind);
%% Save data
est_psf_int_interpolation = psf_est_int.edge_int_str{psf_corr_max_sub_1, psf_corr_max_sub_2}.n_min_edge_int_interpolation.n_min_edge_int;

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

est_psf_str.corr_mat = reshape(psf_corr_mat, psf_array_size);
est_psf_str.corr_max = psf_corr_mat_max;

est_psf_str.best_fit_interpolation = est_psf_int_interpolation;
est_psf_str.best_fit_FWHM = psf_est_int.psf_FWHM{psf_corr_max_sub_1, psf_corr_max_sub_2};
est_psf_str.best_fit_FWHM_x_ind = psf_corr_max_sub_1;
est_psf_str.best_fit_FWHM_z_ind = psf_corr_max_sub_2;
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