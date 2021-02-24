function est_psf_str = fun_radius_estimation_get_best_fit_psf_xy(psf_est_int, ...
    est_r, est_z, est_psf_FWHM, exp_profile_xy, fix_est_dir)

est_psf_str = struct;

[~, vessel_radius_ind] = min(abs(psf_est_int.vessel_radius_list_um - est_r));
[~, vessel_ori_ind] = min(abs(psf_est_int.vessel_ori_z_comp - abs(est_z)));
fprintf('Closest simulated vessels: radius %f um, orientation vector z-component %f\n', ...
    psf_est_int.vessel_radius_list_um(vessel_radius_ind), ...
    sin(psf_est_int.vessel_elevation_agl_rad(vessel_ori_ind)));
if ~isempty(fix_est_dir)
    % Find the closest psf parameters
    switch fix_est_dir
        case {1, 2, [1,2]}
            fix_psf_FWHM = est_psf_FWHM(1);
            assert(~isnan(fix_psf_FWHM) && fix_psf_FWHM == est_psf_FWHM(2));
            [~, simulated_PSF_para_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - fix_psf_FWHM));
            num_psf = numel(psf_est_int.psf_FWHM_z_list);
            check_psf_ind = sub2ind(size(psf_est_int.psf_FWHM), repelem(simulated_PSF_para_ind, num_psf), ...
                1 : num_psf);
            psf_corr_mat = nan(num_psf, 1);       
            fprintf('Fixing the PSF FWHM on xy direction.\n');
        case 3
            fix_psf_FWHM = est_psf_FWHM(fix_est_dir);
            [~, simulated_PSF_para_ind] = min(abs(psf_est_int.psf_FWHM_z_list - fix_psf_FWHM));
            num_xy_psf = numel(psf_est_int.psf_FWHM_xy_list);
            check_psf_ind = sub2ind(size(psf_est_int.psf_FWHM), 1 : num_xy_psf, ...
                repelem(simulated_PSF_para_ind, num_xy_psf));
            psf_corr_mat = nan(num_xy_psf, 1);
            fprintf('Fixing the PSF FWHM on z direction.\n');
    end    
end    
psf_xcorr_max_profile = [];
psf_xcorr_max = 0;
psf_xcorr_max_pxl_ind = 0;
psf_xcorr_max_grid_ind = 0;

num_psf = numel(check_psf_ind);
for iter_psf = 1 : num_psf
    tmp_ind = check_psf_ind(iter_psf);
    tmp_psf_int_r_dist =  psf_est_int.edge_int_str{tmp_ind}.normalized_radial_int_dist{vessel_ori_ind, vessel_radius_ind};
    
    tmp_psf_int_r_dist_interp = griddedInterpolant(tmp_psf_int_r_dist(:, 1), tmp_psf_int_r_dist(:, 2), 'linear', 'nearest');
    tmp_int_profile = tmp_psf_int_r_dist_interp(exp_profile_xy.r_abs);
    
    tmp_int_xcorr = xcorr(exp_profile_xy.int_val_n_smooth, tmp_int_profile, 'normalized');
    [tmp_max_xcorr, tmp_max_ind] = max(tmp_int_xcorr);
    psf_corr_mat(iter_psf) = tmp_max_xcorr;
    
    if tmp_max_xcorr > psf_xcorr_max
        psf_xcorr_max = tmp_max_xcorr;
        psf_xcorr_max_profile = tmp_int_profile;
        psf_xcorr_max_pxl_ind = tmp_max_ind;
        psf_xcorr_max_grid_ind = iter_psf;
    end
end

switch fix_est_dir
    case {1, 2, [1,2]}
        psf_corr_max_sub_1 = simulated_PSF_para_ind;
        psf_corr_max_sub_2 = psf_xcorr_max_grid_ind;        
    case 3
        psf_corr_max_sub_1 = psf_xcorr_max_grid_ind;
        psf_corr_max_sub_2 = simulated_PSF_para_ind;
    otherwise 
        error('Unrecognized');
end

fprintf('PSF with max intensity correlation: (%f, %f, %f).\nMax correlation value: %f\n', ...
    psf_est_int.psf_FWHM{psf_corr_max_sub_1, psf_corr_max_sub_2}, psf_xcorr_max);

est_psf_int_interpolation = psf_est_int.edge_int_str{psf_corr_max_sub_1, psf_corr_max_sub_2}.n_min_edge_int_interpolation.n_min_edge_int;
%% Save data
est_psf_str.corr_mat = psf_corr_mat;
est_psf_str.corr_max = psf_xcorr_max;
est_psf_str.corr_max_ind = psf_xcorr_max_pxl_ind;
est_psf_str.max_corr_int_profile = psf_xcorr_max_profile;
est_psf_str.max_corr_shift_ind = (psf_xcorr_max_pxl_ind - numel(exp_profile_xy.r_abs));
est_psf_str.max_corr_shift_um = est_psf_str.max_corr_shift_ind * abs(exp_profile_xy.r_abs(2) - exp_profile_xy.r_abs(1));

est_psf_str.best_fit_interpolation = est_psf_int_interpolation;
est_psf_str.best_fit_FWHM = psf_est_int.psf_FWHM{psf_corr_max_sub_1, psf_corr_max_sub_2};
est_psf_str.best_fit_FWHM_x_ind = psf_corr_max_sub_1;
est_psf_str.best_fit_FWHM_z_ind = psf_corr_max_sub_2;
%% Visualization
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% scatter(ax_hdl, exp_profile_xy.r, psf_xcorr_max_profile);
% hold(ax_hdl, 'on');
% scatter(ax_hdl, exp_profile_xy.r, exp_profile_xy.int_val_n_smooth);
end