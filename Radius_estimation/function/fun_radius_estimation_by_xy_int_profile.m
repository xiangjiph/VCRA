function vsl_r_est_0 = fun_radius_estimation_by_xy_int_profile(xy_int_prof, ...
    vsl_r_est_0, ori_vec_z, est_psf_int_interpolation, max_num_est)

not_converge_Q = true;
vis_Q = false;
ori_vec_z = abs(ori_vec_z);

scalar_threshold_Q = isscalar(est_psf_int_interpolation) && isnumeric(est_psf_int_interpolation);
if scalar_threshold_Q
    assert(est_psf_int_interpolation >= 0 && est_psf_int_interpolation <= 1);
    est_count = max_num_est - 1;
else
    est_count = 0;
end

while not_converge_Q && est_count < max_num_est
    last_est = vsl_r_est_0;
    est_count = est_count + 1;
    if scalar_threshold_Q
        tmp_int_n_th = est_psf_int_interpolation;
    else
        tmp_int_n_th = est_psf_int_interpolation(ori_vec_z, vsl_r_est_0);
        if tmp_int_n_th <= 0.05 || ~isfinite(tmp_int_n_th)
            vsl_r_est_0 = nan;
            break;
        end
    end
    tmp_solu = fun_find_zeros_by_linear_interpolation(xy_int_prof.r, ...
        xy_int_prof.int_val_n_smooth - tmp_int_n_th);
    tmp_num_solu = numel(tmp_solu);
    switch tmp_num_solu
        case 2
            vsl_r_est_0 = (tmp_solu(2) - tmp_solu(1))/2;
        case {0, 1}
%             warning('Only one side of the vessel edge was founded');
%             vsl_r_est_0 = abs(tmp_solu);
            break;
        otherwise
%         warning('The number of solution is more than 2');
            vsl_r_est_0 = (tmp_solu(end) - tmp_solu(1)) / 2;
    end    
    if abs(last_est - vsl_r_est_0) < 1e-3
        not_converge_Q = false;
    end        
end
if vis_Q
%%
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    plot(ax_hdl, xy_int_prof.r, xy_int_prof.int_val_n_smooth);
end
end