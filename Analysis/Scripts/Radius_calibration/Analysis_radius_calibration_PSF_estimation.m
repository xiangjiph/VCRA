%% Load PSF simulation data
psf_est_int = DataManager.load_data(DataManager.fp_metadata_file('WholeBrain', 'ML_2018_08_15', 'psf_fitting_data'));
% psf_fp = fullfile(DataManager.fp_metadata_file('DKLab', 'Rui_PSF_simulation', 'psf_fitting_data'));
% psf_est_int = load(psf_fp);
psf_FWHM_array = cat(1, psf_est_int.psf_FWHM{:});
initial_est_PSF_FWHM_xy = 0.8;
initial_est_PSF_FWHM_z = 8;
[~, initial_est_PSF_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - initial_est_PSF_FWHM_xy));
[~, initial_est_PSF_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - initial_est_PSF_FWHM_z));
initial_est_PSF_list_ind = sub2ind(size(psf_est_int.psf_FWHM), initial_est_PSF_FWHM_xy_ind, initial_est_PSF_FWHM_z_ind);
num_psf = numel(psf_est_int.psf_FWHM);
for iter_psf = 1 : num_psf
    tmp_psf = psf_est_int.edge_int_str{iter_psf};
    tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
    tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
    tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
    for iter_dist = 1 : tmp_num_dist
        tmp_int_xy_dist = tmp_psf.normalized_radial_int_dist{iter_dist};
        tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
            tmp_int_xy_dist(:, 2), 'linear', 'nearest');
        
        tmp_int_z_dist = tmp_psf.normalized_vertical_radial_int_dist{iter_dist};
        tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
            tmp_int_z_dist(:, 2), 'linear', 'nearest');
    end
    psf_est_int.edge_int_str{iter_psf} = tmp_psf;
end
fprintf('Finish pre-computing the intensity interpolation\n');
%% PSF - radius joint estimation
tile_image = medfilt3(avg_data);

skeleton_cc_sub_um = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, x) .* ...
    vessel_graph.info.voxel_size_um, vessel_graph.link.cc_ind, 'UniformOutput', false);
skeleton_cc_r_um = cellfun(@(x) double(full(vessel_graph.radius(x))), ...
    vessel_graph.link.cc_ind, 'UniformOutput', false);

psf_stat = fun_radius_estimation_get_PSF_stat_in_image_stack(psf_est_int, ...
    tile_image, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um);
%% PSF estimation statistics
% Histogram of FWHM for:
% 1. Small vessel
% 2. Small z-component
tmp_vis_Q = (psf_stat.refine_r_est <= 4 & psf_stat.best_fit_corr > 0.99 & ...
    abs(psf_stat.link_ori_vec(:, 3)) < sqrt(1/2));
% tmp_vis_Q = ~isnan(best_fit_corr);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1] .* 1.5;
ax_hdl_1 = subplot(1,3,1);
tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 1);
histogram(ax_hdl_1, tmp_data, psf_est_int.psf_FWHM_xy_list - 0.05, 'Normalization', 'pdf');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
ax_hdl_2 = subplot(1,3,2);
tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 3);
histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;

ax_hdl_3 = subplot(1,3,3);
histogram(ax_hdl_3, psf_stat.refine_r_est(tmp_vis_Q), 'Normalization', 'pdf');
legend(ax_hdl_3, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(psf_stat.refine_r_est(tmp_vis_Q))), 'Location', 'northeast');
ax_hdl_3.XLabel.String = 'Estimated radius (\mum)';
ax_hdl_3.YLabel.String = 'PDF';
ax_hdl_3.FontSize = 14;

fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation', image_group, sprintf('%s_%s_%s_ROI_%d_PSF_FWHM_estimation_hists.png', ...
    dataset_name, stack, image_group, tmp_roi_id));
fun_print_image_in_several_formats(fig_hdl, fig_fp);