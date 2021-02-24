set_env
DataManager = FileManager;
medfilt_tile_image_Q = true;
simulation_file_fp = fullfile(DataManager.fp_metadata_file('Vessel_radius_calibration', 'In_vivo_image_simulation', 'In_vivo_image_simulation'));
invivo_im_prof_str = DataManager.load_data(simulation_file_fp);
%% Preprocess simulation result
num_parameter_combination = numel(invivo_im_prof_str.edge_int_str);
for iter_p = 1 : num_parameter_combination
    tmp_psf = invivo_im_prof_str.edge_int_str{iter_p};
    tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
    tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
    tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
    
    for iter_dist = 1 : tmp_num_dist
        tmp_int_xy_dist = tmp_psf.max_normalized_radial_int_dist{iter_dist};
        if ~isempty(tmp_int_xy_dist)
            tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
                tmp_int_xy_dist(:, 2), 'linear', 'nearest');
            
            tmp_int_z_dist = tmp_psf.max_normalized_vertical_radial_int_dist{iter_dist};
            tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
                tmp_int_z_dist(:, 2), 'linear', 'nearest');
        end
    end
    invivo_im_prof_str.edge_int_str{iter_p} = tmp_psf;
end
fprintf('Finish pre-computing the intensity interpolation\n');
%% Profile analysis
% Is 50% for capillary good? Assume completely no fluorescence signal from
% RBC? - Plot 2 um, 0 degree,
vis_r_um = 2.5;
vis_z = 0;
[~, vis_r_idx] = min(abs(invivo_im_prof_str.vsl_rds_list_um - vis_r_um));
[~, vis_z_idx] = min(abs(invivo_im_prof_str.vsl_ori_z_comp - vis_z));
vis_plot_profs = cell(size(invivo_im_prof_str.edge_int_str));
for iter_profs = 1 : numel(vis_plot_profs)
    tmp_str = invivo_im_prof_str.edge_int_str{iter_profs};
    tmp_mat = tmp_str.max_normalized_radial_int_dist{vis_z_idx, vis_r_idx};
    vis_plot_profs{iter_profs} = tmp_mat;
end
%%
vis_bc_lyr_t_list = [0.2 : 0.2 : 0.8];

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_vis = 1 : numel(vis_bc_lyr_t_list)
    vis_bc_lyr_t_um = vis_bc_lyr_t_list(iter_vis);
    vis_rbc_int_n = 0.2;
    [~, vis_bc_lyr_idx] = min(abs(invivo_im_prof_str.vsl_bc_lyr_thkn_um_list - ...
        vis_bc_lyr_t_um));
    vis_bc_lyr_t_um = invivo_im_prof_str.vsl_bc_lyr_thkn_um_list(vis_bc_lyr_idx);
    [~, vis_rbc_int_idx] = min(abs(invivo_im_prof_str.rbc_int_n_list - ...
        vis_rbc_int_n));
    vis_rbc_int_n = invivo_im_prof_str.rbc_int_n_list(vis_rbc_int_idx);
    tmp_plot_vecs = vis_plot_profs{vis_rbc_int_idx, vis_bc_lyr_idx};
    
    plot(ax_hdl, tmp_plot_vecs(:, 1), tmp_plot_vecs(:, 2), 'LineWidth', 1);
    hold(ax_hdl, 'on');
end
ax_hdl.XLabel.String = 'Horizontal radial direction (\mum)';
ax_hdl.YLabel.String = 'Normalized intensity';
leg_hdl = legend(ax_hdl, arrayfun(@(x) num2str(x, '%.2f'), vis_bc_lyr_t_list, 'UniformOutput', false));
leg_hdl.Title.String = 'BLT (\mum)';
%% Observations
% 1. Normalized edge intensity can be greater than 0.5 if the RBC tube
% intensity is too weak (e.g. 0). In that case, thresholding the local
% image profile by a single threshold would be a problem - create a hole
% 2. So, it seems that the estimation should be based on correlation of the
% entire intensity profile.