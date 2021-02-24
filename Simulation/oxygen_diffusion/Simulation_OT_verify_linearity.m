diffCoeffO2 = 1.9e-9; % m^2/sec in rodent brain (Clark et al 1978)
% Oxygen solubility
% 1.3e-3 mol(L * atm) in water;
% 9.8214e-4 mol/(L * atm) in small rodent brain (in vivo) Clark et al 1978

% alphaO2 = 1.3e-3 / 760; % mol/(L * mmHg) in water
alphaO2 = 9.8214e-4 / 760;
diffAlpha = diffCoeffO2 * alphaO2;
%% Load capillary network
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
load_skel_ver = '240_cube_re';
grid_info = DataManager.load_grid(dataset_name, stack, '240_cube');
mask_version = '240_cube_recon';
% gpuDevice(2);
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Tissue_pO2_simulation');
%%
vis_cube_label = 31998;
vis_cube_grid_sub = grid_info.bbox_grid_sub_list(vis_cube_label, :);
vis_recon_str = DataManager.load_block_mask(dataset_name, stack, ...
    mask_version, vis_cube_grid_sub(1), vis_cube_grid_sub(2), vis_cube_grid_sub(3));
vis_recon = fun_reconstruct_block_mask(vis_recon_str);
volumeViewer(vis_recon);
%% Test whether the reduction depends linearly on the oxygen metabolism rate
% scan_parameter_list = 78e-6 * linspace(0.5, 2, 10);
scan_parameter_list = 20 : 5 : 60;
num_simu = numel(scan_parameter_list);
% pO2_vsl = 36;
downsample_rate = 2;
simu_result = cell(num_simu, 1);
for iter_simu = 1 : num_simu
    oxygen_metabolism_rate = 78e-6;
    tmp_scan_parameter = scan_parameter_list(iter_simu);
    pO2_result = fun_simulation_OT_solve_ct_diff_itr(vis_recon, tmp_scan_parameter, oxygen_metabolism_rate, downsample_rate, []);
    simu_result{iter_simu} = pO2_result;
end
%%
pO2_n_result_ds2_nn = fun_simulation_OT_solve_ct_diff_itr_n(vis_recon, 2, -300);
pO2_n_result_ds2_nn.mask_dt = bwdist(pO2_n_result_ds2_nn.vessel_mask);
profile on 
pO2_n_result_ds2 = fun_simulation_OT_solve_ct_diff_itr_n(vis_recon, mask_dt, -300);
profile off
profile viewer
pO2_n_result_ds2.mask_dt = bwdist(pO2_n_result_ds2.vessel_mask_rz) .* 2;
pO2_n_result = fun_simulation_OT_solve_ct_diff_itr_n(vis_recon, 1, -300);
pO2_n_result.mask_dt = bwdist(pO2_n_result.vessel_mask);
%% 
internal_offset = 8;
fig_hdl = figure;
tiledlayout(fig_hdl, 'flow');
pO2_binned_str_cell = cell(num_simu, 1);
dt_bin_edge = 2 : 2 : 30;
for iter_simu = 1 : num_simu
    ax_hdl = nexttile;
    tmp_result = simu_result{iter_simu};
    tmp_tissue_mask = tmp_result.dt_array > 0;
    % Exclude boundary
    tmp_tissue_mask(1 : internal_offset, :, :) = false;
    tmp_tissue_mask((end - internal_offset + 1) : end, :, :) = false;
    tmp_tissue_mask(:, 1 : internal_offset, :) = false;
    tmp_tissue_mask(:, (end - internal_offset + 1) : end, :) = false;
    tmp_tissue_mask(:, :, 1 : internal_offset) = false;
    tmp_tissue_mask(:, :, (end - internal_offset + 1) : end) = false;
    %
    tmp_x_data = tmp_result.dt_array(tmp_tissue_mask);
    tmp_y_data = tmp_result.pO2_array(tmp_tissue_mask);
    % Binning
    [pO2_binned_str_cell{iter_simu}, ~, ~] = fun_analysis_get_y_stat_in_x_bin(tmp_x_data, tmp_y_data, dt_bin_edge);
    %
    histogram2(ax_hdl, tmp_x_data, tmp_y_data, 'DisplayStyle', 'tile');
    hold(ax_hdl, 'on');
    plot(ax_hdl, pO2_binned_str_cell{iter_simu}.x_bin_val, pO2_binned_str_cell{iter_simu}.y_mean, 'LineWidth', 2);
end
%% 
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_simu = 1 : num_simu
    plot(ax_hdl, pO2_binned_str_cell{iter_simu}.x_bin_val, (pO2_binned_str_cell{iter_simu}.y_mean - scan_parameter_list(iter_simu)), 'LineWidth', 2);
    hold(ax_hdl, 'on');    
end
leg_str_array = arrayfun(@(x) num2str(x, '%.2e'), scan_parameter_list, 'UniformOutput', false);
leg_hdl = legend(ax_hdl,  leg_str_array, 'Interpreter', 'latex');
leg_hdl.Title.String = '$pO_2(0) (mmHg)$';
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.XLabel.String = '$\rho\;(\mu m)$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$pO_2(\rho) - pO_2(0) \; (mmHg)$';
ax_hdl.FontSize = 14;

fig_fp = fullfile(vis_folder, 'verify_linear_behavior', sprintf(...
    'Linear_response_in_capillary_pO2.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
im_1 = imresize3(pO2_n_result_ds2_nn.pO2_array, size(vis_recon));
im_2 = imresize3(pO2_n_result_ds2.pO2_array, size(vis_recon));
im_3 = pO2_n_result.pO2_array;
% im_1 = pO2_n_result_ds2_nn.pO2_array;
% im_2 = pO2_n_result_ds2.pO2_array;
% im_3 = imresize3(pO2_n_result.pO2_array, 0.5);
%% Visualize the relative error w.r.t. 1um result
diff_hist_bin_edge = -1 : 0.05 : 1;
fig_hdl = figure;
ax_hdl_1 = axes(fig_hdl);
diff_nn = (im_1 - im_3) ./ im_3;
histogram(ax_hdl_1, diff_nn, diff_hist_bin_edge, 'FaceAlpha', 0.3);
ax_hdl_1.YScale = 'log';
hold(ax_hdl_1, 'on');
diff_df = (im_2 - im_3) ./ im_3;
histogram(ax_hdl_1, diff_df, diff_hist_bin_edge, 'FaceAlpha', 0.3);
legend(ax_hdl_1, 'Nearest neighbor', 'Default');
%%
vis_sec = 50;
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
ax_hdl_1 = subplot(1,3,1);
imagesc(ax_hdl_1, im_1(:, :, vis_sec));
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_1.Title.String = '1 \mum resolution';
colorbar(ax_hdl_1);
ax_hdl_2 = subplot(1,3,2);
imagesc(ax_hdl_2, im_2(:, :, vis_sec));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.Title.String = '2 \mum resolution';
colorbar(ax_hdl_2);
ax_hdl_3 = subplot(1,3,3);
imagesc(ax_hdl_3, im_1(:, :, vis_sec) - im_2(:, :, vis_sec));
colorbar(ax_hdl_3);
ax_hdl_3.Title.String = '1 \mu m resolution - 2 \mum resolution';
ax_hdl_3.DataAspectRatio = [1,1,1];
% fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
%     'Oxygen_transportation_simulation', 'Difference_in_voxel_size_3.png');
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Analyze pO2 dependence on distance
[oxy_result_n_str, tmp_x_data, tmp_y_data] = fun_simulation_OT_analyze_result_n(...
    pO2_n_result_ds2.mask_dt, pO2_n_result_ds2.pO2_array, 8, 10);
%% Compare with Krogh model
% Fit the Krogh model
r_cap = 2; %\mum
krogh_coeff = 1/2;
local_median_dt_max = median(vis_recon_str.link.features.nearest_tissue_dt_max, 'omitnan');
% local_median_dt_max = max(tmp_x_data(:));
% Determine correction factor by linear regression without intercept
fit_data_x = tmp_x_data;
% fit_data_x = tmp_y_stat.x_bin_val;
fit_y = tmp_y_data;
% fit_y = tmp_y_stat.y_mean;
fit_x_fun = @(x) (- krogh_coeff * (local_median_dt_max .^ 2 .* log((x + r_cap) ./ r_cap) - ...
    ((x + r_cap).^ 2 - r_cap^2)/2));
fit_x = fit_x_fun(fit_data_x);
lin_fit_hdl = fitlm(fit_x, fit_y, 'Intercept', false);
correction_factor = lin_fit_hdl.Coefficients.Estimate(1);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
histogram2(ax_hdl, tmp_x_data, tmp_y_data, 'DisplayStyle', 'tile');
hold(ax_hdl, 'on');
% med_hdl = plot(ax_hdl, tmp_y_stat.x_bin_val, tmp_y_stat.y_median, 'LineWidth', 2);
mean_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat.x_bin_val, oxy_result_n_str.pO2_stat.y_mean, 'LineWidth', 4);
ptl005_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat.x_bin_val, oxy_result_n_str.pO2_stat.y_prctile(:, 1), 'LineWidth', 4, ...
    'LineStyle', '-.');
ptl095_hdl = plot(ax_hdl, oxy_result_n_str.pO2_stat.x_bin_val, oxy_result_n_str.pO2_stat.y_prctile(:, end), 'LineWidth', 4, ...
    'LineStyle', '-.');
% correction_factor = 1;
plot_krogh_x = 0 : local_median_dt_max;
oxygen_krogh = correction_factor * fit_x_fun(plot_krogh_x);
plt_krogh = plot(ax_hdl, plot_krogh_x, oxygen_krogh, 'LineWidth', 4, 'LineStyle', '-');

oxygen_krogh_0 = fit_x_fun(plot_krogh_x);
plt_krogh_0 = plot(ax_hdl, plot_krogh_x, oxygen_krogh_0, 'LineWidth', 4, 'LineStyle', '-');

dot_dt_lm_hdl = scatter(ax_hdl, oxy_result_n_str.dt_lm.dt_v, oxy_result_n_str.dt_lm.pO2_v, 'g+');
dot_oxy_lm_hdl = scatter(ax_hdl, oxy_result_n_str.pO2_lm.dt_v, oxy_result_n_str.pO2_lm.pO2_v, 'r*');
%
leg_hdl_array = [mean_hdl, ptl005_hdl, ptl095_hdl, plt_krogh, plt_krogh_0, dot_dt_lm_hdl, dot_oxy_lm_hdl];
leg_hdl_string = {'Average', '5% percentile', '95% percentile', ...
    sprintf('Fitting Krogh model:\nCorr. factor: %.2f\nR^2-Adjusted: %.2f', correction_factor, ...
    lin_fit_hdl.Rsquared.Adjusted), 'Krogh model'...
    'Local d maximum', 'Local u_0 minimum'};
legend(ax_hdl, leg_hdl_array, leg_hdl_string, 'Location', 'eastoutside');
ax_hdl.XLabel.String = 'd';
ax_hdl.YLabel.String = 'u_0(d)';
ax_hdl.FontSize = 14;
%%
% record_array_1 = simu_array;
vis_list_ind = find(oxy_lm_dt_v < 20 & oxy_lm_v < 6);
vis_list_sub = oxy_lm_sub(vis_list_ind, :) + internal_dist;
vis_sec = 60;
fig_hdl = figure;
ax_hdl_1 = subplot(1,3,1);
imagesc(ax_hdl_1, simu_array(:, :, vis_sec));
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_2 = subplot(1,3,2);
imagesc(ax_hdl_2, record_array_1(:, :, vis_sec));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_3 = subplot(1,3,3);
imagesc(ax_hdl_3, simu_array(:, :, vis_sec) - record_array_1(:, :, vis_sec));
ax_hdl_3.DataAspectRatio = [1,1,1];