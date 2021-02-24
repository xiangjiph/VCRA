clc;clear;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
DataManager = FileManager;
% psf_est_int = DataManager.load_data(DataManager.fp_metadata_file('WholeBrain', 'ML20190124', 'psf_scan_data_v2'));
% psf_est_int = load('psf_fitting_data_v2.mat');
psf_est_int = DataManager.load_data(DataManager.fp_metadata_file(dataset_name, stack, 'psf_fitting_data'));
im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    stack), 'PSF_simulation_analysis');
%% Find the PSF that are closest to estimation 
psf_FWHM_array = cat(1, psf_est_int.psf_FWHM{:});
initial_est_PSF_FWHM_xy = 0.8;
initial_est_PSF_FWHM_z = 7;
[~, initial_est_PSF_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - initial_est_PSF_FWHM_xy));
[~, initial_est_PSF_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - initial_est_PSF_FWHM_z));
initial_est_PSF_list_ind = sub2ind(size(psf_est_int.psf_FWHM), initial_est_PSF_FWHM_xy_ind, initial_est_PSF_FWHM_z_ind);
num_psf = numel(psf_est_int.psf_FWHM);
%% Plot: Normalized edge intensity vs orientation and radius 
best_fit_psf = DataManager.load_data(DataManager.fp_metadata_file('WholeBrain', ...
    'ML_2018_08_15', 'Estimated_PSF_edge_intensity'));
psf_est_stat = DataManager.load_data(DataManager.fp_metadata_file(dataset_name, ...
    stack, 'PSF_estimation'));
% best_fit_psf = psf_est_int.edge_int_str{initial_est_PSF_list_ind};
[int_x, int_y] = ndgrid(best_fit_psf.vsl_theta_list, best_fit_psf.vsl_r_list);
tmp_int_itp = griddedInterpolant(int_x, int_y, best_fit_psf.normalized_min_edge_int);

tmp_plot_theta_list_deg = [0, 30, 45, 60, 75, 90];
tmp_plot_theta_list = tmp_plot_theta_list_deg / 180 * pi;

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_theta = 1 : numel(tmp_plot_theta_list_deg)
    tmp_theta = tmp_plot_theta_list(iter_theta);
    tmp_r =  best_fit_psf.vsl_r_list(3:end);
    tmp_theta = repelem(tmp_theta, numel(tmp_r), 1);    
    plot(ax_hdl, tmp_r, tmp_int_itp(tmp_theta, tmp_r'),...
        'LineWidth', 2, 'Marker', 'o');
    hold(ax_hdl, 'on');
end
leg_hdl = legend(arrayfun(@(x) num2str(x, '%d'), tmp_plot_theta_list_deg, 'UniformOutput', false), 'Location', 'southwest');
leg_hdl.Title.String = '\theta (\circ)';
ax_hdl.XLabel.String = 'Vessel radius (\mum)';
ax_hdl.YLabel.String = 'Normalized edge intensity';
ax_hdl.YLim(1) = 0;
ax_hdl.YLim(2) = 0.6;
ax_hdl.XScale = 'log';
ax_hdl.FontSize = 12;
ax_hdl.Title.String = sprintf('PSF FWHM (%.2f, %.2f, %.2f) \\mum',  ...
    best_fit_psf.psf_FWHM_um);
box(ax_hdl, 'off');
fig_fp = fullfile(im_save_folder, sprintf('Best_fit_PSF_Normalized_edge_int_vs_vessel_radius.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Intensity profile for different FWHM_xy
vis_r = 2;
vis_theta = pi/inf;
vis_psf_FWHM_xy_list = [0.47, 0.82, 1.06];
vis_psf_FWHM_z_list = 6.8;
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_z = 1 : numel(vis_psf_FWHM_xy_list)
    vis_psf_FWHM_xy = vis_psf_FWHM_xy_list(iter_z);
    vis_psf_FWHM_z_um = vis_psf_FWHM_z_list;
    % Find the closest simulated profile
    [~, tmp_psf_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - vis_psf_FWHM_xy));
    [~, tmp_psf_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - vis_psf_FWHM_z_um));
    tmp_psf_data = psf_est_int.edge_int_str{tmp_psf_xy_ind, tmp_psf_z_ind};
    
    [~, tmp_r_ind] = min(abs(tmp_psf_data.vsl_r_list - vis_r));
    [~, tmp_theta_ind] = min(abs(tmp_psf_data.vsl_theta_list - vis_theta));
    tmp_int_prof = tmp_psf_data.normalized_radial_int_dist{tmp_theta_ind, tmp_r_ind};
    plot(ax_hdl, tmp_int_prof(:, 1), tmp_int_prof(:, 2), 'LineWidth', 1.5);
    hold(ax_hdl, 'on');
end
line_hdl = line(ax_hdl, [vis_r, vis_r], [0, 1.2], 'LineWidth', 2, 'LineStyle', '-.', 'Color', 'k');

ax_hdl.YLim(2) = 1.1;
ax_hdl.XLim(2) = 3;
ax_hdl.XLabel.String = 'Radial distance (\mum)';
ax_hdl.YLabel.String = 'Normalized intensity';
ax_hdl.Title.String = sprintf('Lateral radial intensity profile\n(r = %.1f \\mum, \\theta = %d\\circ, FWHM_z = %.2f\\mum)', ...
    vis_r, round(vis_theta * 180/pi), vis_psf_FWHM_z_um);
leg_hdl = legend(arrayfun(@num2str, vis_psf_FWHM_xy_list, 'UniformOutput', false));
leg_hdl.Title.String = 'FWHM_{xy}(\mum)';
box(ax_hdl, 'off');
ax_hdl.FontSize = 12;

fig_fp = fullfile(im_save_folder, sprintf('Normalized_int_prof_different_FWHM_xy_FWHM_z_%d_nm_r_%d_nm_theta_%d_deg.png', ...
    round(vis_psf_FWHM_z_um * 1000), round(vis_r * 1000), round(vis_theta/pi * 180)));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Intensity profile for different FWHM_z
vis_r = 2;
vis_theta = pi/inf;
vis_psf_FWHM_xy = 0.8;
vis_psf_FWHM_z_list = [1.65 3.06 4.95 8.00];
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for iter_z = 1 : numel(vis_psf_FWHM_z_list)
    vis_psf_FWHM_z_um = vis_psf_FWHM_z_list(iter_z);
    % Find the closest simulated profile
    [~, tmp_psf_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - vis_psf_FWHM_xy));
    [~, tmp_psf_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - vis_psf_FWHM_z_um));
    tmp_psf_data = psf_est_int.edge_int_str{tmp_psf_xy_ind, tmp_psf_z_ind};
    
    [~, tmp_r_ind] = min(abs(tmp_psf_data.vsl_r_list - vis_r));
    [~, tmp_theta_ind] = min(abs(tmp_psf_data.vsl_theta_list - vis_theta));
    tmp_int_prof = tmp_psf_data.normalized_radial_int_dist{tmp_theta_ind, tmp_r_ind};
    plot(ax_hdl, tmp_int_prof(:, 1), tmp_int_prof(:, 2), 'LineWidth', 1.5);
    hold(ax_hdl, 'on');
end
line_hdl = line(ax_hdl, [vis_r, vis_r], [0, 1.2], 'LineWidth', 2, 'LineStyle', '-.', 'Color', 'k');

ax_hdl.YLim(2) = 1.1;
ax_hdl.XLim(2) = 3;
ax_hdl.XLabel.String = 'Radial distance (\mum)';
ax_hdl.YLabel.String = 'Normalized intensity';
leg_hdl = legend(arrayfun(@(x) num2str(x, '%.2f'), vis_psf_FWHM_z_list, 'UniformOutput', false));
leg_hdl.Title.String = 'FWHM_z(\mum)';
ax_hdl.Title.String = sprintf('Lateral radial intensity profile\n(r = %.1f \\mum, \\theta = %d\\circ, FWHM_{xy} = %.2f\\mum)', ...
    vis_r, (vis_theta * 180/pi), vis_psf_FWHM_xy);
box(ax_hdl, 'off');
ax_hdl.FontSize = 12;

fig_fp = fullfile(im_save_folder, sprintf('Normalized_int_prof_different_FWHM_z_FWHM_XY_%d_nm_r_%d_nm_theta_%d_deg.png', ...
    round(vis_psf_FWHM_xy * 1000), round(vis_r * 1000), round(vis_theta/pi * 180)));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Error from incorrect PSF estimation
psf_est_2018 = load('/data/Vessel/WholeBrain/ML_2018_08_15/processed_data/metadata/PSF_estimation.mat');
psf_est_2019 = load('/data/Vessel/WholeBrain/ML20190124/processed_data/metadata/PSF_estimation.mat');

psf_edge_data_2018 = load('/data/Vessel/WholeBrain/ML_2018_08_15/processed_data/metadata/Estimated_PSF_edge_intensity.mat');
psf_edge_data_2019 = load('/data/Vessel/WholeBrain/ML20190124/processed_data/metadata/Estimated_PSF_edge_intensity.mat');

% Compare FWHM_z 7, 8, 9
gt_theta_list = asin(psf_est_int.vessel_ori_z_comp([1, 8, 11, 15, 19, 20, 21]));
num_theta = numel(gt_theta_list);
[abs_error_cell, relative_error_cell] = deal(cell(num_theta, 1));

gt_FWHZ_z = psf_est_2018.avg_PSF_FWHM_z.mean;
gt_FWHZ_xy = psf_est_2018.avg_PSF_FWHM_xy.mean;
% [~, gt_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - gt_FWHZ_z));
% [~, gt_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - gt_FWHZ_xy));
% tmp_gt_data = psf_est_int.edge_int_str{gt_FWHM_xy_ind, gt_FWHM_z_ind};

est_FWHZ_xy = psf_est_2019.avg_PSF_FWHM_xy.mean;
est_FWHZ_z = psf_est_2019.avg_PSF_FWHM_z.mean;
[~, est_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - est_FWHZ_xy));
[~, est_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - est_FWHZ_z));
% tmp_est_data = psf_est_int.edge_int_str{est_FWHM_xy_ind, est_FWHM_z_ind};
%
tmp_gt_data = psf_edge_data_2018;
tmp_est_data = psf_edge_data_2019;

gt_r_list = psf_est_int.vessel_radius_list_um(1:1:end);
for iter_theta = 1 : num_theta
    gt_theta = gt_theta_list(iter_theta);    
    [abs_error, relative_error] = deal(nan(size(gt_r_list)));
    for iter_r = 1 : numel(gt_r_list)
        tmp_r = gt_r_list(iter_r);

        [~, tmp_r_ind] = min(abs(tmp_gt_data.vsl_r_list - tmp_r));
        [~, tmp_theta_ind] = min(abs(tmp_gt_data.vsl_theta_list - gt_theta));
        tmp_gt_int_prof = tmp_gt_data.normalized_radial_int_dist{tmp_theta_ind, tmp_r_ind};
        % Find the first 0
        tmp_first_zero_ind = find(tmp_gt_int_prof(:, 2) == 0, 1, 'first');
        tmp_first_1_ind = find(tmp_gt_int_prof(:, 2) == 1, 1, 'last');
        
        tmp_gt_int_prof = tmp_gt_int_prof(tmp_first_zero_ind : -1 : tmp_first_1_ind, :);
        tmp_int_2_r_itp = griddedInterpolant(tmp_gt_int_prof(:, 2), tmp_gt_int_prof(:, 1));
        
        % Get the normalized edge intensity from the estiamted PSF
        
        tmp_est_n_int = tmp_est_data.n_min_edge_int_interpolation.n_min_edge_int(...
            sin(gt_theta), tmp_r);
        tmp_est_r = tmp_int_2_r_itp(tmp_est_n_int);
        
        abs_error(iter_r) = tmp_est_r - tmp_r;
        relative_error(iter_r) = abs_error(iter_r) / tmp_r;
    end
    abs_error_cell{iter_theta} = abs_error;
    relative_error_cell{iter_theta} = relative_error;
end

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 2];
ax_hdl = subplot(1,2,2);
plot(ax_hdl, repmat(gt_r_list, num_theta, 1)', cat(1, relative_error_cell{:})', 'LineWidth', 2);
ax_hdl.XLabel.String = 'Vessel radius (\mum)';
ax_hdl.YLabel.String = 'Relative error';
leg_hdl = legend(arrayfun(@(x) num2str(x, '%d'), round(gt_theta_list * 180 / pi), 'UniformOutput', false), 'Location', 'southeast');
leg_hdl.Title.String = 'Elevation angle (\circ)';
ax_hdl.Title.String = sprintf('GT PSF FWHM (%.2f, %.2f, %.2f) \\mum\nEST PSF FWHM (%.2f, %.2f, %.2f) \\mum', ...
    tmp_gt_data.psf_FWHM_um, tmp_est_data.psf_FWHM_um);

ax_hdl_2 = subplot(1,2,1);
plot(ax_hdl_2, repmat(gt_r_list, num_theta, 1)', cat(1, abs_error_cell{:})', 'LineWidth', 2);
ax_hdl_2.XLabel.String = 'Vessel radius (\mum)';
ax_hdl_2.YLabel.String = 'Absolute error (\mum)';
leg_hdl = legend(arrayfun(@(x) num2str(x, '%d'), round(gt_theta_list * 180 / pi), 'UniformOutput', false), 'Location', 'southeast');
leg_hdl.Title.String = 'Elevation angle (\circ)';
ax_hdl_2.YLim(1) = -0.1;

ax_hdl_2.Title.String = sprintf('GT PSF FWHM (%.2f, %.2f, %.2f) \\mum\nEST PSF FWHM (%.2f, %.2f, %.2f) \\mum', ...
    tmp_gt_data.psf_FWHM_um, tmp_est_data.psf_FWHM_um);

% fig_fp = fullfile(im_save_folder, sprintf('Radius_estimation_error_from_PSF_est_gt_20180815_est_20190124.png'));
fig_fp = fullfile(im_save_folder, sprintf('Radius_estimation_error_from_PSF_est_gt_20190124_est_20180815.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Fixed Theta, Radius error improvement for vessel of different radius
psf_str = psf_edge_data_2018;
vsl_theta = psf_str.vsl_theta_list(1);
radius_gt_list = [1, 2, 4, 8];
num_gt = length(radius_gt_list);
[error_result, error_initial, est_r_initial, est_r_result] = deal(cell(num_gt,1));
[~, theta_idx] = min(abs(psf_str.vsl_theta_list - vsl_theta));
for tmp_gt_idx = 1 : num_gt
    radius_gt = radius_gt_list(tmp_gt_idx);
    [~, radius_gt_idx] = min(abs(psf_str.vsl_r_list - radius_gt));
    
    radius_est_list = psf_str.vsl_r_list(psf_str.vsl_r_list >= 0.05 * radius_gt & ...
        psf_str.vsl_r_list <= 2 * radius_gt);
    num_est_list = numel(radius_est_list);
    [tmp_error_result, tmp_error_initial, ...
        tmp_est_initial, tmp_est_result] = deal(zeros(num_est_list,1)); 
    int_profile = psf_str.normalized_radial_int_dist{theta_idx, radius_gt_idx};
    
    tmp_first_zero_ind = find(int_profile(:, 2) == 0, 1, 'first');
    tmp_first_1_ind = find(int_profile(:, 2) == 1, 1, 'last');
    
    tmp_gt_int_prof = int_profile(tmp_first_zero_ind : -1 : tmp_first_1_ind, :);
    tmp_int_2_r_itp = griddedInterpolant(tmp_gt_int_prof(:, 2), tmp_gt_int_prof(:, 1));

    int_dist = int_profile(:,2);
    int_r = int_profile(:,1);
    [int_dist, tmp_unique_idx,~] = unique(int_dist);
    for tmp_idx = 1 : num_est_list
        radius_est_0 = radius_est_list(tmp_idx);
        
        int_th_est = psf_str.n_min_edge_int_interpolation.n_min_edge_int(sin(vsl_theta), radius_est_0);
        radius_est_1 = tmp_int_2_r_itp(int_th_est);
        tmp_est_initial(tmp_idx) = radius_est_0;
        tmp_est_result(tmp_idx) = radius_est_1;
        
        tmp_error_result(tmp_idx) = (radius_est_1 - radius_gt) / radius_gt;
        tmp_error_initial(tmp_idx) = (radius_est_0 - radius_gt) / radius_gt;
    end
    error_result{tmp_gt_idx} = tmp_error_result;
    error_initial{tmp_gt_idx} = tmp_error_initial;
    est_r_initial{tmp_gt_idx} = tmp_est_initial;
    est_r_result{tmp_gt_idx} = tmp_est_result;
end

fig_hdl = figure;
ax_hdl = axes(fig_hdl);
for tmp_gt_idx = 1 : num_gt
    plot(ax_hdl, error_initial{tmp_gt_idx}, error_result{tmp_gt_idx}, 'LineWidth', 2);
    hold(ax_hdl, 'on');
end
ax_hdl.FontSize = 13;
ax_hdl.XLabel.String = 'Initial relative error';
ax_hdl.YLabel.String = 'Estimation relative error';
% ax_hdl.Title.String = sprintf('\\theta = %d\\circ', round(180 * vsl_theta / pi));

ax_hdl.Title.String = sprintf('PSF FWHM (%.2f, %.2f, %.2f) \\mum, \\theta = %d\\circ',  ...
    best_fit_psf.psf_FWHM_um, round(180 * vsl_theta / pi));

tmp_lgd = legend(ax_hdl, cellstr(num2str(radius_gt_list', '%.1f')), 'Location', 'southeast');
tmp_lgd.Title.String = 'Radius (\mum)';
grid(ax_hdl, 'on');
% box(ax_hdl, 'off');
ax_hdl.XLim = [-1, 1];
% fig_fp = fullfile(im_save_folder, sprintf('Relative_error_before_vs_after_estimation_theta_%d_deg.png', round(180 * vsl_theta / pi)));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Error from wrong estimation of the tilt angle, fix the radius
% close all;
% plot_data = psf_simulation{1};
% radius_gt = 5;
% image_file_name = sprintf('Relative_error_vs_tilt_angle_error_psf%.3d%.3d_radius%d.png',round(plot_data.psf_size(1)*100), ...
%     round(plot_data.psf_size(3)*100), round(radius_gt*10));
% [~, radius_gt_idx] = min(abs(plot_data.vessel_radius_list - radius_gt));
% theta_gt_list = plot_data.theta_list(1:6:end);
% num_gt = length(theta_gt_list);
% error_r = cell(num_gt,1);
% error_theta = cell(num_gt,1);
% for tmp_gt_idx = 1 : num_gt
%     theta_gt = theta_gt_list(tmp_gt_idx);
%     [~, theta_gt_idx] = min(abs(plot_data.theta_list - theta_gt));
%     theta_est_list = unique(min(plot_data.theta_list(end), max(0, theta_gt + (-45:5:45) * pi/180)));
%     num_est_list = numel(theta_est_list);
%     tmp_error_r_result = zeros(num_est_list,1);
%     tmp_error_theta = zeros(num_est_list,1);
%     for tmp_idx = 1 : num_est_list
%         theta_est = theta_est_list(tmp_idx);
%         int_th_est = interp1(plot_data.theta_list, plot_data.abs_min_edge_int(:,radius_gt_idx), theta_est);
%         int_profile = plot_data.radial_int_distribution{theta_gt_idx, radius_gt_idx};
%         int_dist = int_profile(:,2);
%         int_r = int_profile(:,1);
%         [int_dist, tmp_unique_idx,~] = unique(int_dist);
%         radius_est_1 = interp1(int_dist, int_r(tmp_unique_idx), int_th_est);
%         tmp_error_r_result(tmp_idx) = (radius_est_1 - radius_gt) / radius_gt;
%         tmp_error_theta(tmp_idx) = (theta_est - theta_gt);
%     end
%     error_r{tmp_gt_idx} = tmp_error_r_result;
%     error_theta{tmp_gt_idx} = round(tmp_error_theta * 180/pi);
% end
% 
% tmp_ax = axes();
% for tmp_gt_idx = 1 : num_gt
%     tmp_fig = plot(tmp_ax, error_theta{tmp_gt_idx}, error_r{tmp_gt_idx}, 'LineWidth', 2);
%     hold on
% end
% xlabel(tmp_ax, 'Error in tilt angle/{\circ}');
% ylabel(tmp_ax, 'Relative error in estimated r');
% % title(sprintf('Vessel tilt angle %.2d', theta));
% tmp_lgd = legend(tmp_ax, cellstr(num2str(round(theta_gt_list * 180/pi)' , '%.2d')));
% title(tmp_lgd, 'Tilt angle/{\circ}');
% 
% grid(tmp_ax, 'on')
% [tmp_fig.LineWidth] = deal(1.5);
% tmp_ax.LineWidth = 2;
% tmp_ax.FontSize = 15;
% pbaspect(tmp_ax, [1,1,1]);
% outerpos = tmp_ax.OuterPosition;
% ti = tmp_ax.TightInset;
% tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
%     outerpos(4)-ti(2)-ti(4)];
% print(DataManager.fp_analysis_image(image_file_name, save_folder), '-dpng', '-r300');
