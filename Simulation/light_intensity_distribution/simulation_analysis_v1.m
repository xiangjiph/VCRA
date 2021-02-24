set_env;
save_folder = 'Intensity_simulation';
psf_simulation = cell(4,1);
psf_simulation{1} = load('edge_intensity_020206.mat');
psf_simulation{2} = load('edge_intensity_023023081.mat');
psf_simulation{3} = load('edge_intensity_027027111.mat');
psf_simulation{4} = load('edge_intensity_031031145.mat');
%% Generate interpolation
num_simu = numel(psf_simulation);

psf_est = struct;
psf_est.FWHM_um = nan(3, num_simu);
psf_est.sigma_um = nan(3, num_simu);
psf_est.interpolation = cell(num_simu, 1);

for iter_psf = 1 : numel(psf_simulation)
    tmp_str = psf_simulation{iter_psf};
    
    psf_est.FWHM_um(:, iter_psf) = tmp_str.psf_size .* (2*sqrt(2*log(2)));
    psf_est.sigma_um(:, iter_psf) = tmp_str.psf_size;
    
    tmp_int = struct;
    % x is the radius, y is the z-component in the orientation vector (0 -
    % 1). theta is the elevation angle with respect to the xy plane
    tmp_r = tmp_str.vessel_radius_list;
    tmp_z_comp = sin(tmp_str.theta_list);
    [tmp_int.z_comp, tmp_int.radius] = ndgrid(tmp_z_comp, tmp_r);
    tmp_int.variable_name = {'z_component', 'radius'};
    tmp_int.Interpolation_method = 'linear';
    tmp_int.Extrapolation_method = 'nearest';
    tmp_int.min_edge_int_n = griddedInterpolant(tmp_int.z_comp, tmp_int.radius, tmp_str.normalized_min_edge_int, ...
        tmp_int.Interpolation_method, tmp_int.Extrapolation_method);
    psf_est.interpolation{iter_psf} = tmp_int;    
end
psf_est.filepath = './Radius_estimation/psf_simulation_interpolation.mat';
save(psf_est.filepath, '-struct', 'psf_est');
%% Singl e point spread function
% %%  Nromalized edge intensity vs Radius
% % Assuming the oritentation of vessels are uniformaly distributed in 4pi,
% % the average normalized minimum edge intensity for different size of the
% % vessels are: 
close all;
plot_data = psf_simulation{1};
image_file_name = sprintf('Normalized_edge_intensity_psf%.3d%.3d.png',round(plot_data.psf_size(1)*100), ...
    round(plot_data.psf_size(3)*100));
avg_normalized_min_edge_int = mean(plot_data.normalized_min_edge_int, 1);

tmp_plot_theat_idx_list = [1, 7, 10, 13, 19];
tmp_ax = axes();
for theta_idx = tmp_plot_theat_idx_list
    tmp_fig = plot(tmp_ax, plot_data.vessel_radius_list, plot_data.normalized_min_edge_int(theta_idx,:), 'LineWidth', 2);
    hold on
end
plot(plot_data.vessel_radius_list, avg_normalized_min_edge_int, '-.', 'LineWidth', 2);

% tmp_ldg = legend(cellstr(num2str(theta_list(tmp_plot_theat_idx_list)'*180/pi, '%.2f')));
tmp_ldg = legend(tmp_ax, '0', '30', '45', '60', '90', 'Average');
title(tmp_ldg, 'Tilt Angle/{\circ}');
pbaspect(tmp_ax, [1,1,1]);
tmp_ax.LineWidth = 2;
tmp_ax.FontSize = 15;
grid(tmp_ax, 'on')
xlabel(tmp_ax, 'Radius/\mum');
ylabel(tmp_ax, 'Normalized Minimum Edge Intensity');
% xticks(tmp_ax, 0:3:max(test_radius_est_list));
[tmp_fig.LineWidth] = deal(1.5);
outerpos = tmp_ax.OuterPosition;
ti = tmp_ax.TightInset;
tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
    outerpos(4)-ti(2)-ti(4)];

print(fullfile(save_folder, 'Relative_error_psf_shape_est030320_gt030340.png'), '-dpng', '-r300');


print(DataManager.fp_analysis_image(image_file_name, save_folder), '-dpng', '-r300');
%% Relative error if use average normalized edge intensity
% plot_data = psf_simulation{4};
% image_file_name = sprintf('Relative_error_estimated_by_avg_int_psf%.3d%.3d.png',round(plot_data.psf_size(1)*100), ...
%     round(plot_data.psf_size(3)*100));
% 
% avg_min_edge_int = mean(plot_data.abs_min_edge_int, 1);
% num_radius = length(plot_data.vessel_radius_list);
% num_theta = length(plot_data.theta_list);
% tmp_relative_error = zeros(num_radius,num_theta);
% for tmp_theta_idx = 1 : num_theta
%     for tmp_r_idx = 1 : num_radius
%         int_profile = plot_data.radial_int_distribution{tmp_theta_idx, tmp_r_idx};
%         r_list = int_profile(:,1);
%         int_list = int_profile(:,2);
%         [int_list, tmp_unique_idx, ~] = unique(int_list);
%         radius_with_avg_int = interp1(int_list, r_list(tmp_unique_idx), avg_min_edge_int(tmp_r_idx));
%         tmp_relative_error(tmp_r_idx, tmp_theta_idx) = (radius_with_avg_int - plot_data.vessel_radius_list(tmp_r_idx))/plot_data.vessel_radius_list(tmp_r_idx);
%     end
% end
% tmp_ax = axes();
% tmp_fig = plot(tmp_ax, plot_data.vessel_radius_list, tmp_relative_error(:,1:2:end), 'LineWidth', 2);
% xlabel(tmp_ax, 'Radius/\mum');
% ylabel(tmp_ax, 'Relative error');
% grid(tmp_ax, 'on')
% [tmp_fig.LineWidth] = deal(1.5);
% tmp_ax.LineWidth = 2;
% tmp_ax.FontSize = 15;
% pbaspect(tmp_ax, [1,1,1]);
% outerpos = tmp_ax.OuterPosition;
% ti = tmp_ax.TightInset;
% tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
%     outerpos(4)-ti(2)-ti(4)];
% tmp_ldg = legend(tmp_ax, cellstr(num2str(plot_data.theta_list(1:2:end)'*180/pi, '%.0f')));
% title(tmp_ldg, 'Tilt Angle/{\circ}');
% print(DataManager.fp_analysis_image(image_file_name, save_folder), '-dpng', '-r300');
%% Fixed Theta, Radius error improvement for vessel of different radius
% Conclusion: Worthwhile? 
% close all;
% plot_data = psf_simulation{4};
% theta = pi/2;
% image_file_name = sprintf('Relative_error_improvement_vs_initial_error_psf%.3d%.3d_theta%d.png',round(plot_data.psf_size(1)*100), ...
%     round(plot_data.psf_size(3)*100), round(theta*180/pi));
% radius_gt_list = [1, 2, 3, 4, 5];
% num_gt = length(radius_gt_list);
% error_result = cell(num_gt,1);
% error_initial = cell(num_gt,1);
% for tmp_gt_idx = 1 : num_gt
%     radius_gt = radius_gt_list(tmp_gt_idx);
%     radius_est_list = radius_gt * 0.1 :0.1: min(10, radius_gt*2);
%     num_est_list = numel(radius_est_list);
%     tmp_error_result = zeros(num_est_list,1);
%     tmp_error_initial = zeros(num_est_list,1);
%     for tmp_idx = 1 : num_est_list
%         radius_est = radius_est_list(tmp_idx);
%     %     radius_est = 3.5;
%         radius_gt_idx = find(plot_data.vessel_radius_list >= radius_gt, 1);
%         theta_idx = find(plot_data.theta_list >= theta, 1);
% 
%         int_th_est = interp1(plot_data.vessel_radius_list, plot_data.abs_min_edge_int(theta_idx,:), radius_est);
%         int_profile = plot_data.radial_int_distribution{theta_idx, radius_gt_idx};
%         int_dist = int_profile(:,2);
%         int_r = int_profile(:,1);
%         [int_dist, tmp_unique_idx,~] = unique(int_dist);
% 
%         radius = interp1(int_dist, int_r(tmp_unique_idx), int_th_est);
%         tmp_error_result(tmp_idx) = (radius - radius_gt) / radius_gt;
%         tmp_error_initial(tmp_idx) = (radius_est - radius_gt) / radius_gt;
%     end
%     error_result{tmp_gt_idx} = tmp_error_result;
%     error_initial{tmp_gt_idx} = tmp_error_initial;
% end
% 
% tmp_ax = axes();
% for tmp_gt_idx = 1 : num_gt
%     tmp_fig = plot(tmp_ax, error_initial{tmp_gt_idx}, error_result{tmp_gt_idx}, 'LineWidth', 2);
%     hold on
% end
% xlabel(tmp_ax, 'Initial relative error');
% ylabel(tmp_ax, 'Estimation relative error');
% % title(sprintf('Vessel tilt angle %.2d', theta));
% tmp_lgd = legend(tmp_ax, cellstr(num2str(radius_gt_list', '%.1f')));
% title(tmp_lgd, 'Radius/\mum');
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
%% Edge intensity vs. radius. 
% 1. For different point spread function
% 2. For different theta 
% 3. Question to answer: 
%% Error from wrong estimation of the point spread function
% Conclusion: 
% 1. Among different tilt angle, theta=0 is the most sensitive to the
% expansion of the point spread function ( anisotropic expansion ) -
% Estimated the maximum estimation error for different size of the point
% spread function. 
% 2. For 1.6x lateral expansion, the relative estimation of the vessel
% radius is less than 10% for vessel of radius larger than 1.75um. The
% relative error for vessel of radius 1um is about 25%. 
% 
% simulation_gt = psf_simulation{4};
% simulation_est = psf_simulation{1};
% image_file_name = sprintf('Relative_error_psf_shape_est%.3d%.3d_gt%.3d%.3d.png',round(simulation_est.psf_size(1)*100), round(simulation_est.psf_size(3)*100), ...
%    round(simulation_gt.psf_size(1)*100), round(simulation_gt.psf_size(3)*100));
% test_theta_list = intersect(simulation_gt.theta_list, simulation_est.theta_list);
% test_radius_list = intersect(simulation_gt.vessel_radius_list, simulation_est.vessel_radius_list);
% num_test_theta = numel(test_theta_list);
% num_test_radius = numel(test_radius_list);
% relative_error = zeros(num_test_radius, num_test_theta);
% test_radius_est_list = zeros(num_test_radius, num_test_theta);
% for tmp_idx1 = 1 : num_test_theta
%     theta = test_theta_list(tmp_idx1);
% %     theta = 0;
%     for tmp_idx2 = 1 : num_test_radius
%         radius_gt = test_radius_list(tmp_idx2);
% %         radius_gt = 5;
%         [tmp_1, tmp_2] = ndgrid(simulation_est.theta_list, simulation_est.vessel_radius_list);
%         itp_normalized_min_edge_int = griddedInterpolant(repmat(simulation_est.theta_list', 1, numel(simulation_est.vessel_radius_list)), ...
%             repmat(simulation_est.vessel_radius_list, numel(simulation_est.theta_list),1),...
%             simulation_est.normalized_min_edge_int);
%         % Get the estimated normalized edge intensity from the estimated PSF
%         % profile
%         est_normalized_edge_int = itp_normalized_min_edge_int(theta, radius_gt);
%         % Use the found edge intensity to find the real radius corresponding to the
%         % estimated intensity
%         [~, radius_gt_idx] = min(abs(simulation_gt.vessel_radius_list - radius_gt));
%         [~, theta_idx] = min(abs(simulation_gt.theta_list - theta));
%         radial_profile_gt = simulation_gt.radial_int_distribution{theta_idx, radius_gt_idx};
%         [radial_int, unique_int_idx] = unique(radial_profile_gt(:,2));
%         radius_est = interp1(radial_int, radial_profile_gt(unique_int_idx,1), est_normalized_edge_int);
%         test_radius_est_list(tmp_idx2, tmp_idx1) = radius_est;
%         relative_error(tmp_idx2, tmp_idx1) = (radius_est - radius_gt)/radius_gt;
%     end
% end
% close all;
% tmp_ax = axes();
% tmp_fig = plot(tmp_ax, test_radius_list, relative_error(:, 1:6:end));
% tmp_ax.LineWidth = 2;
% grid(tmp_ax, 'on')
% tmp_ax.FontSize = 15;
% xlabel(tmp_ax, 'Radius/\mum');
% ylabel(tmp_ax, 'Relative error');
% xticks(tmp_ax, 0:3:max(test_radius_est_list));
% pbaspect(tmp_ax, [1,1,1]);
% [tmp_fig.LineWidth] = deal(1.5);
% tmp_legent = legend(tmp_ax, cellstr(num2str(test_theta_list(1:6:end)'*180/pi, '%.1f')));
% title(tmp_legent, 'Tilt angle/{\circ}');
% outerpos = tmp_ax.OuterPosition;
% ti = tmp_ax.TightInset;
% tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
%     outerpos(4)-ti(2)-ti(4)];
% print(DataManager.fp_analysis_image(image_file_name, save_folder), '-dpng', '-r300');


% close all;
% tmp_ax = axes();
% tmp_fig = plot(tmp_ax, test_radius_list, test_radius_est_list(:, 1:6:end));
% tmp_ax.LineWidth = 2;
% grid(tmp_ax, 'on')
% tmp_ax.FontSize = 12;
% xlabel(tmp_ax, 'Radius/\mum');
% ylabel(tmp_ax, 'Estimated radius/\mum');
% pbaspect(tmp_ax, [1,1,1]);
% [tmp_fig.LineWidth] = deal(1.5);
% tmp_legent = legend(tmp_ax, cellstr(num2str(test_theta_list(1:6:end)'*180/pi, '%.1f')));
% title(tmp_legent, 'Tilt angle/{\circ}');
% outerpos = tmp_ax.OuterPosition;
% ti = tmp_ax.TightInset;
% tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
%     outerpos(4)-ti(2)-ti(4)];
% print(fullfile(save_folder, 'Radius_est_vs_gt_psf_shape_est030320_gt060640.png'), '-dpng', '-r300');
%% Error from wrong estimation of the tilt angle, fix the radius
close all;
plot_data = psf_simulation{1};
radius_gt = 5;
image_file_name = sprintf('Relative_error_vs_tilt_angle_error_psf%.3d%.3d_radius%d.png',round(plot_data.psf_size(1)*100), ...
    round(plot_data.psf_size(3)*100), round(radius_gt*10));
[~, radius_gt_idx] = min(abs(plot_data.vessel_radius_list - radius_gt));
theta_gt_list = plot_data.theta_list(1:6:end);
num_gt = length(theta_gt_list);
error_r = cell(num_gt,1);
error_theta = cell(num_gt,1);
for tmp_gt_idx = 1 : num_gt
    theta_gt = theta_gt_list(tmp_gt_idx);
    [~, theta_gt_idx] = min(abs(plot_data.theta_list - theta_gt));
    theta_est_list = unique(min(plot_data.theta_list(end), max(0, theta_gt + (-45:5:45) * pi/180)));
    num_est_list = numel(theta_est_list);
    tmp_error_r_result = zeros(num_est_list,1);
    tmp_error_theta = zeros(num_est_list,1);
    for tmp_idx = 1 : num_est_list
        theta_est = theta_est_list(tmp_idx);
        int_th_est = interp1(plot_data.theta_list, plot_data.abs_min_edge_int(:,radius_gt_idx), theta_est);
        int_profile = plot_data.radial_int_distribution{theta_gt_idx, radius_gt_idx};
        int_dist = int_profile(:,2);
        int_r = int_profile(:,1);
        [int_dist, tmp_unique_idx,~] = unique(int_dist);
        radius = interp1(int_dist, int_r(tmp_unique_idx), int_th_est);
        tmp_error_r_result(tmp_idx) = (radius - radius_gt) / radius_gt;
        tmp_error_theta(tmp_idx) = (theta_est - theta_gt);
    end
    error_r{tmp_gt_idx} = tmp_error_r_result;
    error_theta{tmp_gt_idx} = round(tmp_error_theta * 180/pi);
end

tmp_ax = axes();
for tmp_gt_idx = 1 : num_gt
    tmp_fig = plot(tmp_ax, error_theta{tmp_gt_idx}, error_r{tmp_gt_idx}, 'LineWidth', 2);
    hold on
end
xlabel(tmp_ax, 'Error in tilt angle/{\circ}');
ylabel(tmp_ax, 'Relative error in estimated r');
% title(sprintf('Vessel tilt angle %.2d', theta));
tmp_lgd = legend(tmp_ax, cellstr(num2str(round(theta_gt_list * 180/pi)' , '%.2d')));
title(tmp_lgd, 'Tilt angle/{\circ}');

grid(tmp_ax, 'on')
[tmp_fig.LineWidth] = deal(1.5);
tmp_ax.LineWidth = 2;
tmp_ax.FontSize = 15;
pbaspect(tmp_ax, [1,1,1]);
outerpos = tmp_ax.OuterPosition;
ti = tmp_ax.TightInset;
tmp_ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), ...
    outerpos(4)-ti(2)-ti(4)];
print(DataManager.fp_analysis_image(image_file_name, save_folder), '-dpng', '-r300');

%% Visualization