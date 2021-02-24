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
gpuDevice(2);
%%
vis_cube_label = 31998;
vis_cube_grid_sub = grid_info.bbox_grid_sub_list(vis_cube_label, :);
vis_recon_str = DataManager.load_block_mask(dataset_name, stack, ...
    mask_version, vis_cube_grid_sub(1), vis_cube_grid_sub(2), vis_cube_grid_sub(3));
vis_recon = fun_reconstruct_block_mask(vis_recon_str);
volumeViewer(vis_recon);
%%
oxygen_metabolism_rate = 90e-6 * 1.5;
pO2_vsl = 36;
% pO2_result_nds = fun_simulation_OT_solve_ct_diff_itr(vis_recon, pO2_vsl, oxygen_metabolism_rate, 1, []);
% pO2_result_nds_i20 = fun_simulation_OT_solve_ct_diff_itr(vis_recon, pO2_vsl, oxygen_metabolism_rate, 1, 20);
pO2_result = fun_simulation_OT_solve_ct_diff_itr(vis_recon, pO2_vsl, oxygen_metabolism_rate, 1, []);
% pO2_result_cpO2_i = fun_simulation_OT_solve_ct_diff_itr(vis_recon, pO2_vsl, oxygen_metabolism_rate, 2, 10);
% pO2_result_cpO2_i20 = fun_simulation_OT_solve_ct_diff_itr(vis_recon, pO2_vsl, oxygen_metabolism_rate, 2, 20);
%%
im_1 = pO2_result_nds.pO2_array;
im_2 = imresize3(pO2_result.pO2_array, 2);
%%
vis_sec = 80;
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
fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
    'Oxygen_transportation_simulation', 'Difference_in_voxel_size_3.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Compare with Krogh model
vis_data_str = pO2_result;
% vis_data_str.pO2_array = imresize3(vis_data_str.pO2_array, vis_data_str.downsample_r);
% vis_data_str.dt_array = imresize3(vis_data_str.dt_array, vis_data_str.downsample_r);
internal_dist = 16;
tmp_bbox = [ones(1, 3) .* internal_dist, size(vis_data_str.pO2_array) - internal_dist];
tmp_bbox(4:6) = tmp_bbox(4:6) - tmp_bbox(1:3) + 1;
tmp_mask = crop_bbox3(vis_data_str.vessel_mask, tmp_bbox);
tmp_x_data = crop_bbox3(vis_data_str.dt_array, tmp_bbox);
tmp_y_data = crop_bbox3(vis_data_str.pO2_array, tmp_bbox);
tmp_x_data = tmp_x_data(~tmp_mask);
tmp_y_data = tmp_y_data(~tmp_mask);

%% Finding the local DT maxima
is_dt_lm_Q = fun_array_local_maximum(vis_data_str.dt_array, 5);
is_dt_lm_Q = is_dt_lm_Q & ~vis_data_str.vessel_mask;
dt_lm_ind = find(is_dt_lm_Q);
dt_lm_sub = fun_ind2sub(size(is_oxy_lm_Q), dt_lm_ind);

is_not_near_boundary_Q = all(dt_lm_sub > min_dist_2_boundary_vxl, 2) & ...
    all(dt_lm_sub < target_array_size - min_dist_2_boundary_vxl, 2);
dt_lm_sub = dt_lm_sub(is_not_near_boundary_Q, :);
dt_lm_ind = dt_lm_ind(is_not_near_boundary_Q);

dt_lm_v = vis_data_str.dt_array(dt_lm_ind);
dt_lm_oxy_v = vis_data_str.pO2_array(dt_lm_ind);
%% Find local oxygen concentration minimum
is_oxy_lm_Q = fun_array_local_maximum(- vis_data_str.pO2_array, 5);
is_oxy_lm_Q = is_oxy_lm_Q & ~vis_data_str.vessel_mask;
oxy_lm_ind = find(is_oxy_lm_Q);
oxy_lm_sub = fun_ind2sub(size(is_oxy_lm_Q), oxy_lm_ind);

min_dist_2_boundary_vxl = 10;
is_not_near_boundary_Q = all(oxy_lm_sub > min_dist_2_boundary_vxl, 2) & ...
    all(oxy_lm_sub < target_array_size - min_dist_2_boundary_vxl, 2);
oxy_lm_sub = oxy_lm_sub(is_not_near_boundary_Q, :);
oxy_lm_ind = oxy_lm_ind(is_not_near_boundary_Q);

oxy_lm_v = vis_data_str.dt_array(oxy_lm_ind);
oxy_lm_dt_v = vis_data_str.pO2_array(oxy_lm_ind);

pdist_oxy_2_dt = pdist2(oxy_lm_sub, dt_lm_sub);
[list_idx_1, list_idx_2, min_dist] = fun_find_col_row_co_minimum(pdist_oxy_2_dt);
%%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
histogram2(ax_hdl, tmp_x_data, tmp_y_data, 'DisplayStyle', 'tile');
hold(ax_hdl, 'on');
% med_hdl = plot(ax_hdl, tmp_y_stat.x_bin_val, tmp_y_stat.y_median, 'LineWidth', 2);
mean_hdl = plot(ax_hdl, tmp_y_stat.x_bin_val, tmp_y_stat.y_mean, 'LineWidth', 2);
r_cap = 2; %\mum
krogh_coeff = oxygen_metabolism_rate /(2 * alphaO2 * diffCoeffO2 * 1e12 );
local_median_dt_max = median(vis_recon_str.link.features.nearest_tissue_dt_max, 'omitnan');
% local_median_dt_max = max(vis_recon_str.link.features.nearest_tissue_dt_max);
% local_median_dt_max = 100;
% local_median_dt_max = 24;
% Determine correction factor by linear regression without intercept
[tmp_y_stat, ~, ~] = fun_analysis_get_y_stat_in_x_bin(tmp_x_data, tmp_y_data, 0 : 1 : local_median_dt_max);

% fit_y = pO2_vsl - tmp_y_stat.y_mean;
% fit_x = krogh_coeff * (local_median_dt_max .^ 2 .* log((tmp_y_stat.x_bin_val + r_cap) ./ r_cap) - ((tmp_y_stat.x_bin_val + r_cap).^ 2 - r_cap^2)/2);

fit_y = pO2_vsl - tmp_y_stat.y_mean;
fit_x = krogh_coeff * ((tmp_y_stat.x_bin_val - local_median_dt_max).^2 - (local_median_dt_max)^2);

% fit_y = pO2_vsl - tmp_y_data;
% fit_x = krogh_coeff * (local_median_dt_max .^ 2 .* log((tmp_x_data + r_cap) ./ r_cap) - ((tmp_x_data + r_cap).^ 2 - r_cap^2)/2);

% fit_x = krogh_coeff * (local_median_dt_max .^ 2 .* (tmp_y_stat.x_bin_val) - ((tmp_y_stat.x_bin_val + r_cap).^ 2 - r_cap^2)/2);
lin_fit_hdl = fitlm(fit_x, fit_y, 'Intercept', false);
correction_factor = lin_fit_hdl.Coefficients.Estimate(1);
% correction_factor = 1;

plot_krogh_x = 0 : tmp_y_stat.x_bin_val(end);
oxygen_krogh = pO2_vsl - correction_factor * krogh_coeff * ((tmp_y_stat.x_bin_val - local_median_dt_max).^2 - (local_median_dt_max)^2);
% oxygen_krogh = pO2_vsl -correction_factor *  krogh_coeff * (local_median_dt_max .^ 2 .* log((plot_krogh_x + r_cap) ./ r_cap) - ((plot_krogh_x + r_cap).^ 2 - r_cap^2)/2);
plt_krogh = plot(ax_hdl, plot_krogh_x, oxygen_krogh, 'LineWidth', 4, 'LineStyle', '-');

% ax_hdl.YLim(1) = 0;

dot_dt_lm_hdl = scatter(ax_hdl, dt_lm_v, dt_lm_oxy_v, 'g+');
dot_oxy_lm_hdl = scatter(ax_hdl, oxy_lm_dt_v, oxy_lm_v, 'r*');
legend(ax_hdl, [dot_dt_lm_hdl, dot_oxy_lm_hdl], {'Local DT maximum', 'Local oxygen minimum'});
%% Solve Poisson equation using central finite difference
% oxygen_metabolism_rate = 90e-6; %mol/(L sec)
% eta = (oxygen_metabolism_rate) * (1/(diffCoeffO2 * 1e12)) * (1/alphaO2); %mmHg/um^2
% 
% cap_pO2 = 40;
% tissue_pO2_0 = 10;
% downsample_ratio = 2;
% eta_ds = eta * (downsample_ratio ^2);
% 
% target_array_size = round(size(vis_recon) ./ downsample_ratio);
% vis_recon_rz = imresize3(uint8(vis_recon), target_array_size) > 0;
% vis_recon_rz_dt = bwdist(vis_recon_rz) .* downsample_ratio;
% vessel_ind = find(vis_recon_rz);
% simu_array_size = size(vis_recon_rz);
% % Initialization 
% simu_array = ones(simu_array_size, 'double') * tissue_pO2_0;
% simu_array(vessel_ind) = cap_pO2;
% simu_array = gpuArray(simu_array);
% tic
% num_max_iter = 5000;
% target_accuracy = 0.0005;
% cum_err = zeros(num_max_iter, 1, 'gpuArray');
% max_err = 1;
% iter = 0;
% while iter < num_max_iter && max_err > target_accuracy
%     iter = iter + 1;
%     tmp_array_0 = simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1));
%     tmp_array_new = - eta_ds/6 + ...
%        (simu_array(1 : (end - 2), 2 : (end - 1), 2 : (end - 1)) + ...
%         simu_array(3 : (end    ), 2 : (end - 1), 2 : (end - 1)) + ...
%         simu_array(2 : (end - 1), 1 : (end - 2), 2 : (end - 1)) + ...
%         simu_array(2 : (end - 1), 3 : (end    ), 2 : (end - 1)) + ...
%         simu_array(2 : (end - 1), 2 : (end - 1), 1 : (end - 2)) + ...
%         simu_array(2 : (end - 1), 2 : (end - 1), 3 : (end    ))) / 6;
%     % Update non-boundary voxels
%     simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1)) = tmp_array_new;
% 
%     % Update boundary voxels
%     simu_array(1, :, :) = simu_array(2, :, :);
%     simu_array(end, :, :) = simu_array(end - 1, :, :);
%     simu_array(:, 1, :) = simu_array(:, 2, :);
%     simu_array(:, end, :) = simu_array(:, end - 1, :);
%     simu_array(:, :, 1) = simu_array(:, :, 2);
%     simu_array(:, :, end) = simu_array(:, :, end - 1);
%     % Fix the oxygen partial pressure on the boundary
%     simu_array(vessel_ind) = cap_pO2;
%     % Compute residual error    
%     max_err = max(abs(tmp_array_0 - tmp_array_new), [], 'all');
%     cum_err(iter) = max_err;
%     fprintf('Finish computing iteration %d. Maximum update is %f.\n', ...
%         iter, max_err);
% end
% toc
% simu_array = gather(simu_array);
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