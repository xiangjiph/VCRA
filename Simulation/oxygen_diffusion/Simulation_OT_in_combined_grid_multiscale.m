set_env;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
load_skl_name = '240_cube_rec';
combined_grid_name = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, combined_grid_name);
grid_info = grid_c_info.grid_ori;
output_graph_name = 'OT_simulation';
data_folder_name = 'pO2';
overwrite_Q = false;
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Tissue_pO2_simulation');
% Might be worth to apply a minimum capillary radius for: 
% 1. Downsampling accuracy
% 2. Reflect the real capillary radius in vivo 
%% Input parameters 
min_cap_r = 0; % um
recon_max_error_rate = 0.1;
int_expand_half_length = 80; % before downsampling
downsample_rate = 2;
lm_wz_list_um = [8 : 4 : 20, 30 : 10 : 60, 80 : 20 : 120];

inhomogeneous_term = 1;
gpuDevice_label = 2;
%%
lm_wz_list_pxl = round(lm_wz_list_um ./ downsample_rate);
int_expand_half_length_ds = int_expand_half_length ./ downsample_rate;

gpuDevice(gpuDevice_label);
%% Test input
grid_c_label = 1262;
grid_c_sub = grid_c_info.bbox_grid_sub_list(grid_c_label, :);
grid_c_ind = sub2ind(grid_c_info.grid_size, grid_c_sub(1), ...
    grid_c_sub(2), grid_c_sub(3));
assert(grid_c_info.bbox_xyz_label_array(grid_c_ind) == grid_c_label, 'Incorrect grid ind');
%%
[grid_c_skl_ind, grid_c_r, int_bbox_local_bbox_mmxx, int_bbox_exp_size] = ...
    fun_simulation_OT_load_skeleton_in_grid(grid_c_info, load_skl_name, grid_c_ind, int_expand_half_length);
%% Convert to graph
vessel_graph = fun_skeleton_to_graph(grid_c_skl_ind, int_bbox_exp_size);
vessel_graph.radius = sparse(grid_c_skl_ind, ones(vessel_graph.num.skeleton_voxel, 1), ...
    max(min_cap_r, double(grid_c_r)), vessel_graph.num.block_voxel, 1);
%% Reconstruction
vessel_recon_mask = fun_graph_to_reconstruction_mask(vessel_graph, false, recon_max_error_rate);

vessel_recon_mask_dt = bwdist(vessel_recon_mask);

%% Solve Poisson equation
pO2_n_result_ds2_nn = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    downsample_rate, inhomogeneous_term, 'nearest', true);

pO2_n_result_ds2_ln = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    downsample_rate, inhomogeneous_term, 'linear', true);
pO2_n_result_ds2_ln.pO2_array_us = imresize3(pO2_n_result_ds2_ln.pO2_array, int_bbox_exp_size, ...
    'Method', 'linear');

pO2_n_result_nn = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    1, inhomogeneous_term, 'nearest', true);
pO2_n_result_nn.dist_to_vsl = bwdist(pO2_n_result_nn.vessel_mask);

pO2_n_result_ln = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    1, inhomogeneous_term, 'linear', true);
%% Debug: why two pO2 calculation has systematic error about a few precent? 
fun_error_n = @(x1, x2)2 * (x1 - x2) ./ (x1 + x2);

pO2_n_result_nn_2 = fun_simulation_OT_solve_ct_diff_itr_cvg_in_roi(pO2_n_result_nn.vessel_mask, ...
    pO2_n_result_nn.inhomogeneous_term, pO2_n_result_nn_2.pO2_array, [], true);

pO2_n_result_ln_2 = fun_simulation_OT_solve_ct_diff_itr_cvg_in_roi(pO2_n_result_ln.vessel_mask, ...
    pO2_n_result_ln.inhomogeneous_term, pO2_n_result_ln_2.pO2_array, [], true);

figure;histogram(fun_error_n(pO2_n_result_nn_2.pO2_array,pO2_n_result_nn.pO2_array), 100, 'Normalization', 'pdf');
figure;histogram(fun_error_n(pO2_n_result_ln_2.pO2_array, pO2_n_result_ln.pO2_array), 100, 'Normalization', 'pdf');
figure;histogram(fun_error_n(pO2_n_result_ln_2.pO2_array, pO2_n_result_nn_2.pO2_array), 100, 'Normalization', 'pdf');
figure;histogram(fun_error_n(pO2_n_result_ln.pO2_array, pO2_n_result_nn.pO2_array), 100, 'Normalization', 'pdf');
% Conclusion: 
% 1. NN interpolation approaches convergence from 0+ 
% 2. Linear interpolation approaches convergence from 0-
% 3. Given maximum update 1e-5, linear interpolation is on average 1.76%
% lower than 20000 iteration result; NN interpolation is on average 2.10%
% higher than 20000 iteration result; linear interpolation is about 4.01%
% lower than NN interpolation
% 4. Given target maximum update 1e-6, neither of the method converged
% within 20000 iterations - probably the lower resolution result were not
% converged either. However, linear interpolation is on overage 0.16%
% smaller than nearest neighbor interpolation

%% Analyze DT error due to downsampling
tmp_x = vessel_recon_mask_dt(~vessel_recon_mask);
tmp_x_stat = fun_analysis_get_basic_statistics(tmp_x);
tmp_bin_edge = 0.5 : 1 : 30;

ds_test_cell = cell(2, 0);
ds_test_cell(:, end+1) = {imresize3(single(vessel_recon_mask), 0.5) > 0.5, 'single cubic 0.5'};
ds_test_cell(:, end+1) = {imresize3(single(vessel_recon_mask), 0.5) > 0, 'single cubic 0'};
ds_test_cell(:, end+1) = {imresize3(uint8(vessel_recon_mask), 0.5) > 0.5, 'uint8 cubic 0.5'};

ds_test_cell(:, end+1) = {imresize3(single(vessel_recon_mask), 0.5, 'Method', 'linear') > 0, 'single linear 0'};
ds_test_cell(:, end+1) = {imresize3(single(vessel_recon_mask), 0.5, 'Method', 'linear') > 0.5, 'single linear 0.5'};
ds_test_cell(:, end+1) = {imresize3(uint8(vessel_recon_mask), 0.5, 'Method', 'linear') > 0.5, 'uint8 linear 0.5'};

ds_test_cell(:, end+1) = {imresize3(uint8(vessel_recon_mask), 0.5, 'Method', 'nearest') > 0.5, 'uint8 nearest 0.5'};

fig_hdl = figure;
num_plot = size(ds_test_cell, 2) + 1;
for iter_plot = 1 : num_plot
    ax_hdl = subplot(ceil(num_plot / 2), 2, iter_plot);
    if iter_plot == 1
        tmp_data = tmp_x;
        tmp_title = '1 \mum';
    else
        tmp_mask = ds_test_cell{1, iter_plot - 1};
        tmp_title = ds_test_cell{2, iter_plot - 1};
        tmp_data = bwdist(tmp_mask) .* 2;
        tmp_data = tmp_data(~tmp_mask);
    end
    tmp_stat = fun_analysis_get_basic_statistics(tmp_data);
    tmp_stat_string = fun_analysis_basic_stat_str_to_string(tmp_stat);
    histogram(ax_hdl, tmp_data, tmp_bin_edge, 'Normalization', 'pdf');
    legend(ax_hdl, tmp_stat_string, 'Location', 'northwest');
    ax_hdl.XLabel.Interpreter = 'latex';
    ax_hdl.XLabel.String = '$d(\mathbf{x})$';
    ax_hdl.YLabel.String = 'PDF';
    ax_hdl.Title.String = tmp_title;
end

fig_fp = fullfile(vis_folder, sprintf('%s_%s_%s_%d_dt_pdf_vs_ds_method.png', dataset_name, ...
        stack, combined_grid_name, grid_c_label));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
ds_err = pO2_n_result_ln.pO2_array - pO2_n_result_nn.pO2_array;
ds_err_n = 2 * ds_err ./ (pO2_n_result_ln.pO2_array + pO2_n_result_nn.pO2_array);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, ds_err_n, 100, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';

all((pO2_n_result_ln.pO2_array == 0) == (pO2_n_result_nn.pO2_array == 0), 'all')
%% Histogram of relative error: 1 um vs 2 um voxel resolution
ds_err = pO2_n_result_ds2_ln.pO2_array_us - pO2_n_result_nn.pO2_array;
ds_err_n = 2 * ds_err ./ (pO2_n_result_nn.pO2_array + pO2_n_result_ds2_ln.pO2_array_us);
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, ds_err_n, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
    %%
dt_bin_um = 0 : 2 : 60;
y_in_x_str = fun_analysis_get_y_stat_in_x_bin(pO2_n_result_nn.dist_to_vsl, ...
    abs(ds_err_n), dt_bin_um);
[bin_count, x_bin, y_bin] = histcounts2(pO2_n_result_nn.dist_to_vsl, abs(ds_err_n));
    %% Plot
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    histogram2(ax_hdl, 'XBinEdges', x_bin, 'YBinEdges', y_bin, 'BinCounts', bin_count, 'DisplayStyle', 'tile');
    ax_hdl.ColorScale = 'log';
    ax_hdl.YScale = 'log';
    ax_hdl.YLabel.String = 'Absolute Relative error';
    ax_hdl.XLabel.String = 'Distance to vessel (\mum)';
    cbar_hdl = colorbar(ax_hdl);
    cbar_hdl.Label.String = 'Number of voxels';
    hold(ax_hdl, 'on');
    plt_med_hdl = plot(ax_hdl, y_in_x_str.x_bin_val, y_in_x_str.y_median);
    plt_med_hdl.LineWidth = 1.5;
    plt_25_hdl = plot(ax_hdl, y_in_x_str.x_bin_val, y_in_x_str.y_prctile(:, 3), '-.', 'LineWidth', 1);
    plt_75_hdl = plot(ax_hdl, y_in_x_str.x_bin_val, y_in_x_str.y_prctile(:, 5), '-.', 'LineWidth', 1);
    leg_hdl = legend(ax_hdl, [plt_med_hdl, plt_25_hdl, plt_75_hdl], 'Median', '25%', '75%');
    fig_fp = fullfile(vis_folder, sprintf('%s_%s_%s_%d_ds_err_vs_dist_linear_itp.png', dataset_name, ...
        stack, combined_grid_name, grid_c_label));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Visualize section
vis_sec = 200;
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
ax_hdl_1 = subplot(1,3,1);
imagesc(ax_hdl_1, pO2_n_result_ln.pO2_array(:, :, vis_sec));
ax_hdl_1.DataAspectRatio = [1,1,1];
colorbar(ax_hdl_1);
ax_hdl_1.Colormap(end, :) = [1,0,0];
% ax_hdl_1.Title.String = '2 \mum resolution (linear upsampled)';
ax_hdl_1.XLabel.String = 'X (\mum)';
ax_hdl_1.YLabel.String = 'Y (\mum)';

ax_hdl_2 = subplot(1,3,2);
imagesc(ax_hdl_2, pO2_n_result_nn.pO2_array(:, :, vis_sec));
ax_hdl_2.DataAspectRatio = [1,1,1];
colorbar(ax_hdl_2);
ax_hdl_2.Colormap(end, :) = [1,0,0];
% ax_hdl_2.Title.String = '1 \mum resolution';
ax_hdl_2.XLabel.String = 'X (\mum)';
ax_hdl_2.YLabel.String = 'Y (\mum)';

ax_hdl_3 = subplot(1,3,3);
imagesc(ax_hdl_3, pO2_n_result_nn.pO2_array(:, :, vis_sec) - pO2_n_result_ln.pO2_array(:, :, vis_sec));
ax_hdl_3.DataAspectRatio = [1,1,1];
colorbar(ax_hdl_3);
ax_hdl_3.Colormap(end, :) = [1,0,0];
% ax_hdl_3.Title.String = '1 \mum resolution - 2 \mum resolution';
ax_hdl_3.XLabel.String = 'X (\mum)';
ax_hdl_3.YLabel.String = 'Y (\mum)';
% fig_fp = fullfile(vis_folder, sprintf('%s_%s_%s_%d_diff_1_vs_2_um_vxl_size_sec_%d_linear_itp.png', dataset_name, ...
%     stack, combined_grid_name, grid_c_label, vis_sec));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Local distance transform properties
% local maxima found by larger window size are a subset of local maxima
% found by smaller window size
mask_size_ds = size(pO2_n_result_ds2_ln.vessel_mask);
valid_mask_bbox_mmxx = [repelem(int_expand_half_length_ds, 1, 3) + 1, ...
    mask_size_ds - int_expand_half_length_ds];
valid_extrema_mask = false(mask_size_ds);
valid_extrema_mask(valid_mask_bbox_mmxx(1) : valid_mask_bbox_mmxx(4), ...
    valid_mask_bbox_mmxx(2) : valid_mask_bbox_mmxx(5), ...
    valid_mask_bbox_mmxx(3) : valid_mask_bbox_mmxx(6)) = true;
valid_extrema_mask = valid_extrema_mask & ~pO2_n_result_ds2_ln.vessel_mask;

tic_search_dtlm = tic;
mask_rz_dt = bwdist(pO2_n_result_ds2_ln.vessel_mask) .* downsample_rate;

[dt_max_info, pO2_min_info] = fun_simulation_OT_SA_wdw_sz(mask_rz_dt, pO2_n_result_ds2_ln.pO2_array, ...
    lm_wz_list_pxl, valid_extrema_mask);
fprintf('Finish searching for local extrema. Elapsed time is %f secodns.\n', toc(tic_search_dtlm));
%%
mask_size = size(pO2_n_result_nn.vessel_mask);
valid_mask_bbox_mmxx = [repelem(int_expand_half_length, 1, 3) + 1, ...
    mask_size - int_expand_half_length];
valid_extrema_mask = false(mask_size);
valid_extrema_mask(valid_mask_bbox_mmxx(1) : valid_mask_bbox_mmxx(4), ...
    valid_mask_bbox_mmxx(2) : valid_mask_bbox_mmxx(5), ...
    valid_mask_bbox_mmxx(3) : valid_mask_bbox_mmxx(6)) = true;
valid_extrema_mask = valid_extrema_mask & ~pO2_n_result_nn.vessel_mask;

[dt_max_info_0, pO2_min_info_0] = fun_simulation_OT_SA_wdw_sz(pO2_n_result_nn.dist_to_vsl,...
    pO2_n_result_nn.pO2_array, lm_wz_list_um, valid_extrema_mask);
%%
vis_group = 7;
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax_hdl_1 = subplot(1,2,1);
tmp_data_1 = dt_max_info.dt{vis_group};
tmp_data_2 = dt_max_info_0.dt{vis_group};

tmp_stat_str_1 = sprintf('Voxel size: 2\\mum\n%s', fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data_1)));
tmp_stat_str_2 = sprintf('Voxel size: 1\\mum\n%s', fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data_2)));

histogram(ax_hdl_1, tmp_data_1, 'Normalization', 'count', 'FaceAlpha', 0.5);
hold(ax_hdl_1, 'on');
histogram(ax_hdl_1, tmp_data_2, 'Normalization', 'count', 'FaceAlpha', 0.5);
legend(ax_hdl_1, tmp_stat_str_1, tmp_stat_str_2);
ax_hdl_1.XLabel.String = 'd_{lm} (\mum)';
ax_hdl_1.YLabel.String = 'Counts';
ax_hdl_1.Title.String = sprintf('Window size: %d \\mum', lm_wz_list_um(vis_group));

ax_hdl_2 = subplot(1,2,2);
tmp_data_1 = pO2_min_info.pO2{vis_group};
tmp_data_2 = pO2_min_info_0.pO2{vis_group};
tmp_stat_str_1 = sprintf('Voxel size: 2\\mum\n%s', fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data_1)));
tmp_stat_str_2 = sprintf('Voxel size: 1\\mum\n%s', fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data_2)));
histogram(ax_hdl_2, tmp_data_1, 'Normalization', 'count', 'FaceAlpha', 0.5);
hold(ax_hdl_2, 'on');
histogram(ax_hdl_2, tmp_data_2, 'Normalization', 'count', 'FaceAlpha', 0.5);
legend(ax_hdl_2, tmp_stat_str_1, tmp_stat_str_2, 'Location', 'northwest');
ax_hdl_2.XLabel.String = 'u_{lm}';
ax_hdl_2.YLabel.String = 'Counts';
ax_hdl_2.Title.String = sprintf('Window size: %d \\mum', lm_wz_list_um(vis_group));
fig_fp = fullfile(vis_folder, sprintf('%s_%s_%s_%d_lm_cmp_wz_%d_um_linear_itp.pmg', dataset_name, ...
    stack, combined_grid_name, grid_c_label, lm_wz_list_um(vis_group)));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Visualize
% tmp_x = lm_wz_list * 2;
% tmp_y = cellfun(@(x) mean(x, 'omitnan'), dt_max_info.dt_v);
% fig_hdl = figure;
% ax_hdl_1 = axes(fig_hdl);
% plot(ax_hdl_1, tmp_x, tmp_y, 'LineWidth', 1.5);
% hold(ax_hdl_1, 'on');
% yyaxis(ax_hdl_1, 'right');
% tmp_y = cellfun(@(x) mean(x, 'omitnan'), pO2_min_info.pO2_v);
% plot(ax_hdl_1, tmp_x, tmp_y, 'LineWidth', 1.5);
% %% Analyze each internal grid
% int_bbox_local_bbox_ds = ceil(int_bbox_local_bbox_mmxx ./ downsample_rate);
% int_bbox_grid_label = grid_c_info.internal_subgrid_label{grid_c_ind};
% num_int_cube = numel(int_bbox_grid_label);
% for iter_int_cube = 1 : num_int_cube
%     tmp_cube_label = int_bbox_grid_label(iter_int_cube);
%     tmp_cube_str = fun_grid_get_single_cube_info(grid_info, tmp_cube_label);
%     tmp_local_bbox_mmxx = int_bbox_local_bbox_ds(iter_int_cube, :);
%     tmp_local_bbox_mmll = [tmp_local_bbox_mmxx(1:3), tmp_local_bbox_mmxx(4:6) - tmp_local_bbox_mmxx(1:3) + 1];
%     % Crop the local dt and pO2 array
%     tmp_int_dt = crop_bbox3(mask_rz_dt, tmp_local_bbox_mmll);
%     tmp_int_pO2 = crop_bbox3(pO2_n_result_ds2.pO2_array, tmp_local_bbox_mmll);
%     % Select the voxel outside the vessel mask 
%     tmp_int_tissue_mask = (tmp_int_dt ~= 0);
%     tmp_int_dt = tmp_int_dt(tmp_int_tissue_mask);
%     tmp_int_pO2 = tmp_int_pO2(tmp_int_tissue_mask);
%     %% Regional statistics - seperate stat
%     tmp_cube_str.local_dt_stat = fun_analysis_get_basic_statistics(tmp_int_dt, true);
%     tmp_cube_str.local_dt_stat.prtl2val_itp = fun_analysis_get_prtl_to_value_interpolation(tmp_cube_str.local_dt_stat);
%     tmp_cube_str.local_dt_stat.val2ptrl_itp = fun_analysis_get_value_to_prtl_interpolation(tmp_cube_str.local_dt_stat);
%     
%     tmp_cube_str.local_pO2_stat = fun_analysis_get_basic_statistics(tmp_int_pO2, true);    
%     tmp_cube_str.local_pO2_stat.prtl2val_itp = fun_analysis_get_prtl_to_value_interpolation(tmp_cube_str.local_pO2_stat);
%     tmp_cube_str.local_pO2_stat.val2ptrl_itp = fun_analysis_get_value_to_prtl_interpolation(tmp_cube_str.local_pO2_stat);
%     %% Regional statistics - joint stat
%     tmp_dt_edge = 0.5 : tmp_cube_str.local_dt_stat.max;
%     tmp_cube_str.pO2_stat_in_dt_bin = fun_analysis_get_y_stat_in_x_bin(tmp_int_dt, ...
%         tmp_int_pO2, tmp_dt_edge);
%     %% 2D histogram - count - normalization can be done later if needed
%     [tmp_cube_str.pO2_dt_hist2.count, tmp_cube_str.pO2_dt_hist2.dt_edge, ...
%         tmp_cube_str.pO2_dt_hist2.pO2_edge] = histcounts2(tmp_int_dt, tmp_int_pO2);
%     %% Local DT maxima inside the cube
%     tmp_cube_str.dt_lm = fun_simulation_OT_SA_get_extrema_in_bbox(dt_max_info, ...
%         tmp_local_bbox_mmxx, downsample_rate);
%     %% Local oxygen minimum inside the cube
%     tmp_cube_str.pO2_lm = fun_simulation_OT_SA_get_extrema_in_bbox(pO2_min_info, ...
%         tmp_local_bbox_mmxx, downsample_rate);
%     %% Save result
%     exit_code = DataManager.write_analysis_data_in_grid(tmp_cube_str, ...
%         data_folder_name, dataset_name, stack, tmp_cube_str.grid_version, ...
%         tmp_cube_str.grid_label);
% end
